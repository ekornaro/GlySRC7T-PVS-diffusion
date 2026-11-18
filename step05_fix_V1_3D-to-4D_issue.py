#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step05_fix_V1_3D-to-4D_issue.py
===============================
Warp native-space DKI parametric maps into T1w/PVS space with ANTs.

WHY VECTOR MAPS ARE HANDLED SEPARATELY
--------------------------------------
V1/V2/V3 are 3-component eigenvector fields. antsApplyTransforms neither -e 1
(vector) nor -e 4 (multichannel) produced a 3-component output here, because the
NIfTI header carries a nonzero pixdim4 (9.8, inherited from the DWI series that
DIPY fitted). ITK reads a 4D image with temporal spacing as a TIME SERIES, not a
vector image -- ANTs' own help is explicit that a multichannel image is "a vector
image with only spatial dimensions". The result was a silent dim4=1 output
containing only the x component.

Rather than fight ITK's type inference, vector maps are split into three ordinary
3D scalars with nibabel, warped independently, and recombined. Identical result,
no flags, and the recombined header is written cleanly with pixdim4 = 1.

  split  ->  warp x, warp y, warp z  ->  merge

VERIFICATION IS BUILT IN
------------------------
The silent-truncation failure is exactly the kind that propagates into an
analysis unnoticed, so every output is checked after writing:
  * expected number of components (3 for vectors, 1 for scalars)
  * grid matches the T1w reference
  * for vectors, the median norm inside the brain

NOTE ON THE VECTOR NORM
-----------------------
Warping a direction field as three independent interpolated scalars does NOT
preserve unit length, and ANTs does not reorient the vectors. If the reported
median norm is well below 1.0, the warped V1 is a degraded reference and any
angle computed against it partly measures its own interpolation error rather than
tensor reorientation error. Read that number before trusting the reorientation QC
in pvs_orientation_analysis.py.

USAGE
-----
  python apply_coregistration_of_T1w_and_bsmall.py --maps V1
  python apply_coregistration_of_T1w_and_bsmall.py --maps Dxx Dyy Dzz Dxy Dxz Dyz V1
  python apply_coregistration_of_T1w_and_bsmall.py --maps V1 --subjects sub-149 --dry-run
"""

import os
import glob
import shutil
import argparse
import tempfile
import subprocess
import numpy as np
import nibabel as nib
from pipeline_config import cfg

# ---- CONFIG ----
DKI_ROOT = cfg("data_root")
T1_ROOT = cfg("pvs_root")

SOURCE_SUBDIR = "DATASINK/DKISpace/tensor"      # where the native maps live
OUT_DIRNAME = "DKI_in_T1w_space"        # written here
TRANSFORM_DIRNAME = "DKI_in_T1w_space"  # transforms read from here

ANTS_CMD = "antsApplyTransforms"
INTERP = "Linear"        # never BSpline: it overshoots on non-negative maps
VECTOR_MAPS = ("V1", "V2", "V3")
DEFAULT_MAPS = ["Dxx", "Dyy", "Dzz", "Dxy", "Dxz", "Dyz", "V1"]


# ==========================================================================
# DISCOVERY
# ==========================================================================
def find_source(subj_dir, name):
    p = os.path.join(subj_dir, "pipeline_output", SOURCE_SUBDIR, f"{name}.nii.gz")
    return p if os.path.exists(p) else None


def find_t1(subj_id):
    hits = sorted(glob.glob(os.path.join(
        T1_ROOT, f"**/{subj_id}*T1wbrain_restore.nii.gz"), recursive=True))
    return hits[0] if hits else None


def find_transforms(t1_path, subj_id):
    tdir = os.path.dirname(t1_path).replace("Segmentation", TRANSFORM_DIRNAME)
    warp = os.path.join(tdir, f"{subj_id}_ses-01_T1wbrain_restore_from_b0_1Warp.nii.gz")
    aff = os.path.join(tdir, f"{subj_id}_ses-01_T1wbrain_restore_from_b0_0GenericAffine.mat")
    missing = [os.path.basename(t) for t in (warp, aff) if not os.path.exists(t)]
    return (warp, aff), missing


# ==========================================================================
# WARPING
# ==========================================================================
def _ants(moving, fixed, out, transforms, dry_run=False):
    cmd = [ANTS_CMD, "-d", "3", "-i", moving, "-r", fixed, "-o", out,
           "-t", transforms[0], "-t", transforms[1], "-n", INTERP]
    if dry_run:
        print("    DRY:", " ".join(cmd))
        return True
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"    ANTs failed: {(r.stderr or r.stdout).strip()[:200]}")
        return False
    return True


def warp_scalar(moving, fixed, out, transforms, dry_run=False):
    return _ants(moving, fixed, out, transforms, dry_run)


def warp_vector(moving, fixed, out, transforms, dry_run=False):
    """Split into 3 scalars, warp each, recombine. Avoids ITK type inference."""
    img = nib.load(moving)
    if img.ndim != 4 or img.shape[-1] != 3:
        print(f"    not a 3-component vector: shape {img.shape}")
        return False
    data = np.asarray(img.dataobj, np.float32)

    tmp = tempfile.mkdtemp(prefix="vecwarp_")
    try:
        warped = []
        for k in range(3):
            src = os.path.join(tmp, f"c{k}.nii.gz")
            dst = os.path.join(tmp, f"c{k}_w.nii.gz")
            # 3D header, no temporal spacing -- this is the whole point
            h = img.header.copy()
            h.set_data_shape(data.shape[:3])
            nib.save(nib.Nifti1Image(data[..., k], img.affine, h), src)
            if dry_run:
                print(f"    DRY: component {k}")
                continue
            if not _ants(src, fixed, dst, transforms):
                return False
            warped.append(dst)
        if dry_run:
            return True

        vols = [nib.load(w) for w in warped]
        merged = np.stack([np.asarray(v.dataobj, np.float32) for v in vols], axis=-1)
        ref = vols[0]
        hdr = ref.header.copy()
        hdr.set_data_shape(merged.shape)
        z = ref.header.get_zooms()[:3]
        hdr.set_zooms(tuple(z) + (1.0,))        # pixdim4 = 1, not 9.8
        if os.path.dirname(out):
            os.makedirs(os.path.dirname(out), exist_ok=True)
        nib.save(nib.Nifti1Image(merged, ref.affine, hdr), out)
        return True
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ==========================================================================
# VERIFICATION
# ==========================================================================
def verify(out, fixed, is_vector):
    if not os.path.exists(out):
        return "MISSING output"
    img = nib.load(out)
    ref = nib.load(fixed)
    want = 3 if is_vector else 1
    got = img.shape[-1] if img.ndim == 4 else 1
    if got != want:
        return f"FAIL dim4={got}, expected {want}"
    if img.shape[:3] != ref.shape[:3]:
        return f"FAIL grid {img.shape[:3]} != T1w {ref.shape[:3]}"
    if is_vector:
        v = np.asarray(img.dataobj, np.float32)
        n = np.linalg.norm(v, axis=-1)
        nz = n[n > 1e-6]
        if nz.size == 0:
            return "FAIL all-zero vectors"
        return f"ok dim4=3, median norm {np.median(nz):.3f}"
    return "ok"


# ==========================================================================
# MAIN
# ==========================================================================
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--maps", nargs="+", default=DEFAULT_MAPS)
    ap.add_argument("--subjects", nargs="*", default=None)
    ap.add_argument("--dki-root", default=DKI_ROOT)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    dirs = ([os.path.join(a.dki_root, s) for s in a.subjects] if a.subjects
            else sorted(glob.glob(os.path.join(a.dki_root, "sub-*"))))
    print(f"{len(dirs)} subjects | maps: {a.maps} | interp {INTERP}")
    print(f"vector maps {VECTOR_MAPS} are split/warped/merged\n")

    results = []
    for d in dirs:
        sub = os.path.basename(d.rstrip("/"))
        t1 = find_t1(sub)
        if not t1:
            print(f"{sub}: no T1w -- skipped")
            continue
        transforms, missing = find_transforms(t1, sub)
        if missing:
            print(f"{sub}: missing transform(s) {missing} -- skipped")
            continue
        outdir = os.path.dirname(t1).replace("Segmentation", OUT_DIRNAME)
        os.makedirs(outdir, exist_ok=True)

        for name in a.maps:
            src = find_source(d, name)
            if not src:
                print(f"{sub} {name}: source not found -- skipped")
                continue
            out = os.path.join(outdir, f"{name}_toT1w_Warped.nii.gz")
            if os.path.exists(out) and not a.overwrite:
                print(f"{sub} {name}: exists, not overwritten")
                continue
            is_vec = name.startswith(VECTOR_MAPS)
            fn = warp_vector if is_vec else warp_scalar
            ok = fn(src, t1, out, transforms, a.dry_run)
            if a.dry_run:
                continue
            status = verify(out, t1, is_vec) if ok else "ANTs failed"
            print(f"{sub} {name}: {status}")
            results.append((sub, name, status))

    bad = [r for r in results if not r[2].startswith("ok")]
    print(f"\n{len(results) - len(bad)}/{len(results)} outputs verified")
    if bad:
        print("problems:")
        for s, n, st in bad[:20]:
            print(f"  {s} {n}: {st}")

    vec = [r for r in results if r[1].startswith(VECTOR_MAPS) and "norm" in r[2]]
    if vec:
        norms = [float(r[2].split("norm ")[1]) for r in vec]
        print(f"\nvector median norm across subjects: {np.median(norms):.3f} "
              f"(range {min(norms):.3f}-{max(norms):.3f})")
        if np.median(norms) < 0.9:
            print("  [!] well below 1.0. Component-wise interpolation has shortened")
            print("      the vectors, so the warped V1 is a degraded reference and")
            print("      the reorientation QC angle will partly measure ITS error,")
            print("      not the tensor's. Interpret that angle cautiously.")
        else:
            print("  near 1.0 -- the warped V1 is a usable reference for the QC.")


if __name__ == "__main__":
    main()
