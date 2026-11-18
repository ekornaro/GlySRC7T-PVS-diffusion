#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step06_pvs_orientation_analysis.py (v2)
=======================================
v2 ADDS: AFFINE u-ROTATION CORRECTION (--rotate-u), licensed by the step18/19
audits on the regenerated transforms.

WHY
---
The tensor components Dxx..Dyz were warped into T1w space as scalars, so at
every T1w voxel the tensor D is still EXPRESSED IN THE NATIVE DKI FRAME, while
the PVS axis u is computed from T1w voxel coordinates. u^T D u therefore mixes
frames, and the step18 audit measured the per-subject affine rotation between
them at 0.4-6.8 deg (median ~1.9). Step19 Q1 showed the residual nonlinear-warp
rotation is NOT spatially coincident with PVS (enrichment 0.087 vs 0.100
benchmark) and does not differ by group, so correcting the GLOBAL AFFINE term
alone is sufficient; the warp term remains as uniform, group-fair angular noise
to be stated as a limitation.

The correction rotates u into the tensor's frame (never touching the tensor
volumes, so the det=-1 handedness block never enters): for each subject the
rotation factor R of the ANTs 0GenericAffine (right polar decomposition,
reflection-safe) is composed with the direction-cosine chain

    u_DKIidx = C_dwi^-1 . F . R^(+/-1) . F . C_t1w . u_T1widx

where C_* are the rotation parts of the NIfTI RAS affines, F = diag(-1,-1,1)
converts RAS<->LPS (ITK works in LPS), and the +/-1 selects --rotate-u R or Rt.
C_dwi comes from a native-space map found under --dki-root; if unavailable the
DWI direction cosines are assumed identity (near-axial; a warning is printed).

ON THE R vs Rt DIRECTION
------------------------
ANTs stores the fixed->moving point mapping (fixed = T1w here), which makes R
the analytically correct choice for taking T1w-frame vectors into the native
frame. The Rt mode exists as a SENSITIVITY CHECK: because theta is ~random and
the median affine rotation is ~2 deg, tier results invariant to {R, Rt} are
proof the direction ambiguity is immaterial -- run both once and report that.

MODES (per the no-silent-defaults rule, the default is 'off' and reproduces v1
bit-for-bit; the active mode is echoed in the header and written to the QC csv)
    --rotate-u off   v1 behaviour, u left in the T1w frame
    --rotate-u R     apply the affine rotation (analytically correct direction)
    --rotate-u Rt    apply its transpose (sensitivity)

Original v1 documentation follows.

----------------------------------------------------------------------
Write brain-wide Dpar / Dperp / theta / Linearity maps into
<subj>/revised_DKI_in_T1w_space/ so the tier-extraction script can sample them
exactly like MD, AD, lowb_ADC and every other metric.

THE PROBLEM THIS SOLVES
-----------------------
D_par = u^T D u needs an axis u. u comes from PCA of a PVS cluster's voxels, so
it is defined only INSIDE clusters. That gives a per-cluster table, which answers
"does D_par differ by group?" but NOT "is that difference specific to PVS?" --
and the second question is the one Reviewer #3 is actually disputing.

SOLUTION: NEAREST-CLUSTER AXIS PROPAGATION
------------------------------------------
Every brain voxel inherits the axis of its NEAREST PVS cluster, via
distance_transform_edt(~pvs, return_indices=True). D_par is then defined
brain-wide and means the same thing everywhere: diffusivity along the local
perivascular direction.

That makes the tiers directly comparable. NearPVS voxels are evaluated along the
axis of the PVS they surround; FakeControl voxels along the axis of whichever
cluster is nearest to them. So the contrast

    D_par(PVS) - D_par(FakeControl)

isolates location while holding the KIND of measurement fixed. A generic
white-matter effect appears in both tiers; a perivascular effect does not.

WHY THIS IS THE RIGHT CONTROL
-----------------------------
The alternative -- giving each control ROI its own PCA axis -- is wrong here: a
region-grown blob is not elongated, so its principal axis is essentially
arbitrary and D_par would be measuring noise. Propagation keeps the reference
direction anatomically meaningful everywhere.

CAVEAT TO STATE IN THE PAPER
----------------------------
Far from any PVS the propagated axis is a long extrapolation and the direction
becomes arbitrary. MAX_PROPAGATE_MM caps this; beyond it, voxels are left as NaN.
The default 20 mm comfortably covers the NearPVS shell (~6 mm) and both exclusion
margins (8 and 12 mm), so every tier is inside the valid region.

OUTPUTS (into <subj>/revised_DKI_in_T1w_space/, so the tier script finds them)
------------------------------------------------------------------------------
  Dpar_toT1w_Warped.nii.gz    u^T D u, diffusivity ALONG the local PVS axis
  Dperp_toT1w_Warped.nii.gz   (trace - Dpar)/2, ACROSS it
  Dpar_minus_Dperp_toT1w_Warped.nii.gz   orientation-referenced anisotropy
  theta_toT1w_Warped.nii.gz   angle (deg) between the PVS axis and the tensor's
                              principal eigenvector
  clusterLinearity_toT1w_Warped.nii.gz   PCA lambda1/lambda2 of the source
                              cluster; use to mask out unreliable axes
  distToPVS_toT1w_Warped.nii.gz   mm to nearest PVS, for QC

The _toT1w_Warped suffix is for compatibility with the tier script's exact-match
pattern. These maps were COMPUTED in T1w space, not warped into it.

USAGE
-----
  python pvs_orientation_maps.py --subjects sub-149            # pilot
  python pvs_orientation_maps.py                               # cohort
"""

import os
import glob
import argparse
import numpy as np
import nibabel as nib
from scipy import ndimage
import pandas as pd
from pipeline_config import cfg

T1_ROOT = cfg("pvs_root")
STUDY_XLSX = cfg("study_xlsx")

DKI_SUBDIR = "DKI_in_T1w_space"
MAP_SUFFIX = "_toT1w_Warped.nii.gz"
COMPONENTS = ["Dxx", "Dyy", "Dzz", "Dxy", "Dxz", "Dyz"]

MIN_CLUSTER_SIZE = 3
MAX_PROPAGATE_MM = 20.0      # beyond this the propagated axis is meaningless
CHUNK = 500_000              # voxels per batch in the eigen solve

DKI_ROOT_DEFAULT = cfg("data_root")
NATIVE_REF_SUBPATH = "pipeline_output/DATASINK/DKISpace/tensor/Dxx.nii.gz"
_F_LPS = np.diag([-1.0, -1.0, 1.0])          # RAS <-> LPS


# ==========================================
# u-ROTATION (v2) -- step18's extraction, verbatim where possible
# ==========================================
def _polar_rotation(M):
    """Rotation factor of M = R P, reflection-safe (from step18)."""
    W, s, Vt = np.linalg.svd(M)
    R = W @ Vt
    if np.linalg.det(R) < 0:
        W[:, -1] *= -1
        R = W @ Vt
    return R


def _read_ants_affine(path):
    """3x3 linear part of an ITK/ANTs .mat transform (from step18)."""
    from scipy.io import loadmat
    m = loadmat(path)
    key = next((k for k in m if k.startswith("AffineTransform")), None)
    if key is None:
        raise RuntimeError(f"{path}: no AffineTransform_* variable")
    p = np.asarray(m[key]).ravel()
    if p.size < 12:
        raise RuntimeError(f"{path}: expected >=12 parameters, got {p.size}")
    return p[:9].reshape(3, 3)


def _dircos(affine):
    """Normalized rotation part of a NIfTI RAS affine (columns scaled to unit)."""
    M = np.asarray(affine, float)[:3, :3]
    return M / np.linalg.norm(M, axis=0, keepdims=True)


def build_u_rotation(subj_dir, sub, t1w_affine, mode, dki_root):
    """
    Total rotation taking a T1w-index-frame unit vector into the native DKI
    (bvec/index) frame in which Dxx..Dyz are expressed. Returns (R_total,
    angle_deg, note). mode: 'R' | 'Rt'.
    """
    hits = sorted(glob.glob(os.path.join(subj_dir, "**",
                                         "*from_b0_0GenericAffine.mat"),
                            recursive=True))
    if not hits:
        raise RuntimeError("no *from_b0_0GenericAffine.mat found for --rotate-u")
    A = _read_ants_affine(hits[0])
    R_a = _polar_rotation(A)                          # T1w -> b0, LPS physical
    if mode == "Rt":
        R_a = R_a.T

    C_t1 = _dircos(t1w_affine)
    note = ""
    ref = os.path.join(dki_root, sub, NATIVE_REF_SUBPATH)
    if os.path.exists(ref):
        C_dwi = _dircos(nib.load(ref).affine)
    else:
        C_dwi = np.eye(3)
        note = "no native ref; DWI direction cosines assumed identity"

    R_total = np.linalg.inv(C_dwi) @ _F_LPS @ R_a @ _F_LPS @ C_t1
    R_total = _polar_rotation(R_total)               # renormalize numerics
    ang = float(np.degrees(np.arccos(np.clip((np.trace(R_total) - 1) / 2, -1, 1))))
    return R_total.astype(np.float32), ang, note


def included_subjects(xlsx=STUDY_XLSX):
    """Study inclusion list. Every glob-based script must filter on this."""
    x = pd.ExcelFile(xlsx)
    out = set()
    for sheet in ("HC", "PPCS", "Boxers"):
        t = x.parse(sheet)
        u = t["Will we used it?"].astype(str).str.strip().str.lower()
        out |= {f"sub-{int(n)}" for n in t.loc[u == "yes", "Subj"]}
    return out


def _load_components(map_dir):
    imgs = {}
    for c in COMPONENTS:
        p = os.path.join(map_dir, f"{c}{MAP_SUFFIX}")
        if not os.path.exists(p):
            raise FileNotFoundError(f"missing {c}{MAP_SUFFIX}")
        imgs[c] = nib.load(p)
    return imgs


def process_subject(subj_dir, overwrite=False, max_prop=MAX_PROPAGATE_MM,
                    rotate_u="off", dki_root=DKI_ROOT_DEFAULT):
    sub = os.path.basename(subj_dir.rstrip("/"))
    map_dir = os.path.join(subj_dir, DKI_SUBDIR)
    imgs = _load_components(map_dir)
    ref = imgs["Dxx"]
    shape = ref.shape[:3]

    hits = glob.glob(os.path.join(subj_dir, "**",
                                  "*_ses-01_PVS_postprocessed.nii*"), recursive=True)
    if len(hits) != 1:
        raise RuntimeError(f"PVS mask not unique ({len(hits)} hits)")
    pvs_img = nib.load(hits[0])
    pvs = pvs_img.get_fdata() > 0
    if pvs.shape != shape:
        raise RuntimeError(f"PVS {pvs.shape} != tensor {shape}")
    zooms = tuple(float(z) for z in pvs_img.header.get_zooms()[:3])

    # ---- per-cluster principal axes ----
    labeled, n_feat = ndimage.label(pvs, structure=np.ones((3, 3, 3)))
    axes = np.zeros((n_feat + 1, 3), np.float32)      # index 0 unused
    linear = np.zeros(n_feat + 1, np.float32)
    kept = 0
    for lab, sl in enumerate(ndimage.find_objects(labeled), start=1):
        if sl is None:
            continue
        sub_lab = labeled[sl] == lab
        if int(sub_lab.sum()) < MIN_CLUSTER_SIZE:
            continue
        off = np.array([s.start for s in sl])
        coords = (np.argwhere(sub_lab) + off) * np.array(zooms)
        c = coords - coords.mean(axis=0)
        w, V = np.linalg.eigh((c.T @ c) / max(len(c) - 1, 1))
        axes[lab] = V[:, -1]
        linear[lab] = (w[-1] / w[-2]) if w[-2] > 1e-12 else 999.0
        kept += 1

    # ---- v2: rotate u into the tensor's native frame (never the tensor) ----
    u_rot_deg, u_rot_note = 0.0, "off"
    if rotate_u in ("R", "Rt"):
        R_tot, u_rot_deg, u_rot_note = build_u_rotation(
            subj_dir, sub, ref.affine, rotate_u, dki_root)
        axes = axes @ R_tot.T                        # row vectors: u' = R u
        n_ax = np.linalg.norm(axes, axis=1, keepdims=True)
        axes = np.where(n_ax > 1e-6, axes / np.maximum(n_ax, 1e-12), axes)
        u_rot_note = u_rot_note or rotate_u

    # ---- propagate the nearest cluster's axis to every voxel ----
    dist, inds = ndimage.distance_transform_edt(~pvs, sampling=zooms,
                                                return_indices=True)
    nearest_lab = labeled[tuple(inds)]                # (X,Y,Z) label of nearest PVS
    nearest_lab[pvs] = labeled[pvs]                   # inside PVS: its own label

    valid = (nearest_lab > 0) & (dist <= max_prop) & (linear[nearest_lab] > 0)
    idx = np.flatnonzero(valid.ravel())
    if idx.size == 0:
        raise RuntimeError("no valid voxels after propagation")

    U = axes[nearest_lab.ravel()[idx]]                # (n, 3)
    lin_v = linear[nearest_lab.ravel()[idx]]

    # ---- tensor at those voxels, in chunks ----
    g = {c: np.asarray(imgs[c].dataobj, np.float32).ravel()[idx] for c in COMPONENTS}
    n = idx.size
    D = np.empty((n, 3, 3), np.float32)
    D[:, 0, 0] = g["Dxx"]; D[:, 1, 1] = g["Dyy"]; D[:, 2, 2] = g["Dzz"]
    D[:, 0, 1] = D[:, 1, 0] = g["Dxy"]
    D[:, 0, 2] = D[:, 2, 0] = g["Dxz"]
    D[:, 1, 2] = D[:, 2, 1] = g["Dyz"]
    del g

    d_par = np.einsum("ni,nij,nj->n", U, D, U)
    tr = np.einsum("nii->n", D)
    d_perp = (tr - d_par) / 2.0

    theta = np.empty(n, np.float32)
    for i in range(0, n, CHUNK):
        s = slice(i, min(i + CHUNK, n))
        _, vv = np.linalg.eigh(D[s])
        v1 = vv[..., -1]
        theta[s] = np.degrees(np.arccos(np.clip(
            np.abs(np.einsum("ni,ni->n", v1, U[s])), 0, 1)))
    del D

    # ---- write ----
    os.makedirs(map_dir, exist_ok=True)
    outs = {"Dpar": d_par, "Dperp": d_perp,
            "Dpar_minus_Dperp": d_par - d_perp, "theta": theta,
            "clusterLinearity": lin_v, "distToPVS": dist.ravel()[idx]}
    for name, vec in outs.items():
        p = os.path.join(map_dir, f"{name}{MAP_SUFFIX}")
        if os.path.exists(p) and not overwrite:
            continue
        vol = np.full(int(np.prod(shape)), np.nan, np.float32)
        vol[idx] = vec
        nib.save(nib.Nifti1Image(vol.reshape(shape), ref.affine), p)

    # v2: provenance sidecar -- the maps keep tier-compatible names, so the
    # rotation mode that produced them must be recorded next to them.
    import json, datetime
    with open(os.path.join(map_dir, "pvs_orientation_provenance.json"), "w") as fh:
        json.dump({"script": "step06_pvs_orientation_analysis.py",
                   "rotate_u": rotate_u, "u_rot_deg": round(u_rot_deg, 3),
                   "u_rot_note": u_rot_note,
                   "written": datetime.datetime.now().isoformat(timespec="seconds")},
                  fh, indent=2)

    inside = pvs.ravel()[idx]
    return {
        "Subject": sub, "n_clusters": kept, "n_valid_vox": int(n),
        "rotate_u": rotate_u, "u_rot_deg": round(u_rot_deg, 3),
        "u_rot_note": u_rot_note,
        "frac_brain_covered": float(n / max(int(np.prod(shape)), 1)),
        "Dpar_med_inPVS": float(np.nanmedian(d_par[inside])),
        "Dperp_med_inPVS": float(np.nanmedian(d_perp[inside])),
        "Dpar_med_outPVS": float(np.nanmedian(d_par[~inside])),
        "theta_med_inPVS": float(np.nanmedian(theta[inside])),
        "median_dist_mm": float(np.median(dist.ravel()[idx])),
        "status": "ok",
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=T1_ROOT)
    ap.add_argument("--subjects", nargs="*", default=None)
    ap.add_argument("--out", default="./pvs_orientation_maps_qc")
    ap.add_argument("--max-propagate-mm", type=float, default=MAX_PROPAGATE_MM)
    ap.add_argument("--rotate-u", choices=["off", "R", "Rt"], default="off",
                    help="v2: rotate the PVS axis u into the tensor's native "
                         "frame using the affine rotation. 'off' (default) "
                         "reproduces v1 exactly; 'R' is the analytically "
                         "correct direction; 'Rt' is the sensitivity check.")
    ap.add_argument("--dki-root", default=DKI_ROOT_DEFAULT,
                    help="root holding native tensor maps (for DWI direction "
                         "cosines); only used when --rotate-u is R/Rt")
    ap.add_argument("--no-inclusion-filter", action="store_true")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    dirs = ([os.path.join(a.root, s) for s in a.subjects] if a.subjects
            else sorted(glob.glob(os.path.join(a.root, "sub-*"))))
    if not a.no_inclusion_filter:
        try:
            inc = included_subjects()
            before = len(dirs)
            dirs = [d for d in dirs if os.path.basename(d.rstrip("/")) in inc]
            print(f"study inclusion filter: {before} -> {len(dirs)} subjects")
        except Exception as exc:
            print(f"[!] could not read inclusion list ({exc}); processing all")

    banner = {"off": "u-rotation OFF (v1-identical run)",
              "R":   "u-rotation R  -- axes rotated into the native DKI frame",
              "Rt":  "u-rotation Rt -- SENSITIVITY direction, not the primary"}
    print(f"\n*** {banner[a.rotate_u]} ***")
    print(f"axis propagation capped at {a.max_propagate_mm} mm\n")
    recs = []
    for d in dirs:
        sub = os.path.basename(d.rstrip("/"))
        try:
            r = process_subject(d, a.overwrite, a.max_propagate_mm,
                                a.rotate_u, a.dki_root)
            recs.append(r)
            print(f"  {sub}: {r['n_clusters']} clusters | "
                  f"Dpar in/out PVS {r['Dpar_med_inPVS']:.3e}/"
                  f"{r['Dpar_med_outPVS']:.3e} | "
                  f"{r['frac_brain_covered']*100:.0f}% coverage")
        except Exception as exc:                     # per-subject isolation
            print(f"  {sub}: FAILED -- {exc}")
            recs.append({"Subject": sub, "status": f"FAILED: {exc}"})

    df = pd.DataFrame(recs)
    df.to_csv(os.path.join(a.out, "orientation_maps_qc.csv"), index=False)
    ok = df[df.get("status", "") == "ok"] if "status" in df else df
    if len(ok):
        print(f"\n{len(ok)} subjects written")
        print(f"median Dpar inside PVS : {ok['Dpar_med_inPVS'].median():.3e}")
        print(f"median Dpar outside PVS: {ok['Dpar_med_outPVS'].median():.3e}")
        print(f"median theta inside PVS: {ok['theta_med_inPVS'].median():.1f} deg")
        print("\nNEXT: add to METRICS_OF_INTEREST in the tier-extraction script:")
        print("  'Dpar', 'Dperp', 'Dpar_minus_Dperp', 'theta'")
        print("then read the PVS-minus-FakeControl contrast for Dpar. That is the")
        print("test of whether any D_par effect is PERIVASCULAR rather than global.")
        print("\nMask unreliable axes with clusterLinearity < 2 before interpreting.")
    print(f"\nQC -> {a.out}")


if __name__ == "__main__":
    main()
