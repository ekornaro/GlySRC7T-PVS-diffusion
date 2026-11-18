#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step08_diagnose_block_difference.py
===================================
Find what differs between the two acquisition blocks, given that reorientation
helps block A and hurts block B against the brainstem's own shape axis.

WHAT WE KNOW
------------
  block A (sub <= 145)  mean delta_shape  +0.0069   improved 23/38
  block B (sub >= 149)  mean delta_shape  -0.0166   improved  2/22
  between-block p = 9.6e-06, with a sharp step at the block boundary

Ruled out already: gradient tables byte-identical, b-values and JSON sidecars
matching, same ANTs version, and head positioning (the brainstem shape axis and
the scanner S-I axis coincide to within 0.036 in coherence, corr(delta_shape,
delta_fixed) = 0.975 -- there is almost no tilt variability to exploit).

WHAT THIS CHECKS
----------------
1. DIRECTION COSINES / OBLIQUITY.  This is the leading suspect, and it is a
   limitation of reorient_tensors_to_T1w.py rather than of the data.
   local_jacobian differentiates with respect to voxel index and divides by
   spacing, giving

       G = du/di . diag(1/s) = (du/dx_phys) . D

   where D is the image direction cosine matrix. The true Jacobian needs
   du/dx_phys, so the computed rotation is wrong by a factor of D unless
   D = I. A fixed obliquity would bias every subject equally and cancel; an
   obliquity that DIFFERS BY BLOCK produces exactly the observed pattern.

2. AFFINE ROTATION AXIS, not just magnitude.  Earlier work compared rotation
   angles (which did not differ meaningfully by group). The axis was never
   examined. Two blocks rotating about different axes would interact with any
   orientation reference differently.

3. b0 AVAILABILITY.  The convention check reported "b0 not found" for every
   block B subject. If those files genuinely do not exist, block B went through
   a different preprocessing route, which would be the whole explanation.

4. bvec PROVENANCE.  Eddy-rotated bvecs are motion-corrected per volume per
   subject, so two subjects CANNOT have identical files. The byte-identical
   comparison reported earlier therefore came from the unrotated originals. If
   the DKI fit consumed those, every subject carries a motion-dependent
   orientation error -- and if motion differs by block, so does the error.

USAGE
-----
  python diagnose_block_difference.py --root <ROOT> --dwi-root <DWI_ROOT>
"""

import os
import glob
import argparse
import numpy as np
from pipeline_config import cfg

try:
    import nibabel as nib
    import pandas as pd
    from scipy import stats
    from scipy.io import loadmat
except ImportError as e:
    raise SystemExit(f"missing dependency: {e}")

ROOT_DIR = cfg("pvs_root")
DWI_ROOT = cfg("data_root")
BAD = ("Warp", "Affine", "Warped", "BiasField", "from_b0", "mask", "_seg")
# The b0 ANTs registered, at its known fixed relative path. NOT corrected_b0
# from select_corrected_b0 -- that is a different pipeline stage.
B0_RELPATH = os.path.join("pipeline_output", "DATASINK", "DKISpace", "dwi_proc",
                          "DKI_denoised_gibbs_eddy_b0.nii.gz")
B0_PATTERNS = ["DKI_denoised_gibbs_eddy_b0.nii*", "*_ses-01_b0.nii*", "*_b0.nii*"]


def find_one(d, pats, subj=None):
    for pat in pats:
        for h in sorted(glob.glob(os.path.join(d, "**", pat), recursive=True)):
            if any(t in os.path.basename(h) for t in BAD):
                continue
            if subj and subj not in h:
                continue
            return h
    return None


def direction_cosines(path):
    """Unit direction matrix D, and obliquity in degrees.

    Obliquity is the largest angle between any image axis and the scanner axis
    it is closest to (the AFNI definition). Computed per column, so it stays
    valid for radiological/LAS images: the previous version formed a signed
    permutation P and took arccos((tr(D P^T)-1)/2), which is only a rotation
    angle when det(P) = +1. For an axis-aligned LAS image det(P) = -1 and the
    formula is meaningless.
    """
    A = nib.load(path).affine[:3, :3]
    z = np.sqrt((A ** 2).sum(axis=0))
    D = A / np.where(z > 1e-9, z, 1.0)
    ang = max(float(np.degrees(np.arccos(np.clip(np.abs(D[:, k]).max(), 0, 1))))
              for k in range(3))
    return D, ang


def matrix_angle(D1, D2):
    """Rotation angle between two direction matrices, in degrees.

    Returns NaN if they differ by a reflection (opposite handedness), which is
    a categorical difference rather than an angle and must not be averaged.
    """
    if not (isinstance(D1, np.ndarray) and isinstance(D2, np.ndarray)):
        return np.nan          # missing entries arrive as NaN floats, not None
    R = np.asarray(D1) @ np.asarray(D2).T
    if np.linalg.det(R) < 0:
        return np.nan
    return float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))


def frame_relation(D1, D2):
    """Handedness and the angle of the NEAREST ROTATION between two frames.

    A determinant of -1 means the two voxel-to-world mappings have opposite
    handedness. That is a reflection, and no rotation can represent it -- which
    matters because reorient_tensors_to_T1w.py forces det = +1 in its polar
    decomposition and would silently substitute a proper rotation for it.
    """
    if not (isinstance(D1, np.ndarray) and isinstance(D2, np.ndarray)):
        return np.nan, np.nan
    R = np.asarray(D1) @ np.asarray(D2).T
    det = float(np.linalg.det(R))
    W, sv, Vt = np.linalg.svd(R)          # nearest orthogonal, then force det=+1
    Q = W @ Vt
    if np.linalg.det(Q) < 0:
        W[:, -1] *= -1
        Q = W @ Vt
    ang = float(np.degrees(np.arccos(np.clip((np.trace(Q) - 1) / 2, -1, 1))))
    return det, ang


def affine_axis(path):
    """Rotation angle and unit rotation axis of the ANTs affine."""
    m = loadmat(path)
    k = next((x for x in m if x.startswith("AffineTransform")), None)
    if k is None:
        return np.nan, np.array([np.nan] * 3)
    A = np.asarray(m[k]).ravel()[:9].reshape(3, 3)
    W, s, Vt = np.linalg.svd(A)
    R = W @ Vt
    if np.linalg.det(R) < 0:
        W[:, -1] *= -1
        R = W @ Vt
    ang = float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))
    ax = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    n = np.linalg.norm(ax)
    return ang, (ax / n if n > 1e-9 else np.array([0.0, 0.0, 1.0]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=ROOT_DIR)
    ap.add_argument("--dwi-root", default=DWI_ROOT)
    ap.add_argument("--block-split", type=int, default=149)
    ap.add_argument("--out", default="block_difference_diagnosis.csv")
    a = ap.parse_args()

    subs = sorted(s for s in os.listdir(a.root)
                  if s.startswith("sub-") and os.path.isdir(os.path.join(a.root, s)))
    # scan the DWI root ONCE, not once per subject
    all_bvec = glob.glob(os.path.join(a.dwi_root, "**", "*.bvec"), recursive=True)
    all_rot = glob.glob(os.path.join(a.dwi_root, "**", "*rotated_bvec*"),
                        recursive=True)
    rows = []
    for s in subs:
        sd = os.path.join(a.root, s)
        num = int(s.split("-")[1])
        r = {"Subject": s, "num": num,
             "block": "B" if num >= a.block_split else "A"}

        Dt1 = None
        t1 = find_one(sd, ["*_ses-01_T1wbrain_restore.nii*"])
        if t1:
            Dt1, ob = direction_cosines(t1)
            r["t1w_obliquity_deg"] = ob
            r["t1w_D"] = Dt1

        b0 = os.path.join(a.dwi_root, s, B0_RELPATH)
        if not os.path.exists(b0):
            b0 = find_one(os.path.join(a.dwi_root, s), B0_PATTERNS)
        r["b0_found"] = bool(b0)
        if b0:
            Db0, ob = direction_cosines(b0)
            r["b0_obliquity_deg"] = ob
            r["b0_shape"] = str(nib.load(b0).shape)
            # the quantity that actually matters: how far the DWI frame is
            # rotated relative to the T1w frame the tensors were warped into
            r["b0_vs_t1w_deg"] = matrix_angle(Db0, Dt1)
            det, ang = frame_relation(Db0, Dt1)
            r["b0_vs_t1w_reflected"] = (det < 0) if np.isfinite(det) else np.nan
            r["b0_vs_t1w_rot_deg"] = ang

        aff = find_one(sd, ["*from_b0_0GenericAffine.mat"])
        if aff:
            ang, ax = affine_axis(aff)
            r["affine_deg"] = ang
            r["axis_x"], r["axis_y"], r["axis_z"] = ax
        r["has_rotated_bvecs"] = any(f"/{s}/" in p or s in os.path.basename(p)
                                     for p in all_rot)
        r["has_plain_bvec"] = any(f"/{s}/" in p or s in os.path.basename(p)
                                  for p in all_bvec)
        rows.append(r)

    d = pd.DataFrame(rows)
    A_, B_ = d[d.block == "A"], d[d.block == "B"]
    print(f"subjects: A={len(A_)}  B={len(B_)}\n")

    print("=== 3. b0 AVAILABILITY (explains the 'not found' in the convention check)")
    print(f"  block A: {int(A_.b0_found.sum())}/{len(A_)} have a b0")
    print(f"  block B: {int(B_.b0_found.sum())}/{len(B_)} have a b0")
    if B_.b0_found.sum() == 0 and A_.b0_found.sum() > 0:
        print("  [!] NO block B subject has a b0 in the DWI root. Block B was")
        print("      preprocessed differently -- this is very likely the answer.")

    print("\n=== 1. OBLIQUITY (deg from nearest axis-aligned orientation)")
    for col in ("t1w_obliquity_deg", "b0_obliquity_deg"):
        if col not in d or d[col].isna().all():
            continue
        ga, gb = A_[col].dropna(), B_[col].dropna()
        if len(ga) > 2 and len(gb) > 2:
            p = stats.mannwhitneyu(ga, gb).pvalue
            flag = "  <== DIFFERS BY BLOCK" if p < 0.05 else ""
            print(f"  {col:20s} A {ga.median():6.2f}  B {gb.median():6.2f}"
                  f"   p={p:.2e}{flag}")
    if "b0_vs_t1w_deg" in d and d.b0_vs_t1w_deg.notna().any():
        ga, gb = A_.b0_vs_t1w_deg.dropna(), B_.b0_vs_t1w_deg.dropna()
        print(f"  b0_vs_t1w_deg        A {ga.median():6.2f}  B {gb.median():6.2f}"
              + (f"   p={stats.mannwhitneyu(ga, gb).pvalue:.2e}"
                 if len(ga) > 2 and len(gb) > 2 else ""))
        print("    ^ rotation of the DWI frame relative to the T1w frame; the")
        print("      quantity a block-dependent tensor-frame error would live in")
        if "b0_vs_t1w_rot_deg" in d:
            ga2 = A_.b0_vs_t1w_rot_deg.dropna() if "b0_vs_t1w_rot_deg" in A_ else []
            gb2 = B_.b0_vs_t1w_rot_deg.dropna() if "b0_vs_t1w_rot_deg" in B_ else []
            if len(ga2) > 2 and len(gb2) > 2:
                print(f"  b0_vs_t1w_rot_deg    A {np.median(ga2):6.2f}  "
                      f"B {np.median(gb2):6.2f}   "
                      f"p={stats.mannwhitneyu(ga2, gb2).pvalue:.2e}"
                      "   (nearest rotation, defined even when reflected)")
        if "b0_vs_t1w_reflected" in d:
            rf = d.b0_vs_t1w_reflected.dropna().astype(bool)
            na = int(A_.get("b0_vs_t1w_reflected", pd.Series(dtype=object))
                     .dropna().astype(bool).sum())
            nb = int(B_.get("b0_vs_t1w_reflected", pd.Series(dtype=object))
                     .dropna().astype(bool).sum())
            da = int(A_.b0_vs_t1w_reflected.notna().sum())
            db = int(B_.b0_vs_t1w_reflected.notna().sum())
            print(f"\n  HANDEDNESS: b0 and T1w differ by a REFLECTION in "
                  f"{int(rf.sum())}/{len(rf)} subjects")
            print(f"    block A: {na}/{da}    block B: {nb}/{db}")
            if da and db:
                try:
                    p = stats.fisher_exact([[na, da - na], [nb, db - nb]])[1]
                    print(f"    Fisher exact p={p:.2e}"
                          + ("   <== REFLECTION TRACKS THE BLOCK" if p < 0.05 else ""))
                except Exception:
                    pass
            print("    A reflection cannot be produced by any rotation.")
            print("    reorient_tensors_to_T1w.py forces det=+1 in its polar")
            print("    decomposition, so for reflected subjects it substitutes a")
            print("    proper rotation for an improper one -- wrong, and wrong in")
            print("    a way that depends on which subjects are reflected.")

    if "t1w_D" in d and d.t1w_D.notna().any():
        # tolerant comparison: angle of each subject's D from a reference D.
        # Comparing rounded matrices as STRINGS made near-identical matrices
        # look distinct and fired the warning spuriously.
        ref = d.t1w_D.dropna().iloc[0]
        d["t1w_D_dev_deg"] = [matrix_angle(x, ref) if x is not None else np.nan
                              for x in d.t1w_D]
        ga = d[d.block == "A"].t1w_D_dev_deg.dropna()
        gb = d[d.block == "B"].t1w_D_dev_deg.dropna()
        print(f"\n  T1w direction matrix, deviation from a common reference:")
        print(f"    A max {ga.max():.3f} deg    B max {gb.max():.3f} deg")
        if len(ga) > 2 and len(gb) > 2 and max(ga.max(), gb.max()) > 0.5:
            p = stats.mannwhitneyu(ga, gb).pvalue
            print(f"    between-block p={p:.2e}"
                  + ("   <== DIFFERS BY BLOCK" if p < 0.05 else ""))
        else:
            print("    all subjects share one T1w direction matrix -> the "
                  "direction-cosine term")
            print("    is a constant bias and CANNOT explain the block split.")

    print("\n=== 2. AFFINE ROTATION AXIS (magnitude was already known not to differ)")
    if "axis_z" in d:
        for c in ("axis_x", "axis_y", "axis_z", "affine_deg"):
            ga, gb = A_[c].dropna(), B_[c].dropna()
            if len(ga) > 2 and len(gb) > 2:
                p = stats.mannwhitneyu(ga, gb).pvalue
                flag = "  <== DIFFERS" if p < 0.05 else ""
                print(f"  {c:12s} A {ga.median():+.3f}  B {gb.median():+.3f}"
                      f"   p={p:.2e}{flag}")

    print("\n=== 4. bvec PROVENANCE")
    print(f"  eddy-rotated bvecs present: A {int(A_.has_rotated_bvecs.sum())}"
          f"/{len(A_)}   B {int(B_.has_rotated_bvecs.sum())}/{len(B_)}")
    if not d.has_rotated_bvecs.any():
        print("  [!] no *rotated_bvecs* found anywhere. Confirm which file the")
        print("      DKI fit consumed -- byte-identical bvecs across subjects")
        print("      means the unrotated originals were compared, and possibly")
        print("      used.")

    d.to_csv(a.out, index=False)
    print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()
