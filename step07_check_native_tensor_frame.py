#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step07_check_native_tensor_frame.py
===================================
Determine whether the diffusion tensors are correctly oriented IN NATIVE DWI
SPACE, before any registration or frame composition enters.

WHY THIS TEST AND NOT ANOTHER
-----------------------------
Every previous orientation check composed frames: DWI voxel -> world ->
registration rotation -> T1w voxel. Four confident structural claims were made
on that basis (head positioning, transform convention, direction cosines,
handedness) and each failed its own test. The reason became clear only at the
end: 'identity' -- the raw warped components with NO transform applied -- already
splits by acquisition block (A 0.6515 vs B 0.7220, p = 2.7e-06). A reference
that is block-dependent before anything is applied cannot validate anything
applied afterwards.

This script removes the whole chain. It fits a tensor in native DWI space from
the DWI data and its own gradient table, and asks a single question: does the
principal diffusion direction line up with anatomy, in the space where the
tensor was fitted?

  * coherence high (~0.85-0.95) and EQUAL across blocks
        -> the tensors are correctly oriented at the source. Everything that
           split by block is downstream, in the warping or frame composition,
           and Dpar in T1w space is the only casualty.

  * coherence splits by block
        -> the problem is in the fit itself, i.e. the gradient table or its
           convention, and it is a property of the data rather than of any
           later code.

REGISTRATION IS USED FOR ONE THING ONLY
---------------------------------------
The anatomical reference is the brainstem, which needs SynthSeg, which lives in
T1w space. The mask is warped into DWI space with nearest-neighbour
interpolation. Warping a BINARY MASK carries no orientation information -- a few
degrees of registration error moves the boundary of a large structure and
changes nothing about the conclusion. No tensor is transformed at any point.

THE REFERENCE
-------------
The brainstem's own SHAPE principal axis, computed from the warped mask in DWI
voxel coordinates, with V1 restricted to FA >= FA_MIN. Both the axis and V1 live
in the same native frame, so nothing about head position or image storage
convention enters. Low-FA voxels have an arbitrary V1 and dilute the score
toward 0.500, which is why they are excluded.

BVEC CONVENTION -- A SECONDARY OUTPUT, AND A WEAK ONE
-----------------------------------------------------
Four sign conventions are scored, but do not lean on them. A sign flip MIRRORS
the tensor, and a mirror leaves V1 unchanged whenever V1 lies in the mirror
plane or along its normal. Simulation with a 52-direction scheme:

    V1 exactly along +z      none 1.000  flip_x 1.000  flip_y 1.000  flip_z 1.000
    V1 tilted 20 deg from z  none 0.766  flip_x 0.766  flip_y 1.000  flip_z 1.000
    V1 tilted 45 deg         none 0.000  flip_x 1.000  flip_y 0.000  flip_z 1.000

The brainstem is the WORST possible region for this, because its fibres run
close to +z, and flips also pair up (flip_y and flip_z differ by a 180 deg
rotation, which |V1 . u| cannot separate). Detecting a flip properly needs
obliquely oriented fibres, which this reference does not provide.

This does not matter for the question being asked. The gradient tables are
byte-identical between the two blocks, so any convention error is IDENTICAL in
both and cannot produce a block difference. The flip columns are reported for
completeness; the verdict rests on the block comparison.

Do not expect coherence near 1.0: the SynthSeg brainstem label includes the
pons, whose transverse pontocerebellar fibres run left-right and legitimately
score near zero. What matters is the comparison between conventions and between
blocks, not the absolute level.

USAGE
-----
  python check_native_tensor_frame.py --root <ROOT> --dwi-root <DWI_ROOT>
  python check_native_tensor_frame.py --root <ROOT> --subjects sub-102 sub-149
"""

import os
import glob
import shutil
import argparse
import tempfile
import subprocess
import numpy as np
from pipeline_config import cfg

try:
    import nibabel as nib
    import pandas as pd
    from scipy import stats
except ImportError as e:
    raise SystemExit(f"missing dependency: {e}")

ROOT_DIR = cfg("pvs_root")
DWI_ROOT = cfg("data_root")
DWI_DIR_REL = os.path.join("pipeline_output", "DATASINK", "DKISpace", "dwi_proc")
B0_NAME = "DKI_denoised_gibbs_eddy_b0.nii.gz"

BRAINSTEM_LABEL = 16
FA_MIN = 0.35
B_MAX = 1500.0        # stay in the DTI regime; higher shells need kurtosis
B0_MAX = 50.0
MIN_VOX = 200

FLIPS = {"none": (1, 1, 1), "flip_x": (-1, 1, 1),
         "flip_y": (1, -1, 1), "flip_z": (1, 1, -1)}


def find_one(d, pats):
    for pat in pats:
        hits = sorted(glob.glob(os.path.join(d, "**", pat), recursive=True))
        if hits:
            return hits[0]
    return None


def load_gradients(dwi_dir, n_vol):
    """Locate bval/bvec, tolerating the block-B naming that lacks 'rotated'."""
    bval = find_one(dwi_dir, ["*.bval", "*bvals*"])
    cands = sorted(glob.glob(os.path.join(dwi_dir, "**", "*.bvec"), recursive=True)) \
        + sorted(glob.glob(os.path.join(dwi_dir, "**", "*bvecs*"), recursive=True))
    # prefer an explicitly rotated table, else any table with the right shape
    cands.sort(key=lambda p: (0 if "rotat" in os.path.basename(p).lower() else 1,
                              len(os.path.basename(p))))
    bvec = None
    for c in cands:
        try:
            v = np.loadtxt(c)
        except Exception:
            continue
        if v.ndim == 2 and 3 in v.shape and max(v.shape) == n_vol:
            bvec = c
            break
    if bval is None or bvec is None:
        return None, None, None, None
    b = np.loadtxt(bval).ravel()
    g = np.loadtxt(bvec)
    if g.shape[0] != 3:
        g = g.T
    if len(b) != g.shape[1] or len(b) != n_vol:
        return None, None, None, None
    return b, g, os.path.basename(bval), os.path.basename(bvec)


def fit_dti(S, b, g):
    """Log-linear least-squares DTI at a set of voxels. S is (N, V).

    S0 is estimated as a free parameter (an intercept column) rather than taken
    from b0 volumes. That matters here: the protocol is b = 100/500/1000/2000
    with the b0 stored as a separate file, and block A's 4-D series carries no
    b <= 50 volume at all (52 volumes vs block B's 53). Requiring a b0 made the
    fit return None for every block A subject.

    Model:  ln S = ln S0 - b g^T D g
    """
    use = b <= B_MAX
    if use.sum() < 7:
        return None, None
    Su = np.maximum(S[:, use], 1e-6)
    y = np.log(Su)                                       # (N, Vu)
    gg, bb = g[:, use], b[use]
    X = np.stack([np.ones_like(bb),
                  -bb * gg[0] ** 2, -bb * gg[1] ** 2, -bb * gg[2] ** 2,
                  -2 * bb * gg[0] * gg[1], -2 * bb * gg[0] * gg[2],
                  -2 * bb * gg[1] * gg[2]], axis=1)      # (Vu, 7)
    if np.linalg.matrix_rank(X) < 7:
        return None, None
    coef, *_ = np.linalg.lstsq(X, y.T, rcond=None)       # (7, N); [0] = ln S0
    n = S.shape[0]
    D = np.empty((n, 3, 3))
    D[:, 0, 0], D[:, 1, 1], D[:, 2, 2] = coef[1], coef[2], coef[3]
    D[:, 0, 1] = D[:, 1, 0] = coef[4]
    D[:, 0, 2] = D[:, 2, 0] = coef[5]
    D[:, 1, 2] = D[:, 2, 1] = coef[6]
    w, V = np.linalg.eigh(D)
    md = w.mean(axis=1, keepdims=True)
    num = np.sqrt(((w - md) ** 2).sum(axis=1))
    den = np.sqrt((w ** 2).sum(axis=1))
    fa = np.sqrt(1.5) * num / np.where(den > 1e-12, den, np.nan)
    return V[:, :, -1], fa


def shape_axis(coords):
    c = coords - coords.mean(axis=0)
    w, V = np.linalg.eigh((c.T @ c) / max(len(c) - 1, 1))
    return V[:, -1], (float(w[-1] / w[-2]) if w[-2] > 1e-12 else np.inf)


def warp_mask_to_dwi(seg_p, b0_p, aff_p, iwarp_p, tmp):
    """Brainstem mask -> DWI space, nearest neighbour. Masks carry no
    orientation information, so registration error is immaterial here."""
    seg = nib.load(seg_p)
    m = (np.asarray(seg.dataobj) == BRAINSTEM_LABEL).astype(np.uint8)
    src = os.path.join(tmp, "bs.nii.gz")
    nib.save(nib.Nifti1Image(m, seg.affine), src)
    out = os.path.join(tmp, "bs_dwi.nii.gz")
    cmd = ["antsApplyTransforms", "-d", "3", "-i", src, "-r", b0_p,
           "-o", out, "-n", "NearestNeighbor",
           "-t", f"[{aff_p},1]"]
    if iwarp_p:
        cmd += ["-t", iwarp_p]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        return None, f"antsApplyTransforms unavailable ({e})"
    if r.returncode != 0:
        return None, (r.stderr.strip() or r.stdout.strip() or "failed")[:160]
    return out, None


def process(subj, subj_dir, dwi_root, fa_min):
    dwi_dir = os.path.join(dwi_root, subj, DWI_DIR_REL)
    b0_p = os.path.join(dwi_dir, B0_NAME)
    if not os.path.exists(b0_p):
        print(f"[{subj}] no b0; skipped")
        return None
    dwi_p = find_one(dwi_dir, ["*eddy.nii*", "*denoised_gibbs_eddy.nii*",
                               "*preproc*.nii*", "*.nii.gz"])
    if dwi_p is None or dwi_p == b0_p:
        cands = [p for p in sorted(glob.glob(os.path.join(dwi_dir, "*.nii*")))
                 if len(nib.load(p).shape) == 4]
        dwi_p = cands[0] if cands else None
    if dwi_p is None:
        print(f"[{subj}] no 4D DWI found in {dwi_dir}; skipped")
        return None
    img = nib.load(dwi_p)
    if len(img.shape) != 4:
        print(f"[{subj}] {os.path.basename(dwi_p)} is not 4D; skipped")
        return None
    b, g, bvn, bvcn = load_gradients(os.path.dirname(dwi_p), img.shape[3])
    if b is None:
        print(f"[{subj}] no matching bval/bvec for {img.shape[3]} volumes; skipped")
        return None

    seg_p = find_one(subj_dir, ["*_ses-01_synthseg.nii.gz"])
    aff_p = find_one(subj_dir, ["*from_b0_0GenericAffine.mat"])
    iw = sorted(glob.glob(os.path.join(subj_dir, "**",
                                       "*from_b0_1InverseWarp.nii.gz"),
                          recursive=True))
    if not (seg_p and aff_p):
        print(f"[{subj}] no SynthSeg / affine; skipped")
        return None

    tmp = tempfile.mkdtemp(prefix=f"nat_{subj}_")
    try:
        mp, err = warp_mask_to_dwi(seg_p, b0_p, aff_p, iw[0] if iw else None, tmp)
        if mp is None:
            print(f"[{subj}] mask warp failed: {err}")
            return None
        mask = np.asarray(nib.load(mp).dataobj) > 0
        if mask.sum() < MIN_VOX:
            print(f"[{subj}] brainstem in DWI space too small ({int(mask.sum())})")
            return None
        idx = np.argwhere(mask)
        zooms = np.asarray(img.header.get_zooms()[:3], float)
        u, lin = shape_axis(idx * zooms)

        data = np.asarray(img.dataobj, np.float32)
        S = data[mask]                                   # (N, V)
        row = {"Subject": subj, "n_vox": int(mask.sum()),
               "shape_linearity": round(lin, 2), "n_vol": int(img.shape[3]),
               "bval_file": bvn, "bvec_file": bvcn,
               "dwi_file": os.path.basename(dwi_p)}
        for name, fl in FLIPS.items():
            gf = g * np.array(fl)[:, None]
            v1, fa = fit_dti(S, b, gf)
            if v1 is None:
                row[f"coh_{name}"] = np.nan
                continue
            d = np.abs(v1 @ (u / np.linalg.norm(u)))
            sel = np.isfinite(d) & np.isfinite(fa) & (fa >= fa_min)
            row[f"coh_{name}"] = float(d[sel].mean()) if sel.sum() > 30 else np.nan
            row[f"nfa_{name}"] = int(sel.sum())
        vals = {k: row.get(f"coh_{k}", np.nan) for k in FLIPS}
        if all(not np.isfinite(v) for v in vals.values()):
            row["best_flip"] = None
            print(f"[{subj}] fit failed for every convention "
                  f"(n_vol={img.shape[3]}, b range {b.min():.0f}-{b.max():.0f})")
            return row
        best = max(FLIPS, key=lambda k: (vals[k] if np.isfinite(vals[k]) else -1))
        row["best_flip"] = best
        print(f"[{subj}] " + "  ".join(
            f"{k} {row.get(f'coh_{k}', float('nan')):.3f}" for k in FLIPS)
            + f"  -> {best}   (bvec: {bvcn})", flush=True)
        return row
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=ROOT_DIR)
    ap.add_argument("--dwi-root", default=DWI_ROOT)
    ap.add_argument("--fa-min", type=float, default=FA_MIN)
    ap.add_argument("--block-split", type=int, default=149)
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--out", default="native_tensor_frame_check.csv")
    a = ap.parse_args()

    if shutil.which("antsApplyTransforms") is None:
        print("[!] antsApplyTransforms not on PATH -- needed to bring the "
              "brainstem mask into DWI space (mask only, nearest neighbour).")

    subs = a.subjects or sorted(
        s for s in os.listdir(a.root)
        if s.startswith("sub-") and os.path.isdir(os.path.join(a.root, s)))
    print(f"{len(subs)} subjects | native-space DTI, brainstem shape axis, "
          f"FA >= {a.fa_min}")
    print("no tensor is transformed; the registration moves a binary mask only")
    print("-" * 78)

    rows = [r for r in (process(s, os.path.join(a.root, s), a.dwi_root, a.fa_min)
                        for s in subs) if r]
    if not rows:
        raise SystemExit("no subjects processed")
    d = pd.DataFrame(rows)
    d["num"] = d.Subject.str.extract(r"(\d+)").astype(int)
    d["block"] = np.where(d.num >= a.block_split, "B", "A")

    print("-" * 78)
    print(f"n = {len(d)}   A = {(d.block=='A').sum()}   B = {(d.block=='B').sum()}")
    if "n_vol" in d and d.groupby("block").n_vol.nunique().max() >= 1:
        vc = {b: sorted(g.n_vol.unique().tolist()) for b, g in d.groupby("block")}
        if len({tuple(v) for v in vc.values()}) > 1:
            print(f"\n[!] volume counts differ by block: {vc}")
            print("    Confirm this is only an extra b0 volume and not a")
            print("    difference in the diffusion sampling itself.")
    print("\nCOHERENCE by gradient sign convention (median)")
    print(f"  {'convention':12s} {'all':>8s} {'block A':>9s} {'block B':>9s} "
          f"{'A vs B p':>10s}")
    for k in FLIPS:
        c = f"coh_{k}"
        ga, gb = d[d.block == "A"][c].dropna(), d[d.block == "B"][c].dropna()
        p = (stats.mannwhitneyu(ga, gb).pvalue
             if len(ga) > 2 and len(gb) > 2 else np.nan)
        print(f"  {k:12s} {d[c].median():8.4f} {ga.median():9.4f} "
              f"{gb.median():9.4f} {p:10.2e}")

    print(f"\nbest convention: "
          + ", ".join(f"{k}={v}" for k, v in d.best_flip.value_counts().items()))
    print("  by block:")
    print(pd.crosstab(d.block, d.best_flip).to_string())

    # The verdict rests on the as-is convention, compared across blocks. The
    # gradient tables are identical between blocks, so a convention error is
    # common to both and cannot generate a block difference.
    c = "coh_none"
    ga, gb = d[d.block == "A"][c].dropna(), d[d.block == "B"][c].dropna()
    p = (stats.mannwhitneyu(ga, gb).pvalue
         if len(ga) > 2 and len(gb) > 2 else np.nan)
    print("\nVERDICT  (native-space coherence, gradients as stored)")
    print(f"  block A {ga.median():.4f}   block B {gb.median():.4f}   p={p:.2e}")
    if not np.isfinite(p):
        print("  Too few subjects per block to conclude.")
    elif p < 0.05:
        print("  The tensors DIFFER BY BLOCK in native space, before any warp or")
        print("  frame composition. The problem is at the source -- in the data or")
        print("  the fit -- not in the reorientation code. Note the direction:")
        print("  whichever block scores lower has the less well-oriented tensors.")
        print("  MD/FA/AD/RD remain unaffected: they are rotation-invariant.")
    else:
        print("  No block difference in native space. The tensors are equally")
        print("  well oriented at the source, so everything that split by block")
        print("  arose downstream, in the warp or the frame composition. Only")
        print("  T1w-space Dpar is affected; every rotation-invariant metric")
        print("  (MD, FA, AD, RD) is untouched. That closes the methods loose end.")

    lo = d[[f"coh_{k}" for k in FLIPS]].max(axis=1)
    print(f"\n  reference strength: best-of-conventions coherence, median "
          f"{lo.median():.3f} (random = 0.500)")

    if lo.median() < 0.60:
        print("  [!] below 0.60 -- the reference lacks discriminating power here.")
        print("      Treat the block comparison as indicative, not conclusive.")

    d.to_csv(a.out, index=False)
    print(f"\n-> {a.out}")


if __name__ == "__main__":
    main()
