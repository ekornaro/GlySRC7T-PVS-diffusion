#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step04_fit_data_to_DKI_model.py
===============================
Re-run the SAME standard DKI fit already used for the paper, and additionally
export the diffusion tensor D, its eigenvalues and its eigenvectors -- so the
PVS-orientation analysis uses the DKI diffusion tensor rather than a separate
DTI fit.

WHY THIS IS NOT A CHANGE TO THE PAPER
-------------------------------------
fwdki.fit_standard_dki does exactly:

    model = dki.DiffusionKurtosisModel(gtab, fit_method="WLS")
    fit   = model.fit(data, mask=mask)
    FA, MD, AD, RD = fit.fa, fit.md, fit.ad, fit.rd

The eigen-decomposition was ALREADY computed inside that fit -- fit.evals and
fit.evecs are what fit.fa/md/ad/rd are derived from. It simply was not written to
disk. This script re-runs the identical model, identical fit_method, identical
gradient table, and saves those internals. The scalar maps it produces should be
numerically indistinguishable from the ones already in the manuscript.

That is verified, not assumed: if an existing FA/MD map is found on disk, the
script correlates the new one against it and reports max absolute difference. A
mismatch means something about the inputs differs and you must resolve it BEFORE
using the tensor.

WHY THE DKI TENSOR IS THE RIGHT ONE (a point worth making to the reviewer)
-------------------------------------------------------------------------
The DKI diffusion tensor is the b->0 limit of the signal, estimated with the
non-Gaussian term modelled explicitly. A DTI fit over the same shells absorbs
kurtosis into D and so UNDERESTIMATES diffusivity. For a compartment argued to
contain fast, relatively unhindered fluid, the DKI-derived D is the better
estimate, not a compromise. Sourcing D from DKI also removes any inconsistency
with the rest of the analysis.

ON THE b=0 ASYMMETRY
--------------------
All volumes are used, including the b=0 present in one cohort only. That is
deliberate: the goal is to REPRODUCE the existing maps, and the A/B test showed
uncorrected DKI metrics are insensitive to that volume (r >= 0.94, median change
< 4%) because the ordinary DKI design matrix has a free intercept. The problem
was specific to the free-water stage, which is not run here.

OUTPUTS (per subject, into <subj>/<OUT_SUBDIR>/)
-----------------------------------------------
  Dxx Dyy Dzz Dxy Dxz Dyz   6 unique tensor components  -> warp these
  L1 L2 L3                  eigenvalues, descending
  V1 V2 V3                  eigenvectors (4D, 3 components each)
  FA MD AD RD               recomputed, for verification only

USAGE
-----
  python make_dki_tensor.py --subjects sub-149            # pilot, check the QC
  python make_dki_tensor.py                               # whole cohort
"""

import os
import glob
import argparse
import numpy as np
import nibabel as nib
import pandas as pd

from validate_s0_extrapolation import DEFAULT_ROOT, find_subject_files, load_subject

OUT_SUBDIR = "pipeline_output/DATASINK/DKISpace/tensor"
# where the existing scalar maps live, for the verification step
EXISTING_MAP_DIRS = ["pipeline_output/DATASINK/DKISpace/dki",
                     "pipeline_output/DATASINK/DKISpace/dkifit",
                     "pipeline_output/DATASINK/DKISpace"]
EXISTING_PATTERNS = {"FA": ["dkifit_FA.nii.gz", "FA.nii.gz"],
                     "MD": ["dkifit_MD.nii.gz", "MD.nii.gz"]}

COMPONENT_IDX = {"Dxx": (0, 0), "Dyy": (1, 1), "Dzz": (2, 2),
                 "Dxy": (0, 1), "Dxz": (0, 2), "Dyz": (1, 2)}


def _find_existing(subj_dir, metric):
    for d in EXISTING_MAP_DIRS:
        for pat in EXISTING_PATTERNS.get(metric, []):
            p = os.path.join(subj_dir, d, pat)
            if os.path.exists(p):
                return p
    return None

def _save(path, arr, ref, overwrite):
    if os.path.exists(path) and not overwrite:
        return
    hdr = ref.header.copy()
    hdr.set_data_shape(arr.shape)
    z = ref.header.get_zooms()[:3]
    hdr.set_zooms(tuple(z) + ((1.0,) if arr.ndim == 4 else ()))
    nib.save(nib.Nifti1Image(np.nan_to_num(arr).astype(np.float32), ref.affine, hdr), path)

def process_subject(subj_dir, overwrite=False, save_scalars=True):
    from dipy.reconst import dki
    from dipy.core.gradients import gradient_table

    sub = os.path.basename(subj_dir.rstrip("/"))
    files = find_subject_files(subj_dir)
    img, data, bvals, bvecs, mask, _ = load_subject(files)

    gtab = gradient_table(bvals, bvecs=bvecs, b0_threshold=50.0)
    model = dki.DiffusionKurtosisModel(gtab, fit_method="WLS")
    fit = model.fit(np.nan_to_num(data), mask=mask)

    D = fit.quadratic_form                      # (X, Y, Z, 3, 3)
    evals = fit.evals                           # (X, Y, Z, 3) descending
    evecs = fit.evecs                           # (X, Y, Z, 3, 3), columns

    outdir = os.path.join(subj_dir, OUT_SUBDIR)
    os.makedirs(outdir, exist_ok=True)

    for name, (i, j) in COMPONENT_IDX.items():
        _save(os.path.join(outdir, f"{name}.nii.gz"), D[..., i, j], img, overwrite)
    for k in range(3):
        _save(os.path.join(outdir, f"L{k+1}.nii.gz"), evals[..., k], img, overwrite)
        _save(os.path.join(outdir, f"V{k+1}.nii.gz"), evecs[..., :, k], img, overwrite)

    rec = {"Subject": sub, "n_vol": int(bvals.size),
           "n_b0": int((np.asarray(bvals, float) <= 50).sum()),
           "n_vox_mask": int(mask.sum()), "status": "ok"}

    if save_scalars:
        for name, arr in (("FA", fit.fa), ("MD", fit.md),
                          ("AD", fit.ad), ("RD", fit.rd)):
            _save(os.path.join(outdir, f"{name}_recomputed.nii.gz"), arr, img, overwrite)

    # ---- verification against the maps already used in the manuscript ----
    for metric, new in (("FA", fit.fa), ("MD", fit.md)):
        p = _find_existing(subj_dir, metric)
        if not p:
            rec[f"{metric}_check"] = "no existing map found"
            continue
        old = np.asarray(nib.load(p).dataobj, float)
        if old.shape != new.shape:
            rec[f"{metric}_check"] = f"shape mismatch {old.shape} vs {new.shape}"
            continue
        m = mask & np.isfinite(old) & np.isfinite(new) & (old != 0) & (np.asarray(new) != 0)
        a, b = old[m], np.asarray(new)[m]
        rec[f"{metric}_r"] = float(np.corrcoef(a, b)[0, 1]) if a.size > 10 else np.nan
        rec[f"{metric}_maxabsdiff"] = float(np.max(np.abs(a - b)))
        rec[f"{metric}_median_old"] = float(np.median(a))
        rec[f"{metric}_median_new"] = float(np.median(b))

    # ---- internal consistency: does D reproduce its own eigenvalues? ----
    sm = mask & (evals[..., 0] > 0)
    tr = np.trace(D, axis1=-2, axis2=-1)[sm]
    rec["trace_vs_evalsum_maxdiff"] = float(np.max(np.abs(tr - evals[sm].sum(axis=-1))))
    rec["frac_negative_L3"] = float(np.mean(evals[mask][..., 2] < 0))
    return rec


def run_cohort(root=DEFAULT_ROOT, subjects=None, out="./dki_tensor_qc",
               overwrite=False):
    os.makedirs(out, exist_ok=True)
    dirs = ([os.path.join(root, s) for s in subjects] if subjects
            else sorted(glob.glob(os.path.join(root, "sub-*"))))
    print(f"{len(dirs)} subjects | refitting standard DKI (WLS) to export D")

    recs = []
    for d in dirs:
        sub = os.path.basename(d)
        try:
            r = process_subject(d, overwrite)
            recs.append(r)
            chk = (f"FA r={r['FA_r']:.6f} maxdiff={r['FA_maxabsdiff']:.2e}"
                   if "FA_r" in r else r.get("FA_check", ""))
            print(f"  {sub}: ok  {chk}")
        except Exception as exc:                      # per-subject isolation
            print(f"  {sub}: FAILED -- {exc}")
            recs.append({"Subject": sub, "status": f"FAILED: {exc}"})

    df = pd.DataFrame(recs)
    df.to_csv(os.path.join(out, "dki_tensor_qc.csv"), index=False)

    print("\n=== VERIFICATION: do the recomputed maps match the published ones? ===")
    for metric in ("FA", "MD"):
        c = f"{metric}_r"
        if c in df:
            r = df[c].dropna()
            md_ = df[f"{metric}_maxabsdiff"].dropna()
            print(f"{metric}: r = {r.min():.6f}-{r.max():.6f} | "
                  f"max abs diff = {md_.max():.3e}")
            if len(r) and r.min() < 0.999:
                print(f"  [!] {metric} does NOT reproduce. The inputs differ from")
                print( "      what produced the published maps -- check that")
                print( "      PATTERNS points at the same DWI series. Resolve this")
                print( "      before using the tensor.")
            else:
                print("  reproduces -- the tensor is consistent with the paper.")
        else:
            print(f"{metric}: no existing map found to compare against. Set")
            print( "  EXISTING_MAP_DIRS / EXISTING_PATTERNS to enable verification.")

    if "frac_negative_L3" in df:
        v = df["frac_negative_L3"].dropna()
        print(f"\nnegative L3 fraction: median {v.median()*100:.2f}%, "
              f"max {v.max()*100:.2f}%")
        print("  a few negative eigenvalues are normal at this SNR; a large or")
        print("  group-dependent fraction is not -- check it against group.")
    print(f"\nWritten to {out}")
    print(f"\nNEXT: warp Dxx..Dyz (and V1 for the reorientation QC) into")
    print( "  revised_DKI_in_T1w_space with Linear interpolation, then run")
    print( "  pvs_orientation_analysis.py --analyze --qc-only")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=DEFAULT_ROOT)
    ap.add_argument("--subjects", nargs="*", default=None)
    ap.add_argument("--out", default="./dki_tensor_qc")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    run_cohort(a.root, a.subjects, a.out, a.overwrite)


if __name__ == "__main__":
    main()
