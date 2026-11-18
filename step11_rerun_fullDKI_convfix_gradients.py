#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step11_rerun_fullDKI_convfix_gradients.py (v2 -- convfix gradients)
===================================================================
Re-run the standard DKI fit and export the diffusion tensor, NOW FITTING WITH
THE CONVENTION-CORRECTED GRADIENTS written by step20e
(DKI_corrected_convfix.bvec/.bval: z row negated for block A, identity copy
for block B; evidence = step20d, 41/41 and 24/24 unanimous).

v2 GRADIENT HANDLING -- VERIFIED SUBSTITUTION, NOT A PATH CHANGE
----------------------------------------------------------------
For every subject the script:
  1. requires the _convfix pair AND step20e's bvec_convention_provenance.json
     to exist (refuses otherwise -- fitting some subjects with original
     gradients would silently reintroduce the block mixture);
  2. loads the _convfix pair and CHECKS it against the gradients the loader
     returned: bvals identical, and bvecs identical (action=identity) or
     identical-with-z-negated (action=flip_z), per the provenance json.
     Any other relation -> refuse. This proves the fit input is exactly the
     step20e correction of the same gradient set the pipeline used;
  3. fits with the _convfix gradients and records grad_source/grad_action in
     the QC csv.

WHAT MUST NOT CHANGE, AND IS CHECKED
------------------------------------
A bvec z-flip reflects the fitted tensor in place; eigenvalues are untouched.
FA/MD/AD/RD are therefore expected to match the published maps EXACTLY as
before -- the r >= 0.999 verification below is the guard that the manuscript's
scalar results are unaffected by this refit. What changes is the tensor's
ORIENTATION (components/eigenvectors), which is the entire point.
Run with --overwrite so the previously exported (mis-oriented, block A)
tensor components are actually replaced.

ACCEPTANCE CRITERION FOR THE REFIT
----------------------------------
Re-run step20d on the refit tensors afterwards: every subject, both blocks,
should now win 'none' with block-B-level tract lengths. Only then does the
orientation branch (step15 -> step16 -> theta -> v7) restart.

(The v1 docstring's 'ON THE b=0 ASYMMETRY' section is retired: since the
reprocessing, all three cohorts carry their b=0 and are fit symmetrically.)

WHY THE DKI TENSOR IS THE RIGHT ONE (unchanged from v1)
-------------------------------------------------------------------------
The DKI diffusion tensor is the b->0 limit of the signal, estimated with the
non-Gaussian term modelled explicitly. A DTI fit over the same shells absorbs
kurtosis into D and so UNDERESTIMATES diffusivity. For a compartment argued to
contain fast, relatively unhindered fluid, the DKI-derived D is the better
estimate, not a compromise. Sourcing D from DKI also removes any inconsistency
with the rest of the analysis.

ON THE b=0 ASYMMETRY
--------------------
Retired -- see above.

USAGE
-----
  python step11_rerun_fullDKI_convfix_gradients.py --subjects sub-115 --overwrite   # pilot
  python step11_rerun_fullDKI_convfix_gradients.py --overwrite                      # cohort
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

GRAD_DIR_SUBPATH = "pipeline_output/DATASINK/DKISpace/dwi_proc"
FIX_STEM = "DKI_corrected_convfix"
PROVENANCE_NAME = "bvec_convention_provenance.json"


def load_convfix_gradients(subj_dir, bvals_loaded, bvecs_loaded):
    """
    Load the step20e convention-corrected gradients and PROVE they are the
    correction of the gradient set the loader returned.

    Returns (bvals, bvecs, action). Raises on any inconsistency -- a subject
    must never be silently fit with unverified gradients.
    """
    import json
    d = os.path.join(subj_dir, GRAD_DIR_SUBPATH)
    fb = os.path.join(d, FIX_STEM + ".bvec")
    fv = os.path.join(d, FIX_STEM + ".bval")
    fp = os.path.join(d, PROVENANCE_NAME)
    missing = [os.path.basename(p) for p in (fb, fv, fp) if not os.path.exists(p)]
    if missing:
        raise RuntimeError(f"convfix gradients missing ({missing}); "
                           f"run step20e for this subject first")
    with open(fp) as fh:
        action = json.load(fh).get("action")
    if action not in ("flip_z", "identity"):
        raise RuntimeError(f"unrecognized provenance action '{action}'")

    bvals_fix = np.loadtxt(fv).ravel()
    B = np.loadtxt(fb)
    if B.ndim != 2 or 3 not in B.shape:
        raise RuntimeError(f"convfix bvec has unexpected shape {B.shape}")
    b3 = B if B.shape[0] == 3 else B.T                       # (3, N)

    L = np.asarray(bvecs_loaded, float)
    l3 = L if (L.ndim == 2 and L.shape[0] == 3) else L.T     # (3, N)
    if b3.shape != l3.shape:
        raise RuntimeError(f"convfix bvec {b3.shape} != loaded {l3.shape}")
    if not np.allclose(bvals_fix, np.asarray(bvals_loaded, float), atol=1e-3):
        raise RuntimeError("convfix bvals differ from loaded bvals")

    expect = l3.copy()
    if action == "flip_z":
        expect[2, :] *= -1.0
    if not np.allclose(b3, expect, atol=1e-4):
        raise RuntimeError(
            "convfix bvecs are not the loaded bvecs with the provenance "
            f"action '{action}' applied -- gradient sources are inconsistent; "
            "resolve before fitting")

    bvecs_out = b3 if (L.ndim == 2 and L.shape[0] == 3) else b3.T
    return bvals_fix, bvecs_out, action


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

    # v2: fit with the step20e convention-corrected gradients (verified)
    bvals, bvecs, grad_action = load_convfix_gradients(subj_dir, bvals, bvecs)

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
           "n_vox_mask": int(mask.sum()),
           "grad_source": FIX_STEM, "grad_action": grad_action,
           "status": "ok"}

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
    print(f"{len(dirs)} subjects | refitting standard DKI (WLS) with "
          f"{FIX_STEM} gradients (step20e) to export D")

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
    if "grad_action" in df:
        print("\ngradient actions used: "
              + ", ".join(f"{k}={v}" for k, v in
                          df["grad_action"].value_counts().items()))
    print(f"\nWritten to {out}")
    print("\nNEXT: re-run step20d on these refit tensors -- EVERY subject in "
          "both blocks\n  should now win 'none' at block-B tract lengths. Only "
          "then: step15 (warp\n  Dxx..Dyz + V1, --overwrite), step16 "
          "(--rotate-u off, --overwrite), theta,\n  and the v7 re-extraction "
          "with the orientation metrics.")
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
