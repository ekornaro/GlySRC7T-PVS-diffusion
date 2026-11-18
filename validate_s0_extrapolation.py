#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate_s0_extrapolation.py
============================
Validate no-b0 S0 estimation against ground truth, using the RHI cohort as the
only group that actually has a measured b = 0 volume.

WHY
---
HC and PPCS have 0 x b=0 and 6 x b=100.
RHI          have 1 x b=0 and 6 x b=100.
Stage 2 of the free-water pipeline needs S0 as an ABSOLUTE quantity, so the two
cohorts have been going through two different S0 pathways. RHI is therefore the
only place where the extrapolator can be checked against truth.

WHAT IT DOES
------------
For every subject that HAS a b=0 volume:
  1. Sets the b=0 volume(s) aside as ground truth  ->  S0_measured
  2. Re-estimates S0 from the remaining shells only, with SEVERAL estimators
  3. Reports ratio = S0_estimated / S0_measured, voxelwise, summarised in mask(s)
  4. Bins the ratio by MD decile, because the bias is tissue dependent and the
     whole extension-to-HC/PPCS argument rests on the bias being stable across
     the tissue range that all three groups occupy.

ESTIMATORS REPORTED
-------------------
  dti_low      log-linear DTI intercept, low shells only  (b=100, 500 by default)
  dti_all      log-linear DTI intercept, all non-zero shells
  powder_lin   linear fit to the spherical-mean (powder) signal, low shells
  powder_quad  quadratic (kurtosis-aware) fit to the powder signal, all shells

>>> IMPORTANT <<<
The point of this script is to characterise YOUR pipeline's estimator. Whichever
of the four matches what EvgBrainPipe Stage 1 actually does is the row that
matters; the others are there so you can see whether a different estimator would
have been better. If Stage 1 does something not covered here, import that
function and pass it via `extra_estimators` rather than trusting a lookalike.

USAGE
-----
  python validate_s0_extrapolation.py --root /path/to/processed --out ./s0_qc
  python validate_s0_extrapolation.py --subjects sub-001 sub-002 --save-maps

Importable:
  from validate_s0_extrapolation import validate_subject, run_cohort
"""

import os
import glob
import argparse
import numpy as np
import nibabel as nib
import pandas as pd
from pipeline_config import cfg

# ==========================================================================
# CONFIGURATION  -- confirm these against your tree before the first real run
# ==========================================================================
DEFAULT_ROOT = cfg("data_root")

# Anchored patterns. NO leading wildcards -- that is exactly how `f` matched
# `_mf` last time. Each pattern is relative to the subject directory.
PATTERNS = {
    "dwi":   ["pipeline_output/DATASINK/DKISpace/dwi_proc/DKI_denoised_gibbs_eddy.nii.gz", "pipeline_output/DATASINK/DKISpace/dwi_proc/DKI_denoised_gibbs_eddy.nii.gz"],
    "bval":  ["pipeline_output/DATASINK/DKISpace/dwi_proc/DKI_corrected.bval", "pipeline_output/DATASINK/DKISpace/dwi_proc/DKI_corrected.bval"],
    "bvec":  ["pipeline_output/DATASINK/DKISpace/dwi_proc/DKI_corrected.bvec", "pipeline_output/DATASINK/DKISpace/dwi_proc/DKI_corrected.bvec"],
    "mask":  ["pipeline_output/DATASINK/DKISpace/masks/ANTS_T1_brain_mask.nii.gz", "pipeline_output/DATASINK/DKISpace/masks/ANTS_T1_brain_mask.nii.gz"],
    #"wm":    ["pipeline_output/DATASINK/DKISpace/dwi_proc/*_label-WM_probseg.nii.gz", "pipeline_output/DATASINK/DKISpace/dwi_proc/*_WM*.nii.gz"],   # optional
}

B0_TOL = 10.0        # b <= B0_TOL is treated as b = 0
SHELL_TOL = 40.0     # tolerance for grouping b-values into shells
LOW_SHELLS = (100.0, 500.0)
RIDGE = 1e-8         # per-column ridge, guards singular XtWX
CHUNK = 100_000      # voxels per batch in the WLS solve


# ==========================================================================
# FILE DISCOVERY
# ==========================================================================
def _first_match(subj_dir, patterns, required=True, label=""):
    for p in patterns:
        hits = sorted(glob.glob(os.path.join(subj_dir, p)))
        if len(hits) == 1:
            return hits[0]
        if len(hits) > 1:
            raise RuntimeError(
                f"{label}: pattern '{p}' is ambiguous in {subj_dir}\n  " +
                "\n  ".join(hits))
    if required:
        raise FileNotFoundError(
            f"{label}: nothing matched {patterns} in {subj_dir}")
    return None


def find_subject_files(subj_dir, patterns=None):
    """Resolve the inputs for one subject. Fails loudly, never silently."""
    patterns = patterns or PATTERNS
    out = {}
    for key in ("dwi", "bval", "bvec", "mask"):
        out[key] = _first_match(subj_dir, patterns[key], True, key)
    out["wm"] = _first_match(subj_dir, patterns.get("wm", []), False, "wm")
    return out


# ==========================================================================
# LOADING
# ==========================================================================
def load_subject(files):
    img = nib.load(files["dwi"])
    data = np.asarray(img.dataobj, dtype=np.float64)
    bvals = np.loadtxt(files["bval"]).ravel().astype(float)
    bvecs = np.loadtxt(files["bvec"]).astype(float)
    if bvecs.shape[0] == 3:
        bvecs = bvecs.T
    if not (data.shape[-1] == bvals.size == bvecs.shape[0]):
        raise ValueError(
            f"volume/gradient mismatch: dwi={data.shape[-1]}, "
            f"bval={bvals.size}, bvec={bvecs.shape[0]}")
    # normalise gradient directions, leave b=0 rows as zeros
    n = np.linalg.norm(bvecs, axis=1)
    ok = n > 1e-6
    bvecs[ok] /= n[ok, None]

    mask = np.asarray(nib.load(files["mask"]).dataobj) > 0.5
    wm = None
    if files.get("wm"):
        wm = np.asarray(nib.load(files["wm"]).dataobj) > 0.5
        wm &= mask
    return img, data, bvals, bvecs, mask, wm


def shell_of(bvals, tol=SHELL_TOL):
    """Round b-values onto their nominal shells."""
    uniq = []
    out = np.zeros_like(bvals)
    for b in bvals:
        hit = next((u for u in uniq if abs(u - b) <= tol), None)
        if hit is None:
            uniq.append(b)
            hit = b
        out[bvals == b] = hit
    # snap to the median of each cluster for stability
    for u in set(out):
        out[out == u] = np.median(bvals[out == u])
    return out


# ==========================================================================
# ESTIMATORS
# ==========================================================================
def _design_dti(bvals, bvecs):
    """log S = logS0 - b g^T D g   ->  columns [1, -b gx^2, ...]"""
    b, g = bvals, bvecs
    return np.column_stack([
        np.ones_like(b),
        -b * g[:, 0] ** 2,
        -b * g[:, 1] ** 2,
        -b * g[:, 2] ** 2,
        -2.0 * b * g[:, 0] * g[:, 1],
        -2.0 * b * g[:, 0] * g[:, 2],
        -2.0 * b * g[:, 1] * g[:, 2],
    ])


def _wls_loglinear(X, S, n_iter=2, ridge=RIDGE, chunk=CHUNK):
    """
    Batched weighted log-linear least squares.
      X : (n_vol, n_par)
      S : (n_vox, n_vol) strictly positive signal
    Weights are S_pred^2, the standard variance-stabilising choice for
    log-transformed DWI data. Returns params (n_vox, n_par).
    """
    nvox, nvol = S.shape
    npar = X.shape[1]
    logS = np.log(S)
    pinvX = np.linalg.pinv(X)
    params = logS @ pinvX.T                                   # OLS start

    for _ in range(n_iter):
        out = np.empty_like(params)
        for i in range(0, nvox, chunk):
            sl = slice(i, min(i + chunk, nvox))
            w = np.exp(2.0 * (params[sl] @ X.T))              # (c, n_vol)
            w /= w.max(axis=1, keepdims=True)                 # scale-free
            A = np.einsum("np,cn,nq->cpq", X, w, X, optimize=True)
            A[:, np.arange(npar), np.arange(npar)] += ridge * np.maximum(
                np.einsum("cpp->cp", A), 1e-12)
            rhs = np.einsum("np,cn,cn->cp", X, w, logS[sl], optimize=True)
            out[sl] = np.linalg.solve(A, rhs[..., None])[..., 0]
        params = out
    return params


def _powder(S, shells, keep):
    """Spherical mean per shell. Returns (n_vox, n_shell), (n_shell,)"""
    uniq = np.array(sorted(set(shells[keep])))
    P = np.stack([S[:, keep][:, shells[keep] == u].mean(axis=1) for u in uniq],
                 axis=1)
    return P, uniq


def estimate_s0(S, bvals, bvecs, shells, method):
    """
    S : (n_vox, n_vol), positive. Returns (S0 map (n_vox,), MD map or None).
    """
    nz = shells > B0_TOL
    low = nz & np.isin(shells, np.array(LOW_SHELLS))

    if method == "dti_low":
        keep = low
        if keep.sum() < 8:
            return None, None
        p = _wls_loglinear(_design_dti(bvals[keep], bvecs[keep]), S[:, keep])
        return np.exp(p[:, 0]), (p[:, 1] + p[:, 2] + p[:, 3]) / 3.0

    if method == "dti_all":
        keep = nz
        p = _wls_loglinear(_design_dti(bvals[keep], bvecs[keep]), S[:, keep])
        return np.exp(p[:, 0]), (p[:, 1] + p[:, 2] + p[:, 3]) / 3.0

    if method in ("powder_lin", "powder_quad"):
        keep = low if method == "powder_lin" else nz
        P, ub = _powder(S, shells, keep)
        deg = 1 if method == "powder_lin" else 2
        if ub.size < deg + 1:
            return None, None
        X = np.column_stack([np.ones_like(ub)] +
                            [(-ub) ** k if k == 1 else ub ** k
                             for k in range(1, deg + 1)])
        p = np.log(P) @ np.linalg.pinv(X).T
        return np.exp(p[:, 0]), p[:, 1]

    raise ValueError(f"unknown method '{method}'")


ESTIMATORS = ["dti_low", "dti_all", "powder_lin", "powder_quad"]


# ==========================================================================
# PER-SUBJECT VALIDATION
# ==========================================================================
def _summarise(ratio, name, tag):
    q = np.nanpercentile(ratio, [5, 25, 50, 75, 95])
    return {
        f"{tag}_n": int(np.isfinite(ratio).sum()),
        f"{tag}_median": q[2],
        f"{tag}_iqr": q[3] - q[1],
        f"{tag}_p5": q[0],
        f"{tag}_p95": q[4],
        f"{tag}_pct_within_2pct": float(np.mean(np.abs(ratio - 1) < 0.02) * 100),
    }


def validate_subject(subj_dir, save_maps=False, outdir=None,
                     extra_estimators=None):
    """
    Returns (summary_row, md_bin_rows). Both are lists of dicts.
    `extra_estimators` : {name: callable(S, bvals, bvecs, shells) -> (S0, MD)}
      Use this to plug in EvgBrainPipe's actual Stage-1 S0 function.
    """
    sub = os.path.basename(subj_dir.rstrip("/"))
    files = find_subject_files(subj_dir)
    img, data, bvals, bvecs, mask, wm = load_subject(files)
    shells = shell_of(bvals)

    b0 = shells <= B0_TOL
    if b0.sum() == 0:
        return [{"Subject": sub, "status": "no measured b0 -- skipped",
                 "n_b0": 0, "n_vol": int(bvals.size)}], []

    S = data[mask]                                   # (n_vox, n_vol)
    S = np.clip(S, 1e-6, None)
    S0_meas = S[:, b0].mean(axis=1)
    good = np.isfinite(S0_meas) & (S0_meas > 1e-3)

    row = {"Subject": sub, "status": "ok",
           "n_b0": int(b0.sum()), "n_vol": int(bvals.size),
           "shells": ",".join(str(int(s)) for s in sorted(set(shells))),
           "n_vox_mask": int(mask.sum())}

    md_rows = []
    md_ref = None

    for name in ESTIMATORS + list((extra_estimators or {}).keys()):
        try:
            if extra_estimators and name in extra_estimators:
                S0_est, MD = extra_estimators[name](S, bvals, bvecs, shells)
            else:
                S0_est, MD = estimate_s0(S, bvals, bvecs, shells, name)
        except Exception as exc:                      # per-estimator isolation
            row[f"{name}_error"] = str(exc)[:120]
            continue
        if S0_est is None:
            row[f"{name}_error"] = "insufficient shells"
            continue

        ratio = np.where(good, S0_est / S0_meas, np.nan)
        ratio[~np.isfinite(ratio)] = np.nan
        row.update(_summarise(ratio, name, f"{name}_brain"))

        if wm is not None:
            wm_in_mask = wm[mask]
            row.update(_summarise(np.where(wm_in_mask, ratio, np.nan),
                                  name, f"{name}_wm"))

        if md_ref is None and MD is not None:
            md_ref = MD

        # MD-decile profile: is the bias tissue dependent?
        if md_ref is not None:
            edges = np.nanpercentile(md_ref, np.arange(0, 101, 10))
            edges = np.unique(edges)
            if edges.size > 2:
                idx = np.clip(np.digitize(md_ref, edges[1:-1]), 0,
                              edges.size - 2)
                for k in range(edges.size - 1):
                    sel = idx == k
                    if sel.sum() < 50:
                        continue
                    md_rows.append({
                        "Subject": sub, "estimator": name, "md_decile": k,
                        "md_lo": edges[k], "md_hi": edges[k + 1],
                        "n_vox": int(sel.sum()),
                        "ratio_median": float(np.nanmedian(ratio[sel])),
                    })

        if save_maps and outdir:
            vol = np.full(mask.shape, np.nan, dtype=np.float32)
            vol[mask] = ratio
            os.makedirs(outdir, exist_ok=True)
            nib.save(nib.Nifti1Image(vol, img.affine, img.header),
                     os.path.join(outdir, f"{sub}_S0ratio_{name}.nii.gz"))

    return [row], md_rows


# ==========================================================================
# COHORT
# ==========================================================================
def run_cohort(root=DEFAULT_ROOT, out="./s0_qc", subjects=None,
               save_maps=False, extra_estimators=None):
    os.makedirs(out, exist_ok=True)
    dirs = ([os.path.join(root, s) for s in subjects] if subjects
            else sorted(glob.glob(os.path.join(root, "sub-*"))))
    print(f"{len(dirs)} subject directories")

    rows, md_rows = [], []
    for d in dirs:
        sub = os.path.basename(d)
        try:
            r, m = validate_subject(d, save_maps, out, extra_estimators)
            rows += r
            md_rows += m
            st = r[0].get("status", "")
            key = "dti_low_brain_median"
            extra = f"  ratio(dti_low)={r[0][key]:.4f}" if key in r[0] else ""
            print(f"  {sub}: {st}{extra}")
        except Exception as exc:                      # per-subject isolation
            print(f"  {sub}: FAILED -- {exc}")
            rows.append({"Subject": sub, "status": f"FAILED: {exc}"})

    df = pd.DataFrame(rows)
    md = pd.DataFrame(md_rows)
    df.to_csv(os.path.join(out, "s0_validation_summary.csv"), index=False)
    if not md.empty:
        md.to_csv(os.path.join(out, "s0_validation_md_bins.csv"), index=False)

    ok = df[df.get("status", "") == "ok"] if "status" in df else df
    if len(ok):
        print("\nS0_estimated / S0_measured, across subjects with a b=0:")
        print(f"{'estimator':<14}{'median':>9}{'IQR':>9}{'min':>9}{'max':>9}"
              f"{'%vox within 2%':>17}")
        for name in ESTIMATORS + list((extra_estimators or {}).keys()):
            c = f"{name}_brain_median"
            if c not in ok:
                continue
            v = ok[c].dropna()
            w = ok.get(f"{name}_brain_pct_within_2pct", pd.Series(dtype=float))
            print(f"{name:<14}{v.median():>9.4f}"
                  f"{v.quantile(.75) - v.quantile(.25):>9.4f}"
                  f"{v.min():>9.4f}{v.max():>9.4f}"
                  f"{w.dropna().median() if len(w.dropna()) else np.nan:>17.1f}")
        print("\nRead it like this:")
        print("  median ~ 1.00 and tight  -> extrapolator is unbiased; the S0")
        print("                              pathway is not what broke Stage 2.")
        print("  median off 1.00          -> that is the multiplicative bias HC")
        print("                              and PPCS have been carrying.")
        print("  ratio varies across MD deciles -> the bias is tissue dependent,")
        print("     so it does NOT cancel between groups with different tissue")
        print("     composition. Check s0_validation_md_bins.csv before you")
        print("     argue the bias is group-neutral.")
    print(f"\nWritten to {out}")
    return df, md


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=DEFAULT_ROOT)
    ap.add_argument("--out", default="./s0_qc")
    ap.add_argument("--subjects", nargs="*", default=None)
    ap.add_argument("--save-maps", action="store_true")
    a = ap.parse_args()
    run_cohort(a.root, a.out, a.subjects, a.save_maps)


if __name__ == "__main__":
    main()
