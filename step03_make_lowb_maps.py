#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step03_make_lowb_maps.py
========================
Low-b (b <= 500) diffusion metrics answering Reviewer #3, point 2, using an
S0-FREE per-direction ratio design.

WHY v2
------
v1 fitted a log-linear tensor with an intercept to the 12 low-b volumes
(6 at b~100, 6 at b=500). That is 7 parameters from 12 measurements spanning
only 6 directions at b_max = 500, where b*D ~ 0.5. The result had the right
central tendency (median MD matched the powder ADC to 2%) but essentially zero
voxelwise correlation with it -- the classic signature of a fit dominated by
noise rather than one that is broken.

v2 exploits the fact that the SAME 6 directions are sampled at both shells
(verified: |g_lo . g_hi| = I). Taking the per-direction ratio removes S0
analytically:

    ADC(g) = ln( S_lo(g) / S_hi(g) ) / (b_hi - b_lo)

leaving 6 direction-resolved ADCs and a 6-parameter tensor with NO intercept:

    ADC(g) = g^T D g

  design                      obs / par    condition number
  ------------------------------------------------------------
  v1  log-linear + intercept    12 / 7          4.3e+02
  v2  S0-free ratio              6 / 6          2.6e+00

Two orders of magnitude better conditioned, and S0 -- a quantity this dataset
has already shown to be unreliable -- is eliminated rather than estimated.

The design matrix is voxel-independent, so the whole-brain solve is one
pseudo-inverse. This is faster than v1 as well as better posed.

OUTPUTS (per subject, into <subj>/<OUT_SUBDIR>/)
-----------------------------------------------
  lowb_ADC   powder (spherical-mean) ADC. Direction free, 1 parameter from 12
             volumes -- the most robust low-b estimate. PRIMARY.
  lowb_MD    trace/3 of the ratio tensor. Should track lowb_ADC closely.
  lowb_AD    axial diffusivity at low b -- the direct analogue of the AD result
             that carries the rebuttal, in the encoding regime the reviewer
             specifies.
  lowb_RD    radial diffusivity at low b.
  lowb_FA    anisotropy of the ratio tensor.

QC TO READ
----------
  corr_ADC_vs_MD should now be high. If it is still near zero, the low-b data
  cannot support a tensor at all and the powder ADC alone is what you report --
  which still answers point 2, just without the directional breakdown.

USAGE
-----
  python make_lowb_maps_v2.py --subjects sub-102 sub-149
  python make_lowb_maps_v2.py --overwrite
"""

import os
import glob
import argparse
import numpy as np
import nibabel as nib
import pandas as pd

from validate_s0_extrapolation import DEFAULT_ROOT, find_subject_files, load_subject

B0_THRESH = 50.0
BMAX_LOW = 500.0
SHELL_SNAP = 50.0
OUT_SUBDIR = "pipeline_output/DATASINK/DKISpace/lowb"
MATCH_TOL = 0.98            # |g_lo . g_hi| below this = directions not matched


def _nominal(bvals, snap=SHELL_SNAP):
    """Snap b-values onto nominal shells so 100 vs 101 is not a difference."""
    return np.round(np.asarray(bvals, float) / snap) * snap


def _mean_by_direction(S, bvecs):
    """Average volumes sharing a direction. Returns (n_vox, n_dir), (n_dir, 3)."""
    dirs, cols = [], []
    for i, g in enumerate(bvecs):
        hit = next((k for k, d in enumerate(dirs)
                    if abs(float(np.dot(g, d))) > 0.999), None)
        if hit is None:
            dirs.append(g)
            cols.append([i])
        else:
            cols[hit].append(i)
    D = np.array(dirs)
    M = np.stack([S[:, c].mean(axis=1) for c in cols], axis=1)
    return M, D


def _match(lo_dirs, hi_dirs):
    """Pair each hi direction with its lo counterpart. Sign is irrelevant."""
    C = np.abs(lo_dirs @ hi_dirs.T)
    order = C.argmax(axis=0)
    quality = C[order, np.arange(C.shape[1])]
    return order, quality


def _tensor_from_adc(adc, dirs):
    """Solve ADC(g) = g^T D g. Design is voxel-independent -> one pinv."""
    g = dirs
    M = np.column_stack([g[:, 0] ** 2, g[:, 1] ** 2, g[:, 2] ** 2,
                         2 * g[:, 0] * g[:, 1], 2 * g[:, 0] * g[:, 2],
                         2 * g[:, 1] * g[:, 2]])
    d = adc @ np.linalg.pinv(M).T                      # (n_vox, 6)
    n = adc.shape[0]
    D = np.empty((n, 3, 3))
    D[:, 0, 0], D[:, 1, 1], D[:, 2, 2] = d[:, 0], d[:, 1], d[:, 2]
    D[:, 0, 1] = D[:, 1, 0] = d[:, 3]
    D[:, 0, 2] = D[:, 2, 0] = d[:, 4]
    D[:, 1, 2] = D[:, 2, 1] = d[:, 5]
    ev = np.linalg.eigvalsh(np.nan_to_num(D))[:, ::-1]
    return ev, float(np.linalg.cond(M))


def process_subject(subj_dir, overwrite=False):
    sub = os.path.basename(subj_dir.rstrip("/"))
    files = find_subject_files(subj_dir)
    img, data, bvals, bvecs, mask, _ = load_subject(files)
    nb = _nominal(bvals)

    lo_shell = min(u for u in set(nb) if u > B0_THRESH)
    hi_shell = BMAX_LOW
    sel_lo, sel_hi = nb == lo_shell, nb == hi_shell
    rec = {"Subject": sub, "n_vol_total": int(bvals.size),
           "lowb_shells": f"{int(lo_shell)},{int(hi_shell)}",
           "n_vol_lo": int(sel_lo.sum()), "n_vol_hi": int(sel_hi.sum()),
           "n_b0_excluded": int((nb <= B0_THRESH).sum())}
    if sel_lo.sum() < 6 or sel_hi.sum() < 6:
        rec["status"] = "REFUSED: fewer than 6 volumes in a low shell"
        return rec

    idx = np.flatnonzero(mask.ravel())
    flat = data.reshape(-1, data.shape[-1])[idx]
    S_lo, d_lo = _mean_by_direction(np.clip(flat[:, sel_lo], 1e-6, None),
                                    bvecs[sel_lo])
    S_hi, d_hi = _mean_by_direction(np.clip(flat[:, sel_hi], 1e-6, None),
                                    bvecs[sel_hi])

    order, quality = _match(d_lo, d_hi)
    rec["n_dirs_matched"] = int(d_hi.shape[0])
    rec["min_match_quality"] = float(quality.min())
    if quality.min() < MATCH_TOL:
        rec["status"] = (f"REFUSED: direction matching poor "
                         f"(min |g_lo.g_hi| = {quality.min():.3f}). The two "
                         f"shells do not sample the same directions; report "
                         f"powder ADC only.")
        return rec

    db = hi_shell - lo_shell
    adc_dir = np.log(S_lo[:, order] / S_hi) / db        # (n_vox, n_dir)
    powder = np.log(S_lo.mean(axis=1) / S_hi.mean(axis=1)) / db

    ev, cond = _tensor_from_adc(adc_dir, d_hi)
    rec["design_condition_number"] = cond

    md = ev.mean(axis=1)
    ad = ev[:, 0]
    rd = ev[:, 1:].mean(axis=1)
    num = np.sqrt(1.5 * ((ev - md[:, None]) ** 2).sum(axis=1))
    den = np.sqrt((ev ** 2).sum(axis=1))
    fa = np.divide(num, den, out=np.zeros_like(num), where=den > 0)

    outdir = os.path.join(subj_dir, OUT_SUBDIR)
    os.makedirs(outdir, exist_ok=True)
    for name, vec in (("lowb_ADC", powder), ("lowb_MD", md), ("lowb_AD", ad),
                      ("lowb_RD", rd), ("lowb_FA", fa)):
        p = os.path.join(outdir, f"{name}.nii.gz")
        if os.path.exists(p) and not overwrite:
            continue
        vol = np.zeros(mask.size, np.float32)
        vol[idx] = np.nan_to_num(vec, nan=0.0, posinf=0.0, neginf=0.0)
        nib.save(nib.Nifti1Image(vol.reshape(mask.shape), img.affine), p)

    ok = np.isfinite(powder) & np.isfinite(md)
    rec.update({
        "median_lowb_ADC": float(np.median(powder[ok])),
        "median_lowb_MD": float(np.median(md[ok])),
        "median_lowb_AD": float(np.median(ad[ok])),
        "median_lowb_RD": float(np.median(rd[ok])),
        "corr_ADC_vs_MD": float(np.corrcoef(powder[ok], md[ok])[0, 1]),
        "frac_negative_ADC": float(np.mean(powder[ok] < 0)),
        "frac_negative_MD": float(np.mean(md[ok] < 0)),
        "status": "ok",
    })
    return rec


def run_cohort(root=DEFAULT_ROOT, subjects=None, out="./lowb_qc_v2",
               overwrite=False):
    os.makedirs(out, exist_ok=True)
    dirs = ([os.path.join(root, s) for s in subjects] if subjects
            else sorted(glob.glob(os.path.join(root, "sub-*"))))
    print(f"{len(dirs)} subjects | S0-free low-b maps (b <= {BMAX_LOW:.0f})")

    recs = []
    for d in dirs:
        sub = os.path.basename(d)
        try:
            r = process_subject(d, overwrite)
            recs.append(r)
            print(f"  {sub}: {r['status'][:60]:<60} "
                  f"ADC={r.get('median_lowb_ADC', float('nan')):.2e} "
                  f"r={r.get('corr_ADC_vs_MD', float('nan')):.3f}")
        except Exception as exc:                       # per-subject isolation
            print(f"  {sub}: FAILED -- {exc}")
            recs.append({"Subject": sub, "status": f"FAILED: {exc}"})

    df = pd.DataFrame(recs)
    df.to_csv(os.path.join(out, "lowb_maps_qc_v2.csv"), index=False)

    for col, label in (("lowb_shells", "shells"),
                       ("n_vol_lo", "volumes at low shell"),
                       ("n_vol_hi", "volumes at b=500")):
        if col in df:
            u = df[col].dropna().unique()
            flag = "" if len(u) == 1 else "   [!] NOT identical across subjects"
            print(f"{label}: {list(u)}{flag}")
    if "corr_ADC_vs_MD" in df:
        c = df["corr_ADC_vs_MD"].dropna()
        print(f"\npowder ADC vs ratio-tensor MD: r = {c.median():.3f} "
              f"(range {c.min():.3f}-{c.max():.3f})")
        if c.median() > 0.8:
            print("  -> tensor is well posed. lowb_AD and lowb_RD are usable.")
        else:
            print("  -> still poorly posed. Report lowb_ADC alone; it answers")
            print("     point 2 without the directional breakdown.")
    print(f"\nWritten to {out}")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=DEFAULT_ROOT)
    ap.add_argument("--subjects", nargs="*", default=None)
    ap.add_argument("--out", default="./lowb_qc_v2")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    run_cohort(a.root, a.subjects, a.out, a.overwrite)


if __name__ == "__main__":
    main()
