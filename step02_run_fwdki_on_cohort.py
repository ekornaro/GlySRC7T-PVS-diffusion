#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step02_run_fwdki_on_cohort.py
=============================
Batch-run the free-water-corrected DKI pipeline (fwdki_pvs.run_from_files) over
every sub-* folder in a cohort, writing outputs into each subject's DKISpace.

Layout assumed (per subject):
  <root>/sub-XXX/pipeline_output/DATASINK/DKISpace/dwi_proc/DKI_denoised_gibbs_eddy.nii.gz
  <root>/sub-XXX/pipeline_output/DATASINK/DKISpace/dwi_proc/<gradients>.bval / .bvec
  <root>/sub-XXX/pipeline_output/DATASINK/DKISpace/masks/ANTS_T1_brain_mask.nii.gz
Outputs written to:
  <root>/sub-XXX/pipeline_output/DATASINK/DKISpace/fwe_dki/

DEFENSIVE GRADIENT SELECTION
----------------------------
The bval/bvec basename is NOT assumed. Among all .bval/.bvec pairs in the dwi_proc
folder, the runner picks the pair whose entry count EQUALS the image's 4th dimension
(number of volumes). If the stem-matched name (DKI_denoised_gibbs_eddy.*) exists and
matches, it is preferred. If nothing matches the volume count, the subject is SKIPPED
and logged -- never fit with a mismatched gradient table (the failure mode that
silently produces all-zero maps).

ERROR ISOLATION
---------------
One subject's failure never stops the cohort. Every subject's status (OK / SKIP /
ERROR) and detail is printed live and written to a summary CSV at the end.

USAGE
-----
  python run_cohort.py                         # all sub-* under the default root
  python run_cohort.py --subjects sub-105 sub-106
  python run_cohort.py --skip-existing         # skip subjects that already have f.nii.gz
  python run_cohort.py --root /path/to/cohort
  python run_cohort.py --pvs-mask-rel pipeline_output/DATASINK/DKISpace/masks/PVS_mask.nii.gz
                                               # also write region summaries if a PVS mask exists
"""

import os
import glob
import csv
import argparse
import numpy as np
import nibabel as nib

import fwdki as F
from pipeline_config import cfg

# ---- cohort layout (edit here or override on the CLI) ----------------------- #
ROOT_DEFAULT = cfg("data_root")
DWI_REL  = "pipeline_output/DATASINK/DKISpace/dwi_proc"
DWI_NAME = "DKI_denoised_gibbs_eddy.nii.gz"
MASK_REL = "pipeline_output/DATASINK/DKISpace/masks/ANTS_T1_brain_mask.nii.gz"
OUT_REL  = "pipeline_output/DATASINK/DKISpace/fwe_dki"


def _nvols(nifti_path):
    img = nib.load(nifti_path)
    return int(img.shape[3]) if img.ndim == 4 else 1


def find_gradients(dwi_path, n_vols):
    """Return (bval, bvec, n) for the pair whose bval entry count == n_vols.

    Preference order: (1) stem-matched name (same as the NIfTI), (2) any other
    .bval/.bvec pair in the same folder. Only a pair whose count matches the image
    volume count is accepted -- this is what prevents silent gradient mismatch.
    """
    folder = os.path.dirname(dwi_path)
    stem = dwi_path[:-7] if dwi_path.endswith(".nii.gz") else os.path.splitext(dwi_path)[0]

    candidates = []
    if os.path.exists(stem + ".bval") and os.path.exists(stem + ".bvec"):
        candidates.append((stem + ".bval", stem + ".bvec"))
    for bv in sorted(glob.glob(os.path.join(folder, "*.bval"))):
        bvec = bv[:-5] + ".bvec"
        if os.path.exists(bvec) and (bv, bvec) not in candidates:
            candidates.append((bv, bvec))

    matches = []
    for bval, bvec in candidates:
        try:
            n = len(np.atleast_1d(np.loadtxt(bval)).ravel())
        except Exception:
            continue
        if n == n_vols:
            matches.append((bval, bvec, n))
    if matches:
        return matches[0]            # preference order preserved
    return None, None, None


def process_subject(sub_dir, skip_existing=False, pvs_mask_rel=None):
    sub = os.path.basename(sub_dir.rstrip("/"))
    dwi  = os.path.join(sub_dir, DWI_REL, DWI_NAME)
    mask = os.path.join(sub_dir, MASK_REL)
    out  = os.path.join(sub_dir, OUT_REL)

    if not os.path.exists(dwi):
        return sub, "SKIP", "no DWI file"
    if skip_existing and os.path.exists(os.path.join(out, "f.nii.gz")):
        return sub, "SKIP", "already done"

    try:
        n_vols = _nvols(dwi)
    except Exception as e:
        return sub, "ERROR", f"cannot read DWI: {repr(e)[:100]}"

    bval, bvec, n = find_gradients(dwi, n_vols)
    if bval is None:
        return sub, "SKIP", f"no bval/bvec matching {n_vols} volumes"

    bmask = mask if os.path.exists(mask) else None
    pvs = None
    if pvs_mask_rel:
        cand = os.path.join(sub_dir, pvs_mask_rel)
        pvs = cand if os.path.exists(cand) else None

    try:
        F.run_from_files(dwi, bval, bvec, out, brain_mask=bmask, pvs_mask=pvs)
        detail = f"{os.path.basename(bval)} ({n} vols)"
        if bmask is None:
            detail += " [no brain mask]"
        if pvs is not None:
            detail += " [+PVS summaries]"
        return sub, "OK", detail
    except Exception as e:
        return sub, "ERROR", repr(e)[:150]


def main():
    ap = argparse.ArgumentParser(description="Batch fwdki_pvs over a cohort of sub-* folders.")
    ap.add_argument("--root", default=ROOT_DEFAULT, help="cohort root containing sub-* folders")
    ap.add_argument("--subjects", nargs="*", default=None,
                    help="specific subject IDs (e.g. sub-105 sub-106); default = all sub-*")
    ap.add_argument("--skip-existing", action="store_true",
                    help="skip subjects that already have fwe_dki/f.nii.gz")
    ap.add_argument("--pvs-mask-rel", default=None,
                    help="relative path (from sub dir) to a PVS mask; if present, region "
                         "summaries are written alongside the maps")
    a = ap.parse_args()

    if a.subjects:
        sub_dirs = [os.path.join(a.root, s) for s in a.subjects]
    else:
        sub_dirs = sorted(glob.glob(os.path.join(a.root, "sub-*")))
    if not sub_dirs:
        print(f"No subjects found under {a.root}")
        return

    print(f"Cohort root: {a.root}")
    print(f"Subjects to process: {len(sub_dirs)}")
    print("-" * 70)

    results = []
    for sd in sub_dirs:
        sub, status, msg = process_subject(sd, skip_existing=a.skip_existing,
                                           pvs_mask_rel=a.pvs_mask_rel)
        print(f"[{status:5}] {sub}: {msg}", flush=True)
        results.append((sub, status, msg))

    # summary
    n_ok = sum(1 for _, s, _ in results if s == "OK")
    n_skip = sum(1 for _, s, _ in results if s == "SKIP")
    n_err = sum(1 for _, s, _ in results if s == "ERROR")
    print("-" * 70)
    print(f"Done. OK={n_ok}  SKIP={n_skip}  ERROR={n_err}")

    summary = os.path.join(a.root, "cohort_fwe_summary.csv")
    try:
        with open(summary, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["subject", "status", "detail"])
            w.writerows(results)
        print(f"Summary written: {summary}")
    except Exception as e:
        print(f"(could not write summary CSV: {e})")


if __name__ == "__main__":
    main()
