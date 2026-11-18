#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step15_warp_all_maps_to_T1w.py
==============================
ONE uniform antsApplyTransforms pass for EVERY map of EVERY subject into
<sub>/DKI_in_T1w_space/ -- scalars (12) + tensor components (6) + V1 --
replacing the mixed-history contents revealed by the v8 --check audit
(block A: all 12 scalars missing; block B: MK_mf/AK_mf/RK_mf missing;
everything else of assorted vintage). After this run, every *_toT1w_Warped
file in the tree shares one script, one transform set, one interpolation
policy (Linear; NearestNeighbor never used here), one timestamp window, and
one provenance CSV.

NATIVE SOURCES (from the sub-106 audit; fwe/ and fwe_dki/ are DELIBERATELY
never read -- fwe_dki/MD.nii.gz is the free-water-corrected map hiding under
the standard name)
------------------------------------------------------------------
  DKISpace/dki/dkifit_MD.nii.gz                      -> MD
  DKISpace/dki/dkifit_AD.nii.gz                      -> AD    (same RD, FA)
  DKISpace/dki/MK_*medfiltering_True.nii.gz          -> MK_mf (same AK, RK)
  DKISpace/lowb/lowb_ADC.nii.gz                      -> lowb_ADC (etc x5)
  DKISpace/tensor/Dxx.nii.gz ... Dyz                 -> Dxx..Dyz  (step14_v2)
  DKISpace/tensor/V1.nii.gz                          -> V1 (vector-split warp)

SAFETY / SELF-CHECKS
--------------------
* dry-run by default; --apply to write.
* CALIBRATION GATE: before warping anything, the script warps one tensor
  component for one subject that already has a step15-v1 output on disk and
  requires r >= 0.999 against it. If my transform composition (order,
  direction, reference) differs from the audited v1 machinery, this fails
  and NOTHING is written. Override only with --skip-calibration.
* REPLACEMENT CHECK: any existing *_toT1w_Warped map that gets overwritten
  is first correlated against its replacement; r is logged. For block B's
  9 scalars (natives unchanged) r should be ~1.0 -- a low value there is a
  red flag and is printed loudly.
* mtimes of every native input are recorded in the provenance CSV
  (native_mtime); the console prints the oldest native per map type so a
  stale-native problem is visible immediately.
* V1 is warped with the pixdim4-safe vector-split method (3 scalar warps,
  restack, pixdim4=1), as in step15 v1.

USAGE
-----
  python step15_warp_all_maps_to_T1w.py                    # dry run + calibration
  python step15_warp_all_maps_to_T1w.py --apply
  python step15_warp_all_maps_to_T1w.py --subjects sub-115 --apply
"""

import os
import glob
import argparse
import datetime
import subprocess
import tempfile

import numpy as np
import nibabel as nib
import pandas as pd
from pipeline_config import cfg

# ==========================================================================
# CONFIG
# ==========================================================================
DKI_ROOT = cfg("data_root")
T1_ROOT = cfg("pvs_root")
STUDY_XLSX = cfg("study_xlsx")

DATASINK = "pipeline_output/DATASINK/DKISpace"
OUT_DIRNAME = "DKI_in_T1w_space"
MAP_SUFFIX = "_toT1w_Warped.nii.gz"
INTERP = "Linear"

# warped_name -> (subdir under DATASINK, native filename or glob)
SCALAR_MAPS = {
    "MD":       ("dki", "dkifit_MD.nii.gz"),
    "AD":       ("dki", "dkifit_AD.nii.gz"),
    "RD":       ("dki", "dkifit_RD.nii.gz"),
    "FA":       ("dki", "dkifit_FA.nii.gz"),
    # kurtosis: two DATASINK naming conventions exist (block A long-name,
    # block B 'dkifit_*_medfilt'); both are the medfiltering_True product --
    # verified in block B's working tree (rename_mk/MK_*medfiltering_True).
    "MK_mf":    ("dki", ["MK_*medfiltering_True.nii.gz", "dkifit_MK_medfilt.nii.gz"]),
    "AK_mf":    ("dki", ["AK_*medfiltering_True.nii.gz", "dkifit_AK_medfilt.nii.gz"]),
    "RK_mf":    ("dki", ["RK_*medfiltering_True.nii.gz", "dkifit_RK_medfilt.nii.gz"]),
    "lowb_ADC": ("lowb", "lowb_ADC.nii.gz"),
    "lowb_MD":  ("lowb", "lowb_MD.nii.gz"),
    "lowb_AD":  ("lowb", "lowb_AD.nii.gz"),
    "lowb_RD":  ("lowb", "lowb_RD.nii.gz"),
    "lowb_FA":  ("lowb", "lowb_FA.nii.gz"),
    "Dxx":      ("tensor", "Dxx.nii.gz"),
    "Dyy":      ("tensor", "Dyy.nii.gz"),
    "Dzz":      ("tensor", "Dzz.nii.gz"),
    "Dxy":      ("tensor", "Dxy.nii.gz"),
    "Dxz":      ("tensor", "Dxz.nii.gz"),
    "Dyz":      ("tensor", "Dyz.nii.gz"),
}
VECTOR_MAPS = {"V1": ("tensor", "V1.nii.gz")}

CALIB_MAP = "Dxx"          # map used for the calibration gate
CALIB_R_MIN = 0.999


# ==========================================================================
# DISCOVERY
# ==========================================================================
def included_subjects(xlsx=STUDY_XLSX):
    x = pd.ExcelFile(xlsx)
    subs = []
    for sheet in ("HC", "PPCS", "Boxers"):
        t = x.parse(sheet)
        use = t["Will we used it?"].astype(str).str.strip().str.lower()
        subs += [f"sub-{int(n)}" for n in t.loc[use == "yes", "Subj"]]
    return sorted(set(subs))


def find_native(sub, subdir, pattern):
    """pattern may be a string or a list of alternative patterns; exactly one
    file must match across ALL alternatives (ambiguity is refused)."""
    pats = [pattern] if isinstance(pattern, str) else list(pattern)
    hits = []
    for pat in pats:
        hits += glob.glob(os.path.join(DKI_ROOT, sub, DATASINK, subdir, pat))
    hits = sorted(set(hits))
    if len(hits) == 1:
        return hits[0], None
    if not hits:
        return None, f"missing {subdir}/{'|'.join(pats)}"
    return None, (f"{len(hits)} matches for {subdir}/{'|'.join(pats)} "
                  f"(ambiguous): {[os.path.basename(h) for h in hits]}")


def find_t1_and_transforms(sub):
    t1 = sorted(glob.glob(os.path.join(
        T1_ROOT, sub, "**", f"{sub}_ses-01_T1wbrain_restore.nii.gz"),
        recursive=True))
    aff = sorted(glob.glob(os.path.join(
        T1_ROOT, sub, "**", f"{sub}_ses-01_T1wbrain_restore_from_b0_0GenericAffine.mat"),
        recursive=True))
    warp = sorted(glob.glob(os.path.join(
        T1_ROOT, sub, "**", f"{sub}_ses-01_T1wbrain_restore_from_b0_1Warp.nii.gz"),
        recursive=True))
    return (t1[0] if t1 else None,
            aff[0] if aff else None,
            warp[0] if warp else None)


# ==========================================================================
# WARPING
# ==========================================================================
def ants_warp(native, ref, aff, warp, out, interp=INTERP):
    cmd = ["antsApplyTransforms", "-d", "3", "-i", native, "-r", ref,
           "-o", out, "-n", interp, "-t", warp, "-t", aff]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(out):
        raise RuntimeError(f"antsApplyTransforms failed: "
                           f"{(r.stderr or r.stdout).strip()[:200]}")


def warp_vector(native, ref, aff, warp, out, workdir):
    """pixdim4-safe: split 3 components, warp each as 3D scalar, restack."""
    img = nib.load(native)
    data = np.asarray(img.dataobj, np.float32)
    if data.ndim != 4 or data.shape[-1] != 3:
        raise RuntimeError(f"{native}: expected 4D 3-component, got {data.shape}")
    parts = []
    for k in range(3):
        cp = os.path.join(workdir, f"vc{k}.nii.gz")
        wp = os.path.join(workdir, f"vc{k}_w.nii.gz")
        nib.save(nib.Nifti1Image(data[..., k], img.affine), cp)
        ants_warp(cp, ref, aff, warp, wp)
        parts.append(np.asarray(nib.load(wp).dataobj, np.float32))
    ref_img = nib.load(ref)
    stacked = np.stack(parts, axis=-1)
    out_img = nib.Nifti1Image(stacked, ref_img.affine)
    out_img.header["pixdim"][4] = 1.0
    nib.save(out_img, out)


def corr_with_existing(new_path, old_arr):
    a = np.asarray(nib.load(new_path).dataobj, np.float32).ravel()
    b = old_arr.ravel()
    m = np.isfinite(a) & np.isfinite(b) & ((a != 0) | (b != 0))
    if m.sum() < 1000:
        return np.nan
    return float(np.corrcoef(a[m], b[m])[0, 1])


# ==========================================================================
# CALIBRATION GATE
# ==========================================================================
def calibrate(subs, workdir):
    """Warp CALIB_MAP for the first subject with an existing v1 output and
    require r >= CALIB_R_MIN against it. Returns (subject, r) or raises."""
    for sub in subs:
        old_p = os.path.join(T1_ROOT, sub, OUT_DIRNAME, f"{CALIB_MAP}{MAP_SUFFIX}")
        if not os.path.exists(old_p):
            continue
        native, err = find_native(sub, *SCALAR_MAPS[CALIB_MAP])
        t1, aff, warp = find_t1_and_transforms(sub)
        if err or not all((t1, aff, warp)):
            continue
        tmp = os.path.join(workdir, "calib.nii.gz")
        ants_warp(native, t1, aff, warp, tmp)
        old = np.asarray(nib.load(old_p).dataobj, np.float32)
        r = corr_with_existing(tmp, old)
        return sub, r
    raise RuntimeError("no subject with an existing v1 output found to "
                       "calibrate against (use --skip-calibration only if "
                       "you are certain of the transform composition)")


# ==========================================================================
# MAIN
# ==========================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subjects", nargs="*", default=None)
    ap.add_argument("--apply", action="store_true",
                    help="actually write; default is dry run")
    ap.add_argument("--skip-calibration", action="store_true")
    ap.add_argument("--out-csv", default="./step15_v2_warp_provenance.csv")
    a = ap.parse_args(argv)

    subs = a.subjects or included_subjects()
    print(f"Mode: {'APPLY' if a.apply else 'DRY RUN (nothing will be written)'}")
    print(f"{len(subs)} subjects | {len(SCALAR_MAPS)} scalar + "
          f"{len(VECTOR_MAPS)} vector maps -> <sub>/{OUT_DIRNAME}/  "
          f"[interp {INTERP}]")

    rows = []
    with tempfile.TemporaryDirectory(prefix="warp15v2_") as wd:
        if a.apply and not a.skip_calibration:
            sub_c, r = calibrate(subs, wd)
            verdict = "PASS" if (np.isfinite(r) and r >= CALIB_R_MIN) else "FAIL"
            print(f"\nCALIBRATION [{sub_c}, {CALIB_MAP}]: r = {r:.6f} vs "
                  f"existing v1 output -> {verdict}")
            if verdict == "FAIL":
                raise SystemExit(
                    "calibration failed: this script's transform composition "
                    "does not reproduce the audited v1 outputs. NOTHING was "
                    "written. Investigate before proceeding.")

        for sub in subs:
            t1, aff, warp = find_t1_and_transforms(sub)
            if not all((t1, aff, warp)):
                miss = [n for n, p in (("T1w", t1), ("affine", aff),
                                       ("warp", warp)) if p is None]
                print(f"  {sub}: SKIP -- missing {miss}")
                rows.append({"subject": sub, "map": "", "status":
                             f"SKIP missing {miss}"})
                continue
            out_dir = os.path.join(T1_ROOT, sub, OUT_DIRNAME)
            os.makedirs(out_dir, exist_ok=True) if a.apply else None

            todo = list(SCALAR_MAPS.items()) + list(VECTOR_MAPS.items())
            n_ok = n_prob = 0
            for name, (subdir, pattern) in todo:
                native, err = find_native(sub, subdir, pattern)
                if err:
                    rows.append({"subject": sub, "map": name, "status": err})
                    n_prob += 1
                    continue
                out_p = os.path.join(out_dir, f"{name}{MAP_SUFFIX}")
                rec = {"subject": sub, "map": name, "native": native,
                       "native_mtime": datetime.datetime.fromtimestamp(
                           os.path.getmtime(native)).isoformat(timespec="seconds"),
                       "output": out_p}
                if not a.apply:
                    rec["status"] = ("would REPLACE" if os.path.exists(out_p)
                                     else "would create")
                    rows.append(rec)
                    n_ok += 1
                    continue
                old_arr = (np.asarray(nib.load(out_p).dataobj, np.float32)
                           if os.path.exists(out_p) and name in SCALAR_MAPS
                           else None)
                try:
                    if name in VECTOR_MAPS:
                        warp_vector(native, t1, aff, warp, out_p, wd)
                    else:
                        ants_warp(native, t1, aff, warp, out_p)
                except RuntimeError as e:
                    rec["status"] = f"FAILED: {e}"
                    rows.append(rec)
                    n_prob += 1
                    continue
                if old_arr is not None and old_arr.shape == nib.load(out_p).shape:
                    r = corr_with_existing(out_p, old_arr)
                    rec["r_vs_replaced"] = round(r, 6) if np.isfinite(r) else np.nan
                    if np.isfinite(r) and r < 0.99:
                        print(f"    [!] {sub} {name}: replaced map differs "
                              f"from previous (r={r:.4f}) -- expected for "
                              f"regenerated/refit maps, NOT for unchanged "
                              f"natives")
                rec["status"] = "written"
                rows.append(rec)
                n_ok += 1
            print(f"  {sub}: {n_ok} maps {'planned' if not a.apply else 'written'}"
                  + (f", {n_prob} problems" if n_prob else ""))

    df = pd.DataFrame(rows)
    df.to_csv(a.out_csv, index=False)
    prob = df[df.status.str.startswith(("missing", "FAILED", "SKIP"), na=False)]
    print(f"\n{len(df)} rows -> {a.out_csv} | {len(prob)} problems")
    if "native_mtime" in df and df["native_mtime"].notna().any():
        oldest = (df.dropna(subset=["native_mtime"])
                    .groupby("map")["native_mtime"].min().sort_values())
        print("\noldest native per map (stale-input check):")
        print(oldest.to_string())
    if not a.apply:
        print("\nDry run complete. Re-run with --apply to write "
              "(calibration gate runs first).")


if __name__ == "__main__":
    main()
