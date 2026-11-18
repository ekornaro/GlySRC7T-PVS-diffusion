#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step20_fwe_reliability_gate.py
==============================
The load-bearing check that fwdki.py's docstring calls mandatory and that was
never run: quantify how precisely each free-water-elimination (FWE) metric is
ESTIMATED on this acquisition, in the voxels that matter (PVS) and in matched
white matter, using the exact estimator that produced the maps (fwdki.py is
imported, not re-implemented).

Nothing here uses group labels. This is a precision gate, decided BEFORE any
FWE metric enters the tier analysis, with the pass/fail thresholds stated up
front (below) so inclusion cannot be tuned to the outcome.

THREE TESTS, ONE VERDICT PER METRIC
-----------------------------------
(1) RESIDUAL BOOTSTRAP (the fwdki control). Fit standard DKI on the ROI
    voxels, resample residuals over volumes with replacement, re-run the full
    fwdki pipeline (f -> deconfound -> corrected DKI) on each replicate.
    Reported per metric:
      * vox_sd      median per-voxel bootstrap SD  (precision of a voxel)
      * roimed_sd   SD across replicates of the ROI MEDIAN (precision of the
                    subject-level summary that actually enters the mixed model)
      * amp         noise-amplification factor = vox_sd(corrected metric) /
                    vox_sd(uncorrected counterpart), e.g. MKt vs MK. The
                    deconfound divides by (1 - f), so amp ~ 1/(1 - f) is the
                    expected order; a much larger value means the corrected
                    fit is unstable, not merely rescaled.
(2) DIRECTION SUBSAMPLING. Drop --drop-dirs directions from the b=2000 shell
    (repeat --n-draws times) and refit. SD of the ROI median across draws
    tests the 'too few directions' concern as it applies to THIS two-stage
    estimator (Collier's joint fit was rejected for that reason; this checks
    the method that was actually used).
(3) PHYSICAL SANITY of f: median in lateral ventricles (expect ~1), deep WM
    (expect ~0.05-0.15), and PVS ROI (expect above WM). Failure here means
    the model is not describing the data, whatever the precision says.

VERDICT RULE (per metric, stated a priori)
------------------------------------------
  reliability ratio  R = roimed_sd(PVS) / between_subject_sd(PVS tier)
  where the between-subject SD is taken from the v8 spreadsheet for the
  UNCORRECTED counterpart (MDt -> MD etc.; f -> lowb_ADC as the closest fluid-
  sensitive scalar). Estimation noise then contributes R^2 of the variance
  the model sees:
      PASS  R < 0.20   (estimation noise < 4% of between-subject variance)
      FLAG  0.20-0.50  (usable; report with the ratio; expect attenuation)
      FAIL  R > 0.50   (noise-dominated at the subject level; exclude)
  plus: amp > --amp-max (default 4) -> FAIL regardless; f sanity fail -> all
  FAIL. Without --metrics-xlsx the ratio cannot be formed and only the
  absolute SDs, amp, and sanity are reported.

ROI CONSTRUCTION (native DKI space, 2 mm)
-----------------------------------------
PVS masks live in T1w space (0.8 mm) and PVS are 1-2 voxels wide, so nearest-
neighbour resampling to 2 mm would erase most of them. The mask is warped as
a FLOAT (Linear) and thresholded at --pvs-frac (default 0.25): a native voxel
counts as PVS if >= 25% of it is PVS -- the partial-volume reality the tier
extraction also samples. WM = SynthSeg cerebral WM (2/41) warped NN, minus a
2-voxel dilation of the PVS ROI. Ventricles = SynthSeg 4/43. All warps use
the inverse of the b0->T1w registration (T1w -> native), never touching DWI.

USAGE
-----
  python step20_fwe_reliability_gate.py --subjects sub-115 sub-160 \\
      --metrics-xlsx pvs_tier_metrics_v8.xlsx
  python step20_fwe_reliability_gate.py --subjects sub-105 sub-125 sub-141 \\
      sub-155 sub-160 sub-170 --n-boot 100 --out ./fwe_gate
"""

import os
import glob
import json
import argparse
import subprocess
import tempfile
import importlib.util

import numpy as np
import nibabel as nib
import pandas as pd
from scipy import ndimage
from pipeline_config import cfg

# ==========================================================================
# CONFIG
# ==========================================================================
DKI_ROOT = cfg("data_root")
T1_ROOT = cfg("pvs_root")
FWDKI_PATH = cfg("fwdki_py")

DWI_SUBDIR = "pipeline_output/DATASINK/DKISpace/dwi_proc"
DWI_NAME = "DKI_denoised_gibbs_eddy.nii.gz"
BVEC_NAME = "DKI_corrected_convfix.bvec"       # scalar metrics are convention-
BVAL_NAME = "DKI_corrected_convfix.bval"       # invariant; use the harmonized set
FA_NATIVE = "pipeline_output/DATASINK/DKISpace/dki/dkifit_FA.nii.gz"

WM_LABELS = [2, 41]
VENT_LABELS = [4, 43]

# corrected -> uncorrected counterpart (for amp and for the spreadsheet SD)
COUNTERPART = {"FAt": "FA", "MDt": "MD", "ADt": "AD", "RDt": "RD",
               "MKt": "MK_mf", "AKt": "AK_mf", "RKt": "RK_mf", "f": "lowb_ADC"}
UNCORR_KEY = {"FAt": "FA", "MDt": "MD", "ADt": "AD", "RDt": "RD",
              "MKt": "MK", "AKt": "AK", "RKt": "RK"}          # keys in fwdki output

# a-priori thresholds
R_PASS, R_FAIL = 0.20, 0.50
AMP_MAX = 4.0
F_VENT_MIN, F_WM_MAX = 0.7, 0.25


def load_fwdki(path):
    spec = importlib.util.spec_from_file_location("fwdki", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ==========================================================================
# ROI CONSTRUCTION
# ==========================================================================
def find_transforms(sub):
    hits = sorted(glob.glob(os.path.join(
        T1_ROOT, sub, "**", f"{sub}_ses-01_T1wbrain_restore_from_b0_*"), recursive=True))
    aff = next((h for h in hits if h.endswith("0GenericAffine.mat")), None)
    inv = next((h for h in hits if "1InverseWarp.nii.gz" in os.path.basename(h)), None)
    return aff, inv


def warp_to_native(src, ref, aff, inv, out, interp):
    cmd = ["antsApplyTransforms", "-d", "3", "-i", src, "-r", ref, "-o", out,
           "-n", interp]
    # inverse of the b0 -> T1w registration: the inverted affine is listed
    # first (applied last), as in step07. Runs before the fix listed the
    # inverse warp first.
    cmd += ["-t", f"[{aff},1]"]
    if inv:
        cmd += ["-t", inv]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(out):
        raise RuntimeError(f"antsApplyTransforms failed: {(r.stderr or r.stdout)[:200]}")
    return np.asarray(nib.load(out).dataobj, np.float32)


def build_rois(sub, ref_b0_path, workdir, pvs_frac):
    aff, inv = find_transforms(sub)
    if not aff:
        raise RuntimeError("no b0<->T1w affine found")
    pvs_p = glob.glob(os.path.join(T1_ROOT, sub, "**", "*_ses-01_PVS_postprocessed.nii*"),
                      recursive=True)
    seg_p = glob.glob(os.path.join(T1_ROOT, sub, "**", "*_ses-01_synthseg.nii.gz"),
                      recursive=True)
    if len(pvs_p) != 1 or len(seg_p) != 1:
        raise RuntimeError(f"PVS/SynthSeg not unique ({len(pvs_p)}/{len(seg_p)})")

    # PVS as a partial-volume fraction in native space
    pvs_img = nib.load(pvs_p[0])
    pvs_f = os.path.join(workdir, "pvs_float.nii.gz")
    nib.save(nib.Nifti1Image((np.asarray(pvs_img.dataobj) > 0).astype(np.float32),
                             pvs_img.affine), pvs_f)
    pvs_nat = warp_to_native(pvs_f, ref_b0_path, aff, inv,
                             os.path.join(workdir, "pvs_nat.nii.gz"), "Linear")
    seg_nat = warp_to_native(seg_p[0], ref_b0_path, aff, inv,
                             os.path.join(workdir, "seg_nat.nii.gz"), "NearestNeighbor")

    pvs = pvs_nat >= pvs_frac
    wm = np.isin(seg_nat, WM_LABELS) & ~ndimage.binary_dilation(pvs, iterations=2)
    vent = np.isin(seg_nat, VENT_LABELS)
    return pvs, wm, vent, {"pvs_frac_median": float(np.median(pvs_nat[pvs]))
                           if pvs.any() else np.nan}


def subsample(mask, n_max, rng):
    idx = np.argwhere(mask)
    if idx.shape[0] > n_max:
        idx = idx[rng.choice(idx.shape[0], n_max, replace=False)]
    return idx


# ==========================================================================
# TESTS
# ==========================================================================
def roi_signal(data, idx):
    return np.stack([data[tuple(v)] for v in idx], axis=0)          # (n_vox, N)


def run_on_rows(fw, sig, bvals, bvecs, metrics):
    """Run the full fwdki pipeline on (n_vox, N) signal rows."""
    out = fw.run_subject(sig.reshape(sig.shape[0], 1, 1, -1), bvals, bvecs,
                         brain_mask=None)
    return {m: np.asarray(out[m]).reshape(-1) for m in metrics if m in out}


def bootstrap(fw, sig, bvals, bvecs, metrics, n_boot, rng):
    """Residual bootstrap mirroring fwdki.bootstrap_estimation_sd, but capturing
    corrected AND uncorrected metrics, and the ROI-median statistic."""
    from dipy.reconst import dki
    gtab = fw.build_gtab(bvals, bvecs)
    model = dki.DiffusionKurtosisModel(gtab, fit_method="WLS")
    S0 = fw.estimate_S0(sig, bvals, bvecs)
    fit = model.fit(sig)
    pred = model.predict(fit.model_params, S0=S0)
    resid = sig - pred
    n_vox, N = sig.shape
    samples = {m: np.full((n_boot, n_vox), np.nan) for m in metrics}
    for b in range(n_boot):
        boot = np.clip(pred + resid[:, rng.integers(0, N, size=N)], 0.0, None)
        out = run_on_rows(fw, boot, bvals, bvecs, metrics)
        for m in metrics:
            if m in out:
                samples[m][b] = out[m]
    res = {}
    for m in metrics:
        s = samples[m]
        res[m] = {"vox_sd": float(np.nanmedian(np.nanstd(s, axis=0))),
                  "roimed_sd": float(np.nanstd(np.nanmedian(s, axis=1))),
                  "roimed": float(np.nanmedian(np.nanmedian(s, axis=1)))}
    return res


def direction_subsample(fw, sig, bvals, bvecs, metrics, drop, n_draws, rng):
    hi = np.flatnonzero(bvals >= 1900)
    if len(hi) <= drop + 6:
        return {}
    meds = {m: [] for m in metrics}
    for _ in range(n_draws):
        dropped = rng.choice(hi, drop, replace=False)
        keep = np.ones(len(bvals), bool); keep[dropped] = False
        out = run_on_rows(fw, sig[:, keep], bvals[keep], bvecs[keep], metrics)
        for m in metrics:
            if m in out:
                meds[m].append(np.nanmedian(out[m]))
    return {m: {"dirsub_roimed_sd": float(np.nanstd(v)) if len(v) > 1 else np.nan}
            for m, v in meds.items()}


def between_subject_sd(xlsx, sheet="Case_A_Near_PVS"):
    """Between-subject SD of subject-level PVS-tier medians, per uncorrected metric."""
    a = pd.read_excel(xlsx, sheet_name=sheet)
    out = {}
    for m in set(COUNTERPART.values()):
        col = f"{m}_Original"
        if col in a:
            per_sub = a.groupby("Subject")[col].median()
            out[m] = float(per_sub.std(ddof=1))
    return out


# ==========================================================================
# MAIN
# ==========================================================================
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subjects", nargs="+", required=True)
    ap.add_argument("--fwdki", default=FWDKI_PATH)
    ap.add_argument("--metrics-xlsx", default="",
                    help="v8 spreadsheet; enables the reliability ratio R")
    ap.add_argument("--n-boot", type=int, default=100)
    ap.add_argument("--n-vox", type=int, default=300, help="voxels per ROI")
    ap.add_argument("--pvs-frac", type=float, default=0.25)
    ap.add_argument("--drop-dirs", type=int, default=5)
    ap.add_argument("--n-draws", type=int, default=10)
    ap.add_argument("--amp-max", type=float, default=AMP_MAX)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="./fwe_gate")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    fw = load_fwdki(a.fwdki)
    metrics = fw.PRIMARY_METRICS + fw.UNCORR_METRICS
    bs_sd = between_subject_sd(a.metrics_xlsx) if a.metrics_xlsx else {}

    print(f"FWE reliability gate | fwdki from {a.fwdki}")
    print(f"thresholds (a priori): R PASS<{R_PASS} FAIL>{R_FAIL} | amp<={a.amp_max} "
          f"| f(vent)>={F_VENT_MIN} f(WM)<={F_WM_MAX}")
    print(f"n_boot={a.n_boot}  n_vox/ROI={a.n_vox}  dir-subsample drop {a.drop_dirs} "
          f"x {a.n_draws}\n")

    rows = []
    with tempfile.TemporaryDirectory(prefix="fwegate_") as wd:
        for sub in a.subjects:
            sd_ = os.path.join(DKI_ROOT, sub, DWI_SUBDIR)
            dwi = os.path.join(sd_, DWI_NAME)
            bvec_p, bval_p = os.path.join(sd_, BVEC_NAME), os.path.join(sd_, BVAL_NAME)
            if not all(os.path.exists(p) for p in (dwi, bvec_p, bval_p)):
                print(f"[{sub}] SKIP -- DWI/gradients missing"); continue
            img = nib.load(dwi)
            data = np.asarray(img.dataobj, np.float32)
            bvals = np.loadtxt(bval_p).ravel()
            bvecs = np.loadtxt(bvec_p)
            if bvecs.shape[0] == 3:
                bvecs = bvecs.T
            if not (bvals <= 50).any():
                print(f"[{sub}] SKIP -- no b=0 volume: this gate is for the "
                      f"harmonized (b0-present) data only"); continue
            b0 = os.path.join(wd, f"{sub}_b0.nii.gz")
            nib.save(nib.Nifti1Image(data[..., int(np.argmin(bvals))], img.affine), b0)
            try:
                pvs, wm, vent, info = build_rois(sub, b0, wd, a.pvs_frac)
            except RuntimeError as e:
                print(f"[{sub}] SKIP -- {e}"); continue
            n_pvs, n_wm = int(pvs.sum()), int(wm.sum())
            print(f"[{sub}] ROIs: PVS {n_pvs} vox (median PVS fraction "
                  f"{info['pvs_frac_median']:.2f}), WM {n_wm}, ventricles {int(vent.sum())}")
            if n_pvs < 30:
                print(f"    [!] too few PVS voxels at pvs_frac={a.pvs_frac}; lower it"); continue

            # ---- (3) physical sanity of f on the full-ROI fit ----
            san = {}
            for name, mask in (("vent", vent), ("wm", wm), ("pvs", pvs)):
                idx = subsample(mask, a.n_vox, rng)
                out = run_on_rows(fw, roi_signal(data, idx), bvals, bvecs, ["f"])
                san[name] = float(np.nanmedian(out["f"]))
            f_ok = san["vent"] >= F_VENT_MIN and san["wm"] <= F_WM_MAX and san["pvs"] > san["wm"]
            print(f"    f sanity: ventricles {san['vent']:.2f} | WM {san['wm']:.3f} | "
                  f"PVS {san['pvs']:.3f}  -> {'OK' if f_ok else 'FAIL'}")

            # ---- (1)+(2) on PVS and WM ROIs ----
            for roi_name, mask in (("PVS", pvs), ("WM", wm)):
                idx = subsample(mask, a.n_vox, rng)
                sig = roi_signal(data, idx)
                print(f"    bootstrap on {roi_name} ({len(idx)} vox, {a.n_boot} boots) ...",
                      flush=True)
                bs = bootstrap(fw, sig, bvals, bvecs, metrics, a.n_boot, rng)
                ds = direction_subsample(fw, sig, bvals, bvecs, fw.PRIMARY_METRICS,
                                         a.drop_dirs, a.n_draws, rng)
                for m in fw.PRIMARY_METRICS:
                    r = {"Subject": sub, "ROI": roi_name, "metric": m,
                         "roimed": bs[m]["roimed"], "vox_sd": bs[m]["vox_sd"],
                         "roimed_sd": bs[m]["roimed_sd"],
                         "dirsub_roimed_sd": ds.get(m, {}).get("dirsub_roimed_sd", np.nan),
                         "f_sanity_ok": f_ok}
                    uk = UNCORR_KEY.get(m)
                    r["amp"] = (bs[m]["vox_sd"] / bs[uk]["vox_sd"]
                                if uk and bs[uk]["vox_sd"] > 0 else np.nan)
                    cp = COUNTERPART[m]
                    r["between_subject_sd"] = bs_sd.get(cp, np.nan)
                    r["R"] = (bs[m]["roimed_sd"] / r["between_subject_sd"]
                              if np.isfinite(r["between_subject_sd"]) and r["between_subject_sd"] > 0
                              else np.nan)
                    rows.append(r)

    df = pd.DataFrame(rows)
    if df.empty:
        raise SystemExit("no results")
    df.to_csv(os.path.join(a.out, "fwe_reliability_gate.csv"), index=False)

    # ---- verdict: PVS ROI, aggregated across subjects (median) ----
    print("\n" + "=" * 78)
    print("VERDICT per metric (PVS ROI, medians across subjects)")
    print("=" * 78)
    pv = df[df.ROI == "PVS"].groupby("metric").agg(
        roimed=("roimed", "median"), vox_sd=("vox_sd", "median"),
        roimed_sd=("roimed_sd", "median"), dirsub=("dirsub_roimed_sd", "median"),
        amp=("amp", "median"), R=("R", "median"), f_ok=("f_sanity_ok", "all"))
    verdicts = {}
    for m, r in pv.iterrows():
        if not r.f_ok:
            v = "FAIL (f sanity)"
        elif np.isfinite(r.amp) and r.amp > a.amp_max:
            v = f"FAIL (amp {r.amp:.1f}x)"
        elif np.isfinite(r.R):
            v = "PASS" if r.R < R_PASS else ("FLAG" if r.R <= R_FAIL else "FAIL")
        else:
            v = "n/a (no --metrics-xlsx)"
        verdicts[m] = v
        print(f"  {m:5s} median {r.roimed:9.3e}  vox_sd {r.vox_sd:8.2e}  "
              f"roimed_sd {r.roimed_sd:8.2e}  dirsub_sd {r.dirsub:8.2e}  "
              f"amp {r.amp:5.2f}  R {r.R:5.2f}  -> {v}")
    with open(os.path.join(a.out, "fwe_gate_verdict.json"), "w") as fh:
        json.dump({"verdicts": verdicts, "thresholds": {"R_pass": R_PASS,
                   "R_fail": R_FAIL, "amp_max": a.amp_max}, "subjects": a.subjects,
                   "n_boot": a.n_boot}, fh, indent=2)
    print(f"\n-> {a.out}/fwe_reliability_gate.csv, fwe_gate_verdict.json")
    print("\nInclusion rule: PASS metrics enter the tier analysis; FLAG metrics "
          "enter with the ratio reported;\nFAIL metrics stay excluded with the "
          "number as the stated reason. Decide from this table, not from the "
          "group results.")


if __name__ == "__main__":
    main()
