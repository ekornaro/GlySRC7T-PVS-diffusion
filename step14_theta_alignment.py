#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step14_theta_alignment.py
=========================
Per-cluster theta statistics from the CURRENT (v3, header-corrected) maps,
testing ALIGNMENT against the uniform null -- replacing the v1 script, whose
input CSV predates both frame corrections and is no longer generated.

WHAT CHANGED FROM v1
--------------------
* INPUT: per-cluster values are sampled directly from the step16 v3 outputs
  (theta_toT1w_Warped, clusterLinearity_toT1w_Warped + the PVS mask), so the
  numbers can never silently come from a stale intermediate. The maps'
  provenance sidecar is read and echoed; a run on non-header maps warns.
* FRAMING: v1 confirmed uniformity; v2 quantifies alignment. Tests are
  one-sided against the uniform-null references, per subject and per group.
* CORRECTED NULL: under orientation independence P(theta < t) = 1 - cos(t),
  so the uniform MEDIAN is 60.0 deg (not 57.3 -- that is the uniform MEAN,
  1 radian). v1's prose compared medians to 57.3; v2 uses:
      median theta        -> 60.0 deg
      mean   theta        -> 57.30 deg
      frac(theta < 30)    -> 1 - cos(30) = 0.1340
      frac(theta > 60)    -> cos(60)     = 0.5000
      cos^2(theta)        -> E[cos^2] = 1/3 (a clean scalar effect size)

OUTPUTS
-------
  <out>/theta_per_cluster.csv     Subject, PVS_ID, n_vox, theta_med, Linearity
  <out>/theta_per_subject.csv     per-subject summaries + group + u_rel_det
  console: group table, alignment tests vs uniform, group KW,
           det=+1 vs det=-1 consistency check (should be null after v3)

USAGE
-----
  python step14_theta_alignment.py
  python step14_theta_alignment.py --min-linearity 2.0 --out ./theta_alignment
"""

import os
import glob
import json
import argparse

import numpy as np
import pandas as pd
import nibabel as nib
from scipy import ndimage
from scipy.stats import kruskal, mannwhitneyu, wilcoxon
from pipeline_config import cfg

T1_ROOT = cfg("pvs_root")
STUDY_XLSX = cfg("study_xlsx")
QC_CSV = "./pvs_orientation_maps_qc/orientation_maps_qc.csv"

DKI_SUBDIR = "DKI_in_T1w_space"
MAP_SUFFIX = "_toT1w_Warped.nii.gz"
MIN_CLUSTER_SIZE = 3

# uniform-null references (see docstring)
UNIF_MEDIAN = 60.0
UNIF_MEAN = 57.2958
UNIF_FRAC30 = 1.0 - np.cos(np.radians(30))     # 0.1340
UNIF_FRAC60 = np.cos(np.radians(60))           # 0.5000
UNIF_COS2 = 1.0 / 3.0


def group_map(xlsx=STUDY_XLSX):
    x = pd.ExcelFile(xlsx)
    out = {}
    for sheet, label in [("HC", "HC"), ("PPCS", "PPCS"), ("Boxers", "RHI")]:
        t = x.parse(sheet)
        use = t["Will we used it?"].astype(str).str.strip().str.lower()
        for n in t.loc[use == "yes", "Subj"].astype(int):
            out[f"sub-{n}"] = label
    return out


def per_cluster_rows(subj_dir, sub, min_size=MIN_CLUSTER_SIZE):
    md = os.path.join(subj_dir, DKI_SUBDIR)
    th_p = os.path.join(md, f"theta{MAP_SUFFIX}")
    li_p = os.path.join(md, f"clusterLinearity{MAP_SUFFIX}")
    if not (os.path.exists(th_p) and os.path.exists(li_p)):
        raise RuntimeError("theta/clusterLinearity maps missing (run step16 v3)")

    prov_p = os.path.join(md, "pvs_orientation_provenance.json")
    prov = {}
    if os.path.exists(prov_p):
        with open(prov_p) as fh:
            prov = json.load(fh)

    hits = glob.glob(os.path.join(subj_dir, "**",
                                  "*_ses-01_PVS_postprocessed.nii*"), recursive=True)
    if len(hits) != 1:
        raise RuntimeError(f"PVS mask not unique ({len(hits)})")
    pvs = nib.load(hits[0]).get_fdata() > 0
    theta = np.asarray(nib.load(th_p).dataobj, np.float32)
    lin = np.asarray(nib.load(li_p).dataobj, np.float32)
    if theta.shape != pvs.shape:
        raise RuntimeError("map/mask grid mismatch")

    labeled, n_feat = ndimage.label(pvs, structure=np.ones((3, 3, 3)))
    rows = []
    for lab, sl in enumerate(ndimage.find_objects(labeled), start=1):
        if sl is None:
            continue
        m = labeled[sl] == lab
        n = int(m.sum())
        if n < min_size:
            continue
        tvals = theta[sl][m]
        lvals = lin[sl][m]
        tvals = tvals[np.isfinite(tvals)]
        if tvals.size == 0:
            continue
        rows.append({"Subject": sub, "PVS_ID": lab, "n_vox": n,
                     "theta_med": float(np.median(tvals)),
                     "Linearity": float(np.nanmedian(lvals))})
    return rows, prov


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=T1_ROOT)
    ap.add_argument("--study-xlsx", default=STUDY_XLSX,
                    help="study inclusion spreadsheet (group mapping)")
    ap.add_argument("--qc-csv", default=QC_CSV,
                    help="step16 v3 QC csv (for u_rel_det per subject)")
    ap.add_argument("--min-linearity", type=float, default=2.0)
    ap.add_argument("--min-cluster-size", type=int, default=MIN_CLUSTER_SIZE)
    ap.add_argument("--out", default="./theta_alignment")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)

    groups = group_map(a.study_xlsx)
    dirs = sorted(glob.glob(os.path.join(a.root, "sub-*")))

    all_rows, prov_modes = [], set()
    for d in dirs:
        sub = os.path.basename(d.rstrip("/"))
        if sub not in groups:
            continue
        try:
            rows, prov = per_cluster_rows(d, sub, a.min_cluster_size)
        except RuntimeError as e:
            print(f"  {sub}: SKIP -- {e}")
            continue
        all_rows.extend(rows)
        if prov:
            prov_modes.add(prov.get("rotate_u", "?"))

    print(f"map provenance rotate_u modes seen: {sorted(prov_modes) or 'none recorded'}")
    if prov_modes - {"header"}:
        print("  [!] WARNING: some maps were not produced with --rotate-u header. "
              "Re-run step16 v3 before trusting these numbers.")

    d = pd.DataFrame(all_rows)
    d["Group"] = d.Subject.map(groups)
    d.to_csv(os.path.join(a.out, "theta_per_cluster.csv"), index=False)
    n_all = len(d)
    d = d[d.Linearity >= a.min_linearity]
    print(f"clusters: {n_all} total, {len(d)} with Linearity >= {a.min_linearity}")

    per_sub = (d.groupby(["Subject", "Group"])
                 .agg(theta_med=("theta_med", "median"),
                      theta_mean=("theta_med", "mean"),
                      frac_aligned=("theta_med", lambda x: float((x < 30).mean())),
                      frac_perp=("theta_med", lambda x: float((x > 60).mean())),
                      cos2_med=("theta_med",
                                lambda x: float(np.median(np.cos(np.radians(x))**2))),
                      n_clusters=("PVS_ID", "size"))
                 .reset_index())

    # attach u_rel_det from the step16 QC
    if os.path.exists(a.qc_csv):
        qc = pd.read_csv(a.qc_csv)
        if "u_rel_det" in qc:
            per_sub = per_sub.merge(qc[["Subject", "u_rel_det"]],
                                    on="Subject", how="left")
    per_sub.to_csv(os.path.join(a.out, "theta_per_subject.csv"), index=False)

    print("\n=== per-group medians (uniform-null reference in brackets) ===")
    t = per_sub.groupby("Group")[["theta_med", "frac_aligned",
                                  "frac_perp", "cos2_med", "n_clusters"]].median()
    print(t.round(3).to_string())
    print(f"  [uniform: theta_med {UNIF_MEDIAN:.1f} | frac_aligned "
          f"{UNIF_FRAC30:.3f} | frac_perp {UNIF_FRAC60:.3f} | cos^2 {UNIF_COS2:.3f}]")

    print("\n=== ALIGNMENT vs the uniform null (across subjects) ===")
    for col, ref, direction in (("theta_med", UNIF_MEDIAN, "less"),
                                ("frac_aligned", UNIF_FRAC30, "greater"),
                                ("frac_perp", UNIF_FRAC60, "less"),
                                ("cos2_med", UNIF_COS2, "greater")):
        v = per_sub[col].dropna()
        try:
            p = wilcoxon(v - ref, alternative=direction).pvalue
        except ValueError:
            p = np.nan
        n_side = int(((v < ref) if direction == "less" else (v > ref)).sum())
        print(f"  {col:13s} median {v.median():7.3f} vs {ref:6.3f}  "
              f"({n_side}/{len(v)} subjects on the aligned side)  "
              f"Wilcoxon one-sided p = {p:.2e}")

    print("\n=== group comparisons (Kruskal-Wallis across subjects) ===")
    for col in ["theta_med", "frac_aligned", "frac_perp", "cos2_med"]:
        arrs = [g[col].dropna().values for _, g in per_sub.groupby("Group")]
        print(f"  {col:13s} KW p = {kruskal(*arrs).pvalue:.4f}")

    if "u_rel_det" in per_sub and per_sub["u_rel_det"].notna().any():
        r = per_sub[per_sub.u_rel_det < 0]["theta_med"]
        n = per_sub[per_sub.u_rel_det > 0]["theta_med"]
        if len(r) and len(n):
            p = mannwhitneyu(r, n).pvalue
            print(f"\n=== reflection consistency (must be null after v3) ===")
            print(f"  det=+1 {n.median():.1f} deg (n={len(n)}) vs "
                  f"det=-1 {r.median():.1f} deg (n={len(r)}), MW p = {p:.3f}"
                  + ("" if p > 0.05 else "   <-- [!] investigate"))

    print(f"\n-> {a.out}/theta_per_cluster.csv, theta_per_subject.csv")


if __name__ == "__main__":
    main()
