#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step25_make_supplementary_tables.py  (v1, 2026-10-17)

Builds the Supplementary Tables of the revised Brain Communications manuscript
from the statistics workbook written by step23_plot_all_markers_across_cohorts_and_ROIs.py
(stats_all_corrected_v13.xlsx) and, optionally, from the mixed-effects workbook
(MixedEffects_v12_3TIERS_cov_10k.xlsx, the primary run of step18_test_group_by_tier_interaction.py).

Tables produced
  Table S2  Within-PVS diffusion/kurtosis metrics and relative volumes, per region
            (family A: 13 measures x 3 regions): group descriptives, Kruskal-Wallis
            H / p / q, and the three pairwise Mann-Whitney contrasts with
            U / p / q / Cliff's delta [95% bootstrap CI].
  Table S3  PVS-weighted metrics (family E: 7 x 3), same layout.
  Table S4  Within-subject contrasts PVS minus reference tissue (adjacent WM 3-6 mm,
            matched WM >= 6 mm), from the mixed-effects workbook: per-group median
            difference, Wilcoxon p / q, RHI-vs-HC q (and adjusted p if present).
            Only written when --mixed-effects is given.
  Table S5  Group sizes per region after the technical exclusions, and the number of
            excluded participant x region cells by reason (value_exclusions sheet).
  Table S6  SynthSeg structure volumes relative to eTIV (family F), if present.

Descriptives are computed from the per-metric spreadsheets listed in the workbook's
'file_selection' sheet (the same files the statistics were computed from), applying
the same inclusion list (study xlsx, or the 'group_sizes' sheet when the study xlsx
is not reachable) and the same value exclusions. Median [IQR] by default
(--descriptive mean gives mean +/- SD).

Outputs (in --out): Supplementary_Tables_generated.docx, one CSV per table,
and Supplementary_Tables_generated.md.

Usage
  python step25_make_supplementary_tables.py \
      --stats  .../stats_families_AEF/based_on_.../stats_all_corrected_v13.xlsx \
      --mixed-effects MixedEffects_v12_3TIERS_cov_10k.xlsx \
      --out supplementary_tables
  Options: --data-root (root of the 'spreadsheets/' tree; default: found by walking up
  from --stats), --study-xlsx (overrides run_info), --descriptive median|mean,
  --list-excluded-ids (adds the pseudonymous IDs of excluded cells to Table S5).
"""
import argparse
import os
import re
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd

try:
    from docx import Document
    from docx.enum.section import WD_ORIENT
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.shared import Pt, Cm
except ImportError:                      # pragma: no cover
    sys.exit("python-docx is required:  pip install python-docx")

# ----------------------------------------------------------------------------
# naming
# ----------------------------------------------------------------------------
REGION_NAME = OrderedDict([("cerebral_white_matter", "Centrum semiovale"),
                           ("basal_ganglia", "Basal ganglia"),
                           ("midbrain", "Midbrain")])
GROUPS = [("hc", "HC"), ("ppcs", "PPCS"), ("boxers", "RHI")]
CONTRASTS = [("hc_vs_ppcs_u", "PPCS vs HC", -1),      # sign: workbook delta is first-vs-second
             ("hc_vs_boxers_u", "RHI vs HC", -1),
             ("ppcs_vs_boxers_u", "RHI vs PPCS", -1)]
# The workbook computes Cliff's delta as (group named first) - (group named second),
# i.e. hc_vs_boxers = HC relative to RHI. The figures and the manuscript report the
# contrast the other way round (RHI relative to HC), so the sign is flipped here.
# Set --keep-workbook-sign to disable the flip.

METRIC_NAME = {"MD": "MD", "AD": "AD", "RD": "RD", "FA": "FA",
               "MK_mf": "MK", "AK_mf": "AK", "RK_mf": "RK",
               "lowb_MD": "Low-b MD", "lowb_AD": "Low-b AD",
               "lowb_RD": "Low-b RD", "lowb_FA": "Low-b FA", "VOL": "Volume"}
MEASURE_ORDER = ["ROI volume / eTIV", "PVS volume / ROI volume",
                 "MD", "AD", "RD", "FA", "MK", "AK", "RK",
                 "Low-b MD", "Low-b AD", "Low-b RD", "Low-b FA"]
COLUMN_LABEL = {"diffusion_within_pvs": "",
                "PVS-Weighted-diffusion": "PVS-weighted ",
                "relative-tissue-size-in-mm": "ROI volume / eTIV",
                "relative-pvs-size-within-tissue-in-mm": "PVS volume / ROI volume",
                "tissue-size-in-mm": "ROI volume (mm3)",
                "pvs-size-within-tissue-in-mm": "PVS volume (mm3)"}
DIFFUSIVITIES = {"MD", "AD", "RD", "lowb_MD", "lowb_AD", "lowb_RD"}


def measure_label(metric, column):
    if column in ("relative-tissue-size-in-mm", "relative-pvs-size-within-tissue-in-mm",
                  "tissue-size-in-mm", "pvs-size-within-tissue-in-mm"):
        return COLUMN_LABEL[column]
    return COLUMN_LABEL.get(column, column + " ") + METRIC_NAME.get(metric, metric)


def unit_and_scale(metric, column):
    """Return (unit string for the table, multiplier applied to raw values)."""
    if column.startswith("relative"):
        return ("%", 100.0)
    if column in ("tissue-size-in-mm", "pvs-size-within-tissue-in-mm"):
        return ("mm3", 1.0)
    if metric in DIFFUSIVITIES:
        return ("10⁻³ mm²/s", 1e3)
    return ("", 1.0)


def norm_id(v):
    s = str(v).strip()
    if re.fullmatch(r"\d+(\.0+)?", s):
        return f"sub-{int(float(s))}"
    return s


# ----------------------------------------------------------------------------
# formatting helpers
# ----------------------------------------------------------------------------
def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "n/a"
    if p < 0.001:
        return "< 0.001"
    return f"{p:.3f}"


def fmt_num(x, nd=3):
    return "n/a" if x is None or not np.isfinite(x) else f"{x:.{nd}f}"


def descriptive(values, how, nd):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return "n/a"
    if how == "mean":
        return f"{v.mean():.{nd}f} ± {v.std(ddof=1):.{nd}f}" if v.size > 1 else f"{v.mean():.{nd}f}"
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    return f"{med:.{nd}f} [{q1:.{nd}f}, {q3:.{nd}f}]"


def decimals_for(metric, column):
    if column.startswith("relative"):
        return 2 if "pvs" in column else 2
    if column in ("tissue-size-in-mm", "pvs-size-within-tissue-in-mm"):
        return 0
    return 3


# ----------------------------------------------------------------------------
# loading
# ----------------------------------------------------------------------------
def find_data_root(stats_path, explicit):
    if explicit:
        if not os.path.isdir(os.path.join(explicit, "spreadsheets")):
            print(f"[!] --data-root {explicit} has no 'spreadsheets/' folder; "
                  "relative paths of file_selection may not resolve")
        return explicit
    d = os.path.dirname(os.path.abspath(stats_path))
    while True:
        if os.path.isdir(os.path.join(d, "spreadsheets")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            sys.exit("could not find a 'spreadsheets/' folder above the stats workbook; "
                     "give --data-root")
        d = parent


def load_inclusion(study_xlsx):
    sheets = {"HC": "hc", "PPCS": "ppcs", "Boxers": "boxers"}
    xls = pd.ExcelFile(study_xlsx)
    inc = {}
    for sheet, grp in sheets.items():
        if sheet not in xls.sheet_names:
            continue
        df = xls.parse(sheet)
        subj = "Subj" if "Subj" in df.columns else df.columns[0]
        col = "Will we used it?"
        if col not in df.columns:
            raise SystemExit(f"'{col}' missing in sheet {sheet} of {study_xlsx}")
        keep = df[df[col].astype(str).str.strip().str.lower() == "yes"][subj]
        inc[grp] = {norm_id(v) for v in keep.dropna()
                    if str(v).strip().lower() not in ("total:", "")}
    if set(inc) != {"hc", "ppcs", "boxers"}:
        raise SystemExit(f"study xlsx gave groups {sorted(inc)}; need hc, ppcs, boxers")
    return inc


def removed_from_group_sizes(gs):
    """{file basename -> {group -> set(removed ids)}} from the group_sizes sheet."""
    out = {}
    for _, r in gs.iterrows():
        rem = r.get("removed")
        ids = set()
        if isinstance(rem, str) and rem.strip():
            ids = {norm_id(x) for x in rem.split(",") if x.strip()}
        out.setdefault(os.path.normpath(str(r["file"])), {})[r["group"]] = ids
    return out


def read_metric_file(path, inc, removed_by_file, exclude, region, map_name, kind):
    dfs = pd.read_excel(path, sheet_name=["hc", "ppcs", "boxers"])
    key = os.path.normpath(path)
    for grp, df in dfs.items():
        ids = df.iloc[:, 0].map(norm_id)
        if inc is not None:
            keep = ids.isin(inc[grp])
        else:
            rem = None
            for k, v in removed_by_file.items():
                if os.path.basename(k) == os.path.basename(key) and \
                        os.path.basename(os.path.dirname(os.path.dirname(k))) == region:
                    rem = v.get(grp, set())
                    break
            keep = ~ids.isin(rem) if rem is not None else pd.Series(True, index=ids.index)
        df = df[keep.values]
        if kind == "diffusion" and map_name == "PVSROI" and exclude.get(region):
            df = df[~df.iloc[:, 0].map(norm_id).isin(exclude[region]).values]
        dfs[grp] = df
    return dfs


# ----------------------------------------------------------------------------
# tables
# ----------------------------------------------------------------------------
STRUCTURE_NAME = {
    "3rd_ventricle": "Third ventricle", "4th_ventricle": "Fourth ventricle", "CSF": "Extra-ventricular CSF",
    "lateral_ventricle": "Lateral ventricles", "inferior_lateral_ventricle": "Inferior lateral ventricles",
    "cerebral_white_matter": "Cerebral white matter", "cortex": "Cerebral cortex (all)",
    "cerebellum_cortex": "Cerebellar cortex", "cerebellum_white_matter": "Cerebellar white matter",
    "brain-stem": "Brainstem", "ventral_DC": "Ventral diencephalon", "thalamus": "Thalamus",
    "caudate": "Caudate", "putamen": "Putamen", "pallidum": "Pallidum", "hippocampus": "Hippocampus",
    "amygdala": "Amygdala", "accumbens_area": "Nucleus accumbens",
    "bankssts": "Banks of superior temporal sulcus", "caudalanteriorcingulate": "Caudal anterior cingulate",
    "caudalmiddlefrontal": "Caudal middle frontal", "cuneus": "Cuneus", "entorhinal": "Entorhinal",
    "frontalpole": "Frontal pole", "fusiform": "Fusiform", "inferiorparietal": "Inferior parietal",
    "inferiortemporal": "Inferior temporal", "insula": "Insula", "isthmuscingulate": "Isthmus cingulate",
    "lateraloccipital": "Lateral occipital", "lateralorbitofrontal": "Lateral orbitofrontal", "lingual": "Lingual",
    "medialorbitofrontal": "Medial orbitofrontal", "middletemporal": "Middle temporal", "paracentral": "Paracentral",
    "parahippocampal": "Parahippocampal", "parsopercularis": "Pars opercularis", "parsorbitalis": "Pars orbitalis",
    "parstriangularis": "Pars triangularis", "pericalcarine": "Pericalcarine", "postcentral": "Postcentral",
    "posteriorcingulate": "Posterior cingulate", "precentral": "Precentral", "precuneus": "Precuneus",
    "rostralanteriorcingulate": "Rostral anterior cingulate", "rostralmiddlefrontal": "Rostral middle frontal",
    "superiorfrontal": "Superior frontal", "superiorparietal": "Superior parietal",
    "superiortemporal": "Superior temporal", "supramarginal": "Supramarginal", "temporalpole": "Temporal pole",
    "transversetemporal": "Transverse temporal",
}


def build_family_table(xls, family, data_root, inc, removed_by_file, exclude, how, flip):
    fam = xls.parse("fdr_families")
    fam = fam[fam.family == family]
    if fam.empty:
        return None
    kw = xls.parse("Kruskal-Wallis").set_index("Case")
    pw = {k: xls.parse(k).set_index("Case") for k, _, _ in CONTRASTS}
    sel = xls.parse("file_selection")
    sel = sel[sel.used.astype(str).str.lower().isin(["true", "1"])]

    def segmentation_of(region, case):
        # Case = <region>_<segmentation>_<metric>_...; the segmentation tells
        # PVSROI (study ROI) and synthseg (whole structure) apart for the same region
        for seg in sel.segmentation.dropna().unique():
            if str(case).startswith(f"{region}_{seg}_"):
                return seg
        return None

    def file_for(region, metric, column, seg):
        tree = "volumetry" if column.startswith(("relative", "tissue", "pvs")) else "diffusion"
        m = sel[(sel.region == region) & (sel.metric == metric) & (sel.tree == tree)]
        if seg is not None:
            m = m[m.segmentation == seg]
        if m.empty:
            return None
        return os.path.join(data_root, str(m.iloc[0]["file"]))

    rows = []
    cache = {}
    for _, f in fam.iterrows():
        region, metric, column, case = f.region, f.metric, f.column, f.Case
        seg = segmentation_of(region, case)
        # REGION_NAME names the PVSROI study ROIs (cerebral_white_matter = centrum
        # semiovale); a synthseg row of the same region is the whole structure
        if region not in REGION_NAME or seg == "synthseg":
            region_label = STRUCTURE_NAME.get(str(region), str(region).replace("_", " "))
        else:
            region_label = REGION_NAME[region]
        path = file_for(region, metric, column, seg)
        unit, scale = unit_and_scale(metric, column)
        nd = decimals_for(metric, column)
        desc = {}
        if path and os.path.exists(path):
            kind = "volumetry" if "/volumetry/" in path.replace("\\", "/") else "diffusion"
            map_name = os.path.basename(os.path.dirname(path))
            if path not in cache:
                cache[path] = read_metric_file(path, inc, removed_by_file, exclude,
                                               region, map_name, kind)
            dfs = cache[path]
            for g, G in GROUPS:
                if column in dfs[g].columns:
                    desc[G] = descriptive(dfs[g][column].values * scale, how, nd)
                else:
                    desc[G] = "n/a"
        else:
            for _, G in GROUPS:
                desc[G] = "[[file not found]]"
        if case not in kw.index:
            print(f"[!] {case}: not in Kruskal-Wallis sheet, skipped")
            continue
        k = kw.loc[case]
        row = OrderedDict()
        row["Region"] = region_label
        row["Measure"] = measure_label(metric, column)
        row["Unit"] = unit
        row["n (HC/PPCS/RHI)"] = f"{int(k.n_hc)}/{int(k.n_ppcs)}/{int(k.n_boxers)}"
        for _, G in GROUPS:
            row[G] = desc.get(G, "n/a")
        row["Kruskal–Wallis H"] = fmt_num(float(k.H_stat), 2)
        row["KW p"] = fmt_p(float(k.p_value))
        row["KW q"] = fmt_p(float(k.p_value_FDR))
        for key, label, _ in CONTRASTS:
            if case in pw[key].index:
                r = pw[key].loc[case]
                sgn = -1.0 if flip else 1.0
                d, lo, hi = (sgn * float(r.delta_cl), sgn * float(r.delta_ci_lo),
                             sgn * float(r.delta_ci_hi))
                if flip:
                    lo, hi = hi, lo
                row[f"{label}: U"] = fmt_num(float(r.u_stat), 0)
                row[f"{label}: p"] = fmt_p(float(r.p_value))
                row[f"{label}: q"] = fmt_p(float(r.p_value_FDR))
                row[f"{label}: Δ [95% CI]"] = f"{d:+.2f} [{lo:+.2f}, {hi:+.2f}]"
            else:
                for suf in ("U", "p", "q", "Δ [95% CI]"):
                    row[f"{label}: {suf}"] = "n/a"
        rows.append(row)
    df = pd.DataFrame(rows)
    # order
    reg_order = {v: i for i, v in enumerate(REGION_NAME.values())}
    meas_order = {v: i for i, v in enumerate(MEASURE_ORDER)}
    df["_r"] = df.Region.map(lambda r: reg_order.get(r, 99))
    df["_m"] = df.Measure.map(lambda m: meas_order.get(m.replace("PVS-weighted ", ""), 99))
    df = df.sort_values(["_r", "_m", "Measure"]).drop(columns=["_r", "_m"]).reset_index(drop=True)
    return df


def build_exclusion_table(xls, list_ids):
    ve = xls.parse("value_exclusions")
    kw = xls.parse("Kruskal-Wallis")
    rows = []
    for region, rname in REGION_NAME.items():
        sub = ve[ve.region == region]
        k = kw[(kw.ROI == region) & (kw.Column == "diffusion_within_pvs")]
        n = k.iloc[0] if not k.empty else None
        for g, G in GROUPS:
            s = sub[sub.group == g]
            exc = s[s.excluded.astype(int) == 1]
            clauses = []
            for _, e in exc.iterrows():
                parts = [x.strip() for x in re.split(r";\s*", str(e.reasons)) if x.strip()]
                parts = [re.sub(r"^(\d+) PVS voxels < (\d+)$", lambda m: f"PVS mask {m.group(1)} voxels (< {m.group(2)})", p)
                         for p in parts]
                def _val(m):
                    return f"{m.group(1)} = {float(m.group(2)) * 1e3:.2f} × 10⁻³ mm²/s (> 3)"
                parts = [re.sub(r"^(\w+) = ([0-9.e\-+]+) mm2/s outside \[0, 3\.0e-03\]$", _val, p) for p in parts]
                parts = [p.replace("lowb_", "low-b ") for p in parts]
                clause = "; ".join(parts)
                if list_ids:
                    clause = f"{e.subject}: {clause}"
                clauses.append(clause)
            reason_txt = " | ".join(clauses) if clauses else "—"
            rows.append(OrderedDict([
                ("Region", rname), ("Group", G),
                ("Participants with PVS in region", int(len(s))),
                ("Excluded cells", int(len(exc))),
                ("Reason(s)", reason_txt),
                ("n tested (within-PVS metrics)",
                 "n/a" if n is None else int(n[f"n_{g}"])),
                ("Median PVS voxels (all)", "n/a" if s.empty else f"{s.pvs_voxels.median():.0f}"),
                ("Min PVS voxels (kept)",
                 "n/a" if s[s.excluded.astype(int) == 0].empty
                 else f"{s[s.excluded.astype(int) == 0].pvs_voxels.min():.0f}"),
            ]))
    return pd.DataFrame(rows)


KEEP_BETA = False


def build_within_subject_table(me_path):
    xls = pd.ExcelFile(me_path)
    sheet = "within_subject" if "within_subject" in xls.sheet_names else \
        next((s for s in xls.sheet_names if "within" in s.lower()), None)
    if sheet is None:
        print(f"[!] {me_path}: no within_subject sheet (sheets: {xls.sheet_names})")
        return None, None
    W = xls.parse(sheet)
    ri = xls.parse("run_info") if "run_info" in xls.sheet_names else None
    cols = list(W.columns)
    print(f"mixed-effects workbook sheet '{sheet}' columns: {cols}")
    ref_col = "Ref_tier" if "Ref_tier" in cols else next((c for c in cols if "ref" in c.lower()), None)
    met_col = "Metric" if "Metric" in cols else cols[0]
    adj_cols = [c for c in cols if re.search(r"adj", c, re.I) and re.search(r"RHI", c)]
    rows = []
    ref_name = {"NearPVS": "Adjacent WM (3–6 mm)", "Fake_m6": "Matched WM (≥6 mm)",
                "Fake_m10": "Matched WM (≥10 mm)"}
    pretty = {"MD": "MD", "AD": "AD", "RD": "RD", "FA": "FA", "MK_mf": "MK", "AK_mf": "AK",
              "RK_mf": "RK", "Dpar": "D∥", "Dperp": "D⊥", "Dpar_minus_Dperp": "D∥ − D⊥",
              "lowb_MD": "Low-b MD", "lowb_AD": "Low-b AD", "lowb_RD": "Low-b RD",
              "lowb_FA": "Low-b FA", "lowb_ADC": "Low-b ADC"}
    order = ["MD", "AD", "RD", "FA", "MK_mf", "AK_mf", "RK_mf", "Dpar", "Dperp",
             "Dpar_minus_Dperp", "lowb_MD", "lowb_AD", "lowb_RD", "lowb_FA"]
    diffus = {"MD", "AD", "RD", "Dpar", "Dperp", "Dpar_minus_Dperp", "lowb_MD", "lowb_AD",
              "lowb_RD", "lowb_ADC"}
    for ref in ["NearPVS", "Fake_m6", "Fake_m10"]:
        sub = W[W[ref_col] == ref] if ref_col else W
        if sub.empty:
            continue
        for met in order:
            r = sub[sub[met_col] == met]
            if r.empty:
                continue
            r = r.iloc[0]
            sc = 1e3 if met in diffus else 1.0
            row = OrderedDict([("Reference tissue", ref_name.get(ref, ref)),
                               ("Metric", pretty.get(met, met)),
                               ("Unit", "10⁻³ mm²/s" if met in diffus else "")])
            for G in ["HC", "PPCS", "RHI"]:
                d = r.get(f"delta_{G}", np.nan)
                p = r.get(f"p_{G}_vs_0", np.nan)
                q = r.get(f"q_{G}_vs_0", np.nan)
                row[f"{G}: median PVS − ref"] = fmt_num(float(d) * sc, 3) if np.isfinite(float(d)) else "n/a"
                row[f"{G}: Wilcoxon p"] = fmt_p(float(p)) if np.isfinite(float(p)) else "n/a"
                row[f"{G}: q"] = fmt_p(float(q)) if np.isfinite(float(q)) else "n/a"
            for c in ["p_RHI_vs_HC", "q_RHI_vs_HC"]:
                if c in cols:
                    row["RHI vs HC: " + c.split("_")[0]] = fmt_p(float(r[c]))
            for c in adj_cols:
                v = r[c]
                if c.startswith("p_"):
                    row["RHI vs HC: adjusted p"] = fmt_p(float(v)) if np.isfinite(float(v)) else "n/a"
                elif c.startswith("beta_"):
                    if KEEP_BETA:
                        row["RHI vs HC: adjusted β (workbook units)"] = fmt_num(float(v), 3) if np.isfinite(float(v)) else "n/a"
                else:
                    row[f"RHI vs HC: {c}"] = fmt_p(float(v)) if isinstance(v, (int, float, np.floating)) and np.isfinite(float(v)) else str(v)
            rows.append(row)
    note = ""
    if ri is not None and len(ri.columns) >= 2:
        kv = dict(zip(ri.iloc[:, 0].astype(str), ri.iloc[:, 1].astype(str)))
        note = "; ".join(f"{k} = {v}" for k, v in kv.items()
                         if any(t in k.lower() for t in ("covar", "perm", "n_subj", "within", "fdr", "metric")))
    return pd.DataFrame(rows), note


def build_interaction_table(me_path):
    """Group × tissue-tier interaction tests from the 'interaction_tests' sheet of the mixed-effects
    workbook (one model per metric over all tiers; participant as random effect / cluster)."""
    xls = pd.ExcelFile(me_path)
    if "interaction_tests" not in xls.sheet_names:
        print(f"[!] {me_path}: no interaction_tests sheet (sheets: {xls.sheet_names})")
        return None, ""
    T = xls.parse("interaction_tests")
    ri = xls.parse("run_info") if "run_info" in xls.sheet_names else None
    pretty = {"MD": "MD", "AD": "AD", "RD": "RD", "FA": "FA", "MK_mf": "MK", "AK_mf": "AK",
              "RK_mf": "RK", "Dpar": "D∥", "Dperp": "D⊥", "Dpar_minus_Dperp": "D∥ − D⊥",
              "lowb_MD": "Low-b MD", "lowb_AD": "Low-b AD", "lowb_RD": "Low-b RD", "lowb_FA": "Low-b FA"}
    order = ["MD", "AD", "RD", "FA", "MK_mf", "AK_mf", "RK_mf", "Dpar", "Dperp",
             "Dpar_minus_Dperp", "lowb_MD", "lowb_AD", "lowb_RD", "lowb_FA"]
    method_name = {"OLS+cluster-robust": "OLS, cluster-robust SE (Wald F)", "MixedLM(ML)": "Mixed model, ML (LR χ²)"}
    tier_name = {"PVS": "PVS", "NearPVS": "adjacent WM (3–6 mm)", "NearPVSfull": "adjacent WM (full shell)",
                 "Fake_m6": "matched WM (≥6 mm)", "Fake_m10": "matched WM (≥10 mm)", "WholeWM": "whole WM"}
    def P(v):
        return fmt_p(float(v)) if v is not None and np.isfinite(float(v)) else "n/a"
    rows = []
    for met in order:
        r = T[T["Metric"] == met]
        if r.empty:
            continue
        r = r.iloc[0]
        row = OrderedDict([("Metric", pretty.get(met, met)),
                           ("n", int(r["n_subjects"]) if "n_subjects" in T.columns else ""),
                           ("Model", method_name.get(str(r.get("method", "")), str(r.get("method", "")))),
                           ("Group × tier: statistic (LR χ² or Wald F)", fmt_num(float(r["stat"]), 2)),
                           ("df", int(r["df"])),
                           ("Group × tier: p", P(r.get("p_group_x_tier"))),
                           ("Group × tier: q", P(r.get("q_group_x_tier"))),
                           ("Group × tier: permutation p", P(r.get("p_perm_allgroups"))),
                           ("Group × tier: permutation q", P(r.get("q_perm_allgroups"))),
                           ("RHI vs HC × tier: df", int(r["df_HCvRHI"]) if "df_HCvRHI" in T.columns else ""),
                           ("RHI vs HC × tier: permutation p", P(r.get("p_HCvRHI_perm"))),
                           ("RHI vs HC × tier: permutation q", P(r.get("q_HCvRHI_perm")))])
        rows.append(row)
    df = pd.DataFrame(rows)
    note = ""
    if ri is not None and len(ri.columns) >= 2:
        kv = dict(zip(ri.iloc[:, 0].astype(str), ri.iloc[:, 1].astype(str)))
        tiers = [tier_name.get(t.strip(), t.strip()) for t in str(kv.get("tiers", "")).split(",") if t.strip()]
        tier_codes = [t.strip() for t in str(kv.get("tiers", "")).split(",") if t.strip()]
        # the full-shell tier is missing for D∥, D⊥ and D∥ − D⊥, so their df is
        # smaller only in workbooks that include that tier
        full_shell = ("; D∥, D⊥ and D∥ − D⊥ lack the full-shell tier, hence their smaller df"
                      if "NearPVSfull" in tier_codes else "")
        cov_name = {"icv_ml": "intracranial volume",
                    "wm_absldj_med": "the median absolute log-Jacobian determinant of the non-linear "
                                     "warp within white matter",
                    "scale_x": "the affine left–right scaling factor"}
        covs = [cov_name.get(c.strip(), c.strip()) for c in str(kv.get("covariates_used", "")).split(",")
                if c.strip() and c.strip().lower() != "nan"]
        cov_txt = (", ".join(covs[:-1]) + " and " + covs[-1]) if len(covs) > 1 else (covs[0] if covs else "none")
        note = (f"Tissue tiers entered in the models ({len(tiers)}): {', '.join(tiers)}{full_shell}. "
                f"Covariates: {cov_txt}; n = {kv.get('n_subjects', '?')}.")
    if not df.empty:
        q_perm = pd.to_numeric(T["q_perm_allgroups"], errors="coerce")
        q_wald = pd.to_numeric(T["q_group_x_tier"], errors="coerce")
        i1, i2 = q_perm.idxmin(), q_wald.idxmin()
        print(f"interaction: smallest permutation q = {q_perm.min():.3f} ({pretty.get(T.loc[i1, 'Metric'])}); "
              f"smallest Wald q = {q_wald.min():.3f} ({pretty.get(T.loc[i2, 'Metric'])})")
    return df, note


# ----------------------------------------------------------------------------
# docx writer
# ----------------------------------------------------------------------------
def add_table(doc, title, note, df, font_pt=7):
    doc.add_paragraph().add_run(title).bold = True
    if df is None or df.empty:
        doc.add_paragraph("[[no data]]")
        return
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Table Grid"
    for i, c in enumerate(df.columns):
        cell = t.rows[0].cells[i]
        cell.text = ""
        run = cell.paragraphs[0].add_run(str(c))
        run.bold = True
        run.font.size = Pt(font_pt)
    for _, r in df.iterrows():
        cells = t.add_row().cells
        for i, c in enumerate(df.columns):
            cells[i].text = ""
            run = cells[i].paragraphs[0].add_run(str(r[c]))
            run.font.size = Pt(font_pt)
    if note:
        p = doc.add_paragraph()
        run = p.add_run(note)
        run.font.size = Pt(font_pt + 1)
        run.italic = True
    doc.add_paragraph()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stats", required=True, help="stats_all_corrected_v13.xlsx")
    ap.add_argument("--mixed-effects", default=None, help="MixedEffects_v11_ALLTIERS_cov_10k.xlsx")
    ap.add_argument("--data-root", default=None)
    ap.add_argument("--study-xlsx", default=None)
    ap.add_argument("--descriptive", choices=["median", "mean"], default="median")
    ap.add_argument("--keep-workbook-sign", action="store_true",
                    help="report Cliff's delta as stored (first-named group minus second) "
                         "instead of RHI-vs-HC / RHI-vs-PPCS / PPCS-vs-HC")
    ap.add_argument("--list-excluded-ids", action="store_true")
    ap.add_argument("--keep-beta", action="store_true",
                    help="keep the adjusted-β column of Table 4 (its units in the workbook are not those of the deltas)")
    ap.add_argument("--out", default="supplementary_tables")
    a = ap.parse_args()
    global KEEP_BETA
    KEEP_BETA = a.keep_beta

    os.makedirs(a.out, exist_ok=True)
    xls = pd.ExcelFile(a.stats)
    ri = xls.parse("run_info")
    kv = dict(zip(ri.key.astype(str), ri.value.astype(str)))
    data_root = find_data_root(a.stats, a.data_root)
    print(f"data root: {data_root}")

    study = a.study_xlsx or kv.get("study_xlsx")
    inc = None
    if study and os.path.exists(study):
        inc = load_inclusion(study)
        print(f"inclusion from {study}: " + ", ".join(f"{g}={len(v)}" for g, v in inc.items()))
    else:
        print(f"[i] study xlsx not reachable ({study}); using the 'group_sizes' sheet "
              "to reproduce the subject selection")
    removed_by_file = removed_from_group_sizes(xls.parse("group_sizes"))

    ve = xls.parse("value_exclusions")
    exclude = {}
    for _, r in ve[ve.excluded.astype(int) == 1].iterrows():
        exclude.setdefault(r.region, set()).add(norm_id(r.subject))
    print("value exclusions: " + (", ".join(f"{k}: {len(v)}" for k, v in exclude.items()) or "none"))

    flip = not a.keep_workbook_sign
    how = a.descriptive
    desc_txt = ("median [interquartile range]" if how == "median" else "mean ± SD")

    tA = build_family_table(xls, "A", data_root, inc, removed_by_file, exclude, how, flip)
    tE = build_family_table(xls, "E", data_root, inc, removed_by_file, exclude, how, flip)
    tF = build_family_table(xls, "F", data_root, inc, removed_by_file, exclude, how, flip)
    tX = build_exclusion_table(xls, a.list_excluded_ids)
    tW, wnote, tI, inote = None, "", None, ""
    if a.mixed_effects:
        tW, wnote = build_within_subject_table(a.mixed_effects)
        tI, inote = build_interaction_table(a.mixed_effects)

    common_note = (f"Values are {desc_txt} per group; diffusivities in 10⁻³ mm²/s, relative volumes in %. "
                   "Kruskal–Wallis H test across the three groups; q = Benjamini–Hochberg-corrected p within the "
                   "family of tests ({fam}). Pairwise Mann–Whitney U tests; q corrected over the three contrasts "
                   "of all tests in the family together. Cliff's Δ is positive when the first-named group is "
                   "larger; 95% CI from a percentile bootstrap ({ci}). n = participants with valid PVS values "
                   "in the region after the technical exclusions (min. PVS voxels = {mv}; diffusivity range "
                   "{md}). HC, healthy controls; PPCS, athletes with persistent post-concussion symptoms; RHI, "
                   "elite boxers exposed to repetitive head impacts; MD/AD/RD, mean/axial/radial diffusivity; "
                   "FA, fractional anisotropy; MK/AK/RK, mean/axial/radial kurtosis; low-b, tensor fit to "
                   "b ≤ 500 s/mm²; ROI, region of interest; eTIV, estimated total intracranial volume.")
    ci_txt = kv.get("ci", "percentile bootstrap")
    noteA = common_note.format(fam=f"family A, {kv.get('family_sizes_kruskal_wallis', '')}",
                               ci=ci_txt, mv=kv.get("min_pvs_voxels", "?"),
                               md=kv.get("max_diffusivity_mm2_s", "?"))
    noteE = common_note.format(fam="family E, PVS-weighted metrics (median within PVS × PVS volume / ROI volume)",
                               ci=ci_txt, mv=kv.get("min_pvs_voxels", "?"),
                               md=kv.get("max_diffusivity_mm2_s", "?"))
    noteF = (f"Structure volumes from the SynthSeg parcellation of the T1-weighted image (cortical regions are "
             "bilateral Desikan–Killiany parcels), expressed as a "
             f"percentage of eTIV. Values are {desc_txt} per group. Kruskal–Wallis H test across the three groups; "
             "q = Benjamini–Hochberg-corrected p within this family of 52 omnibus tests; pairwise Mann–Whitney U "
             "tests, q corrected over the three contrasts of all 52 structures together. Cliff's Δ is positive when "
             f"the first-named group is larger; 95% CI from a percentile bootstrap ({ci_txt}). n = all included "
             "participants (no technical exclusions apply to volumes). HC, healthy controls; PPCS, athletes with "
             "persistent post-concussion symptoms; RHI, elite boxers exposed to repetitive head impacts; CSF, "
             "cerebrospinal fluid; eTIV, estimated total intracranial volume.")
    noteX = ("A participant × region cell was excluded from all within-PVS diffusion metrics when the PVS mask "
             f"held fewer than {kv.get('min_pvs_voxels', '?')} voxels or any fitted diffusivity/FA was outside "
             f"the plausible range ({kv.get('max_diffusivity_mm2_s', '?')}). Volumetric measures were never "
             f"excluded. Workbook summary: {kv.get('excluded_cells', '')}.")
    noteW = ("Within-subject contrast = median within PVS minus median within the participant's own reference "
             "tissue, projected onto the PVS axis for D∥ and D⊥. Wilcoxon signed-rank test against zero per "
             "group; q = Benjamini–Hochberg across the metrics tested against the same reference tissue. "
             "RHI vs HC: difference of the contrast between groups (q unadjusted; adjusted columns, where "
             "present, include the covariates listed). " + wnote)

    tables = [("Supplementary Table 2. Within-PVS diffusion and kurtosis metrics and relative volumes by region "
               "(primary family of tests).", noteA, tA, "TableS2_family_A.csv"),
              ("Supplementary Table 3. PVS-weighted diffusion and kurtosis metrics by region.", noteE, tE,
               "TableS3_family_E.csv")]
    noteI = ("One model per metric with the per-participant median value as outcome; fixed effects: group (HC, "
             "PPCS, RHI), tissue tier (listed below) and their interaction; "
             "participant as random effect (mixed model fitted by maximum likelihood, likelihood-ratio test against "
             "the model without interaction) or, where the mixed model did not converge, as cluster for robust "
             "standard errors (OLS, Wald F test). Group × tier: joint test of all interaction terms; permutation p "
             "from 10,000 relabellings of group between participants, keeping each participant's tiers together; RHI vs HC × tier: the same restricted to the HC–RHI contrast. q = Benjamini–Hochberg "
             "across the 14 metrics. A significant interaction would indicate that the PVS-versus-tissue "
             "difference depends on group. " + inote)
    if tW is not None:
        tables.append(("Supplementary Table 4. Within-subject contrasts: PVS minus reference tissue.", noteW, tW,
                       "TableS4_within_subject.csv"))
    if tI is not None and not tI.empty:
        tables.append(("Supplementary Table 4 (continued). Group × tissue-tier interaction tests from the "
                       "mixed-effects models.", noteI, tI, "TableS4b_interaction.csv"))
    tables.append(("Supplementary Table 5. Group sizes per region and technical exclusions.", noteX, tX,
                   "TableS5_exclusions.csv"))
    if tF is not None and not tF.empty:
        tables.append(("Supplementary Table 6. SynthSeg structure volumes relative to eTIV.", noteF, tF,
                       "TableS6_family_F.csv"))

    doc = Document()
    sec = doc.sections[0]
    sec.orientation = WD_ORIENT.LANDSCAPE
    sec.page_width, sec.page_height = sec.page_height, sec.page_width
    for side in ("left_margin", "right_margin", "top_margin", "bottom_margin"):
        setattr(sec, side, Cm(1.5))
    style = doc.styles["Normal"]
    style.font.name = "Arial"
    style.font.size = Pt(9)
    h = doc.add_paragraph()
    h.add_run("Supplementary Tables (generated)").bold = True
    doc.add_paragraph(f"Source workbook: {os.path.abspath(a.stats)} (script {kv.get('script', '?')}, "
                      f"{kv.get('date', '?')}; included {kv.get('included', '?')}; "
                      f"pairwise FDR {kv.get('pairwise_fdr', '?')}; low-b family {kv.get('lowb_family', '?')})."
                      + (f" Mixed-effects workbook: {os.path.abspath(a.mixed_effects)}." if a.mixed_effects else ""))

    md_lines = ["# Supplementary Tables (generated)", ""]
    for title, note, df, csv in tables:
        add_table(doc, title, note, df)
        if df is not None:
            df.to_csv(os.path.join(a.out, csv), index=False)
            md_lines += [f"**{title}**", "", df.to_markdown(index=False), "", note, ""]
    docx_path = os.path.join(a.out, "Supplementary_Tables_generated.docx")
    doc.save(docx_path)
    with open(os.path.join(a.out, "Supplementary_Tables_generated.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(md_lines))
    print(f"wrote {docx_path}")
    for title, _, df, csv in tables:
        print(f"  {csv}: {0 if df is None else len(df)} rows")


if __name__ == "__main__":
    main()
