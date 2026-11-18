#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step19_plot_pvs_tier_figures.py
===============================
v16-3 (2026-10-16): the word "control" no longer appears in any figure text.
In a study with a healthy-control group it reads as "the HC group". The
matched reference tissue is "matched white matter" everywhere (as in Figure 2);
the tissue available for matching is "eligible white matter". Fig6 title:
"Matched white matter: comparability across groups"; its first panel:
"Eligible white matter", y-axis "Eligible WM / WM volume". Fig3: y-axis
"Eligible white matter (relative to WM volume)", title "Eligible white matter
by exclusion margin". A Kruskal-Wallis p below 0.001 prints as "p < 0.001"
instead of "p = 0.000". Statistics and all other figures unchanged.

v16-2 (2026-10-05): schematic E, left drawing: only the adjacent white matter
(the 3-6 mm ring between the two dashed circles) is shaded and carries the
fibre lines; the PVS with its 0-3 mm margin and everything beyond 6 mm stay
blank, so the drawing shows where the adjacent-WM fibre angle is measured.
The note under it reads "shaded ring". Statistics unchanged; everything else
as v16-1.

v16-1 (2026-10-04): the right-hand schematic compared the matched-WM angle with
"60 deg", the median angle of ONE random direction. The matched-WM angle is the
mean over the control draws of each cluster, whose random reference is ~58 deg
(and cos^2 1/3). The label now reports cos^2 against 0.33, which needs no such
caveat. Everything else as v16.

CHANGELOG v16 (2026-10-04)
---------------------------
 * --fdr-families single is now the default, matching mixed-effects v11 (low-b
   and full-range fits are one family). The run still stops if the workbook
   used a different scheme.
 * --exclude-metrics (default lowb_ADC): not plotted and not in any FDR family.

CHANGELOG v15 (2026-10-04)
---------------------------
 * Reads the v9 extraction by default (pvs_tier_metrics_v9.xlsx),
   which adds the fibre angle measured in the adjacent and matched white matter.
 * Schematic (E): both drawings use MEASURED angles (median over PVS with a
   well-defined axis, SE <= 5 deg), computed at run time: left = angle between
   the PVS axis and the fibres of the adjacent WM; right = the same axis vs the
   fibres of the matched WM (>= 6 mm). Without the v9 columns it falls back to
   v14 (angle at the PVS; matched WM 'not preferred').
 * Alignment panels plot each participant's MEAN cos^2 (v12-v14 plotted the
   median, whose random-direction reference is 0.25, not the 1/3 drawn).
 * --alignment-panel precision|tier for panel G: 'precision' (default, = v14)
   shows alignment at the PVS by precision of the axis; 'tier' shows alignment
   in PVS, adjacent WM and matched WM side by side (axes with SE <= 5 deg).

CHANGELOG v14 (2026-10-04)
---------------------------
 * FA removed from the main figure (row A: MD, D||, D-perp, MK -> panels A-D);
   the schematic is E, D|| - D-perp by tier F, PVS-fibre alignment G.
 * Schematic (E) redrawn. Left: the PVS with its axis u and the white-matter
   bundle running past it, crossing u at the MEASURED angle theta (median
   angle between the PVS axis and the principal fibre direction at the PVS,
   for well-defined axes, computed from the data at run time), with the
   adjacent band (3-6 mm) marked: the same bundle. Right, zoomed out: the
   small PVS and, >= 6 mm away, the matched white matter with its own bundle.
   No angle is measured there (the extraction stores theta at the PVS only);
   the drawing says what the data show: D|| ~ D-perp, as for fibres with no
   preferred orientation relative to u.

CHANGELOG v13 (2026-10-04) -- main Figure 2 simplified
-------------------------------------------------------
 Row A (A-E: MD, D||, D-perp, MK, FA): only PVS and adjacent WM (3-6 mm);
   no statistics text; significance as stars on a light line joining each
   group's PVS violin to its adjacent-WM violin (q from the workbook, PVS vs
   adjacent, within group). y-axis starts at the data minimum.
 Row B: F = schematic (PVS cluster with its axis, adjacent fibres running along
   it, matched WM with fibres at another angle, and what D|| - D-perp means);
   G = D|| - D-perp by tier (PVS, adjacent, matched) with stars on lines for
   PVS vs adjacent and PVS vs matched; H = PVS-fibre alignment by axis
   precision. G and H sit to the right; F on the left.
 Matched-WM panels for the row-A metrics remain in the supplementary figures.

CHANGELOG v12 (2026-10-04) -- main Figure 2 split into two messages
--------------------------------------------------------------------
 Row A (five panels: MD, D-parallel, D-perp, MK, FA): "PVS are a faster,
   less directional, more heterogeneous medium than their surroundings".
 Row B: B1 = along-axis anisotropy (D|| - D-perp) by tissue tier, showing it is
   shared with the adjacent white matter; B2 = alignment of each PVS axis with
   the local fibre direction (cos^2 theta) by precision of the axis estimate,
   showing that PVS run along fibres (random = 0.33). Together: "the along-axis
   anisotropy is inherited from the fibres PVS run along".
 q-lines as before (mixed-effects workbook, PVS vs adjacent and vs matched WM).
 --main-metrics now sets the row-A metrics (default MD Dpar Dperp MK_mf FA).
 Output: Fig_main_PVS_distinct.png/.pdf (same name as before).

v11-1 (2026-10-03): the control-pool figure (Fig6) gets a bold y-axis title on
every panel ("Control pool size / WM volume" or "Clusters (%)"); panel titles
name the quantity and the margin. Everything else as v11.

CHANGELOG v11 (2026-10-03) -- display only, statistics unchanged
------------------------------------------------------------------
 * READABLE NAMES EVERYWHERE. No internal code appears in any figure text:
     tiers    : NearPVS -> "Adjacent WM (3-6 mm)", NearPVSfull -> "Adjacent WM
                (3-6 mm, full shell)", Fake m6 / m10 -> "Matched WM (>=6 mm)" /
                "(>=10 mm)" (the same terms in every figure)
     metrics  : MK_mf -> MK, Dpar -> D-parallel, Dpar_minus_Dperp -> D|| - D-perp,
                lowb_* -> "low-b ...", with the full name in titles
     QC panels: pool_frac_m6 -> "Control pool size (>=6 mm)", fakeFail ->
                "Unmatched clusters (%)", relax_region -> "Depth-relaxed
                matches (%)", relax_any -> "Matched in another region (%)"
 * Control-pool figure (Fig6): relaxations and failures are PERCENTAGES of each
   subject's clusters (Subject_QC n_clusters), not counts, so high-burden
   subjects are not penalised for having more PVS; one row per margin; panels
   whose values are all zero are dropped and named in a footnote ("no region
   relaxation was needed at 6 mm"); footnote states the Kruskal-Wallis
   p-values are uncorrected.
 * Diffusivities plotted in 10^-3 mm^2/s in every figure (Fig1/2/5 were raw).
 * Fonts: all text larger; axis titles bold and larger (FS_LAB 22); tick labels
   16; annotations 15; the main and supplementary figures scaled accordingly.
Figures + statistics from pvs_tier_metrics_v9.xlsx, with the
within-subject statistics taken from the v10 mixed-effects workbook so that
figures and text carry identical numbers.

CHANGELOG v10 (2026-09-29)
--------------------------
 * FDR FAMILIES, matching mixed-effects v10. --fdr-families split (default):
   the depth-bin FDR (Fig2) runs separately for the main metrics (MD, AD, RD,
   FA, MK, AK, RK, Dpar, Dperp, Dpar_minus_Dperp) and for the lowb_* metrics
   (the low-b fit, reported as the answer to Reviewer 3, major concern 2).
   'single' = one family (v9).
 * The workbook's own family setting (run_info 'fdr_families'; absent = v9 =
   single) must equal --fdr-families, otherwise the run stops: the main and
   supplementary figures print the workbook's q-values, and the figures must
   not mix two correction schemes.
 * Statistics file: Figure_statistics_v10.xlsx (run_info records the families).

CHANGELOG v9 (2026-09-29)
-------------------------
NEW FIGURES
 * Fig_main_PVS_distinct (new main Figure 2): "PVS are distinct from their
   surroundings, in every group". Tier on x (PVS, adjacent WM 3-6 mm, matched
   WM >= 6 mm), group as hue, one panel per metric (default Dpar-Dperp, Dpar,
   Dperp, MK, FA, MD). Under each panel: within-group PVS-vs-adjacent q-values
   read from the v9 workbook (Ref_tier = NearPVS). PNG + PDF.
 * Fig_suppl_ref-<tier> (regenerated Suppl. Figs 10/11): paired PVS vs
   reference per group, one figure per reference tier (default NearPVS and
   Fake_m6), subject lines joining the pairs, median paired difference in %
   with stars from the v9 workbook, and the RHI-vs-HC test of the PVS EXCESS
   (q_RHI_vs_HC), which is the pairwise version of the interaction. PNG + PDF.
FIXES
 * ORIENTATION FILTER. v8 kept clusters with Axis_SE_deg <= 10 OR NaN. The
   keep-NaN rule was copied from the Linearity filter, where NaN meant
   'perfectly linear'; for Axis_SE_deg NaN means a 3-voxel cluster with NO axis
   estimate -- the worst clusters, not the best (4,699 of 22,255). v8 also
   applied the filter while the mixed-effects script applies none, so Dpar /
   Dperp figures and statistics came from different cluster sets. v9: no
   filter by default (= the v9 mixed-effects workbook); --axis-se-max DEG
   applies one and then EXCLUDES NaN. What was applied is printed and stored.
 * Fig2 depth bins: a subject enters a bin only with >= --min-clusters-per-bin
   clusters (default 5). The deep bin (0.8-1.0) holds 58 clusters in total,
   median 1 per subject, 14-24% with non-physical diffusivities (CSF partial
   volume at the ventricles); under the rule it drops out instead of
   producing FA > 1 and negative RD. FDR now runs across ALL metrics x bins
   (one family), not within each metric, so a lone star means something.
 * Fig1: the '(x)' ratio to the PVS difference is printed only when the PVS
   difference is at least 0.25 pooled SD (FA printed -100x before); the two
   annotation lines no longer overlap.
 * Fig4: NearPVS added as a reference, tiers ordered by distance (m6 before
   m10), shorter y-label. With --mixed-xlsx, the stars come from the v9
   workbook, and the per-group median deltas are cross-checked against it;
   a mismatch is printed loudly and the figure falls back to its own tests.
 * Fig3 y-label: pool_frac sums the region pools over WM voxels and can exceed
   1, so it is labelled a relative pool size, not a fraction.
 * Stale end-of-run text removed (it still said orientation was skipped and
   quoted v6 p-values). Statistics file is now Figure_statistics_v9.xlsx with
   a run_info sheet.

v7 CHANGES (post frame-audit)
-----------------------------
* ORIENTATION UNSUPPRESSED. v6 set SUPPRESS_ORIENTATION = True under the
  mid-forensics belief that the block-dependent gradient convention made
  Dpar/Dperp incomparable across groups. That diagnosis is superseded: the
  convention was corrected AT THE SOURCE (step20d determination 41/41+24/24,
  step20e verified bvec fix, step14_v2 refit, step20d post-fix acceptance
  65/65 'none'), and the remaining b0-vs-T1w storage REFLECTION is handled
  inside the v8 extraction itself (header relation on u, det recorded per
  subject, step16 v3 convergence 44.1 vs 42.6 deg, p = 0.905). Orientation
  metrics are therefore group-comparable and PLOTTED BY DEFAULT.
* Inputs: v8 spreadsheet + the 2026-08-28 study xlsx (sub-127/130 excluded).

WHAT CHANGED FROM THE v4 PLOTTING SCRIPT
----------------------------------------
MECHANICAL
 * Control tiers are named by MARGIN IN MILLIMETRES (_m6, _m10), not by dilation
   iterations (_d10, _d15). The suffix regex accepts both so v4 sheets still
   parse, but the script now PRINTS which it read -- an earlier set of figures
   was generated from the v4 spreadsheet without anyone noticing.
 * Subject_QC exposes pool_frac_m*, not validWM_frac_d*.

SUBSTANTIVE
 * Fig1 annotates the RAW group difference and the tier SD alongside Cliff's
   delta. In v4 the apparent PVS -> FakeControl gradient came almost entirely
   from the control tiers having 2-4x the between-subject SD while the raw
   difference stayed flat out to 11 mm. A standardised effect size alone cannot
   show that, and it is what made the gradient look real.
 * Fig5 (new) plots raw difference against Cliff's delta across tiers. READ IT
   BEFORE Fig1.
 * Fig6 (new) tests whether the control pool is comparable across groups
   (pool_frac and relaxation counts). Both differ by group in this cohort, so m6
   is primary and m10 is a sensitivity check.
 * Fig2 stratifies by NORMALIZED DEPTH, not thickness. On a 0.80 mm grid the
   EDT thickness has a hard floor at 2*voxel = 1.597 mm which 99.1% of clusters
   sit at exactly; T1_thin was unreachable by construction and T3_thick held 33
   of 33,333 clusters. It never tested partial-volume dependence.
 * Orientation metrics (Dpar/Dperp/Dpar_minus_Dperp) are SKIPPED by default. The
   diffusion gradient sign convention differs between the two acquisition
   blocks, so orientation-dependent quantities cannot be compared across groups.
   Rotation-invariant metrics (MD, FA, AD, RD, kurtosis) are unaffected. Set
   SUPPRESS_ORIENTATION = False to plot them anyway.
 * Cluster-axis quality uses Axis_SE_deg (jackknife angular SE) rather than
   Linearity, since it measures directly what the filter is for.

WHAT IT PRODUCES
----------------
FIGURE 1 - spatial-specificity gradient, per metric. Tier on x, group as hue,
    with Cliff's delta, the raw difference and the tier SD printed beneath.
FIGURE 2 - stratification by cortical depth, per metric.
FIGURE 3 - control-pool fraction by exclusion margin and group.
FIGURE 4 - PVS minus that subject's OWN control ROIs. The key panel: a
    within-subject difference cancels head size, global diffusivity, scanner
    state and registration deformation, and it tests directly whether the
    PVS-vs-WM excess differs by group -- which Fig1 never does.
FIGURE 5 - raw difference vs standardised effect across tiers.
FIGURE 6 - control-pool comparability across groups.
"""

import os
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import kruskal, mannwhitneyu, wilcoxon
from pipeline_config import cfg

# ==========================================
# CONFIGURATION
# ==========================================
RESULTS_FILE = os.path.join(cfg("workbooks_dir"),  # v15: v9 extraction
                            "pvs_tier_metrics_v9.xlsx")
GROUPS_FILE = cfg("study_xlsx")

# The spreadsheet names the boxer cohort "Boxers"; the paper calls it RHI.
SHEET_TO_GROUP = {"HC": "HC", "PPCS": "PPCS", "Boxers": "RHI", "RHI": "RHI"}
INCLUSION_COL = "Will we used it?"   # same filter the extraction applied
OUTPUT_DIR = cfg("tier_figures_dir")

GROUPS_ORDER = ["HC", "PPCS", "RHI"]
PALETTE = {"HC": "#BDB76B", "PPCS": "#B8860B", "RHI": "#800000"}
BLACK_PAL = {g: "k" for g in ["HC", "PPCS", "RHI"]}
# v3: the partial-volume axis is THICKNESS (inscribed diameter), not volume.
# PVS are tubes -- volume grows by LENGTH, so it is nearly uncorrelated with the
# per-voxel fluid fraction. Only cross-sectional diameter sets that.
# THICKNESS STRATIFICATION IS REMOVED. On a 0.80 mm grid the EDT-based thickness
# has a hard floor at 2*voxel = 1.597 mm, which 99.1% of clusters sit at exactly;
# T1_thin was unreachable by construction and T3_thick held 33 of 33,333
# clusters. It never tested partial-volume dependence. v6 stratifies by
# NORMALIZED DEPTH (0 = cortical interface, 1 = ventricle), which is continuous,
# populated in all five bins, and is the axis along which the control pool
# actually varies with PVS burden.
BIN_COL = "Depth_Bin"
BIN_ORDER = [0, 1, 2, 3, 4]
BIN_LABEL = {0: "superficial\n(0$-$0.2)", 1: "0.2$-$0.4", 2: "mid\n(0.4$-$0.6)",
             3: "0.6$-$0.8", 4: "deep\n(0.8$-$1.0)"}

# Dpar/Dperp use a per-cluster PCA axis. For blob-like clusters that axis is
# noise, so those rows are dropped. Perfectly linear clusters give lambda2 = 0
# -> Linearity = inf -> NaN after rounding, and those are the MOST tube-like
# clusters, so NaN must be KEPT, not discarded.
MIN_LINEARITY = 2.0
# v9: set from --axis-se-max. None = no filter, which is what the v9
# mixed-effects script uses, so figures and workbook see the same clusters.
AXIS_SE_MAX = None
ORIENTATION_METRICS = {"Dpar", "Dperp", "Dpar_minus_Dperp"}
# v7: orientation metrics are group-comparable again. The block-dependent
# gradient convention was fixed at the source (steps 20d/20e/14_v2; post-fix
# acceptance 65/65 'none'), and the b0-vs-T1w storage reflection is handled
# inside the v8 extraction (header relation on u; u_rel_det in Subject_QC;
# convergence p = 0.905). The v6 suppression rationale is retained here only
# as history -- do not resurrect it without re-running the step20d audit.
SUPPRESS_ORIENTATION = False

FS_TITLE, FS_LAB, FS_TICK, FS_ANNOT = 22, 22, 16, 15

# v6 names control tiers by MARGIN IN MILLIMETRES (_m6, _m10) rather than by
# dilation iterations (_d10, _d15). Both accepted so old sheets still plot.
SUFFIX_RE = re.compile(r"_(Original|NearPVS|NearPVSfull|FakeControl_[dm]\d+)$")
FAKE_PREFIX = "FakeControl_"

# ---------- v11: the only names that may appear in figure text ----------
TIER_NAME = {"PVS": "PVS",
             "NearPVS": "Adjacent WM\n(3$-$6 mm)",
             "NearPVSfull": "Adjacent WM\n(3$-$6 mm, full shell)",
             "Fake m6": "Matched WM\n($\\geq$6 mm)",
             "Fake m10": "Matched WM\n($\\geq$10 mm)"}
TIER_NAME_1L = {k: v.replace("\n", " ") for k, v in TIER_NAME.items()}
METRIC_SHORT = {"MD": "MD", "AD": "AD", "RD": "RD", "FA": "FA",
                "MK_mf": "MK", "AK_mf": "AK", "RK_mf": "RK",
                "Dpar": "D$_\\parallel$", "Dperp": "D$_\\perp$",
                "Dpar_minus_Dperp": "D$_\\parallel$ $-$ D$_\\perp$",
                "lowb_MD": "low-b MD", "lowb_AD": "low-b AD", "lowb_RD": "low-b RD",
                "lowb_FA": "low-b FA", "lowb_ADC": "low-b ADC"}
METRIC_FULL = {"MD": "Mean diffusivity (MD)", "AD": "Axial diffusivity (AD)",
               "RD": "Radial diffusivity (RD)", "FA": "Fractional anisotropy (FA)",
               "MK_mf": "Mean kurtosis (MK)", "AK_mf": "Axial kurtosis (AK)",
               "RK_mf": "Radial kurtosis (RK)",
               "Dpar": "Diffusivity along the PVS axis (D$_\\parallel$)",
               "Dperp": "Diffusivity across the PVS axis (D$_\\perp$)",
               "Dpar_minus_Dperp": "Along-axis anisotropy (D$_\\parallel$ $-$ D$_\\perp$)",
               "lowb_MD": "Low-b mean diffusivity", "lowb_AD": "Low-b axial diffusivity",
               "lowb_RD": "Low-b radial diffusivity", "lowb_FA": "Low-b fractional anisotropy",
               "lowb_ADC": "Low-b ADC"}


def mshort(m):
    return METRIC_SHORT.get(m, m)


def mfull(m):
    return METRIC_FULL.get(m, m)


def tname(t, one_line=False):
    return (TIER_NAME_1L if one_line else TIER_NAME).get(t, t)


def margin_name(tag):
    """'m6' / 'm10' / 'Fake_m6' -> '>=6 mm' (matched-WM margin)."""
    mm = re.search(r"(\d+)", str(tag)).group(1)
    return f"$\\geq${mm} mm"


# ==========================================
# STATS
# ==========================================
def cliffs_delta(x, y):
    """Cliff's delta for x vs y (positive => x tends to exceed y). Returns (delta, p)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if x.size < 2 or y.size < 2:
        return np.nan, np.nan
    try:
        U, p = mannwhitneyu(x, y, alternative="two-sided")
    except ValueError:
        return np.nan, np.nan
    return (2.0 * U / (x.size * y.size)) - 1.0, p


def fdr(pvals):
    """Benjamini-Hochberg q-values."""
    p = np.asarray(pvals, float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    if ok.sum() == 0:
        return q
    pv = p[ok]; n = pv.size
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


def stars(q):
    if not np.isfinite(q):
        return "n/a"
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else "ns"


# ==========================================
# LOADING
# ==========================================
def load_groups(path):
    """Subject -> group, applying the SAME inclusion filter as the extraction.

    Two failure modes this guards against, both of which happened:
      * the sheet is called "Boxers", not "RHI" -- the old code printed a warning
        and silently continued with two groups;
      * without the inclusion column, subjects excluded from the study would be
        plotted, so the figures' N would not match the spreadsheet's.
    """
    m = {}
    xls = pd.ExcelFile(path)
    for sheet, g in SHEET_TO_GROUP.items():
        if sheet not in xls.sheet_names:
            continue
        df = pd.read_excel(xls, sheet_name=sheet)
        subj_col = "Subj" if "Subj" in df.columns else df.columns[0]
        if INCLUSION_COL in df.columns:
            use = df[INCLUSION_COL].astype(str).str.strip().str.lower()
            df = df[use == "yes"]
        else:
            print(f"  [!] '{INCLUSION_COL}' missing in sheet '{sheet}'; "
                  f"taking ALL rows -- N may not match the extraction")
        for v in df[subj_col].astype(str).str.strip():
            if v.lower() in ("nan", "total:", ""):
                continue
            key = v if v.startswith("sub-") else f"sub-{int(float(v))}"
            m[key] = g
    if not m:
        raise SystemExit(f"no groups loaded. Sheets present: {xls.sheet_names}\n"
                         f"Expected one of {list(SHEET_TO_GROUP)}.")
    print("  group sizes:", pd.Series(m).value_counts().to_dict())
    return m


def parse_columns(df):
    """{base_metric: {tier: column}} -- strips KNOWN suffixes, so bases containing
    underscores (MK_mf, RKt_mf) survive intact."""
    out = {}
    for c in df.columns:
        mt = SUFFIX_RE.search(c)
        if not mt:
            continue
        base, tier = c[:mt.start()], mt.group(1)
        out.setdefault(base, {})[tier] = c
    return out


def subject_level(df, col, by_bin=False):
    """Median across PVS clusters, per subject (optionally per thickness bin)."""
    keys = ["Subject", "Group"] + ([BIN_COL] if by_bin else [])
    return df.groupby(keys)[col].median().reset_index()


def apply_linearity_filter(df, metric):
    """Orientation metrics: optional axis-precision filter (v9).

    AXIS_SE_MAX None -> no filter (matches the v9 mixed-effects workbook).
    Otherwise keep Axis_SE_deg <= AXIS_SE_MAX and DROP NaN: a NaN axis SE is a
    3-voxel cluster whose axis could not be estimated at all. (v8 kept NaN,
    a rule copied from the Linearity filter where NaN meant 'perfectly linear'.)
    """
    if metric not in ORIENTATION_METRICS or AXIS_SE_MAX is None:
        return df
    if "Axis_SE_deg" not in df.columns:
        raise SystemExit("--axis-se-max given but Axis_SE_deg column missing")
    se = pd.to_numeric(df["Axis_SE_deg"], errors="coerce")
    return df[se <= AXIS_SE_MAX]


# ==========================================
# FIGURES
# ==========================================
def _fmt_axis(ax, ylabel, title, xlabel=""):
    ax.set_ylabel(ylabel, fontsize=FS_LAB, fontweight="bold")
    ax.set_xlabel(xlabel, fontsize=FS_LAB, fontweight="bold")
    ax.set_title(title, fontsize=FS_TITLE, fontweight="bold", pad=14)
    ax.tick_params(labelsize=FS_TICK)
    sns.despine(ax=ax, top=True, right=True)


def fig_tiers(long_df, metric, tiers, outdir, stats_rows):
    """Figure 1: group separation as a function of distance from the PVS."""
    d = long_df[long_df["Metric"] == metric].copy()
    if d.empty:
        return
    present = [t for t in tiers if t in d["Tier"].unique()]
    sc, unit = _scale(metric)                      # v11: 10^-3 mm^2/s
    d["Value"] = d["Value"] * sc
    fig, ax = plt.subplots(figsize=(3.0 * len(present) + 3, 7.5))
    sns.violinplot(data=d, x="Tier", y="Value", hue="Group", order=present,
                   hue_order=GROUPS_ORDER, palette=PALETTE, inner="quartile",
                   linewidth=1, cut=0, ax=ax)
    sns.stripplot(data=d, x="Tier", y="Value", hue="Group", order=present,
                  hue_order=GROUPS_ORDER, dodge=True, palette=BLACK_PAL, alpha=.45,
                  size=2.5, ax=ax, legend=False)

    # RHI vs HC delta under each tier
    y0, y1 = ax.get_ylim()
    span = y1 - y0
    ax.set_ylim(y0 - 0.34 * span, y1)          # v9: room for 4 text lines
    # A STANDARDISED effect size alone is misleading here. In v4 the apparent
    # PVS -> FakeControl gradient came almost entirely from the control tiers
    # having 2-4x the between-subject SD, while the RAW RHI-HC difference stayed
    # flat out to 11 mm. Both are annotated so the reader can see which is
    # moving.
    ann = []
    for i, t in enumerate(present):
        a = d[(d.Tier == t) & (d.Group == "RHI")]["Value"]
        b = d[(d.Tier == t) & (d.Group == "HC")]["Value"]
        dl, p = cliffs_delta(a, b)
        raw = float(np.nanmedian(a) - np.nanmedian(b))
        na, nb = a.dropna(), b.dropna()
        sd = (np.sqrt(((len(na) - 1) * na.var(ddof=1)
                       + (len(nb) - 1) * nb.var(ddof=1))
                      / max(len(na) + len(nb) - 2, 1))
              if len(na) > 1 and len(nb) > 1 else np.nan)
        ann.append((i, dl, p, t, raw, sd))
    qs = fdr([x[2] for x in ann])
    raw0 = ann[0][4] if ann and np.isfinite(ann[0][4]) and ann[0][4] != 0 else np.nan
    sd0 = ann[0][5] if ann else np.nan
    # v9: a ratio to a near-zero PVS difference is meaningless (FA printed
    # -100x). Only print it when the PVS difference is >= 0.25 pooled SD.
    ratio_ok = (np.isfinite(raw0) and np.isfinite(sd0) and sd0 > 0
                and abs(raw0) >= 0.25 * sd0)
    for (i, dl, p, t, raw, sd), q in zip(ann, qs):
        rel = raw / raw0 if ratio_ok else np.nan
        ax.text(i, y0 - 0.03 * span,
                f"RHI vs HC\n$\\Delta$={dl:+.2f} {stars(q)}",
                ha="center", va="top", fontsize=FS_ANNOT, fontweight="bold",
                color="#800000" if np.isfinite(q) and q < 0.05 else "0.35")
        ax.text(i, y0 - 0.17 * span,
                f"difference {raw:+.2f}" + (f"  ({rel:.2f}$\\times$)"
                                            if np.isfinite(rel) else "")
                + f"\nSD {sd:.2f}",
                ha="center", va="top", fontsize=FS_ANNOT - 2, color="0.25")
        stats_rows.append(["Fig1_tier", metric, t, "RHI vs HC", dl, p, q,
                           raw / sc, sd / sc, rel])

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:3], labels[:3], title="Group", fontsize=FS_TICK,
              title_fontsize=FS_TICK, frameon=True, loc="upper right")
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([tname(t) for t in present])
    _fmt_axis(ax, f"{mshort(metric)}{unit}",
              f"{mfull(metric)}: RHI vs HC by tissue tier", xlabel="Tissue tier")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"Fig1_tiers_{metric}.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig)


def fdr_family_of(metric, mode):
    """v10: same rule as the mixed-effects v10 script."""
    if mode == "single":
        return "all"
    return "lowb" if str(metric).startswith("lowb_") else "main"


def depthbin_subject_table(df, col, min_clusters):
    """Subject x depth-bin medians, keeping a (subject, bin) cell only if it
    holds >= min_clusters clusters. The deep bin has a median of 1 cluster per
    subject and many non-physical fits; a single cluster's value is not a
    subject-level estimate."""
    keys = ["Subject", "Group", BIN_COL]
    g = df.groupby(keys)[col]
    sl = g.median().to_frame(col).join(g.size().rename("n_clusters")).reset_index()
    return sl[(sl["n_clusters"] >= min_clusters) & sl[BIN_COL].isin(BIN_ORDER)]


def depthbin_tests(sl, col, min_n):
    """RHI vs HC per bin -> {bin: (delta, p)} (NaN where n is insufficient)."""
    out = {}
    for b in BIN_ORDER:
        a = sl[(sl[BIN_COL] == b) & (sl.Group == "RHI")][col]
        c = sl[(sl[BIN_COL] == b) & (sl.Group == "HC")][col]
        out[b] = (np.nan, np.nan) if (a.size < min_n or c.size < min_n) \
            else cliffs_delta(a, c)
    return out


def fig_thicknessbins(sl, metric, col, tests, qmap, outdir, stats_rows):
    """Figure 2: RHI vs HC by normalised depth (v9: subject-bin cells need
    >= --min-clusters-per-bin clusters; q-values are FDR across ALL metrics x
    bins, passed in as qmap)."""
    shown = [b for b in BIN_ORDER if (sl[BIN_COL] == b).any()]
    if len(shown) < 2:
        print(f"  [!] {metric}: fewer than 2 populated depth bins, Fig2 skipped")
        return
    sc, unit = _scale(metric)                      # v11: 10^-3 mm^2/s
    sl = sl.copy(); sl[col] = sl[col] * sc
    fig, ax = plt.subplots(figsize=(2.6 * len(shown) + 2.5, 7.5))
    sns.violinplot(data=sl, x=BIN_COL, y=col, hue="Group", order=shown,
                   hue_order=GROUPS_ORDER, palette=PALETTE, inner="quartile",
                   linewidth=1, cut=0, ax=ax)
    sns.stripplot(data=sl, x=BIN_COL, y=col, hue="Group", order=shown,
                  hue_order=GROUPS_ORDER, dodge=True, palette=BLACK_PAL, alpha=.45,
                  size=3, ax=ax, legend=False)
    y0, y1 = ax.get_ylim(); span = y1 - y0
    ax.set_ylim(y0 - 0.26 * span, y1)
    for i, b in enumerate(shown):
        dl, p = tests[b]
        q = qmap.get((metric, b), np.nan)
        txt = "insufficient n" if not np.isfinite(dl) else \
              f"RHI vs HC\n$\\Delta$={dl:+.2f} {stars(q)}"
        ax.text(i, y0 - 0.06 * span, txt, ha="center", va="top",
                fontsize=FS_ANNOT, fontweight="bold",
                color="#800000" if np.isfinite(q) and q < 0.05 else "0.35")
        stats_rows.append(["Fig2_depthbin", metric, b, "RHI vs HC", dl, p, q,
                           np.nan, np.nan, np.nan])
    ax.set_xticks(range(len(shown)))
    ax.set_xticklabels([BIN_LABEL[b] for b in shown])
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:3], labels[:3], title="Group", fontsize=FS_TICK,
              title_fontsize=FS_TICK, frameon=True, loc="upper right")
    _fmt_axis(ax, f"{mshort(metric)} within PVS{unit}",
              f"{mfull(metric)}: RHI vs HC within PVS by cortical depth",
              xlabel="Normalised depth (0 = cortex, 1 = ventricle)")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"Fig2_depthbins_{metric}.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig)


def fig_raw_vs_standardised(long_df, metric, tiers, outdir):
    """Figure 5: the diagnostic that overturned the v4 gradient.

    Left axis: RAW RHI-HC difference per tier, in the metric's own units.
    Right axis: the same contrast as Cliff's delta.

    If the standardised effect decays while the raw difference does not, the
    'spatial gradient' is a denominator artefact -- the control tiers are simply
    noisier -- and no claim of PVS specificity can rest on it.
    """
    d = long_df[long_df["Metric"] == metric]
    present = [t for t in tiers if t in d["Tier"].unique()]
    if len(present) < 2:
        return
    raws, dls, sds = [], [], []
    for t in present:
        a = d[(d.Tier == t) & (d.Group == "RHI")]["Value"].dropna()
        b = d[(d.Tier == t) & (d.Group == "HC")]["Value"].dropna()
        if len(a) < 2 or len(b) < 2:
            raws.append(np.nan); dls.append(np.nan); sds.append(np.nan); continue
        raws.append(float(np.median(a) - np.median(b)))
        dls.append(cliffs_delta(a, b)[0])
        sds.append(float(np.sqrt(((len(a) - 1) * a.var(ddof=1)
                                  + (len(b) - 1) * b.var(ddof=1))
                                 / (len(a) + len(b) - 2))))
    x = np.arange(len(present))
    sc, unit = _scale(metric)                      # v11: 10^-3 mm^2/s
    raws = [r * sc for r in raws]; sds = [v * sc for v in sds]
    fig, ax = plt.subplots(figsize=(3.0 * len(present) + 3, 7))
    ax.plot(x, raws, "o-", color="#800000", lw=2.5, ms=10, label="RHI $-$ HC difference")
    ax.axhline(0, color="0.6", lw=1, ls=":")
    ax.set_ylabel(f"RHI $-$ HC difference, {mshort(metric)}{unit}", fontsize=FS_LAB,
                  fontweight="bold", color="#800000")
    ax.set_xlabel("Tissue tier", fontsize=FS_LAB, fontweight="bold")
    ax.tick_params(axis="y", labelcolor="#800000", labelsize=FS_TICK)
    ax.set_xticks(x)
    ax.set_xticklabels([tname(t) for t in present], fontsize=FS_TICK)
    ax2 = ax.twinx()
    ax2.plot(x, dls, "s--", color="#1f4e79", lw=2.5, ms=10,
             label="effect size (Cliff's $\\Delta$)")
    ax2.set_ylabel("Effect size (Cliff's $\\Delta$)", fontsize=FS_LAB,
                   fontweight="bold", color="#1f4e79")
    ax2.tick_params(axis="y", labelcolor="#1f4e79", labelsize=FS_TICK)
    for xi, sd in zip(x, sds):
        ax2.annotate(f"SD\n{sd:.2f}", (xi, ax2.get_ylim()[0]),
                     ha="center", va="bottom", fontsize=FS_ANNOT - 2,
                     color="0.35")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=FS_TICK, loc="best", frameon=True)
    ax.set_title(f"{mfull(metric)}: difference and effect size by tissue tier",
                 fontsize=FS_TITLE, fontweight="bold")
    sns.despine(right=False)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"Fig5_raw_vs_d_{metric}.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig)


def fig_match_quality(qc, outdir, stats_rows):
    """Figure 6: is the matched control tissue comparable across groups?

    v11: relaxations and failures are shown as PERCENTAGES of each subject's
    clusters (Subject_QC n_clusters); one row per exclusion margin; panels
    whose values are all zero are dropped and named in the footnote.
    """
    if "Group" not in qc.columns:
        return
    margins = sorted({re.search(r"m(\d+)$", c).group(1) for c in qc.columns
                      if re.match(r"(pool_frac|fakeFail|relax_region|relax_any)_m\d+$", c)},
                     key=int)
    if not margins:
        return
    n_cl = qc["n_clusters"] if "n_clusters" in qc.columns else None
    panels = [("pool_frac", "Eligible white matter", False),
              ("relax_region", "Depth-relaxed matches", True),
              ("relax_any", "Matched in another region", True),
              ("fakeFail", "Unmatched clusters", True)]
    ylab = {False: "Eligible WM / WM volume", True: "Clusters (%)"}
    rows_data, omitted = [], []
    for mm in margins:
        row = []
        for key, title, as_pct in panels:
            c = f"{key}_m{mm}"
            if c not in qc.columns or (as_pct and n_cl is None):
                row.append(None); continue
            v = qc[c].astype(float)
            if as_pct:
                v = 100.0 * v / n_cl.replace(0, np.nan)
            if np.nanmax(np.abs(v.values)) == 0:
                omitted.append(f"{title.lower()} at {margin_name(mm)}")
                row.append(("zero", title)); continue
            row.append((c, title, v, as_pct))
        rows_data.append((mm, row))
    ncol = len(panels)                       # fixed column per panel type
    fig, axes = plt.subplots(len(rows_data), ncol, figsize=(5.2 * ncol, 5.8 * len(rows_data)),
                             squeeze=False)
    for ri, (mm, row) in enumerate(rows_data):
        for ci, cell in enumerate(row):
            ax = axes[ri, ci]
            if cell is None:
                ax.axis("off"); continue
            if cell[0] == "zero":
                ax.axis("off")
                ax.text(0.5, 0.5, f"{cell[1]}\nmargin {margin_name(mm)}\n\nall zero:\nnot needed",
                        transform=ax.transAxes, ha="center", va="center",
                        fontsize=FS_ANNOT, color="0.4", fontweight="bold")
                continue
            c, title, v, as_pct = cell
            dd = pd.DataFrame({"Group": qc["Group"], "v": v.values})
            sns.violinplot(data=dd, x="Group", y="v", hue="Group", order=GROUPS_ORDER,
                           hue_order=GROUPS_ORDER, palette=PALETTE, legend=False,
                           inner="quartile", cut=0, ax=ax)
            sns.stripplot(data=dd, x="Group", y="v", order=GROUPS_ORDER,
                          color="k", alpha=.5, size=3.5, ax=ax)
            gs = [dd.loc[dd.Group == g, "v"].dropna() for g in GROUPS_ORDER]
            gs = [g for g in gs if len(g) > 2]
            p = np.nan
            if len(gs) > 1 and np.ptp(np.concatenate(gs)) > 0:
                p = float(kruskal(*gs).pvalue)
            ptxt = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"       # v16-3
            ax.set_title(f"{title}\nmargin {margin_name(mm)},  {ptxt}",
                         fontsize=FS_ANNOT + 1, fontweight="bold",
                         color="#800000" if np.isfinite(p) and p < 0.05 else "0.25")
            ax.set_xlabel("")
            ax.set_ylabel(ylab[as_pct], fontsize=FS_LAB - 2, fontweight="bold")
            ax.tick_params(labelsize=FS_TICK)
            sns.despine(ax=ax)
            stats_rows.append(["Fig6_match", c, "", "KW across groups (per cent of clusters)"
                               if "pool" not in c else "KW across groups",
                               np.nan, p, np.nan, np.nan, np.nan, np.nan])
    fig.suptitle("Matched white matter: comparability across groups",
                 fontsize=FS_TITLE, fontweight="bold")
    foot = "Kruskal-Wallis p-values are uncorrected."
    if omitted:
        foot += "  Not shown (all values zero): " + "; ".join(omitted) + "."
    fig.text(0.01, 0.005, foot, fontsize=FS_ANNOT - 1, color="0.3", ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.03, 1, 0.96))
    fig.savefig(os.path.join(outdir, "Fig6_match_quality.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig)


def fig_wm_availability(qc, outdir, stats_rows):
    """Figure 3: is the exclusion margin a group-dependent confound?"""
    fracs = [c for c in qc.columns if c.startswith("pool_frac_m")]
    if not fracs:
        return
    m = qc.melt(id_vars=["Subject", "Group"], value_vars=fracs,
                var_name="Margin", value_name="ValidWM")
    m["Margin"] = m["Margin"].str.replace("pool_frac_", "", regex=False)
    morder = sorted(m["Margin"].unique(), key=lambda t: int(re.search(r"(\d+)", t).group(1)))
    fig, ax = plt.subplots(figsize=(4 + 3.0 * len(fracs), 7))
    sns.violinplot(data=m, x="Margin", y="ValidWM", hue="Group", order=morder,
                   hue_order=GROUPS_ORDER, palette=PALETTE, inner="quartile",
                   linewidth=1, cut=0, ax=ax)
    sns.stripplot(data=m, x="Margin", y="ValidWM", hue="Group", order=morder,
                  hue_order=GROUPS_ORDER, dodge=True, palette=BLACK_PAL, alpha=.5,
                  size=3, ax=ax, legend=False)
    y0, y1 = ax.get_ylim(); span = y1 - y0
    ax.set_ylim(y0 - 0.26 * span, y1)
    for i, mar in enumerate(morder):
        groups = [m[(m.Margin == mar) & (m.Group == g)]["ValidWM"].dropna()
                  for g in GROUPS_ORDER]
        groups = [g for g in groups if g.size > 1]
        p = kruskal(*groups)[1] if len(groups) >= 2 else np.nan
        ax.text(i, y0 - 0.06 * span,
                f"Kruskal-Wallis\np={p:.3f}" if np.isfinite(p) else "n/a",
                ha="center", va="top", fontsize=FS_ANNOT, fontweight="bold",
                color="#800000" if np.isfinite(p) and p < 0.05 else "0.35")
        stats_rows.append(["Fig3_pool_frac", "pool_frac", mar,
                           "Kruskal across groups", np.nan, p, np.nan, np.nan, np.nan, np.nan])
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:3], labels[:3], title="Group", fontsize=FS_TICK,
              title_fontsize=FS_TICK, frameon=False, loc="upper left",
              bbox_to_anchor=(1.01, 1.0))
    ax.set_xticks(range(len(morder)))
    ax.set_xticklabels([margin_name(x) for x in morder])
    _fmt_axis(ax, "Eligible white matter\n(relative to WM volume)",
              "Eligible white matter by exclusion margin",
              xlabel="Minimum distance from any PVS")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "Fig3_validWM_availability.png"), dpi=300,
                bbox_inches="tight")
    plt.close(fig)


# ---- v9 helpers shared by Fig4 and the paired supplementary figures --------
DIFFUSIVITY = {"MD", "AD", "RD", "Dpar", "Dperp", "Dpar_minus_Dperp",
               "lowb_AD", "lowb_MD", "lowb_RD", "lowb_ADC"}
PRETTY = {"lowb_MD": "Low-b mean diffusivity", "lowb_AD": "Low-b axial diffusivity",
          "lowb_RD": "Low-b radial diffusivity", "lowb_FA": "Low-b fractional anisotropy",
          "lowb_ADC": "Low-b ADC",
          "Dpar_minus_Dperp": "D$_\\parallel$ $-$ D$_\\perp$",
          "Dpar": "D$_\\parallel$ (along PVS axis)",
          "Dperp": "D$_\\perp$ (across PVS axis)",
          "MK_mf": "Mean kurtosis", "AK_mf": "Axial kurtosis",
          "RK_mf": "Radial kurtosis", "FA": "Fractional anisotropy",
          "MD": "Mean diffusivity", "AD": "Axial diffusivity",
          "RD": "Radial diffusivity"}
REF_PRETTY = {"NearPVS": "adjacent white matter (3$-$6 mm, size-matched)",
              "Fake_m6": "matched white matter ($\\geq$6 mm)",
              "Fake_m10": "matched white matter ($\\geq$10 mm)"}
TIER_TICK = TIER_NAME                     # v11: one naming table for every figure
CROSSCHECK = []          # (metric, ref, group, here, workbook) mismatches


def _scale(metric):
    """Diffusivities are plotted in 10^-3 mm^2/s."""
    return (1e3, " ($\\times10^{-3}$ mm$^2$/s)") if metric in DIFFUSIVITY else (1.0, "")


def paired_table(near, fake, metric, cols_near, cols_fake, ref):
    """Per-subject PVS, reference and delta = PVS - reference.
    ref: 'NearPVS' or 'Fake_m6' / 'Fake_m10' (the mixed-effects workbook's Ref_tier names)."""
    if "Original" not in cols_near.get(metric, {}):
        return None
    n_f = apply_linearity_filter(near, metric)
    pvs = subject_level(n_f, cols_near[metric]["Original"])
    pvs.columns = ["Subject", "Group", "PVS"]
    if ref == "NearPVS":
        if "NearPVS" not in cols_near[metric]:
            return None
        ctl = subject_level(n_f, cols_near[metric]["NearPVS"])
    else:
        col = cols_fake.get(metric, {}).get("FakeControl_" + ref.split("_")[-1])
        if col is None:
            return None
        ctl = subject_level(apply_linearity_filter(fake, metric), col)
    ctl.columns = ["Subject", "Group", "Ref"]
    m = pvs.merge(ctl, on=["Subject", "Group"], how="inner").dropna()
    m["Delta"] = m["PVS"] - m["Ref"]
    return m


def workbook_row(W, metric, ref):
    if W is None:
        return None
    r = W[(W.Metric == metric) & (W.Ref_tier == ref)]
    return None if r.empty else r.iloc[0]


def crosscheck(m, row, metric, ref):
    """Per-group median delta here vs the workbook. True if they agree."""
    ok = True
    for g in GROUPS_ORDER:
        here = float(m.loc[m.Group == g, "Delta"].median())
        wb = float(row[f"delta_{g}"])
        if np.isfinite(here) and np.isfinite(wb) and \
                not np.isclose(here, wb, rtol=1e-6, atol=1e-12):
            CROSSCHECK.append((metric, ref, g, here, wb)); ok = False
    return ok


def fig_within_subject_contrast(near, fake, metric, cols_near, cols_fake,
                                W, outdir, stats_rows):
    """Figure 4: PVS minus each subject's OWN reference tissue (v9).

    References: adjacent WM (NearPVS, size-matched 3-6 mm) and the matched
    distant controls, in order of distance. For the orientation metrics
    NearPVS is the fair reference: the PVS axis transplanted into distant
    tissue lands in different fibre geometry. With a mixed-effects workbook the RHI-vs-HC
    stars are the workbook's q_RHI_vs_HC (FDR across metrics within a
    reference tier), after checking that the per-group median deltas match.
    """
    fakes = sorted((t for t in cols_fake.get(metric, {}) if t.startswith(FAKE_PREFIX)),
                   key=lambda t: int(re.search(r"(\d+)", t).group(1)))
    refs = ["NearPVS"] + ["Fake_" + t.split("_")[-1] for t in fakes]
    recs, meta = [], []
    for ref in refs:
        m = paired_table(near, fake, metric, cols_near, cols_fake, ref)
        if m is None or m.empty:
            continue
        lab = ("PVS $-$ Adjacent WM\n(3$-$6 mm)" if ref == "NearPVS"
               else f"PVS $-$ Matched WM\n({margin_name(ref)})")
        m["Tier"] = lab
        recs.append(m); meta.append((ref, lab, m))
    if not recs:
        return
    cd = pd.concat(recs, ignore_index=True)
    order = [lab for _, lab, _ in meta]
    sc, unit = _scale(metric)
    cd["Plot"] = cd["Delta"] * sc

    fig, ax = plt.subplots(figsize=(3.6 * len(order) + 3, 7.5))
    sns.violinplot(data=cd, x="Tier", y="Plot", hue="Group", order=order,
                   hue_order=GROUPS_ORDER, palette=PALETTE, inner="quartile",
                   linewidth=1, cut=0, ax=ax)
    sns.stripplot(data=cd, x="Tier", y="Plot", hue="Group", order=order,
                  hue_order=GROUPS_ORDER, dodge=True, palette=BLACK_PAL,
                  alpha=.5, size=3, ax=ax, legend=False)
    ax.axhline(0, color="k", lw=1.2, ls="--", zorder=0)
    y0, y1 = ax.get_ylim(); span = y1 - y0
    ax.set_ylim(y0 - 0.20 * span, y1)

    own = []
    for i, (ref, lab, m) in enumerate(meta):
        a = m.loc[m.Group == "RHI", "Delta"]; b = m.loc[m.Group == "HC", "Delta"]
        dl, p = cliffs_delta(a, b)
        own.append((i, ref, lab, m, dl, p))
    own_q = fdr([x[5] for x in own])
    for (i, ref, lab, m, dl, p), q_own in zip(own, own_q):
        row = workbook_row(W, metric, ref)
        use_wb = row is not None and crosscheck(m, row, metric, ref)
        q = float(row["q_RHI_vs_HC"]) if use_wb else q_own
        src = "mixed-effects workbook" if use_wb else "figure-level FDR"
        ax.text(i, y0 - 0.07 * span,
                f"RHI vs HC\n$\\Delta$={dl:+.2f} {stars(q)}",
                ha="center", va="top", fontsize=FS_ANNOT, fontweight="bold",
                color="#800000" if np.isfinite(q) and q < 0.05 else "0.35")
        stats_rows.append(["Fig4_within_subject", metric, lab.replace("\n", " "), f"RHI vs HC ({src})",
                           dl, p, q, np.nan, np.nan, np.nan])
        for g in GROUPS_ORDER:
            v = m.loc[m.Group == g, "Delta"]
            qg = float(row[f"q_{g}_vs_0"]) if use_wb else np.nan
            try:
                p1 = float(wilcoxon(v).pvalue) if v.size >= 3 and v.abs().sum() else np.nan
            except ValueError:
                p1 = np.nan
            stats_rows.append(["Fig4_within_subject", metric, f"{lab.replace(chr(10), ' ')} | {g}",
                               f"delta vs 0, Wilcoxon ({src})", float(v.median()),
                               p1, qg, np.nan, np.nan, np.nan])

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:3], labels[:3], title="Group", fontsize=FS_TICK,
              title_fontsize=FS_TICK, frameon=True, loc="upper right")
    _fmt_axis(ax, f"PVS $-$ reference tissue, {mshort(metric)}{unit}",
              f"{mfull(metric)}\nPVS minus each subject's own reference tissue",
              xlabel="Within-subject contrast")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"Fig4_within_subject_{metric}.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)


def _save_both(fig, outdir, stem):
    fig.savefig(os.path.join(outdir, stem + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(outdir, stem + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def _q_line(row, prefix):
    parts = []
    for g in GROUPS_ORDER:
        q = float(row[f"q_{g}_vs_0"]) if row is not None else np.nan
        parts.append(f"{g} {stars(q)}")
    return prefix + "   ".join(parts)


def _tier_panel(ax, long_df, met, tiers, W, letter, stats_rows, zero_line=False):
    """One tier panel (violins by group, q-lines from the workbook)."""
    d = long_df[(long_df.Metric == met) & long_df.Tier.isin(tiers)].copy()
    present = [t for t in tiers if t in set(d.Tier)]
    sc, unit = _scale(met)
    d["Plot"] = d["Value"] * sc
    sns.violinplot(data=d, x="Tier", y="Plot", hue="Group", order=present,
                   hue_order=GROUPS_ORDER, palette=PALETTE, inner="quartile",
                   linewidth=1, cut=0, density_norm="width", ax=ax, legend=False)
    sns.stripplot(data=d, x="Tier", y="Plot", hue="Group", order=present,
                  hue_order=GROUPS_ORDER, dodge=True, palette=BLACK_PAL,
                  alpha=.4, size=2.2, ax=ax, legend=False)
    if zero_line:
        ax.axhline(0, color="k", lw=1, ls="--", zorder=0)
    y0, y1 = ax.get_ylim(); span = y1 - y0
    ax.set_ylim(y0 - 0.30 * span, y1)
    r_near = workbook_row(W, met, "NearPVS")
    r_far = workbook_row(W, met, "Fake_m6")
    ax.text(0.5, 0.115, _q_line(r_near, "PVS vs adjacent WM:   "),
            transform=ax.transAxes, ha="center", va="center",
            fontsize=12, fontweight="bold", color="0.15")
    ax.text(0.5, 0.045, _q_line(r_far, "PVS vs matched WM:   "),
            transform=ax.transAxes, ha="center", va="center",
            fontsize=11, color="0.35")
    for ref, r in (("NearPVS", r_near), ("Fake_m6", r_far)):
        for g in GROUPS_ORDER:
            stats_rows.append(["Main_Fig2", met, f"PVS - {ref} | {g}",
                               "delta vs 0 (mixed-effects workbook)",
                               np.nan if r is None else float(r[f"delta_{g}"]),
                               np.nan if r is None else float(r[f"p_{g}_vs_0"]),
                               np.nan if r is None else float(r[f"q_{g}_vs_0"]),
                               np.nan, np.nan, np.nan])
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([TIER_TICK.get(t, t) for t in present], fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel(f"{mshort(met)}{unit}", fontsize=15, fontweight="bold")
    ax.tick_params(axis="y", labelsize=12)
    ax.set_title(letter, loc="left", fontsize=22, fontweight="bold")
    sns.despine(ax=ax)


def _alignment_panel(ax, near, letter, stats_rows):
    """B2: alignment of the PVS axis with the fibre direction (cos^2 theta) by
    precision of the axis estimate; per-subject medians, violins by group."""
    if not {"cos2_theta", "Axis_SE_deg", "Group"} <= set(near.columns):
        ax.axis("off"); return
    d = near[["Subject", "Group", "cos2_theta", "Axis_SE_deg"]].dropna().copy()
    edges = [0, 5, 10, 20, 180]
    labels = ["$\\leq$5°", "5$-$10°", "10$-$20°", ">20°"]
    d["bin"] = pd.cut(d["Axis_SE_deg"], edges, labels=labels, include_lowest=True)
    # per-subject MEAN: the random-direction reference (1/3) is the expected
    # mean of cos^2; the median of cos^2 for random directions is 0.25
    sl = d.groupby(["Subject", "Group", "bin"], observed=True)["cos2_theta"].mean().reset_index()
    cnt = d.groupby("bin", observed=True).size()
    sns.violinplot(data=sl, x="bin", y="cos2_theta", hue="Group", order=labels,
                   hue_order=GROUPS_ORDER, palette=PALETTE, inner="quartile",
                   linewidth=1, cut=0, ax=ax, legend=False)
    sns.stripplot(data=sl, x="bin", y="cos2_theta", hue="Group", order=labels,
                  hue_order=GROUPS_ORDER, dodge=True, palette=BLACK_PAL,
                  alpha=.4, size=2.2, ax=ax, legend=False)
    ax.axhline(1 / 3, color="k", lw=1.2, ls="--", zorder=0)
    ax.text(0.99, 1 / 3, " random\n direction", transform=ax.get_yaxis_transform(),
            ha="left", va="center", fontsize=11, color="0.3")
    for i, lab in enumerate(labels):
        m = float(d.loc[d["bin"] == lab, "cos2_theta"].mean()) if lab in cnt.index else np.nan
        n = int(cnt.get(lab, 0))
        stats_rows.append(["Main_Fig2_alignment", "cos2_theta", lab,
                           "cohort mean over clusters", m, np.nan, np.nan, n, np.nan, np.nan])
    ax.set_xlabel("Precision of the PVS axis (angular SE)", fontsize=14, fontweight="bold")
    ax.set_ylabel("Alignment with fibre direction\n(cos$^2\\theta$, PVS axis vs V1)",
                  fontsize=15, fontweight="bold")
    ax.tick_params(labelsize=12)
    ax.set_title(letter, loc="left", fontsize=22, fontweight="bold")
    sns.despine(ax=ax)


def _dodge_offsets(n_hue=3, width=0.8):
    """x offsets seaborn uses for dodged violins."""
    each = width / n_hue
    return [(-width / 2 + each / 2) + i * each for i in range(n_hue)]


def _tier_panel_stars(ax, long_df, met, tiers, refs, W, letter, stats_rows,
                      zero_line=False, ylabel_full=True):
    """v13: violins by group for the tiers given; no text, only stars on a
    light line joining each group's PVS violin to its reference violin.
    refs: list of (tier name in long_df, workbook Ref_tier name)."""
    d = long_df[(long_df.Metric == met) & long_df.Tier.isin(tiers)].copy()
    present = [t for t in tiers if t in set(d.Tier)]
    sc, unit = _scale(met)
    d["Plot"] = d["Value"] * sc
    sns.violinplot(data=d, x="Tier", y="Plot", hue="Group", order=present,
                   hue_order=GROUPS_ORDER, palette=PALETTE, inner="quartile",
                   linewidth=1, cut=0, density_norm="width", ax=ax, legend=False)
    sns.stripplot(data=d, x="Tier", y="Plot", hue="Group", order=present,
                  hue_order=GROUPS_ORDER, dodge=True, palette=BLACK_PAL,
                  alpha=.4, size=2.4, ax=ax, legend=False)
    if zero_line:
        ax.axhline(0, color="k", lw=1, ls="--", zorder=0)
    lo, hi = float(d.Plot.min()), float(d.Plot.max())
    rng = hi - lo
    offs = _dodge_offsets()
    # stars on light lines: one level per (reference, group)
    step = 0.07 * rng
    level = hi + 0.05 * rng
    for ref_tier, ref_key in refs:
        if ref_tier not in present:
            continue
        j = present.index(ref_tier)
        row = workbook_row(W, met, ref_key)
        for gi, g in enumerate(GROUPS_ORDER):
            q = float(row[f"q_{g}_vs_0"]) if row is not None else np.nan
            x0, x1 = 0 + offs[gi], j + offs[gi]
            ax.plot([x0, x1], [level, level], color=PALETTE[g], lw=1.4, alpha=0.9,
                    solid_capstyle="round", zorder=1)
            ax.plot([x0, x0], [level - 0.012 * rng, level], color=PALETTE[g], lw=1.4)
            ax.plot([x1, x1], [level - 0.012 * rng, level], color=PALETTE[g], lw=1.4)
            st = stars(q) or "ns"
            ax.text((x0 + x1) / 2, level + 0.004 * rng, st, ha="center", va="bottom",
                    fontsize=14 if st != "ns" else 10, fontweight="bold",
                    color=PALETTE[g] if st != "ns" else "0.45")
            stats_rows.append(["Main_Fig2", met, f"PVS - {ref_key} | {g}",
                               "delta vs 0 (mixed-effects workbook)",
                               np.nan if row is None else float(row[f"delta_{g}"]),
                               np.nan if row is None else float(row[f"p_{g}_vs_0"]),
                               q, np.nan, np.nan, np.nan])
            level += step
    ax.set_ylim(lo - 0.03 * rng, level + 0.02 * rng)          # tight to the data
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([TIER_TICK.get(t, t) for t in present], fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel(f"{mshort(met)}{unit}", fontsize=15, fontweight="bold")
    ax.tick_params(axis="y", labelsize=12)
    ax.set_title(letter, loc="left", fontsize=22, fontweight="bold")
    sns.despine(ax=ax)


def _measured_theta_tier(df, col, se_max=5.0):
    """Median angle in column `col` over clusters with Axis_SE_deg <= se_max."""
    if df is None or col not in df.columns or "Axis_SE_deg" not in df.columns:
        return np.nan, 0
    t = df.loc[df["Axis_SE_deg"] <= se_max, col].dropna()
    return (float(t.median()), int(t.size)) if t.size else (np.nan, 0)


def _measured_theta(near, se_max=5.0):
    """Median angle (deg) between the PVS axis and the principal fibre direction
    at the PVS, over clusters with a well-defined axis (Axis_SE_deg <= se_max)."""
    if near is None or not {"Theta_deg", "Axis_SE_deg"} <= set(near.columns):
        return np.nan, 0
    t = near.loc[near["Axis_SE_deg"] <= se_max, "Theta_deg"].dropna()
    return (float(t.median()), int(t.size)) if t.size else (np.nan, 0)


def _bundle(ax, cx, cy, angle_deg, half_len, half_width, n=7, color="#a89a6c", lw=1.8,
            z=1, clip=None):
    """Parallel fibre lines through (cx, cy) along angle_deg."""
    ca, sa = np.cos(np.radians(angle_deg)), np.sin(np.radians(angle_deg))
    for k in np.linspace(-half_width, half_width, n):
        x0, y0 = cx - half_len * ca - k * sa, cy - half_len * sa + k * ca
        x1, y1 = cx + half_len * ca - k * sa, cy + half_len * sa + k * ca
        ln, = ax.plot([x0, x1], [y0, y1], color=color, lw=lw, zorder=z,
                      solid_capstyle="round")
        if clip is not None:
            ln.set_clip_path(clip)


def _schematic_panel(ax, letter, near=None, fake=None):
    """E: what is measured.
    Left  - the PVS (axis u) and the white-matter bundle running past it at the
            measured angle theta; the adjacent band (3-6 mm) lies in the same bundle.
    Right - zoomed out: the small PVS and, >= 6 mm away, matched white matter with
            its own bundle (no preferred angle to u: D|| ~ D-perp)."""
    from matplotlib.patches import Ellipse, Circle, Arc, Annulus
    th_pvs, n_th = _measured_theta(near)
    th_adj, _ = _measured_theta_tier(near, "NearPVS_Theta_deg")
    th_far, _ = _measured_theta_tier(fake, "FakeControl_m6_Theta_deg")
    v9 = np.isfinite(th_adj) and np.isfinite(th_far)
    th = th_adj if v9 else th_pvs                      # angle drawn on the left
    th_draw = th if np.isfinite(th) else 40.0
    ax.set_xlim(0, 11); ax.set_ylim(0, 7.0); ax.set_aspect("equal"); ax.axis("off")
    navy, maroon, tan, band = "#2b4c7e", "#800000", "#a89a6c", "#f3ecd9"

    # ---------------- left: PVS and the bundle running past it ----------------
    cx, cy, u_ang = 2.8, 3.9, 20.0
    r_in, r_out = 0.95, 1.95                            # adjacent WM: 3-6 mm from the PVS
    # v16-2: only the adjacent-WM ring is shaded and carries the fibre lines
    adj = Annulus((cx, cy), r_out, r_out - r_in, facecolor=band, edgecolor="none", zorder=0)
    ax.add_patch(adj)
    _bundle(ax, cx, cy, u_ang + th_draw, 2.6, 1.9, n=9, color=tan, z=1, clip=adj)
    ring_out = Circle((cx, cy), r_out, facecolor="none", edgecolor="0.45", lw=1.1, ls=(0, (4, 3)), zorder=2)
    ring_in = Circle((cx, cy), r_in, facecolor="none", edgecolor="0.45", lw=1.1, ls=(0, (4, 3)), zorder=2)
    ax.add_patch(ring_out); ax.add_patch(ring_in)
    pvs = Ellipse((cx, cy), 1.75, 0.38, angle=u_ang, facecolor="#8fb3d9", edgecolor=navy,
                  lw=1.6, zorder=4)
    ax.add_patch(pvs)
    ua, ub = np.cos(np.radians(u_ang)), np.sin(np.radians(u_ang))
    ax.annotate("", xy=(cx + 1.6 * ua, cy + 1.6 * ub), xytext=(cx - 1.6 * ua, cy - 1.6 * ub),
                arrowprops=dict(arrowstyle="<->", color=navy, lw=1.8), zorder=5)
    fa, fb = np.cos(np.radians(u_ang + th_draw)), np.sin(np.radians(u_ang + th_draw))
    ax.annotate("", xy=(cx + 1.6 * fa, cy + 1.6 * fb), xytext=(cx, cy),
                arrowprops=dict(arrowstyle="-|>", color=maroon, lw=1.6), zorder=5)
    ax.add_patch(Arc((cx, cy), 1.5, 1.5, angle=0, theta1=u_ang, theta2=u_ang + th_draw,
                     color=maroon, lw=1.6, zorder=5))
    mid = np.radians(u_ang + th_draw / 2)
    ax.text(cx + 0.98 * np.cos(mid), cy + 0.98 * np.sin(mid), r"$\theta$", fontsize=15,
            color=maroon, fontweight="bold", ha="center", va="center", zorder=6)
    ax.text(cx + 1.65 * ua + 0.05, cy + 1.65 * ub - 0.05, "$u$", fontsize=15, color=navy,
            fontweight="bold", ha="left", va="center")
    ax.text(cx + 1.7 * fa, cy + 1.7 * fb + 0.08, "fibre", fontsize=12, color=maroon,
            fontweight="bold", ha="center", va="bottom")
    ax.text(cx - 2.35, 6.75, "At the PVS", fontsize=14, fontweight="bold", ha="left", va="top")
    if v9:
        th_txt = rf"$\theta$ = {th_adj:.0f}° (adjacent WM, measured)"
    elif np.isfinite(th):
        th_txt = rf"$\theta$ = {th:.0f}° (at the PVS, measured)"
    else:
        th_txt = r"$\theta$ (not available)"
    ax.text(cx, 1.15, th_txt, fontsize=12.5, color=maroon, fontweight="bold", ha="center")
    ax.text(cx, 0.72, "shaded ring: adjacent WM (3$-$6 mm), same bundle", fontsize=11,
            color="0.35", ha="center")

    # ---------------- right: zoomed out, matched WM >= 6 mm away ----------------
    px, py = 6.55, 3.9                                   # small PVS
    pvs2 = Ellipse((px, py), 0.62, 0.15, angle=u_ang, facecolor="#8fb3d9", edgecolor=navy,
                   lw=1.3, zorder=4)
    ax.add_patch(pvs2)
    ax.annotate("", xy=(px + 0.55 * ua, py + 0.55 * ub), xytext=(px - 0.55 * ua, py - 0.55 * ub),
                arrowprops=dict(arrowstyle="<->", color=navy, lw=1.4), zorder=5)
    ax.text(px - 0.05, py + 0.35, "$u$", fontsize=13, color=navy, fontweight="bold", ha="center")
    mx, my = 9.25, 3.9                                   # matched WM
    disc2 = Circle((mx, my), 1.15, facecolor=band, edgecolor="0.45", lw=1.1, zorder=0)
    ax.add_patch(disc2)
    far_ang = u_ang + (th_far if v9 else 78.0)
    _bundle(ax, mx, my, far_ang, 1.5, 1.05, n=6, color=tan, z=1, clip=disc2)
    ax.annotate("", xy=(mx + 0.75 * ua, my + 0.75 * ub), xytext=(mx - 0.75 * ua, my - 0.75 * ub),
                arrowprops=dict(arrowstyle="<->", color=navy, lw=1.4, ls="--"), zorder=5)
    ga, gb = np.cos(np.radians(far_ang)), np.sin(np.radians(far_ang))
    ax.annotate("", xy=(mx + 0.95 * ga, my + 0.95 * gb), xytext=(mx, my),
                arrowprops=dict(arrowstyle="-|>", color=maroon, lw=1.4), zorder=5)
    ax.add_patch(Arc((mx, my), 0.9, 0.9, angle=0, theta1=u_ang, theta2=far_ang,
                     color=maroon, lw=1.4, zorder=5))
    mid2 = np.radians(u_ang + (far_ang - u_ang) / 2)
    ax.text(mx + 0.62 * np.cos(mid2), my + 0.62 * np.sin(mid2), r"$\theta$", fontsize=13,
            color=maroon, fontweight="bold", ha="center", va="center", zorder=6)
    ax.annotate("", xy=(mx - 1.2, my), xytext=(px + 0.45, my),
                arrowprops=dict(arrowstyle="<->", color="0.45", lw=1.1))
    ax.text((px + mx - 0.75) / 2, my - 0.32, "$\\geq$6 mm", fontsize=12, color="0.4", ha="center")
    ax.text(5.9, 6.75, "Zoomed out: matched WM", fontsize=14, fontweight="bold", ha="left", va="top")
    if v9:
        c2_far = np.nan
        if fake is not None and "FakeControl_m6_cos2_theta" in fake.columns \
                and near is not None and "Axis_SE_deg" in near.columns:
            mm_ = fake[["Subject", "PVS_ID", "FakeControl_m6_cos2_theta"]].merge(
                near[["Subject", "PVS_ID", "Axis_SE_deg"]], on=["Subject", "PVS_ID"], how="left")
            c2_far = float(mm_.loc[mm_["Axis_SE_deg"] <= 5, "FakeControl_m6_cos2_theta"].mean())
        ax.text(mx, 1.15, rf"$\theta$ = {th_far:.0f}° (matched WM, measured)",
                fontsize=12.5, color=maroon, fontweight="bold", ha="center")
        ax.text(mx, 0.72, (f"= random orientation (cos²θ {c2_far:.2f} vs 0.33)"
                           if np.isfinite(c2_far) else "random orientation: cos²θ 0.33"),
                fontsize=11, color="0.35", ha="center")
    else:
        ax.text(mx, 1.15, r"$\theta$ not preferred: D$_\parallel$ $\approx$ D$_\perp$",
                fontsize=12.5, color=maroon, fontweight="bold", ha="center")
        ax.text(mx, 0.72, "same $u$ projected on another bundle", fontsize=11,
                color="0.35", ha="center")
    ax.plot([5.45, 5.45], [0.5, 6.6], color="0.75", lw=1)       # divider
    ax.set_title(letter, loc="left", fontsize=22, fontweight="bold")


def _alignment_tier_panel(ax, near, fake, letter, stats_rows, se_max=5.0):
    """G ('tier'): cos^2 between the PVS axis and the fibre direction in PVS,
    adjacent WM and matched WM (>= 6 mm), axes with SE <= se_max; per-subject
    medians, violins by group, random level 1/3."""
    need = {"cos2_theta", "NearPVS_cos2_theta", "Axis_SE_deg", "Group", "Subject", "PVS_ID"}
    if near is None or fake is None or not need <= set(near.columns) \
            or "FakeControl_m6_cos2_theta" not in fake.columns:
        return False
    a = near.loc[near["Axis_SE_deg"] <= se_max,
                 ["Subject", "Group", "PVS_ID", "cos2_theta", "NearPVS_cos2_theta"]]
    b = fake[["Subject", "PVS_ID", "FakeControl_m6_cos2_theta"]]
    d = a.merge(b, on=["Subject", "PVS_ID"], how="left")
    long = d.melt(id_vars=["Subject", "Group"],
                  value_vars=["cos2_theta", "NearPVS_cos2_theta", "FakeControl_m6_cos2_theta"],
                  var_name="Tier", value_name="c2").dropna()
    tmap = {"cos2_theta": "PVS", "NearPVS_cos2_theta": "NearPVS",
            "FakeControl_m6_cos2_theta": "Fake m6"}
    long["Tier"] = long["Tier"].map(tmap)
    sl = long.groupby(["Subject", "Group", "Tier"])["c2"].mean().reset_index()   # mean: see above
    order = ["PVS", "NearPVS", "Fake m6"]
    sns.violinplot(data=sl, x="Tier", y="c2", hue="Group", order=order,
                   hue_order=GROUPS_ORDER, palette=PALETTE, inner="quartile",
                   linewidth=1, cut=0, density_norm="width", ax=ax, legend=False)
    sns.stripplot(data=sl, x="Tier", y="c2", hue="Group", order=order,
                  hue_order=GROUPS_ORDER, dodge=True, palette=BLACK_PAL,
                  alpha=.4, size=2.4, ax=ax, legend=False)
    ax.axhline(1 / 3, color="k", lw=1.2, ls="--", zorder=0)
    ax.text(0.99, 1 / 3, " random\n direction", transform=ax.get_yaxis_transform(),
            ha="left", va="center", fontsize=11, color="0.3")
    for t in order:
        v = long.loc[long.Tier == t, "c2"]
        stats_rows.append(["Main_Fig2_alignment_tier", "cos2_theta", t,
                           f"cohort mean, axis SE <= {se_max:g} deg", float(v.mean()),
                           np.nan, np.nan, int(v.size), np.nan, np.nan])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([TIER_TICK.get(t, t) for t in order], fontsize=12)
    ax.set_xlabel("")
    ax.set_ylabel("Alignment with fibre direction\n(cos$^2\\theta$, PVS axis vs V1)",
                  fontsize=15, fontweight="bold")
    ax.tick_params(labelsize=12)
    ax.set_title(letter, loc="left", fontsize=22, fontweight="bold")
    sns.despine(ax=ax)
    return True


def fig_main_two_messages(long_df, near, metrics_a, tiers, W, outdir, stats_rows,
                          fake=None, alignment="precision"):
    """Main Figure 2 (v14). Row A (A-D): MD, D||, D-perp, MK, PVS vs adjacent WM,
    stars on lines. Row B: E schematic, F D|| - D-perp by tier, G PVS-fibre
    alignment."""
    metrics_a = [m for m in metrics_a if m in set(long_df.Metric)]
    if not metrics_a or "Dpar_minus_Dperp" not in set(long_df.Metric):
        print("  [!] main figure: required metrics missing"); return
    import matplotlib.gridspec as gridspec
    nA = len(metrics_a)
    fig = plt.figure(figsize=(23, 14))
    gs = gridspec.GridSpec(2, 12, figure=fig, hspace=0.5, wspace=1.4, height_ratios=[1, 1.05])
    letters = [chr(65 + i) for i in range(nA + 3)]
    tiers_a = [t for t in tiers if t in ("PVS", "NearPVS")]
    width = 12 // nA
    for k, met in enumerate(metrics_a):
        ax = fig.add_subplot(gs[0, width * k:width * k + width])
        _tier_panel_stars(ax, long_df, met, tiers_a, [("NearPVS", "NearPVS")], W,
                          letters[k], stats_rows)
    axE = fig.add_subplot(gs[1, 0:5])
    _schematic_panel(axE, letters[nA], near, fake)
    axF = fig.add_subplot(gs[1, 5:8])
    _tier_panel_stars(axF, long_df, "Dpar_minus_Dperp", ["PVS", "NearPVS", "Fake m6"],
                      [("NearPVS", "NearPVS"), ("Fake m6", "Fake_m6")], W,
                      letters[nA + 1], stats_rows, zero_line=True)
    axG = fig.add_subplot(gs[1, 8:12])
    if alignment == "tier" and _alignment_tier_panel(axG, near, fake, letters[nA + 2], stats_rows):
        pass
    else:
        if alignment == "tier":
            print("  [i] --alignment-panel tier needs the v9 extraction; using 'precision'")
        _alignment_panel(axG, near, letters[nA + 2], stats_rows)
    th, n_th = _measured_theta(near)
    stats_rows.append(["Main_Fig2_schematic", "Theta_deg", "Axis_SE_deg <= 5",
                       "median angle PVS axis vs fibre direction at the PVS", th,
                       np.nan, np.nan, n_th, np.nan, np.nan])
    handles = [plt.Rectangle((0, 0), 1, 1, color=PALETTE[g]) for g in GROUPS_ORDER]
    fig.legend(handles, GROUPS_ORDER, loc="upper right", ncol=3, fontsize=15,
               frameon=False, bbox_to_anchor=(0.995, 0.995))
    fig.text(0.01, 0.975, f"A$-${letters[nA - 1]}   PVS vs adjacent white matter: PVS are distinct structures characterized by faster and "
             "more heterogeneous diffusion",
             fontsize=18, fontweight="bold", ha="left", va="top")
    fig.text(0.01, 0.49, f"{letters[nA]}$-${letters[nA + 2]}   The along-axis anisotropy of "
             "PVS is shared with the adjacent white matter: PVS run along the adjacent fibres",
             fontsize=18, fontweight="bold", ha="left", va="top")
    fig.subplots_adjust(top=0.90, bottom=0.06, left=0.05, right=0.985)
    _save_both(fig, outdir, "Fig_main_PVS_distinct")


def fig_main_distinct(long_df, metrics, tiers, W, outdir, stats_rows):
    """NEW main Figure 2: PVS are distinct from their surroundings, in every
    group. Tier on x, group as hue. Under each panel: within-group q of
    PVS vs adjacent WM and PVS vs matched distant WM, from the mixed-effects workbook
    (Wilcoxon on per-subject deltas, FDR across metrics within a reference)."""
    metrics = [m for m in metrics if m in set(long_df.Metric)]
    if not metrics:
        print("  [!] main figure: none of the requested metrics found"); return
    ncol = 3; nrow = int(np.ceil(len(metrics) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.8 * ncol, 6.2 * nrow), squeeze=False)
    for k, (ax, met) in enumerate(zip(axes.flat, metrics)):
        d = long_df[(long_df.Metric == met) & long_df.Tier.isin(tiers)].copy()
        present = [t for t in tiers if t in set(d.Tier)]
        sc, unit = _scale(met)
        d["Plot"] = d["Value"] * sc
        sns.violinplot(data=d, x="Tier", y="Plot", hue="Group", order=present,
                       hue_order=GROUPS_ORDER, palette=PALETTE, inner="quartile",
                       linewidth=1, cut=0, ax=ax, legend=False)
        sns.stripplot(data=d, x="Tier", y="Plot", hue="Group", order=present,
                      hue_order=GROUPS_ORDER, dodge=True, palette=BLACK_PAL,
                      alpha=.4, size=2.2, ax=ax, legend=False)
        if met == "Dpar_minus_Dperp":
            ax.axhline(0, color="k", lw=1, ls="--", zorder=0)
        y0, y1 = ax.get_ylim(); span = y1 - y0
        ax.set_ylim(y0 - 0.30 * span, y1)
        r_near = workbook_row(W, met, "NearPVS")
        r_far = workbook_row(W, met, "Fake_m6")
        ax.text(0.5, 0.115, _q_line(r_near, "PVS vs adjacent WM:   "),
                transform=ax.transAxes, ha="center", va="center",
                fontsize=13, fontweight="bold", color="0.15")
        ax.text(0.5, 0.045, _q_line(r_far, "PVS vs matched WM ($\\geq$6 mm):   "),
                transform=ax.transAxes, ha="center", va="center",
                fontsize=12, color="0.35")
        for ref, r in (("NearPVS", r_near), ("Fake_m6", r_far)):
            for g in GROUPS_ORDER:
                stats_rows.append(["Main_Fig2", met, f"PVS - {ref} | {g}",
                                   "delta vs 0 (mixed-effects workbook)",
                                   np.nan if r is None else float(r[f"delta_{g}"]),
                                   np.nan if r is None else float(r[f"p_{g}_vs_0"]),
                                   np.nan if r is None else float(r[f"q_{g}_vs_0"]),
                                   np.nan, np.nan, np.nan])
        ax.set_xticks(range(len(present)))
        ax.set_xticklabels([TIER_TICK.get(t, t) for t in present], fontsize=14)
        ax.set_xlabel(""); ax.set_ylabel(f"{PRETTY.get(met, met)}{unit}", fontsize=16,
                                         fontweight="bold")
        ax.tick_params(axis="y", labelsize=13)
        ax.set_title(chr(65 + k), loc="left", fontsize=22, fontweight="bold")
        sns.despine(ax=ax)
    for ax in list(axes.flat)[len(metrics):]:
        ax.axis("off")
    handles = [plt.Rectangle((0, 0), 1, 1, color=PALETTE[g]) for g in GROUPS_ORDER]
    fig.legend(handles, GROUPS_ORDER, loc="upper right", ncol=3, fontsize=15,
               frameon=False, bbox_to_anchor=(0.99, 1.0))
    fig.suptitle("PVS differ from the surrounding white matter",
                 fontsize=20, fontweight="bold", x=0.02, ha="left", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save_both(fig, outdir, "Fig_main_PVS_distinct")


def fig_suppl_reference(near, fake, metrics, cols_near, cols_fake, ref, W,
                        outdir, stats_rows):
    """Regenerated Suppl. Fig 10/11: paired PVS vs one reference tier, per
    group, subject lines joining the pairs. Above each pair: median paired
    difference as % of the reference median, stars = workbook q_<G>_vs_0.
    Panel subtitle: RHI vs HC on the PVS excess (workbook q_RHI_vs_HC)."""
    metrics = [m for m in metrics if "Original" in cols_near.get(m, {})]
    ncol = 3; nrow = int(np.ceil(len(metrics) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.0 * ncol, 5.6 * nrow), squeeze=False)
    reflab = "Adjacent WM" if ref == "NearPVS" else f"Matched WM\n{margin_name(ref)}"
    for k, (ax, met) in enumerate(zip(axes.flat, metrics)):
        m = paired_table(near, fake, met, cols_near, cols_fake, ref)
        row = workbook_row(W, met, ref)
        if m is None or m.empty or row is None:
            ax.axis("off"); continue
        if not crosscheck(m, row, met, ref):
            ax.text(0.5, 0.5, "numbers differ from workbook\n(see console)",
                    transform=ax.transAxes, ha="center", color="red")
        sc, unit = _scale(met)
        pos, cols, vals = [], [], []
        for gi, g in enumerate(GROUPS_ORDER):
            mg = m[m.Group == g]
            for j, (c, alpha) in enumerate((("PVS", 1.0), ("Ref", 0.45))):
                x = gi * 2.6 + j
                parts = ax.violinplot(mg[c] * sc, positions=[x], widths=0.8,
                                      showextrema=False)
                for b in parts["bodies"]:
                    b.set_facecolor(PALETTE[g]); b.set_edgecolor("0.2")
                    b.set_alpha(alpha)
                ax.hlines(np.median(mg[c] * sc), x - 0.25, x + 0.25, color="k", lw=1.5)
            for _, r in mg.iterrows():
                ax.plot([gi * 2.6, gi * 2.6 + 1], [r.PVS * sc, r.Ref * sc],
                        color="0.3", lw=0.6, alpha=0.35)
            ax.scatter(np.full(len(mg), gi * 2.6), mg.PVS * sc, s=6, color="k", alpha=.5, zorder=3)
            ax.scatter(np.full(len(mg), gi * 2.6 + 1), mg.Ref * sc, s=6, color="k", alpha=.5, zorder=3)
        y0, y1 = ax.get_ylim(); span = y1 - y0
        ax.set_ylim(y0, y1 + 0.16 * span)
        for gi, g in enumerate(GROUPS_ORDER):
            mg = m[m.Group == g]
            pct = 100 * mg.Delta.median() / abs(mg.Ref.median()) \
                if mg.Ref.median() != 0 else np.nan
            q = float(row[f"q_{g}_vs_0"])
            colr = "0.4" if not (np.isfinite(q) and q < 0.05) else \
                   ("#1b7a4a" if mg.Delta.median() > 0 else "#5b2a86")
            ax.text(gi * 2.6 + 0.5, y1 + 0.06 * span, f"{pct:+.1f}% {stars(q)}",
                    ha="center", fontsize=14, fontweight="bold", color=colr)
            stats_rows.append([f"Suppl_ref_{ref}", met, g, "PVS - ref, % of ref (mixed-effects workbook q)",
                               float(mg.Delta.median()), float(row[f"p_{g}_vs_0"]), q,
                               pct, np.nan, np.nan])
        ax.set_xticks([gi * 2.6 + j for gi in range(3) for j in (0, 1)])
        ax.set_xticklabels([f"{g}\n{'PVS' if j == 0 else reflab}"
                            for g in GROUPS_ORDER for j in (0, 1)], fontsize=12)
        ax.set_ylabel(f"{PRETTY.get(met, met)}{unit}", fontsize=15, fontweight="bold")
        ax.tick_params(axis="y", labelsize=13)
        qx = float(row["q_RHI_vs_HC"])
        ax.set_title(f"{chr(65 + k)}   {PRETTY.get(met, met)}\n"
                     f"PVS excess, RHI vs HC: q = {qx:.2f}",
                     loc="left", fontsize=15, fontweight="bold")
        sns.despine(ax=ax)
    for ax in list(axes.flat)[len(metrics):]:
        ax.axis("off")
    fig.suptitle(f"PVS vs {REF_PRETTY.get(ref, ref)}: within-subject contrasts",
                 fontsize=20, fontweight="bold", x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    _save_both(fig, outdir, f"Fig_suppl_ref-{ref}")


# ==========================================
# MAIN
# ==========================================
def main():
    global AXIS_SE_MAX
    ap = argparse.ArgumentParser(description="v9 figures (main Fig 2, suppl., sanity checks).")
    ap.add_argument("--results", default=RESULTS_FILE)
    ap.add_argument("--groups", default=GROUPS_FILE)
    ap.add_argument("--out", default=OUTPUT_DIR)
    ap.add_argument("--metrics", nargs="*", default=None,
                    help="restrict the per-metric figures to these base metrics")
    ap.add_argument("--min-n", type=int, default=5,
                    help="min subjects per group per depth bin to test (default 5)")
    ap.add_argument("--min-clusters-per-bin", type=int, default=5,
                    help="v9: a subject enters a depth bin only with >= this many "
                         "clusters there (default 5)")
    ap.add_argument("--axis-se-max", type=float, default=None, metavar="DEG",
                    help="v9: axis-precision filter for Dpar/Dperp (clusters with "
                         "Axis_SE_deg <= DEG; NaN excluded). Default: none, which "
                         "matches the v9 mixed-effects workbook.")
    ap.add_argument("--mixed-xlsx", default=None,
                    help="MixedEffects_v10 workbook (v9 accepted with --fdr-families single). Required for the main and "
                         "supplementary figures; used for the Fig4 stars.")
    ap.add_argument("--exclude-metrics", nargs="*", default=["lowb_ADC"],
                    help="v16: metrics left out of every figure and FDR family "
                         "(default lowb_ADC, a duplicate of lowb_MD)")
    ap.add_argument("--fdr-families", choices=["split", "single"], default="single",
                    help="v10: main and lowb_* metrics as separate FDR families "
                         "(default) or one family (v9). Must match the workbook.")
    ap.add_argument("--main-metrics", nargs="+",
                    default=["MD", "Dpar", "Dperp", "MK_mf"],
                    help="row-A metrics of the main figure, in order (v12)")
    ap.add_argument("--alignment-panel", choices=["precision", "tier"], default="precision",
                    help="v15: panel G = alignment at the PVS by axis precision (default) "
                         "or alignment in PVS / adjacent / matched WM ('tier', needs v9)")
    ap.add_argument("--main-tiers", nargs="+", default=["PVS", "NearPVS", "Fake m6"],
                    help="tiers shown in the main figure (long-table names)")
    ap.add_argument("--suppl-refs", nargs="+", default=["NearPVS", "Fake_m6"],
                    help="one paired supplementary figure per reference tier")
    ap.add_argument("--suppl-metrics", nargs="+",
                    default=["MD", "AD", "RD", "MK_mf", "AK_mf", "RK_mf", "FA", "Dpar", "Dperp"])
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    AXIS_SE_MAX = a.axis_se_max

    gmap = load_groups(a.groups)
    print(f"Groups loaded: {len(gmap)} subjects")
    print(f"results file: {a.results}")
    near = pd.read_excel(a.results, sheet_name="Case_A_Near_PVS")
    fake = pd.read_excel(a.results, sheet_name="Case_B_Fake_PVS")
    qc = pd.read_excel(a.results, sheet_name="Subject_QC")
    for d in (near, fake, qc):
        d["Subject"] = d["Subject"].astype(str).str.strip()
        d["Group"] = d["Subject"].map(gmap)
    near = near.dropna(subset=["Group"]); fake = fake.dropna(subset=["Group"])
    qc = qc.dropna(subset=["Group"])
    print(f"Rows with a group: Case_A {len(near)}, Case_B {len(fake)}, QC {len(qc)}")

    if AXIS_SE_MAX is None:
        print("orientation filter: NONE (all clusters; matches the v9 mixed-effects workbook)")
    else:
        kept = int((pd.to_numeric(near["Axis_SE_deg"], errors="coerce") <= AXIS_SE_MAX).sum())
        print(f"orientation filter: Axis_SE_deg <= {AXIS_SE_MAX} deg, NaN excluded -> "
              f"{kept}/{len(near)} clusters for Dpar/Dperp")

    W = None
    if a.mixed_xlsx:
        W = pd.read_excel(a.mixed_xlsx, sheet_name="within_subject")
        # v10: the workbook's FDR families must match this run's
        try:
            ri = pd.read_excel(a.mixed_xlsx, sheet_name="run_info")
            ri = dict(zip(ri.iloc[:, 0].astype(str), ri.iloc[:, 1]))
        except ValueError:
            ri = {}
        wb_fam = str(ri.get("fdr_families", "single"))
        if wb_fam != a.fdr_families:
            raise SystemExit(f"--fdr-families {a.fdr_families}, but {a.mixed_xlsx} was "
                             f"corrected with '{wb_fam}' families. Regenerate the workbook "
                             f"with mixed-effects v10 (--fdr-families {a.fdr_families}) or "
                             f"pass --fdr-families {wb_fam}.")
        print(f"FDR families: {a.fdr_families} (workbook agrees)")
        print(f"mixed-effects workbook: {a.mixed_xlsx}  "
              f"(reference tiers: {sorted(W.Ref_tier.unique())})")

    cols_near, cols_fake = parse_columns(near), parse_columns(fake)
    for m in (a.exclude_metrics or []):                 # v16
        cols_near.pop(m, None); cols_fake.pop(m, None)
    if a.exclude_metrics:
        print(f"excluded metrics: {sorted(a.exclude_metrics)}")
    bases = sorted(set(cols_near) | set(cols_fake))
    all_bases = list(bases)
    if a.metrics:
        bases = [b for b in bases if b in a.metrics]
    print(f"Metrics found ({len(bases)}): {bases}")
    tiers_found = sorted({t for v in cols_fake.values() for t in v} |
                         {t for v in cols_near.values() for t in v})
    print(f"Tiers found: {tiers_found}")

    # long table: one row per subject x metric x tier (all metrics, so the main
    # figure is independent of --metrics)
    recs = []
    for base in all_bases:
        near_b = apply_linearity_filter(near, base)
        fake_b = apply_linearity_filter(fake, base)
        if base in cols_near and "Original" in cols_near[base]:
            s = subject_level(near_b, cols_near[base]["Original"])
            s.columns = ["Subject", "Group", "Value"]; s["Metric"], s["Tier"] = base, "PVS"
            recs.append(s)
        if base in cols_near and "NearPVS" in cols_near[base]:
            s = subject_level(near_b, cols_near[base]["NearPVS"])
            s.columns = ["Subject", "Group", "Value"]; s["Metric"], s["Tier"] = base, "NearPVS"
            recs.append(s)
        for tier, col in cols_fake.get(base, {}).items():
            if tier.startswith(FAKE_PREFIX):
                s = subject_level(fake_b, col)
                s.columns = ["Subject", "Group", "Value"]
                s["Metric"] = base; s["Tier"] = "Fake " + tier.split("_")[-1]
                recs.append(s)
    long_df = pd.concat(recs, ignore_index=True)
    def _margin(t):
        m = re.search(r"(\d+)", t)
        return int(m.group(1)) if m else 0
    tier_order = ["PVS", "NearPVS"] + sorted(
        (t for t in long_df["Tier"].unique() if t.startswith("Fake")), key=_margin)

    # Fig2 first pass: one FDR family across ALL metrics x depth bins
    depth, pool = {}, []
    for base in bases:
        if base in cols_near and "Original" in cols_near[base] and BIN_COL in near.columns:
            col = cols_near[base]["Original"]
            sl = depthbin_subject_table(apply_linearity_filter(near, base), col,
                                        a.min_clusters_per_bin)
            tests = depthbin_tests(sl, col, a.min_n)
            depth[base] = (sl, col, tests)
            pool += [((base, b), dp[1]) for b, dp in tests.items()]
    # v10: one FDR family per metric family (main / lowb), or one overall
    qmap, fam_n = {}, {}
    for fam in sorted({fdr_family_of(k[0], a.fdr_families) for k, _ in pool}):
        sub = [(k, p) for k, p in pool if fdr_family_of(k[0], a.fdr_families) == fam]
        for (k, _), q in zip(sub, fdr([p for _, p in sub])):
            qmap[k] = q
        fam_n[fam] = int(np.isfinite([p for _, p in sub]).sum())
    n_tests = sum(fam_n.values())
    print(f"depth bins: >= {a.min_clusters_per_bin} clusters per subject-bin; FDR "
          + ", ".join(f"{f}: {n} tests" for f, n in fam_n.items())
          + " (metrics x bins within each family)")

    stats_rows = []
    for base in bases:
        fig_tiers(long_df, base, tier_order, a.out, stats_rows)
        fig_raw_vs_standardised(long_df, base, tier_order, a.out)
        if base in depth:
            sl, col, tests = depth[base]
            fig_thicknessbins(sl, base, col, tests, qmap, a.out, stats_rows)
        fig_within_subject_contrast(near, fake, base, cols_near, cols_fake,
                                    W, a.out, stats_rows)
        print(f"  figures done: {base}")
    fig_wm_availability(qc, a.out, stats_rows)
    fig_match_quality(qc, a.out, stats_rows)

    made_main = False
    if W is None:
        print("\n  [skip] main + supplementary figures: they take their statistics from "
              "the mixed-effects workbook -- pass --mixed-xlsx")
    elif AXIS_SE_MAX is not None:
        print("\n  [skip] main + supplementary figures: the mixed-effects workbook used NO axis "
              "filter, so its q-values do not describe filtered clusters. Rerun "
              "without --axis-se-max (or refit the workbook with the same filter).")
    else:
        fig_main_two_messages(long_df, near, a.main_metrics, a.main_tiers, W, a.out, stats_rows,
                              fake=fake, alignment=a.alignment_panel)
        for ref in a.suppl_refs:
            if ref not in set(W.Ref_tier):
                raise SystemExit(f"--suppl-refs {ref}: not in the workbook "
                                 f"({sorted(W.Ref_tier.unique())})")
            fig_suppl_reference(near, fake, a.suppl_metrics, cols_near, cols_fake,
                                ref, W, a.out, stats_rows)
        made_main = True

    if CROSSCHECK:
        print("\n  [!] FIGURE vs WORKBOOK MISMATCH (median per-subject delta):")
        for met, ref, g, here, wb in CROSSCHECK:
            print(f"      {met:17s} {ref:8s} {g:4s} here {here:+.4e}  workbook {wb:+.4e}")
        print("      Fig4 fell back to its own tests for these; do not use the "
              "workbook q-values with these figures.")
    elif W is not None:
        print("\n  cross-check: per-group median deltas match the mixed-effects workbook for "
              "every metric and reference plotted")

    st = pd.DataFrame(stats_rows, columns=["Figure", "Metric", "Level",
                                           "Contrast", "cliffs_delta_or_median", "p",
                                           "q_fdr", "raw_diff_or_pct", "pooled_sd",
                                           "raw_rel_to_PVS"])
    run_info = pd.DataFrame([{
        "script": os.path.basename(__file__), "results": a.results,
        "groups": a.groups, "mixed_xlsx": a.mixed_xlsx or "",
        "axis_se_max": "none" if AXIS_SE_MAX is None else AXIS_SE_MAX,
        "min_clusters_per_bin": a.min_clusters_per_bin, "min_n": a.min_n,
        "fdr_families": a.fdr_families,
        "depthbin_fdr_family": "; ".join(f"{f}: {n} tests" for f, n in fam_n.items()),
        "n_subjects": {g: int((qc.Group == g).sum()) for g in GROUPS_ORDER},
        "workbook_crosscheck": "ok" if not CROSSCHECK else f"{len(CROSSCHECK)} mismatches",
        "main_and_suppl_made": made_main}]).T.reset_index()
    run_info.columns = ["key", "value"]
    xlsx = os.path.join(a.out, "Figure_statistics_v16-3.xlsx")
    with pd.ExcelWriter(xlsx) as w:
        run_info.to_excel(w, sheet_name="run_info", index=False)
        st.to_excel(w, sheet_name="all_tests", index=False)
        for fam in st["Figure"].unique():
            st[st.Figure == fam].to_excel(w, sheet_name=fam[:31], index=False)
    print(f"\nFigures -> {a.out}")
    print(f"Statistics -> {xlsx}")
    print("\nWHAT EACH FIGURE ANSWERS")
    print("  Fig_main_PVS_distinct : PVS vs adjacent WM and matched WM, per group (main Fig 2)")
    print("  Fig_suppl_ref-*       : the same contrasts, paired, all metrics (suppl.)")
    print("  Fig1 / Fig5           : RHI-HC by tier, raw vs standardised (not PVS-specific)")
    print("  Fig2                  : RHI-HC by depth (sensitivity)")
    print("  Fig3 / Fig6           : control-pool comparability (m6 primary, m10 sensitivity)")
    print("  Fig4                  : PVS minus own reference, per group")


if __name__ == "__main__":
    main()
