#!/usr/bin/env python3
"""
step23_plot_all_markers_across_cohorts_and_ROIs.py
(v12: all_stats_and_plot_on_PVS_saving_statistics_v12_families.py)

CHANGELOG v13 (2026-10-15) -- absolute volumes plotted on request, never tested
---------------------------------------------------------------------------------
 * --also-plot-absolute: the absolute volumes in mm3 are plotted as well:
     synthseg/VOL.xlsx  tissue-size-in-mm   (every structure, eTIV included)
     PVSROI/VOL.xlsx    tissue-size-in-mm, pvs-size-within-tissue-in-mm
   They belong to NO FDR family (pseudo-family 'ABS' in the output): the
   family sizes, q-values, stars, heatmaps and the completeness check are
   exactly as without the flag. Their labels show Cliff's delta with its CI
   (bold when |delta| >= --bold-min-delta) and never stars or colour, since
   no q exists for them. Raw p-values are in the Excel sheets with
   fdr_family = 'ABS', in_FDR_family = 0 and p_value_FDR empty.
   Voxel counts are not plotted (duplicates of the mm3 values).
   Rationale (v6): absolute volumes mostly measure head size, so the relative
   volumes carry the atrophy / PVS-burden tests; the absolute ones are shown
   for reference only.

CHANGELOG v12 (2026-10-05) -- display only, statistics unchanged
------------------------------------------------------------------
 * Two independent cues on each pairwise effect-size label:
     colour  = significance: --significant-colour (blue) when q < 0.05 within
               the family (the star criterion), black otherwise -- as v11;
     weight  = magnitude: bold when |delta| >= --bold-min-delta (default 0.474,
               a large effect by the Romano et al. 2006 convention), regular
               otherwise. Applied to the unrounded delta, so a label printed
               as 0.47 is bold only if |delta| >= 0.474.
   v11 drew every label in bold.

CHANGELOG v11 (2026-10-05) -- display only, statistics unchanged
------------------------------------------------------------------
 * EFFECT-SIZE LABELS IN ONE COLOUR, ONLY WHEN SIGNIFICANT. A pairwise label
   ('RHI > HC (delta=+0.45 [+0.11, +0.75]) **') is drawn in --significant-colour
   (default #1f5aa6, blue) when its q within the family is < 0.05 -- the same
   criterion as the stars, so colour and stars always agree -- and in black
   otherwise. v8-v10 coloured every label by the group showing the
   pathological sign, i.e. up to three colours per plot whatever q was.
   pathology_direction() is kept but no longer used for colour.

CHANGELOG v10 (2026-10-04) -- technical exclusions of within-PVS values
-------------------------------------------------------------------------
 Same change as v8-1 (which was built on v8), here applied to v9.
 * TWO EXCLUSION CRITERIA, fixed before the statistics and applied to every
   subject and every PVS region alike. A subject's within-PVS diffusion in a
   region (every column of diffusion/<region>/PVSROI/*.xlsx: the diffusion
   tests of family A, including the low-b fit, and family E; family B with
   --lowb-family separate) is left out when
     1. its PVS mask in that region has fewer than --min-pvs-voxels voxels
        (pvs-size-within-tissue-in-voxels in volumetry/<region>/PVSROI/VOL.xlsx):
        a median over a handful of voxels is not a reliable estimate;
     2. it is physically impossible: any within-PVS diffusivity (MD, AD, RD of
        the full-range AND the low-b fit) outside [0, --max-diffusivity] mm2/s
        (free water at 37 C is about 3.0e-3), or FA / low-b FA outside [0, 1]
        (only possible with a negative eigenvalue, i.e. a failed fit).
   The whole subject x region is excluded, not only the offending metric: all
   metrics of a region come from the same voxels, so if one is impossible the
   others cannot be trusted either, and all metrics of a region keep the same
   subjects. The low-b and full-range values are checked whichever families
   are run, so the exclusions never depend on --families or --lowb-family.
 * Volumetry (PVS volume / ROI, ROI volume / eTIV, SynthSeg) is never
   affected: a small PVS volume is a valid measurement.
 * Both flags are REQUIRED (no silent default). --min-pvs-voxels 0 and
   --max-diffusivity off switch the criteria off and give v9's statistics.
 * Every included subject x PVS region is listed with its voxel count, the
   checked values and the reason for any exclusion in value_exclusions.csv
   and the Excel sheet value_exclusions. --list-tests writes it too, before
   any statistic is computed.
 * Every result row records the group sizes actually tested (n_hc, n_ppcs,
   n_boxers): they can now differ between regions.

CHANGELOG v9 (2026-10-04)
-------------------------
 * LOW-B AND FULL-RANGE FITS ARE ONE FAMILY (--lowb-family merged, default):
   family A now holds, per PVS region, the 7 DKI metrics within PVS, PVS
   volume / ROI, ROI volume / eTIV AND the 4 low-b metrics within PVS
   (MD AD RD FA) -> 13 per region (39 with three regions). Low-b ADC is not
   used (duplicate of low-b MD). --lowb-family separate reproduces v8
   (family B on its own).
 * Heatmap A gets a fifth block, "Low-b fit within PVS".
 * Violin y-labels use readable metric names ("MK within PVS", "Low-b MD
   within PVS"), never the file stems (MK_mf, lowb_MD).


CHANGELOG v8 (2026-10-01) -- display only, statistics unchanged
------------------------------------------------------------------
 * EFFECT-SIZE COLOUR = the group that shows the pathological sign:
   maroon (RHI), darkgoldenrod (PPCS), darkkhaki (HC). For metrics where an
   INCREASE is pathological (MD AD RD MK AK RK Dperp, their low-b and
   PVS-weighted versions, PVS volume, ventricle/CSF volume) that is the group
   with the higher value; where a DECREASE is pathological (FA, Dpar, tissue /
   ROI volume = atrophy) the group with the lower value. Black only when
   delta = 0.00.
   Two bugs in the old colour logic are fixed:
     - contrasts involving HC were black whenever HC carried the sign
       (e.g. 'PPCS < HC' for AD), instead of darkkhaki;
     - the direction was found by substring search in the Case ID, and every
       '..._diffusion_within_pvs' test contains 'pvs', which is on the
       'increase' list -- so FA (and low-b FA) WITHIN PVS was coloured as if
       an increase were pathological, i.e. inverted.
   The direction now comes from the metric and the column, not the Case ID.
 * LABELS show only the effect size, its 95% CI and the stars:
     'RHI > HC (delta=+0.27 [-0.09, +0.60]) *'
   p and q are no longer printed on the plots (they remain in the Excel
   output); stars are still from q within the family.

CHANGELOG v7 (2026-10-01)
-------------------------
 * PAIRWISE FDR ACROSS THE THREE CONTRASTS. v6 corrected HC-PPCS, HC-RHI and
   PPCS-RHI separately (27 tests each in family A), so the 5% FDR did not
   hold over the union of the three. --pairwise-fdr union (default) pools
   the three contrasts WITHIN each family (family A: 3 x 27 = 81 tests),
   separately for the Mann-Whitney and the t-test sheets.
   --pairwise-fdr per-contrast reproduces v6.
 * Omnibus tests (Kruskal-Wallis, ANOVA) are unchanged: FDR within each
   family. They are the primary inference (the heatmaps); pairwise tests show
   which groups drive an omnibus result.
 * Tests whose p-value could not be computed are left out of every family
   (never counted as p = 1).
 * Plot labels name the scope: "q(A, 3 contrasts)" (union) or "q(A)".

CHANGELOG v6 (2026-10-01)
-------------------------
 The paper compares PVS diffusion with PVS volume and with atrophy, so the
 competing measures share ONE family (same correction for all of them):
   A  group discrimination in the PVS regions, all from PVSROI, per region:
        diffusion_within_pvs of MD AD RD FA MK AK RK            (7)
        PVS volume / ROI volume   (PVSROI/VOL.xlsx)              (1)
        ROI volume / eTIV         (PVSROI/VOL.xlsx)              (1)
      -> 9 x regions (27 with three regions; = submitted Figure 2 without the
         PVS-weighted columns and without the absolute volumes)
   B  low-b, diffusion_within_pvs of lowb_MD lowb_AD lowb_RD lowb_FA x regions
   E  PVS-weighted, PVS-Weighted-diffusion of the 7 metrics of A x regions
   F  anatomical atrophy, synthseg/VOL.xlsx tissue volume / eTIV, every
      SynthSeg structure (eTIV itself excluded: its ratio is 1)
 v5's family C is gone: PVS and ROI volume moved into A, SynthSeg into F.
 Volumes are RELATIVE only: absolute volumes mostly measure head size, and
 voxel / mm3 versions duplicate each other.
 The whole-brain PVS / eTIV (eTIV/PVSROI/VOL.xlsx) is not a PVS region and
 belongs to no family.
 NEW: one Kruskal-Wallis heatmap per family (rows = regions or structures,
 cells = q within the family and H), styled like submitted Figure 2,
 written as PNG + PDF next to the statistics.

CHANGELOG v5 (2026-09-30)
-------------------------
 * FDR FAMILIES. Each test belongs to exactly one family, defined by the
   question it answers; q-values are computed WITHIN a family only (and, as
   before, separately for each pairwise comparison and for ANOVA / KW):
     A  PVS diffusion, primary   : diffusion_within_pvs  of MD AD RD FA MK AK RK
                                   x every PVSROI region          (7 x regions)
     B  low-b (Reviewer 3, R3.2) : diffusion_within_pvs  of lowb_MD lowb_AD
                                   lowb_RD lowb_FA x regions      (4 x regions)
     E  PVS-weighted, secondary  : PVS-Weighted-diffusion of the 7 metrics of A
                                   x regions                      (7 x regions)
     C  volumetry                : PVS volume (PVSROI/VOL.xlsx) and tissue
                                   volume (synthseg/VOL.xlsx), --volume-units
   (Family D, PVS vs adjacent vs matched WM, is the tier analysis in the
   mixed-effects script, not here.)
   --families selects which are run (default A B E). Tests outside the chosen
   families are neither computed nor plotted: Dpar, Dperp, Dpar_minus_Dperp,
   lowb_ADC, diffusion_within_tissue, Dex-*, PVS-weighted low-b, scPVSROI,
   CSF_fast, synthseg diffusion.
 * Each family's size is checked: every metric must be present for every
   region, otherwise the run stops (--allow-incomplete-families to continue),
   so 'A = 7 x 3 = 21' is guaranteed, not assumed. Sizes and members are
   printed and saved (sheet fdr_families); every result row has fdr_family.
 * The plot label shows the family: "q(A) = 0.012".
 * Replaces v4's --family / --diffusion-columns / --also-plot-excluded, which
   defined a single family; the family definitions above are now the only
   selection mechanism.

CHANGELOG v4 (2026-09-30)
-------------------------
 * Diffusion columns tested by default: diffusion_within_pvs and
   PVS-Weighted-diffusion only. Dex-ratio, Dex-tissue and
   diffusion_within_tissue are no longer read: they are not computed, not in
   the FDR family, and not plotted -- also not with --also-plot-excluded
   (v3 plotted every column of the excluded spreadsheets).
   --diffusion-columns can still name any of the five columns explicitly.

CHANGELOG v3 (2026-09-29)
-------------------------
 * WHICH SPREADSHEETS ARE TESTED (the FDR family is built from these only):
     diffusion  tree, <tissue>/PVSROI/<metric>.xlsx  -> diffusion columns
     volumetry  tree, <tissue>/PVSROI/VOL.xlsx       -> PVS-volume columns only
     volumetry  tree, <tissue>/synthseg/VOL.xlsx     -> tissue-volume columns only
   Everything else (scPVSROI, CSF_fast, synthseg diffusion, the tissue columns
   of PVSROI/VOL.xlsx, the PVS columns of synthseg/VOL.xlsx) is not tested and
   not in the family. The folder name must be EXACTLY 'PVSROI' (scPVSROI is a
   different segmentation). Every file is listed with its role in
   file_selection.csv. --also-plot-excluded draws the others with p only.
 * The DIFFUSION tree is now scanned (v1/v2 read only the volumetry tree), and
   its tissue / segmentation / metric are parsed from the path. Before, a
   diffusion file got ROI = map = metric = '_', so every tissue and
   segmentation of one metric shared the same Case ID.
 * --volume-units mm (default) keeps one unit: voxel counts are the mm3
   values / 0.512, i.e. duplicate tests. Relative (normalised) columns are
   always kept. --diffusion-columns restricts the diffusion columns
   (default: all five in the spreadsheet).
 * Columns with no variance across the included subjects are dropped and
   listed (e.g. relative-tissue-size of eTIV, which is 1 by construction).
 * --family now defaults to 'all' = every SELECTED test; a family file still
   narrows it further.
 * Plot fixes: the y-label comes from the column (not the Case ID, which
   contains 'PVSROI' and mislabelled tissue/eTIV as 'to ROI'), diffusion files
   get '<metric> within PVS (mm2/s)' etc.; the eTIV y-range and the hidden
   delta label now apply only to the eTIV volume itself, not to the
   whole-brain PVS / eTIV columns in eTIV/PVSROI/VOL.xlsx, which v2 plotted on
   a 0.5e6-2.5e6 axis (data invisible).

CHANGELOG v2 (2026-09-29)
-------------------------
 * LABELS: each pairwise bracket now reads
       RHI > HC (delta = +0.45 [0.11, 0.75])
       p = 0.015, q = 0.12
   i.e. Cliff's delta with its 95% bootstrap CI, the raw Mann-Whitney p and
   the FDR q, all as numbers. Stars are still drawn from q ONLY, exactly as in
   the submitted figures, so the significance criterion does not change.
   Tests outside the FDR family show "p = ... (not in FDR family)" and no stars.
 * CONFIDENCE INTERVALS: percentile bootstrap of Cliff's delta, subjects
   resampled within each group (--n-boot, default 10000, fixed --seed).
   Saved as delta_ci_lo / delta_ci_hi in the Excel output.
 * FDR FAMILY IS EXPLICIT. v1 pooled every numeric column of every .xlsx under
   the input folder, so the family (and how harsh q was) depended on what
   happened to be in the directory, including duplicates such as volume in mm3
   AND in voxels. v2 requires --family: either a text file listing the tests
   (one Case ID or glob pattern per line) or the literal 'all' (= v1
   behaviour). Run --list-tests first: it writes every candidate Case ID and
   the group sizes WITHOUT computing any statistic, so the family can be chosen
   blind to the results. Family size per comparison is printed and saved.
 * INCLUSION FROM THE STUDY SPREADSHEET. v1 dropped a hard-coded set
   {101, 127, 130, 143, 146, 153}, which disagrees with the 2026-08-28
   inclusion list (that list also drops 114, 128, 137, 140, 147, 148, 171,
   172). v2 keeps exactly the subjects marked 'Yes' in --study-xlsx, checks
   that every data sheet holds each included subject in the right group, and
   stops with an error if one is missing (--allow-missing-subjects to
   downgrade to a warning) or sits in the wrong group sheet.
 * run_info and fdr_family sheets in the Excel output; --out-suffix so a v2
   run never overwrites the v1 plots.

Original description
--------------------
Extended process_group_violinplots.py with FDR-Corrected Plotting

Structure:
1. Pass 1: Scan all files, calculate stats (ANOVA, KW, U-tests), store in memory.
2. Correction: Apply Benjamini-Hochberg FDR to all stored p-values.
3. Pass 2: Re-scan files, generate plots using the CORRECTED Q-values for significance stars.
"""
import os
import re
import argparse
import datetime
import zlib
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, f_oneway, ttest_ind
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import ScalarFormatter
import numpy as np
from statsmodels.stats.multitest import multipletests
from pipeline_config import cfg

ROIs_of_increase_interest = ["Dperp","MD","RD","AD","MK","AK","RK","3rd_ventricle","4th_ventricle","CSF","lateral_ventricle","inferior_lateral_ventricle","pvs"]
#ROIs_of_increase_interest = ["MD","RD","AD","MK","AK","RK","3rd_ventricle","4th_ventricle","CSF","lateral_ventricle","inferior_lateral_ventricle","pvs"]
ROIs_of_decrease_interest = ["Dpar","FA","cortex","cerebral_white_matter","basal_ganglia","brain-stem","thalamus","caudate","putamen","pallidum","hippocampus","amygdala","midbrain","cerebellum_cortex","cerebellum_white_matter"]
#ROIs_of_decrease_interest = ["FA","cortex","cerebral_white_matter","basal_ganglia","brain-stem","thalamus","caudate","putamen","pallidum","hippocampus","amygdala","midbrain","cerebellum_cortex","cerebellum_white_matter"]

# --- SHARED HELPERS ---
def _to_float_or_nan(v):
    if v is None:
        return np.nan
    try:
        return float(v)
    except Exception:
        return np.nan

def cliffs_delta_from_u(x, y, alternative_now='two-sided'):
    # Ensure numpy arrays
    x = np.asarray(x)
    y = np.asarray(y)
    
    # Filter for finite values only
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    
    n1 = len(x)
    n2 = len(y)

    # SAFETY CHECK: If either group is empty, we cannot calculate stats
    if n1 == 0 or n2 == 0:
        return np.nan, np.nan, np.nan

    # Calculate stats
    try:
        U, p = mannwhitneyu(x, y, alternative=alternative_now)
        delta = (2 * U / (n1 * n2)) - 1
        return U, p, delta
    except ValueError:
        # Catch specific scipy errors (e.g. all values identical)
        return np.nan, np.nan, np.nan

# ---------------- v2 configuration (set in main) ----------------
INCLUDED = None          # {'hc': set(ids), 'ppcs': set(ids), 'boxers': set(ids)}
ALLOW_MISSING = False
N_BOOT = 10000
SEED = 0
_CHECKED_FILES = set()
N_REPORT = []            # (file, group, n_in_sheet, n_kept, removed, missing)
STUDY_SHEETS = {'HC': 'hc', 'PPCS': 'ppcs', 'Boxers': 'boxers', 'RHI': 'boxers'}

# ---------------- v6: FDR families ----------------
TISSUE_VOL_COLS = ("tissue-size-in-voxels", "tissue-size-in-mm", "relative-tissue-size-in-mm")
PVS_VOL_COLS = ("pvs-size-within-tissue-in-voxels", "pvs-size-within-tissue-in-mm",
                "relative-pvs-size-within-tissue-in-mm")
DIFFUSIVITY_UNITS = {"MD", "AD", "RD", "Dpar", "Dperp", "Dpar_minus_Dperp",
                     "lowb_ADC", "lowb_MD", "lowb_AD", "lowb_RD"}
PRIMARY_METRICS = ("MD", "AD", "RD", "FA", "MK_mf", "AK_mf", "RK_mf")
LOWB_METRICS = ("lowb_MD", "lowb_AD", "lowb_RD", "lowb_FA")
A_VOLUME_COLS = ("relative-tissue-size-in-mm", "relative-pvs-size-within-tissue-in-mm")
ABS_VOLUME_COLS = ("tissue-size-in-mm", "pvs-size-within-tissue-in-mm")   # v13
PLOT_ABSOLUTE = False                    # v13: set in main (--also-plot-absolute)
FAMILY_DESC = {
    "A": "group discrimination in the PVS regions (PVSROI): diffusion_within_pvs of "
         + ", ".join(PRIMARY_METRICS) + "; PVS volume / ROI volume; ROI volume / eTIV"
         + "; with --lowb-family merged also diffusion_within_pvs of " + ", ".join(LOWB_METRICS),
    "B": "low-b (Reviewer 3, R3.2): diffusion_within_pvs of " + ", ".join(LOWB_METRICS),
    "E": "PVS-weighted (secondary): PVS-Weighted-diffusion of " + ", ".join(PRIMARY_METRICS),
    "F": "anatomical atrophy: SynthSeg tissue volume / eTIV, every structure except eTIV",
}
ACTIVE_FAMILIES = ("A", "E", "F")        # set in main
LOWB_FAMILY = "merged"                   # v9: 'merged' (low-b in A) or 'separate' (B)
PAIRWISE_FDR = "union"                   # v7: 'union' or 'per-contrast'; set in main
U_KEYS = ("hc_vs_ppcs_u", "hc_vs_boxers_u", "ppcs_vs_boxers_u")
T_KEYS = ("hc_vs_ppcs_t", "hc_vs_boxers_t", "ppcs_vs_boxers_t")
A_REGIONS = set()                        # PVS regions = diffusion/PVSROI folders; set in main
DROPPED_CONST = []                       # (file, column) with no variance

# ---------------- v10: technical exclusions ----------------
PVS_VOXEL_COL = "pvs-size-within-tissue-in-voxels"
PLAUSIBILITY_DIFFUSIVITY = ("MD", "AD", "RD", "lowb_MD", "lowb_AD", "lowb_RD")   # mm2/s
PLAUSIBILITY_FA = ("FA", "lowb_FA")                                                # [0, 1]
EXCLUDE = {}                             # {PVS region: set(subject ids)}; set in main
PRETTY_GROUP = {'hc': 'HC', 'ppcs': 'PPCS', 'boxers': 'RHI'}
SIG_COLOR = "#1f5aa6"                    # v11: colour of significant effect-size labels; set in main
LARGE_DELTA = 0.474                      # v12: |delta| drawn in bold (Romano et al. 2006); set in main


def _read_group_sheets(path, why):
    """The three group sheets of one spreadsheet (named explicitly: without
    sheet_name, read_excel returns only the first sheet)."""
    if not os.path.isfile(path):
        raise SystemExit(f"{why}: file not found: {path}")
    try:
        return pd.read_excel(path, sheet_name=['hc', 'ppcs', 'boxers'])
    except Exception as e:
        raise SystemExit(f"{why}: cannot read sheets hc/ppcs/boxers of {path}: {e}")


def _values_by_subject(path, col, why):
    """{subject: value} for the INCLUDED subjects of their own group sheet."""
    out = {}
    for grp, df in _read_group_sheets(path, why).items():
        if col not in df.columns:
            raise SystemExit(f"{why}: {path}, sheet '{grp}' has no column '{col}'")
        for s, v in zip(df.iloc[:, 0].map(norm_id), df[col]):
            if s in INCLUDED[grp]:
                out[s] = _to_float_or_nan(v)
    return out


def _subject_key(s):
    d = re.sub(r"\D", "", str(s))
    return (int(d) if d else 0, str(s))


def build_value_exclusions(diff_root, vol_root, regions, min_vox, max_diff):
    """v10: one row per included subject x PVS region with its PVS voxel
    count, the checked within-PVS values, 'excluded' and the reasons.
    Uses data quality only -- no group statistic is computed or looked at.
    Returns (table, {region: set(excluded subjects)})."""
    group_of = {s: g for g, ids in INCLUDED.items() for s in ids}
    subjects = sorted(group_of, key=_subject_key)
    rows, excl = [], {}
    if min_vox == 0 and max_diff is None:
        return pd.DataFrame(rows), excl
    checked = PLAUSIBILITY_DIFFUSIVITY + PLAUSIBILITY_FA if max_diff is not None else ()
    for region in sorted(regions, key=_region_key):
        vox = {}
        if min_vox > 0:
            p = os.path.join(vol_root, region, "PVSROI", "VOL.xlsx")
            vox = _values_by_subject(p, PVS_VOXEL_COL, "PVS voxel count (criterion 1)")
        vals = {m: _values_by_subject(os.path.join(diff_root, region, "PVSROI", m + ".xlsx"),
                                      "diffusion_within_pvs", "plausibility (criterion 2)")
                for m in checked}
        for s in subjects:
            reasons, notes = [], []
            nv = vox.get(s, np.nan)
            if min_vox > 0:
                if s not in vox:
                    msg = f"{region}: {s} ({group_of[s]}) is not in the PVS voxel-count sheet"
                    if not ALLOW_MISSING:
                        raise SystemExit(msg + "  (use --allow-missing-subjects to continue)")
                    print("  [!] " + msg + " -- criterion 1 not checked")
                    notes.append("voxel count missing: criterion 1 not checked")
                elif not np.isfinite(nv):
                    print(f"  [!] {region}: {s} has an empty PVS voxel count -- criterion 1 not checked")
                    notes.append("voxel count empty: criterion 1 not checked")
                elif nv < min_vox:
                    reasons.append(f"{nv:.0f} PVS voxels < {min_vox}")
            for m in PLAUSIBILITY_DIFFUSIVITY if max_diff is not None else ():
                v = vals[m].get(s, np.nan)
                if np.isfinite(v) and not (0 <= v <= max_diff):
                    reasons.append(f"{m} = {v:.2e} mm2/s outside [0, {max_diff:.1e}]")
            for m in PLAUSIBILITY_FA if max_diff is not None else ():
                v = vals[m].get(s, np.nan)
                if np.isfinite(v) and not (0 <= v <= 1):
                    reasons.append(f"{m} = {v:.3f} outside [0, 1]")
            row = {"region": region, "group": group_of[s], "subject": s, "pvs_voxels": nv}
            row.update({m: vals[m].get(s, np.nan) for m in checked})
            row.update({"excluded": int(bool(reasons)), "reasons": "; ".join(reasons),
                        "notes": "; ".join(notes)})
            rows.append(row)
            if reasons:
                excl.setdefault(region, set()).add(s)
    return pd.DataFrame(rows), excl


def print_value_exclusions(tab, min_vox, max_diff):
    print("\nTechnical exclusions (within-PVS diffusion only; volumetry untouched):")
    print("  criterion 1: " + (f"fewer than {min_vox} PVS voxels in the region" if min_vox
                               else "off"))
    print("  criterion 2: " + (f"MD/AD/RD (full-range and low-b) outside [0, {max_diff:.1e}] "
                               "mm2/s, or FA/low-b FA outside [0, 1]" if max_diff is not None
                               else "off"))
    if tab.empty:
        print("  nothing checked")
        return
    ex = tab[tab.excluded == 1]
    for _, r in ex.iterrows():
        print(f"  {r.region:<24s} {PRETTY_GROUP[r.group]:<5s} {r.subject:<9s} {r.reasons}")
    per = ", ".join(f"{PRETTY_GROUP[g]}={int((ex.group == g).sum())}"
                    for g in ('hc', 'ppcs', 'boxers'))
    print(f"  excluded: {len(ex)} of {len(tab)} subject x region cells ({per})")


def is_absolute_volume(kind, seg, metric, col):
    """v13: an absolute volume in mm3 that --also-plot-absolute plots (no family)."""
    if kind != "volumetry" or metric != "VOL":
        return False
    if seg == "synthseg":
        return col == "tissue-size-in-mm"
    if seg == "PVSROI":
        return col in ABS_VOLUME_COLS
    return False


def family_of(kind, seg, metric, col, region=None):
    """Family letter of one test, or None if it belongs to no family."""
    if kind == "diffusion" and seg == "PVSROI":
        if col == "diffusion_within_pvs" and metric in PRIMARY_METRICS:
            return "A"
        if col == "diffusion_within_pvs" and metric in LOWB_METRICS:
            return "A" if LOWB_FAMILY == "merged" else "B"
        if col == "PVS-Weighted-diffusion" and metric in PRIMARY_METRICS:
            return "E"
        return None
    if kind == "volumetry" and seg == "PVSROI" and metric == "VOL":
        # only the PVS regions themselves (not eTIV / whole brain)
        return "A" if (col in A_VOLUME_COLS and region in A_REGIONS) else None
    if kind == "volumetry" and seg == "synthseg" and metric == "VOL":
        return "F" if (col == "relative-tissue-size-in-mm" and region != "eTIV") else None
    return None


GROUP_COLOR = {'hc': 'darkkhaki', 'ppcs': 'darkgoldenrod', 'boxers': 'maroon'}
INCREASE_IS_PATHOLOGICAL = {"MD", "AD", "RD", "ADC", "MK_mf", "AK_mf", "RK_mf", "Dperp"}
DECREASE_IS_PATHOLOGICAL = {"FA", "Dpar"}


def pathology_direction(kind, metric, col, region):
    """'increase' / 'decrease' (which change is the pathological sign) or None.
    Decided from the metric and the column -- never from the Case ID."""
    if kind == "diffusion":
        m = metric[5:] if metric.startswith("lowb_") else metric
        if m in INCREASE_IS_PATHOLOGICAL:
            return "increase"
        if m in DECREASE_IS_PATHOLOGICAL:
            return "decrease"
        return None
    if col.startswith(("pvs-size", "relative-pvs")):
        return "increase"                                   # PVS enlargement
    if col.startswith(("tissue-size", "relative-tissue")):
        r = str(region).lower()
        return "increase" if ("ventricle" in r or "csf" in r) else "decrease"
    return None


def file_could_hold_family(kind, seg, metric, region=None):
    """Cheap pre-filter on the path (no file opened)."""
    if kind == "diffusion":
        return seg == "PVSROI" and (metric in PRIMARY_METRICS or metric in LOWB_METRICS)
    if metric != "VOL":
        return False
    if seg == "PVSROI":
        return ("A" in ACTIVE_FAMILIES and region in A_REGIONS) or PLOT_ABSOLUTE
    if seg == "synthseg":
        return ("F" in ACTIVE_FAMILIES and region != "eTIV") or PLOT_ABSOLUTE
    return False


# ---------------- v6: per-family Kruskal-Wallis heatmaps ----------------
REGION_ORDER = ("centrum", "cerebral", "basal", "midbrain")
# PVSROI folder name -> name used in the paper (matched on lower case, '_' or ' ')
REGION_LABEL = {"cerebral white matter": "Centrum semiovale"}
COL_LABEL = {"MD": "MD", "AD": "AD", "RD": "RD", "FA": "FA", "MK_mf": "MK",
             "AK_mf": "AK", "RK_mf": "RK", "lowb_MD": "low-b MD", "lowb_AD": "low-b AD",
             "lowb_RD": "low-b RD", "lowb_FA": "low-b FA"}


def _pretty_region(r):
    name = str(r).replace("_", " ").replace("-", " ").strip()
    return REGION_LABEL.get(name.lower(), name.capitalize())


def _region_key(r):
    rl = str(r).lower()
    for i, k in enumerate(REGION_ORDER):
        if k in rl:
            return (i, rl)
    return (len(REGION_ORDER), rl)


def _heat_columns(fam, rows):
    """(column key, label, block) in display order for one family."""
    if fam == "A":
        cols = [("VOL|relative-tissue-size-in-mm", "ROI volume / eTIV", "ROI\nvolumetry"),
                ("VOL|relative-pvs-size-within-tissue-in-mm", "PVS volume / ROI", "PVS\nvolumetry")]
        cols += [(f"{m}|diffusion_within_pvs", f"{COL_LABEL[m]} within PVS", "Diffusivity\nwithin PVS")
                 for m in ("MD", "AD", "RD", "FA")]
        cols += [(f"{m}|diffusion_within_pvs", f"{COL_LABEL[m]} within PVS",
                  "Diffusion kurtosis\nwithin PVS") for m in ("MK_mf", "AK_mf", "RK_mf")]
        cols += [(f"{m}|diffusion_within_pvs", f"{COL_LABEL[m]} within PVS",
                  "Low-b fit\nwithin PVS") for m in LOWB_METRICS]          # v9
    elif fam == "B":
        cols = [(f"{m}|diffusion_within_pvs", f"{COL_LABEL[m]} within PVS", "Low-b fit")
                for m in LOWB_METRICS]
    elif fam == "E":
        cols = [(f"{m}|PVS-Weighted-diffusion", f"PVS-weighted {COL_LABEL[m]}", "PVS-weighted")
                for m in PRIMARY_METRICS]
    else:
        cols = [("VOL|relative-tissue-size-in-mm", "Volume / eTIV", "SynthSeg")]
    present = {f"{r['Metric']}|{r['Column']}" for r in rows}
    return [c for c in cols if c[0] in present]


def plot_family_heatmap(kw_rows, fam, out_dir, vmax=0.2):
    """Kruskal-Wallis heatmap of one family: cells show q (within the family)
    and H; colour = q on [0, vmax] like submitted Figure 2."""
    rows = [r for r in kw_rows if r.get("fdr_family") == fam]
    if not rows:
        return None
    cols = _heat_columns(fam, rows)
    regions = sorted({r["ROI"] for r in rows}, key=_region_key)
    Q = np.full((len(regions), len(cols)), np.nan)
    H = np.full_like(Q, np.nan)
    for r in rows:
        key = f"{r['Metric']}|{r['Column']}"
        ci = [c[0] for c in cols].index(key) if key in [c[0] for c in cols] else None
        if ci is None:
            continue
        ri = regions.index(r["ROI"])
        Q[ri, ci] = r.get("p_value_FDR", np.nan)
        H[ri, ci] = r.get("H_stat", np.nan) if r.get("H_stat") is not None else np.nan
    transpose = fam == "F"                      # many structures: one column, rows = structures
    data, hh = (Q, H)
    ylabels = [_pretty_region(r) for r in regions]
    xlabels = [c[1] for c in cols]
    fw = max(6.0, 1.15 * len(cols) + 3.5) if not transpose else 6.5
    fh = max(3.2, 0.9 * len(regions) + 2.6) if not transpose else max(4.0, 0.42 * len(regions) + 2.0)
    fig, ax = plt.subplots(figsize=(fw, fh))
    im = ax.imshow(np.ma.masked_invalid(data), cmap="YlGnBu_r", vmin=0, vmax=vmax, aspect="auto")
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if not np.isfinite(data[i, j]):
                continue
            txt = f"{data[i, j]:.3f}" + (f"\n(H={hh[i, j]:.1f})" if np.isfinite(hh[i, j]) else "")
            ax.text(j, i, txt, ha="center", va="center", fontsize=9,
                    color="white" if data[i, j] < 0.6 * vmax else "0.15")
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels, rotation=40, ha="right", fontsize=10)
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels, fontsize=11)
    if not transpose:                           # block separators, as in the submitted figure
        blocks = [c[2] for c in cols]
        for j in range(1, len(blocks)):
            if blocks[j] != blocks[j - 1]:
                ax.axvline(j - 0.5, color="k", lw=2)
        starts = [0] + [j for j in range(1, len(blocks)) if blocks[j] != blocks[j - 1]]
        ends = starts[1:] + [len(blocks)]
        for s0, e0 in zip(starts, ends):
            ax.text((s0 + e0 - 1) / 2, -0.6, blocks[s0], ha="center", va="bottom",
                    fontsize=9.5, fontweight="bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("q-value (FDR within family %s)" % fam, fontsize=10)
    fig.suptitle(f"Family {fam}: Kruskal-Wallis across HC, PPCS, RHI (cell: q and H)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.90 if not transpose else 0.95))
    stem = os.path.join(out_dir, f"heatmap_family_{fam}_KruskalWallis")
    fig.savefig(stem + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(stem + ".pdf", bbox_inches="tight")
    plt.close(fig)
    return stem + ".png"


def norm_id(v):
    s = str(v).strip()
    if re.fullmatch(r"\d+(\.0+)?", s):
        return f"sub-{int(float(s))}"
    return s


def load_inclusion(study_xlsx, col="Will we used it?"):
    """{data_sheet_group: set(ids)} of subjects marked 'Yes' in the study xlsx."""
    xls = pd.ExcelFile(study_xlsx)
    inc = {}
    for sheet, grp in STUDY_SHEETS.items():
        if sheet not in xls.sheet_names:
            continue
        df = xls.parse(sheet)
        subj = "Subj" if "Subj" in df.columns else df.columns[0]
        if col not in df.columns:
            raise SystemExit(f"'{col}' missing in sheet {sheet} of {study_xlsx}")
        keep = df[df[col].astype(str).str.strip().str.lower() == "yes"][subj]
        inc.setdefault(grp, set()).update(
            norm_id(v) for v in keep.dropna() if str(v).strip().lower() not in ("total:", ""))
    if set(inc) != {"hc", "ppcs", "boxers"}:
        raise SystemExit(f"study xlsx gave groups {sorted(inc)}; need hc, ppcs, boxers")
    return inc


def cliffs_delta_ci(x, y, n_boot, rng, alpha=0.05, chunk=1000):
    """Percentile bootstrap CI of Cliff's delta (subjects resampled within
    group). delta = mean(sign(x_i - y_j)), identical to 2U/(n1 n2) - 1."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    x = x[np.isfinite(x)]; y = y[np.isfinite(y)]
    if x.size < 3 or y.size < 3 or n_boot <= 0:
        return np.nan, np.nan
    out = []
    for start in range(0, n_boot, chunk):
        b = min(chunk, n_boot - start)
        xb = x[rng.integers(0, x.size, (b, x.size))]
        yb = y[rng.integers(0, y.size, (b, y.size))]
        out.append(np.sign(xb[:, :, None] - yb[:, None, :]).mean(axis=(1, 2)))
    d = np.concatenate(out)
    return float(np.percentile(d, 100 * alpha / 2)), float(np.percentile(d, 100 * (1 - alpha / 2)))


def fmt_p(p, name="p"):
    if p is None or not np.isfinite(p):
        return f"{name} = n/a"
    if p < 0.001:
        return f"{name} < 0.001"
    return f"{name} = {p:.3f}" if p < 0.1 else f"{name} = {p:.2f}"


def annotate_significance(ax, x1, x2, y, pval, line_height=0.1,
                          fontsize_asterisks=24, delta=None, delta_color='k',
                          fontsize_delta=None):
    if pval is None or (isinstance(pval, float) and np.isnan(pval)):
        return

    # Determine stars
    if pval < 0.001:   stars = '***'
    elif pval < 0.01:  stars = '**'
    elif pval < 0.05:  stars = '*'
    else:              stars = '' 

    if fontsize_delta is None:
        fontsize_delta = max(10, int(fontsize_asterisks * 0.55))

    y_line = y + line_height
    ax.plot([x1, x2], [y_line, y_line], lw=1.5, c='k')

    if stars:
        ax.text((x1 + x2) * 0.5, y_line,
                stars, ha='center', va='bottom',
                fontsize=fontsize_asterisks, fontweight='bold', color='k')

    if delta is not None and np.isfinite(delta):
        # Determine Arrow Direction (Assuming Calc = Left - Right)
        # If delta is NEGATIVE, Right Group is BIGGER -> Increase -> UP Arrow
        if delta < 0:
            arrow = "↑"  
        # If delta is POSITIVE, Left Group is BIGGER -> Decrease -> DOWN Arrow
        elif delta > 0:
            arrow = "↓"
        else:
            arrow = ""

        # Construct label
        label = f"Δ={delta:+.2f} {arrow}"
        
        ax.text((x1 + x2) * 0.5, y_line - line_height * 0.35,
                label,
                ha='center', va='top',
                fontsize=fontsize_delta, fontweight='bold', color=delta_color)

def annotate_significance_by_chatgpt(ax, x1, x2, y, pval, line_height=0.1,
                          fontsize_asterisks=24,
                          delta=None, delta_possible_color=['k'],
                          fontsize_delta=None,
                          g1=None, g2=None,
                          show_comp_label=True,
                          case_now=None,
                          show_delta_text=True,
                          p_raw=None, ci=None, in_family=True, family=None,
                          direction=None):
    # v7: the q label names the family and, under union control, the contrasts
    """
    Draws a significance bar. 
    Updates: Only draws the bracket if there is actually something (Stars or Text) to show.
    """
    # --- Determine Stars ---
    stars = ''
    if pval is not None and not (isinstance(pval, float) and np.isnan(pval)):
        if pval < 0.001:   stars = '***'
        elif pval < 0.01:  stars = '**'
        elif pval < 0.05:  stars = '*'
    
    is_delta_valid = (delta is not None and np.isfinite(delta))
    
    # --- FIX START: STRICTER CHECK ---
    # We only proceed if there are stars OR if we are effectively showing the delta text.
    # If there are no stars AND (delta is invalid OR we've hidden the text), we stop.
    should_show_delta = (is_delta_valid and show_delta_text)
    
    if not stars and not should_show_delta:
        return
    # --- FIX END ---

    if fontsize_delta is None:
        fontsize_delta = max(10, int(fontsize_asterisks * 0.55))

    # ---------- Enforce canonical comparison direction ----------
    swapped = False
    if (g1 is not None) and (g2 is not None):
        a, b = g1, g2
        pair = {a, b}
        def want_first(pair_set):
            if pair_set == {'boxers', 'hc'}:   return 'boxers'
            if pair_set == {'boxers', 'ppcs'}: return 'boxers'
            if pair_set == {'ppcs', 'hc'}:     return 'ppcs'
            return None 
        wf = want_first(pair)
        if wf is not None and a != wf:
            x1, x2 = x2, x1
            g1, g2 = g2, g1
            swapped = True

    if swapped and is_delta_valid:
        delta = -delta
    if swapped and ci is not None:                      # v2: CI follows delta
        ci = (-ci[1], -ci[0])

    # ---------- Colour (v11): one colour, only when significant ----------
    # significant = the star criterion (q within the family < 0.05), so the
    # colour never disagrees with the stars; 'direction' is no longer used here
    delta_color = SIG_COLOR if stars else 'k'
    # ---------- Weight (v12): bold only for a large effect ----------
    # magnitude, independent of significance; on the unrounded delta
    delta_weight = 'bold' if (is_delta_valid and abs(delta) >= LARGE_DELTA) else 'normal'

    # ---------- Draw Bracket ----------
    y_line = y + line_height
    ax.plot([x1, x2], [y_line, y_line], lw=1.5, c='k')

    # ---------- Prepare Label ----------
    label = ""
    pretty = {'hc': 'HC', 'ppcs': 'PPCS', 'boxers': 'RHI'}   # v2: was defined
                                                            # only inside an if
    # We use the same flag here to ensure safety
    if should_show_delta: 
        if delta > 0: arrow = "↑"
        elif delta < 0: arrow = "↓"
        else: arrow = ""

        comp = ""
        if show_comp_label and (g1 is not None) and (g2 is not None):
            if delta > 0: comp = f"{pretty.get(g1,g1)} > {pretty.get(g2,g2)} ("
            elif delta < 0: comp = f"{pretty.get(g1,g1)} < {pretty.get(g2,g2)} ("
            else: comp = f"{pretty.get(g1,g1)} == {pretty.get(g2,g2)} ("

        if abs(delta) < 0.005:
            delta_str = "0.00"
            if "==" not in comp: comp = f"{pretty.get(g1,g1)} == {pretty.get(g2,g2)} ("
        else:
            delta_str = f"{delta:+.2f}"

        ci_str = ""
        if ci is not None and np.all(np.isfinite(ci)):
            ci_str = f" [{ci[0]:+.2f}, {ci[1]:+.2f}]"
        # v8: effect size, CI and stars only (p and q are in the Excel output)
        label = f"{comp}Δ={delta_str}{ci_str})" + (f" {stars}" if stars else "")

    # ---------- PLOT TEXT & STARS ----------
    
    text_y_pos = y_line + (line_height * 1.5) 
    
    if label:
        ax.text((x1 + x2) * 0.5, text_y_pos,
                label,
                ha='center', va='bottom',
                fontsize=fontsize_delta, fontweight=delta_weight, color=delta_color)

    if stars and not label:        # v2: with a label the stars are inline
        # If text is hidden (empty label), star_offset is small
        star_offset = (line_height * 9.0) if label else (line_height * 1.5)
        stars_y_pos = text_y_pos + star_offset

        ax.text((x1 + x2) * 0.5, stars_y_pos,
                stars, ha='center', va='bottom',
                fontsize=fontsize_asterisks, fontweight='bold', color='k')
 
def get_data_and_metadata(case, xlsx_path):
    """
    Helper to read data and determine metadata (ROI, metric, etc.) 
    Used by both the Stats pass and the Plotting pass to ensure consistency.
    """
    # v3: both trees have the layout <tissue>/<segmentation>/<metric>.xlsx,
    # so tissue, segmentation and metric are parsed the same way for both
    # (v2 gave diffusion files ROI = map = metric = '_').
    parts = xlsx_path.split("/")
    ROI = parts[-3]
    map_name = parts[-2]
    metric = parts[-1].split(".xlsx")[0]
    base = ROI
    kind = "volumetry" if case == "PVS-volumetry" else "diffusion"

    try:
        dfs = pd.read_excel(xlsx_path, sheet_name=['hc', 'ppcs', 'boxers'])
    except Exception as e:
        print(f"[!] could not read {xlsx_path}: {e}")
        return None, None, None

    # v2: keep exactly the subjects marked 'Yes' in the study spreadsheet
    # (v1 dropped a hard-coded set that disagreed with the inclusion list)
    first = xlsx_path not in _CHECKED_FILES
    all_inc = set().union(*INCLUDED.values())
    for grp, df in dfs.items():
        ids = df.iloc[:, 0].map(norm_id)
        wrong = sorted(set(ids) & (all_inc - INCLUDED[grp]))
        if wrong:
            raise SystemExit(f"{xlsx_path}: sheet '{grp}' contains subjects that the "
                             f"study xlsx assigns to ANOTHER group: {wrong}")
        keep = ids.isin(INCLUDED[grp])
        removed = sorted(set(ids[~keep]))
        missing = sorted(INCLUDED[grp] - set(ids[keep]))
        dfs[grp] = df[keep.values]
        if first:
            N_REPORT.append((xlsx_path, grp, len(df), int(keep.sum()),
                             ",".join(removed), ",".join(missing)))
            if missing:
                msg = (f"{xlsx_path}: sheet '{grp}' lacks included subjects {missing}")
                if not ALLOW_MISSING:
                    raise SystemExit(msg + "  (use --allow-missing-subjects to continue)")
                print("  [!] " + msg)
    _CHECKED_FILES.add(xlsx_path)

    # v10: technical exclusions (subject x region), within-PVS diffusion only
    if kind == "diffusion" and map_name == "PVSROI" and EXCLUDE.get(ROI):
        for grp, df in dfs.items():
            drop = df.iloc[:, 0].map(norm_id).isin(EXCLUDE[ROI])
            dfs[grp] = df[~drop.values]

    # find numeric columns
    common = set(dfs['hc'].columns) & set(dfs['ppcs'].columns) & set(dfs['boxers'].columns)
    numeric_cols = [c for c in dfs['hc'].columns[1:]
                    if c in common and pd.api.types.is_numeric_dtype(dfs['hc'][c])]
    
    # v5: keep only the columns that belong to an active family, then drop
    # columns with no variance across the included subjects (cannot be tested)
    col_family = {c: family_of(kind, map_name, metric, c, ROI) for c in numeric_cols}
    if PLOT_ABSOLUTE:                       # v13: plotted, in no family
        for c in numeric_cols:
            if col_family[c] is None and is_absolute_volume(kind, map_name, metric, c):
                col_family[c] = "ABS"
    numeric_cols = [c for c in numeric_cols
                    if col_family[c] in ACTIVE_FAMILIES or col_family[c] == "ABS"]
    keep_cols = []
    for c in numeric_cols:
        vals = pd.concat([dfs[g][c] for g in ('hc', 'ppcs', 'boxers')]).dropna()
        if vals.nunique() <= 1:
            if first:
                DROPPED_CONST.append((xlsx_path, c))
            continue
        keep_cols.append(c)
    numeric_cols = keep_cols
    metadata = {
        'ROI': ROI, 'map': map_name, 'metric': metric, 'kind': kind,
        'col_family': col_family, 'base': base, 'numeric_cols': numeric_cols
    }
    return dfs, metadata, numeric_cols

# --- PASS 1: STATS CALCULATION ---

def compute_stats_only(case, xlsx_path, results):
    """
    Reads file, calculates all stats, appends to global 'results' dict.
    Does NOT plot.
    """
    dfs, meta, numeric_cols = get_data_and_metadata(case, xlsx_path)
    if dfs is None: return

    group_order = ['hc', 'ppcs', 'boxers']

    for col in numeric_cols:
        # Construct the unique Case ID exactly as before
        case_now = "%s_%s_%s_%s_%s" % (meta['ROI'], meta['map'], meta['metric'], meta['base'], col)
        
        data = {g: dfs[g][col].dropna().values for g in group_order}
        if any(len(v) == 0 for v in data.values()):
            continue

        # 1) ANOVA
        try:
            F_stat, p_f = f_oneway(data['hc'], data['ppcs'], data['boxers'])
            df1 = len(group_order)-1
            df2 = sum(len(v) for v in data.values()) - len(group_order)
            F_str = f"F({df1},{df2}) = {F_stat:.2f}"
        except:
            F_stat, p_f, F_str = None, None, ''

        # Store ANOVA
        row_base = {'Metric': meta['metric'], 'Case': case_now, 'Map': meta['map'],
                    'ROI': meta['ROI'], 'Column': col,
                    'fdr_family': meta['col_family'][col],
                    # v10: group sizes actually tested (exclusions can differ by region)
                    'n_hc': len(data['hc']), 'n_ppcs': len(data['ppcs']),
                    'n_boxers': len(data['boxers'])}
        
        row = row_base.copy()
        row.update({'F_stat': F_stat, 'p_value': p_f, 'F': F_str})
        results['ANOVA'].append(row)

        # 2) Kruskal-Wallis
        try:
            H_stat, p_kw = kruskal(data['hc'], data['ppcs'], data['boxers'])
        except:
            H_stat, p_kw = None, None

        row = row_base.copy()
        row.update({'H_stat': H_stat, 'p_value': p_kw})
        results['Kruskal-Wallis'].append(row)

        # 3) Pairwise
        pairs = [('hc', 'ppcs'), ('hc', 'boxers'), ('ppcs', 'boxers')]
        for g1, g2 in pairs:
            # Mann-Whitney U & Cliff's Delta
            # We no longer use a bare 'try/except' here to avoid hiding bugs.
            # The safety check is now inside the function.
            ustat, p_u, delta_cl = cliffs_delta_from_u(data[g1], data[g2], alternative_now='two-sided')

            # Optional: Debug print if values are missing
            if np.isnan(p_u) or p_u is None:
                # This helps you see in the terminal which comparisons are failing
                print(f"[Warn] No data for {case_now}: {g1} vs {g2}")

            # v2: bootstrap CI; seed derived from the case so it is reproducible
            # regardless of file order
            rng = np.random.default_rng(
                [SEED, zlib.crc32(f"{case_now}|{g1}|{g2}".encode())])
            lo, hi = cliffs_delta_ci(data[g1], data[g2], N_BOOT, rng)
            row_u = row_base.copy()
            row_u.update({'u_stat': ustat, 'p_value': p_u, 'delta_cl': delta_cl,
                          'delta_ci_lo': lo, 'delta_ci_hi': hi,
                          'n1': len(data[g1]), 'n2': len(data[g2])})
            results[f"{g1}_vs_{g2}_u"].append(row_u)           

            # T-test (Apply similar safety logic if needed)
            try:
                # ttest_ind also fails on empty data
                if len(data[g1]) < 2 or len(data[g2]) < 2:
                    tstat, p_t, t_str = np.nan, np.nan, ''
                else:
                    tstat, p_t = ttest_ind(data[g1], data[g2], equal_var=False)
                    df_t = len(data[g1]) + len(data[g2]) - 2
                    t_str = f"t({df_t}) = {tstat:.2f}"
            except Exception as e:
                tstat, p_t, t_str = np.nan, np.nan, ''

            row_t = row_base.copy()
            row_t.update({'t_stat': tstat, 'p_value': p_t, 't': t_str})
            results[f"{g1}_vs_{g2}_t"].append(row_t)

# --- PASS 2: PLOTTING ---

def generate_plot_only(case, xlsx_path, output_dir, q_value_lookup, delta_lookup,
                       fontsize_title=14, fontsize_rest=12, fig_size1=6, fig_size2=4,
                       fontsize_asterisks=24, row_lookup=None):
    """
    Reads file, generates plot. 
    Uses 'q_value_lookup' to decide stars (FDR Corrected) instead of raw p-values.
    """
    dfs, meta, numeric_cols = get_data_and_metadata(case, xlsx_path)
    if dfs is None: return

    group_order = ['hc', 'ppcs', 'boxers']
    group_colors = {'hc': 'darkkhaki', 'ppcs': 'darkgoldenrod', 'boxers': 'maroon'}

    for col in numeric_cols:
        case_now = "%s_%s_%s_%s_%s" % (meta['ROI'], meta['map'], meta['metric'], meta['base'], col)
        
        data = {g: dfs[g][col].dropna().values for g in group_order}
        if any(len(v) == 0 for v in data.values()):
            continue

        # Build plot
        plot_df = pd.concat([pd.DataFrame({'group': g, 'value': data[g]})
                             for g in group_order], ignore_index=True)
        fig, ax = plt.subplots(figsize=(fig_size1, fig_size2))
        sns.violinplot(x='group', y='value', data=plot_df,
                       order=group_order,
                       palette=[group_colors[g] for g in group_order],
                       inner=None, ax=ax)
        sns.swarmplot(x='group', y='value', data=plot_df,
                      order=group_order,
                      color='yellow', size=4, alpha=0.7, ax=ax)

        # Annotate Pairwise using Q-VALUES
        y_max = plot_df['value'].max()
        # Increased from 0.2 to 0.35 to accommodate stacked text/stars
        gap = (plot_df['value'].max() - plot_df['value'].min()) * 0.38   # v8: one-line labels
        pairs = [('hc', 'ppcs'), ('hc', 'boxers'), ('ppcs', 'boxers')]
        
        for i, (g1, g2) in enumerate(pairs):
            # RETRIEVE THE CORRECTED Q-VALUE
            # We look into the dictionary we built in main()
            # Key structure: q_lookup[case_now][f"{g1}_vs_{g2}_u"]
            
            q_val = None
            # --- q-value lookup (Mann–Whitney corrected) ---
            try:
                q_val = q_value_lookup[case_now][f"{g1}_vs_{g2}_u"]
            except KeyError:
                q_val = 1.0

            # --- delta lookup (raw Cliff's delta) ---
            delta_val = np.nan
            try:
                    delta_val = delta_lookup[case_now][f"{g1}_vs_{g2}_u"]
            except KeyError:
                pass
            delta_val = _to_float_or_nan(delta_val)

            # v2: raw p, CI and family membership for the label
            r = (row_lookup or {}).get(case_now, {}).get(f"{g1}_vs_{g2}_u", {})
            p_raw = _to_float_or_nan(r.get('p_value'))
            ci = (_to_float_or_nan(r.get('delta_ci_lo')), _to_float_or_nan(r.get('delta_ci_hi')))
            in_family = bool(r.get('in_FDR_family', 1))
            fam_letter = r.get('fdr_family', None)

            # Optional: color delta by direction
            delta_txt_color = 'k'
            delta_possible_txt_color = ['k']        # v8: colour now set in the annotator
            if np.isfinite(delta_val):
                if i == 0:  # hc vs ppcs
                        delta_possible_txt_color = ['k','darkkhaki','darkgoldenrod']
                elif i == 1:  # hc vs boxers
                        delta_possible_txt_color = ['k','darkkhaki','maroon']
                else:        # ppcs vs boxers
                        delta_possible_txt_color = ['k','darkgoldenrod','maroon']


            # ... existing code ...
            
            # --- MODIFIED BLOCK START ---
            if q_val is not None:
                y = y_max + gap * (i + 1)
                
                # Check if this is an eTIV case
                # v3: only the eTIV VOLUME itself (not PVS / eTIV) is treated as eTIV
                is_etiv = ((meta['ROI'] == "eTIV" or "EstimatedTotalIntraCranial" in meta['base'])
                           and col.startswith("tissue-size"))

                annotate_significance_by_chatgpt(
                ax,
                group_order.index(g1),
                group_order.index(g2),
                y=y,
                pval=q_val,
                
                # IMPORTANT: Keep this at 0.05 to prevent bracket overlap!
                line_height=gap * 0.05,  
                
                fontsize_asterisks=fontsize_asterisks,
                delta=delta_val,
                delta_possible_color=delta_possible_txt_color,
                g1=g1, g2=g2,
                show_comp_label=True,
                case_now=case_now,
                
                # NEW: Hide text if it is eTIV, otherwise show it
                show_delta_text=not is_etiv,
                p_raw=p_raw, ci=ci, in_family=in_family, family=fam_letter,
                direction=pathology_direction(meta['kind'], meta['metric'], col, meta['ROI'])
                )
            # --- MODIFIED BLOCK END ---

        # v2: keep all three 2-line labels inside the axes
        y_lo = ax.get_ylim()[0]
        ax.set_ylim(y_lo, y_max + gap * 3.7)

        # --- Formatting (Same as original) ---
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0)) 
        ax.yaxis.set_major_formatter(formatter)
        ax.yaxis.get_offset_text().set_fontsize(fontsize_rest)
        ax.yaxis.get_offset_text().set_fontweight('bold')

        if (meta['ROI'] == "eTIV" or "EstimatedTotalIntraCranial" in meta['base']) \
                and col == "tissue-size-in-mm":            # v3: the eTIV volume only
            ax.set_ylim(0.5e6, 2.5e6)
        
        ax.set_xlabel('')
        
        # Y-Label Logic (v3: from the column and the file kind, not the Case ID)
        unit = " (mm$^2$/s)" if meta['metric'] in DIFFUSIVITY_UNITS else ""
        if meta['kind'] == "volumetry":
            what = "PVS" if col.startswith(("pvs", "relative-pvs")) else \
                   ("eTIV" if meta['ROI'] == "eTIV" else "Tissue")
            if col.startswith("relative-pvs"):
                y_label = "PVS volume / eTIV" if meta['ROI'] == "eTIV" \
                          else "PVS volume / ROI volume"
            elif col.startswith("relative-tissue"):
                y_label = "Tissue volume / eTIV"
            elif col.endswith("voxels"):
                y_label = f"{what} volume (#voxels)"
            else:
                y_label = f"{what} volume (mm$^3$)"
        else:
            m = COL_LABEL.get(meta['metric'], meta['metric'])        # v9: readable
            m = m[0].upper() + m[1:] if m.startswith("low-b") else m
            y_label = {"diffusion_within_pvs": f"{m} within PVS{unit}",
                       "diffusion_within_tissue": f"{m} in ROI{unit}",
                       "PVS-Weighted-diffusion": f"PVS-weighted {m}{unit}",
                       "Dex-tissue": f"Dex-tissue ({m})",
                       "Dex-ratio": f"Dex-ratio ({m})"}.get(col, f"{m} {col}")

        if "Normalised" in y_label or " / " in y_label: fontsize_ylabel = fontsize_rest - 5
        else: fontsize_ylabel = fontsize_rest
        if len(y_label.replace("$", "").replace("^", "")) > 22:      # v9: long labels fit
            fontsize_ylabel = min(fontsize_ylabel, fontsize_rest - 8)
        
        ax.set_ylabel(y_label, fontsize=fontsize_ylabel, fontweight='bold')
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(['HC', 'PPCS', 'RHI'], fontsize=fontsize_rest, fontweight='bold')
        ax.tick_params(axis='y', labelsize=fontsize_rest)
        ax.tick_params(axis='both', which='major', length=8, width=4, labelsize=fontsize_rest, labelcolor='black', pad=6)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False) 
        ax.spines['bottom'].set_linewidth(4)
        ax.spines['left'].set_linewidth(4)
        fig.tight_layout()
        for lbl in ax.get_yticklabels(): lbl.set_fontweight('bold')

        # Save
        out_png = os.path.join(output_dir, f"{case_now}.png")
        fig.savefig(out_png, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"→ saved plot {out_png}")


# --- MAIN ---

def _max_diffusivity(v):
    """v10: upper diffusivity limit in mm2/s, or 'off'."""
    if str(v).strip().lower() == "off":
        return None
    try:
        x = float(v)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{v}' is neither a number nor 'off'")
    if not (1e-4 <= x <= 1e-2):
        raise argparse.ArgumentTypeError(
            f"{x:g} mm2/s is implausible as a limit -- the unit is mm2/s "
            "(free water at 37 C is about 3.0e-3)")
    return x


def _min_voxels(v):
    try:
        n = int(v)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{v}' is not an integer")
    if n < 0:
        raise argparse.ArgumentTypeError("must be >= 0 (0 = criterion off)")
    return n


def _colour(v):
    """v11: any matplotlib colour (name or hex)."""
    import matplotlib.colors as mcolors
    if not mcolors.is_color_like(v):
        raise argparse.ArgumentTypeError(f"'{v}' is not a matplotlib colour")
    return v


def _abs_delta(v):
    """v12: a |Cliff's delta| threshold in (0, 1]."""
    try:
        x = float(v)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{v}' is not a number")
    if not (0 < x <= 1):
        raise argparse.ArgumentTypeError("Cliff's |delta| lies in (0, 1]")
    return x


def main():
    global INCLUDED, ALLOW_MISSING, N_BOOT, SEED, ACTIVE_FAMILIES, A_REGIONS, PAIRWISE_FDR
    global EXCLUDE, SIG_COLOR, LARGE_DELTA, PLOT_ABSOLUTE
    ap = argparse.ArgumentParser(description="PVS group statistics and violin plots "
                                             "(v6: FDR families A, B, E, F).")
    ap.add_argument("--root", default=cfg("results_dir"))
    ap.add_argument("--etiv-strategy", default="tissues-PVS_VOL-DIFF_nifti-based_and_SPM12-eTIV-based")
    ap.add_argument("--pvs-strategy", default="conservative")
    ap.add_argument("--basedon", default="based_on_basic_info_created_by_Evgenios_updated_on_"
                    "2026-08-28_PVS_segmentation_xlsx")
    ap.add_argument("--study-xlsx", default=cfg("study_xlsx"),
                    help="inclusion list: subjects marked 'Yes' are analysed")
    ap.add_argument("--families", nargs="+", choices=sorted(FAMILY_DESC),
                    default=["A", "E", "F"],
                    help="families to compute and plot (default: all four)")
    ap.add_argument("--lowb-family", choices=["merged", "separate"], default="merged",
                    help="v9: low-b metrics in family A (default) or as family B (v8)")
    ap.add_argument("--pairwise-fdr", choices=["union", "per-contrast"], default="union",
                    help="v7: 'union' (default) = FDR over the three pairwise contrasts "
                         "together within each family; 'per-contrast' = v6")
    ap.add_argument("--list-tests", action="store_true",
                    help="list every test per family and the group sizes, compute NO "
                         "statistics, and exit")
    ap.add_argument("--allow-incomplete-families", action="store_true",
                    help="continue when a family lacks a metric for some region")
    ap.add_argument("--n-boot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--allow-missing-subjects", action="store_true")
    ap.add_argument("--min-pvs-voxels", type=_min_voxels, required=True,
                    help="v10, criterion 1 (required): exclude a subject's within-PVS "
                         "diffusion in a region whose PVS mask has fewer voxels than this "
                         "(column '%s' of volumetry/<region>/PVSROI/VOL.xlsx); 0 = off"
                         % PVS_VOXEL_COL)
    ap.add_argument("--max-diffusivity", type=_max_diffusivity, required=True,
                    help="v10, criterion 2 (required): upper limit in mm2/s for the within-PVS "
                         "MD, AD, RD (full-range and low-b); also excludes FA / low-b FA "
                         "outside [0, 1]. E.g. 3.0e-3 (free water at 37 C); 'off' = off")
    ap.add_argument("--significant-colour", type=_colour, default=SIG_COLOR,
                    help="v11: colour of the effect-size label of a pairwise contrast with "
                         "q < 0.05 (same criterion as the stars); others are black")
    ap.add_argument("--bold-min-delta", type=_abs_delta, default=LARGE_DELTA,
                    help="v12: effect-size labels with |Cliff's delta| >= this are bold, others "
                         "regular (default 0.474 = large effect, Romano et al. 2006)")
    ap.add_argument("--also-plot-absolute", action="store_true",
                    help="v13: also plot the absolute volumes in mm3 (SynthSeg structures incl. "
                         "eTIV, PVS-region and PVS volumes) with delta and CI only; they join "
                         "no FDR family, so every test and q-value is unchanged")
    ap.add_argument("--out-suffix", default="_v13",
                    help="appended to the plot output folder so earlier plots are kept")
    a = ap.parse_args()
    SIG_COLOR = a.significant_colour
    PLOT_ABSOLUTE = a.also_plot_absolute
    LARGE_DELTA = a.bold_min_delta
    ALLOW_MISSING, N_BOOT, SEED = a.allow_missing_subjects, a.n_boot, a.seed
    ACTIVE_FAMILIES = tuple(sorted(set(a.families)))
    global LOWB_FAMILY
    LOWB_FAMILY = a.lowb_family
    if LOWB_FAMILY == "merged" and "B" in ACTIVE_FAMILIES:
        raise SystemExit("--families B is empty with --lowb-family merged (low-b is in A)")
    PAIRWISE_FDR = a.pairwise_fdr

    fontsize_title = 44
    fontsize_rest = 40
    asterisks_size = 40
    fig_size1 = 12
    fig_size2 = 9
    root = a.root
    eTIV_strategy, PVS_strategy = a.etiv_strategy, a.pvs_strategy
    basedon = a.basedon
    trees = {"PVS-diffusion": "diffusion", "PVS-volumetry": "volumetry"}

    INCLUDED = load_inclusion(a.study_xlsx)
    print(f"study xlsx: {a.study_xlsx}")
    print("  included per group: " + ", ".join(f"{g}={len(v)}" for g, v in INCLUDED.items()))
    print("FDR families run: " + ", ".join(ACTIVE_FAMILIES))
    for f in ACTIVE_FAMILIES:
        print(f"  {f}: {FAMILY_DESC[f]}")

    # 1. collect the spreadsheets that can hold a test of an active family
    seg_root = os.path.join(root, "plots_and_stats/%s/%s_PVS_segmentation%s"
                            % (eTIV_strategy, PVS_strategy, a.out_suffix))
    files_to_process, seen = [], []
    # the PVS regions are the region folders that hold diffusion/PVSROI spreadsheets
    diff_root = os.path.join(root, "spreadsheets", eTIV_strategy, PVS_strategy, "diffusion", basedon)
    if {"A", "B", "E"} & set(ACTIVE_FAMILIES):
        if not os.path.isdir(diff_root):
            raise SystemExit(f"diffusion tree not found: {diff_root}")
        A_REGIONS = {d for d in os.listdir(diff_root)
                     if os.path.isdir(os.path.join(diff_root, d, "PVSROI"))}
        print("PVS regions (diffusion/PVSROI): " + ", ".join(sorted(A_REGIONS, key=_region_key)))
    # v10: technical exclusions, from data quality only (before any statistic)
    excl_df = pd.DataFrame()
    if {"A", "B", "E"} & set(ACTIVE_FAMILIES):
        vol_root = os.path.join(root, "spreadsheets", eTIV_strategy, PVS_strategy,
                                "volumetry", basedon)
        excl_df, EXCLUDE = build_value_exclusions(diff_root, vol_root, A_REGIONS,
                                                  a.min_pvs_voxels, a.max_diffusivity)
        print_value_exclusions(excl_df, a.min_pvs_voxels, a.max_diffusivity)
    for case, kind in trees.items():
        if kind == "volumetry" and not ({"A", "F"} & set(ACTIVE_FAMILIES)):
            continue
        if kind == "diffusion" and not ({"A", "B", "E"} & set(ACTIVE_FAMILIES)):
            continue
        in_root = os.path.join(root, "spreadsheets", eTIV_strategy, PVS_strategy, kind, basedon)
        if not os.path.isdir(in_root):
            raise SystemExit(f"{kind} tree not found: {in_root}")
        print(f"Scanning: {in_root}")
        for dp, _, files in sorted(os.walk(in_root)):
            rel = os.path.relpath(dp, in_root)
            parts = rel.split(os.sep)
            for fn in sorted(files):
                if not fn.lower().endswith(".xlsx") or len(parts) != 2:
                    continue
                metric = os.path.splitext(fn)[0]
                ok = file_could_hold_family(kind, parts[1], metric, parts[0])
                seen.append({"file": os.path.relpath(os.path.join(dp, fn), root), "tree": kind,
                             "region": parts[0], "segmentation": parts[1], "metric": metric,
                             "used": ok})
                if ok:
                    files_to_process.append((case, os.path.join(dp, fn),
                                             os.path.join(seg_root, case, basedon, rel)))
    seen_df = pd.DataFrame(seen)
    if not files_to_process:
        raise SystemExit("no spreadsheet matches the active families")
    final_stats_dir = os.path.join(seg_root, "stats_families_" + "".join(ACTIVE_FAMILIES), basedon)
    os.makedirs(final_stats_dir, exist_ok=True)
    seen_df.to_csv(os.path.join(final_stats_dir, "file_selection.csv"), index=False)
    excl_df.to_csv(os.path.join(final_stats_dir, "value_exclusions.csv"), index=False)

    # 2. enumerate the tests per family (opens the files, no statistics) and
    #    check that every family is complete: each metric x each region
    members = []
    for case_type, fpath, _ in files_to_process:
        dfs, meta, cols = get_data_and_metadata(case_type, fpath)
        if dfs is None:
            continue
        for col in cols:
            members.append({"family": meta['col_family'][col], "region": meta['ROI'],
                            "metric": meta['metric'], "column": col,
                            "Case": "%s_%s_%s_%s_%s" % (meta['ROI'], meta['map'], meta['metric'],
                                                        meta['base'], col)})
    M = pd.DataFrame(members)
    if M.empty:
        raise SystemExit("no tests found for the active families")
    problems = []
    for f in ACTIVE_FAMILIES:
        Mf = M[M.family == f]
        if Mf.empty:
            problems.append(f"family {f}: no tests found")
            continue
        if f == "F":
            continue                       # one test per SynthSeg structure, no fixed set
        expected = {"A": [(m, "diffusion_within_pvs") for m in PRIMARY_METRICS]
                         + [("VOL", c) for c in A_VOLUME_COLS]
                         + ([(m, "diffusion_within_pvs") for m in LOWB_METRICS]
                            if LOWB_FAMILY == "merged" else []),
                    "B": [(m, "diffusion_within_pvs") for m in LOWB_METRICS],
                    "E": [(m, "PVS-Weighted-diffusion") for m in PRIMARY_METRICS]}[f]
        regions = sorted(A_REGIONS, key=_region_key)
        have = set(zip(Mf.region, Mf.metric, Mf.column))
        missing = [(r, m, c) for r in regions for m, c in expected if (r, m, c) not in have]
        if missing:
            problems.append(f"family {f}: missing {missing}")
    print("\nFDR family sizes (tests per pairwise comparison):")
    for f in ACTIVE_FAMILIES:
        Mf = M[M.family == f]
        regs = sorted(Mf.region.unique())
        if f == "F":
            print(f"  {f}: {len(Mf):3d} tests = 1 per SynthSeg structure")
        else:
            nper = Mf.groupby("region").size()
            print(f"  {f}: {len(Mf):3d} tests = {int(nper.max())} per region x {len(regs)} regions "
                  f"({', '.join(sorted(regs, key=_region_key))})")
    if PLOT_ABSOLUTE:
        print(f"  absolute volumes plotted, in no family: {int((M.family == 'ABS').sum())} "
              "(fdr_family 'ABS'; p only, no q, no stars)")
    if DROPPED_CONST:
        print(f"  dropped (no variance): {len(DROPPED_CONST)} columns")
    if problems:
        msg = "INCOMPLETE FAMILIES:\n  " + "\n  ".join(problems)
        if not a.allow_incomplete_families:
            raise SystemExit(msg + "\n(--allow-incomplete-families to continue)")
        print("  [!] " + msg)
    M.to_csv(os.path.join(final_stats_dir, "fdr_family_members.csv"), index=False)
    pd.DataFrame(N_REPORT, columns=["file", "group", "n_rows_in_sheet", "n_kept",
                                    "removed", "missing"]).to_csv(
        os.path.join(final_stats_dir, "group_sizes_per_file.csv"), index=False)
    if a.list_tests:
        print(f"\nlisted -> {os.path.join(final_stats_dir, 'fdr_family_members.csv')}")
        print(f"exclusions -> {os.path.join(final_stats_dir, 'value_exclusions.csv')}")
        return

    for _, _, odir in files_to_process:
        os.makedirs(odir, exist_ok=True)

    # --- PASS 1: statistics ---
    print(f"\n--- PASS 1: statistics for {len(files_to_process)} files "
          f"(bootstrap {N_BOOT}, seed {SEED}) ---")
    results = {
        'ANOVA': [], 'Kruskal-Wallis': [],
        'hc_vs_ppcs_t': [], 'hc_vs_boxers_t': [], 'ppcs_vs_boxers_t': [],
        'hc_vs_ppcs_u': [], 'hc_vs_boxers_u': [], 'ppcs_vs_boxers_u': []
    }
    for case_type, fpath, odir in files_to_process:
        compute_stats_only(case_type, fpath, results)

    # --- FDR (v7) ---
    #   omnibus sheets (KW, ANOVA): within each family
    #   pairwise sheets: within each family, over the three contrasts together
    #   ('union', default) or per contrast ('per-contrast' = v6)
    print(f"\n--- FDR (Benjamini-Hochberg) within each family; pairwise: {PAIRWISE_FDR} ---")
    for rows in results.values():
        for r in rows:
            r['in_FDR_family'] = 0 if r.get('fdr_family') == "ABS" else 1   # v13
            r['p_value_FDR'] = np.nan
            r['Significant_FDR'] = 0
            r['fdr_scope'] = ""

    def _valid(r):
        p = r.get('p_value')
        return p is not None and np.isfinite(p)

    def _bh(cells, scope):
        """cells: list of (sheet, index); BH over their p-values."""
        if not cells:
            return
        reject, q, _, _ = multipletests([results[k][i]['p_value'] for k, i in cells],
                                        alpha=0.05, method='fdr_bh')
        for (k, i), qq, sig in zip(cells, q, reject):
            results[k][i]['p_value_FDR'] = qq
            results[k][i]['Significant_FDR'] = int(sig)
            results[k][i]['fdr_scope'] = scope

    if PAIRWISE_FDR == "union":
        groups = [("Kruskal-Wallis",), ("ANOVA",), U_KEYS, T_KEYS]
    else:
        groups = [("Kruskal-Wallis",), ("ANOVA",)] + [(k,) for k in U_KEYS + T_KEYS]
    fam_sizes = {}
    for f in ACTIVE_FAMILIES:
        for keys in groups:
            cells = [(k, i) for k in keys for i, r in enumerate(results[k])
                     if r.get('fdr_family') == f and _valid(r)]
            scope = f"{f}: " + ("+".join(keys) if len(keys) > 1 else keys[0])
            _bh(cells, scope)
            fam_sizes[(keys, f)] = len(cells)
    for f in ACTIVE_FAMILIES:
        if PAIRWISE_FDR == "union":
            print(f"  {f}: Kruskal-Wallis={fam_sizes[(('Kruskal-Wallis',), f)]}, "
                  f"Mann-Whitney (3 contrasts)={fam_sizes[(U_KEYS, f)]}, "
                  f"t-test (3 contrasts)={fam_sizes[(T_KEYS, f)]}")
        else:
            print(f"  {f}: Kruskal-Wallis={fam_sizes[(('Kruskal-Wallis',), f)]}, " +
                  ", ".join(f"{k}={fam_sizes[((k,), f)]}" for k in U_KEYS))

    # --- v6: one Kruskal-Wallis heatmap per family ---
    for f in ACTIVE_FAMILIES:
        hm = plot_family_heatmap(results['Kruskal-Wallis'], f, final_stats_dir)
        if hm:
            print(f"  heatmap {f} -> {hm}")

    # --- lookups for plotting ---
    pairwise_keys = ['hc_vs_ppcs_u', 'hc_vs_boxers_u', 'ppcs_vs_boxers_u',
                     'hc_vs_ppcs_t', 'hc_vs_boxers_t', 'ppcs_vs_boxers_t']
    q_lookup, row_lookup, delta_lookup = {}, {}, {}
    for key in pairwise_keys:
        for row in results[key]:
            q_lookup.setdefault(row['Case'], {})[key] = row.get('p_value_FDR', np.nan)
            row_lookup.setdefault(row['Case'], {})[key] = row
    for key in ['hc_vs_ppcs_u', 'hc_vs_boxers_u', 'ppcs_vs_boxers_u']:
        for row in results[key]:
            delta_lookup.setdefault(row['Case'], {})[key] = row.get('delta_cl', np.nan)

    # --- PASS 2: plots ---
    print("\n--- PASS 2: plots (stars from q within the family; p, q and CI printed) ---")
    for case_type, fpath, odir in files_to_process:
        generate_plot_only(case_type, fpath, odir, q_lookup, delta_lookup,
                           fontsize_title, fontsize_rest, fig_size1, fig_size2,
                           fontsize_asterisks=asterisks_size, row_lookup=row_lookup)

    # --- Excel ---
    excel_path = os.path.join(final_stats_dir, "stats_all_corrected_v13.xlsx")
    run_info = pd.DataFrame([
        ("script", os.path.basename(__file__)),
        ("date", datetime.datetime.now().isoformat(timespec="seconds")),
        ("study_xlsx", a.study_xlsx),
        ("included", ", ".join(f"{g}={len(v)}" for g, v in INCLUDED.items())),
        ("families_run", ", ".join(ACTIVE_FAMILIES)),
        *[(f"family_{f}", FAMILY_DESC[f]) for f in ACTIVE_FAMILIES],
        ("pairwise_fdr", PAIRWISE_FDR),
        ("family_sizes_kruskal_wallis",
         ", ".join(f"{f}={fam_sizes[(('Kruskal-Wallis',), f)]}" for f in ACTIVE_FAMILIES)),
        ("family_sizes_mann_whitney",
         ", ".join(f"{f}={fam_sizes[(U_KEYS, f)] if PAIRWISE_FDR == 'union' else sum(fam_sizes[((k,), f)] for k in U_KEYS)}"
                   for f in ACTIVE_FAMILIES)
         + (" (pooled over the 3 contrasts)" if PAIRWISE_FDR == "union" else " (corrected per contrast)")),
        ("fdr_method", "Benjamini-Hochberg within each family; omnibus tests per family; "
                       + ("pairwise tests over the three contrasts together"
                          if PAIRWISE_FDR == "union" else "pairwise tests per contrast")),
        ("lowb_family", LOWB_FAMILY),
        ("volumes", "relative only (ROI / eTIV, PVS / ROI, structure / eTIV)"
                    + ("; absolute mm3 volumes plotted without any test (fdr_family ABS)"
                       if PLOT_ABSOLUTE else "")),
        ("pvs_regions", ", ".join(sorted(A_REGIONS, key=_region_key))),
        ("ci", f"percentile bootstrap of Cliff's delta, {N_BOOT} resamples, seed {SEED}"),
        ("stars", "from q within the family: * <0.05, ** <0.01, *** <0.001"),
        ("label_colour", f"{SIG_COLOR} when q < 0.05 (same criterion as the stars), else black"),
        ("label_weight", f"bold when |delta| >= {LARGE_DELTA} (large effect, Romano et al. 2006), "
                         "else regular; independent of significance"),
        ("allow_incomplete_families", a.allow_incomplete_families),
        ("min_pvs_voxels", a.min_pvs_voxels if a.min_pvs_voxels else "off"),
        ("max_diffusivity_mm2_s", "off" if a.max_diffusivity is None
         else f"{a.max_diffusivity:.2e} (MD/AD/RD full-range and low-b); FA, low-b FA in [0, 1]"),
        ("exclusion_scope", "subject x PVS region; all within-PVS diffusion columns "
                            "(diffusion/<region>/PVSROI); volumetry never excluded"),
        ("excluded_cells", "0" if excl_df.empty else
         f"{int(excl_df.excluded.sum())} of {len(excl_df)} subject x region cells (" +
         ", ".join(f"{g}={int(((excl_df.group == g) & (excl_df.excluded == 1)).sum())}"
                   for g in ('hc', 'ppcs', 'boxers')) + ")"),
        ("allow_missing_subjects", ALLOW_MISSING)], columns=["key", "value"])
    with pd.ExcelWriter(excel_path) as writer:
        run_info.to_excel(writer, sheet_name="run_info", index=False)
        M.to_excel(writer, sheet_name="fdr_families", index=False)
        seen_df.to_excel(writer, sheet_name="file_selection", index=False)
        pd.DataFrame(N_REPORT, columns=["file", "group", "n_rows_in_sheet", "n_kept",
                                        "removed", "missing"]).to_excel(
            writer, sheet_name="group_sizes", index=False)
        pd.DataFrame(DROPPED_CONST, columns=["file", "column"]).to_excel(
            writer, sheet_name="dropped_no_variance", index=False)
        excl_df.to_excel(writer, sheet_name="value_exclusions", index=False)
        for sheet, rows in results.items():
            pd.DataFrame(rows).to_excel(writer, sheet_name=sheet, index=False)
    print(f"\n✎ wrote all stats to {excel_path}")


if __name__=="__main__":
    main()