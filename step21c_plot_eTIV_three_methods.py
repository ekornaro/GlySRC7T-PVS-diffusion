#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step21c_plot_eTIV_three_methods.py
==================================
Supplementary figure: estimated total intracranial volume (eTIV) of HC, PPCS
and RHI by three methods, one violin panel per method, with group statistics.

Sources (per subject), exactly as in the spreadsheet-creation scripts
----------------------------------------------------------------------
  SPM12 tissue summation   TIV_spm12_results.xlsx, sheets HC / PPCS / Boxers,
                           column TIV_ml (ml -> mm3 x 1000)
  SynthStrip (FreeSurfer)  <seg-root>/<sub>/Segmentation/<sub>_ses-01_ICV_synthstrip.nii.gz
                           non-zero voxels x voxel volume (0.512 mm3)
                           [create_csv_..._synthsegICV_2026-01-21_median.py]
  FreeSurfer recon-all     <seg-root>/<sub>/FreeSurfer/synthseg/volumetry_scores.csv,
                           first row, second column (mm3)
                           [create_csv_..._FreesurferICV_2025-11-21_median.py]

Cohort: subjects marked 'Yes' in --study-xlsx (sheets HC, PPCS, Boxers), the
same inclusion list as the statistics scripts. A subject missing from any
source stops the run (--allow-missing-subjects to drop it with a warning).
Every value is checked to lie in --plausible (default 0.8e6-2.5e6 mm3) so a
unit mix-up (ml, voxels) is caught before plotting.

Statistics per method: Kruskal-Wallis across the three groups; Mann-Whitney U
per pair with Cliff's delta; Benjamini-Hochberg FDR over the three pairs
(--correction fdr, default) or none. Brackets with stars (* <0.05, ** <0.01,
*** <0.001) are drawn only for pairs significant after correction.

Outputs: <out>.png, <out>.pdf, <out>_stats.csv, <out>_values.csv (one row per
subject with the three eTIVs; the console prints their agreement).

Usage
-----
  python step21c_plot_eTIV_three_methods.py --out eTIV_three_methods.png
  (defaults point at the GlySRC7T paths; override with --study-xlsx,
   --spm-xlsx, --seg-root)
"""
import argparse
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import ScalarFormatter
from scipy.stats import kruskal, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests
from pipeline_config import cfg

GROUPS = ("HC", "PPCS", "RHI")
STUDY_SHEETS = {"HC": "HC", "PPCS": "PPCS", "Boxers": "RHI"}
COLORS = {"HC": "darkkhaki", "PPCS": "darkgoldenrod", "RHI": "maroon"}
PAIRS = (("HC", "PPCS"), ("HC", "RHI"), ("PPCS", "RHI"))
VOXEL_MM3 = 0.512                     # 0.8 mm isotropic, as in the creation scripts
METHODS = (                           # (key, panel title)
    ("spm12", "eTIV by $\\it{Tissue\\ Summation}$\n(SPM12)"),
    ("synthstrip", "eTIV by $\\it{mri\\_synthstrip}$\n(FreeSurfer SynthStrip)"),
    ("freesurfer", "eTIV by $\\it{recon}$-$\\it{all}$\n(FreeSurfer)"),
)
DEFAULT_STUDY = cfg("study_xlsx")
DEFAULT_SPM = cfg("tiv_spm12_xlsx")
DEFAULT_SEG_ROOT = cfg("pvs_root")


# ---------------------------------------------------------------- cohort ----
def norm_id(v):
    s = str(v).strip()
    if re.fullmatch(r"\d+(\.0+)?", s):
        return f"sub-{int(float(s))}"
    m = re.search(r"(sub-[A-Za-z0-9]+)", s)
    return m.group(1) if m else s


def load_inclusion(study_xlsx, col="Will we used it?"):
    """{group label: [sub-ids]} of subjects marked 'Yes'."""
    if not os.path.isfile(study_xlsx):
        raise SystemExit(f"study xlsx not found: {study_xlsx}")
    book = pd.ExcelFile(study_xlsx)
    inc = {}
    for sheet, grp in STUDY_SHEETS.items():
        if sheet not in book.sheet_names:
            raise SystemExit(f"{study_xlsx}: sheet '{sheet}' missing (has {book.sheet_names})")
        d = book.parse(sheet)
        subj = "Subj" if "Subj" in d.columns else d.columns[0]
        if col not in d.columns:
            raise SystemExit(f"{study_xlsx}, sheet {sheet}: no column '{col}'")
        keep = d[col].astype(str).str.strip().str.lower() == "yes"
        ids = [norm_id(v) for v in d.loc[keep, subj].dropna()
               if str(v).strip().lower() not in ("total:", "")]
        inc[grp] = sorted(set(ids), key=lambda s: int(re.sub(r"\D", "", s) or 0))
    seen = {}
    for g, ids in inc.items():
        for s in ids:
            if s in seen:
                raise SystemExit(f"{s} is listed under both {seen[s]} and {g}")
            seen[s] = g
    return inc


# ---------------------------------------------------------------- sources ---
def etiv_spm12(spm_xlsx):
    """{sub-id: eTIV mm3} from TIV_spm12_results.xlsx (all group sheets)."""
    if not os.path.isfile(spm_xlsx):
        raise SystemExit(f"SPM12 xlsx not found: {spm_xlsx}")
    book = pd.ExcelFile(spm_xlsx)
    out = {}
    for sheet in STUDY_SHEETS:
        if sheet not in book.sheet_names:
            raise SystemExit(f"{spm_xlsx}: sheet '{sheet}' missing (has {book.sheet_names})")
        d = book.parse(sheet)
        for c in ("subject", "TIV_ml"):
            if c not in d.columns:
                raise SystemExit(f"{spm_xlsx}, sheet {sheet}: no column '{c}'")
        d = d[d["subject"].notna()]
        for s, v in zip(d["subject"].map(norm_id), pd.to_numeric(d["TIV_ml"], errors="coerce")):
            if s in out:
                raise SystemExit(f"{spm_xlsx}: {s} appears twice")
            out[s] = float(v) * 1000.0 if np.isfinite(v) else np.nan
    return out


VOXEL_REPORT = []                     # (subject, zooms, voxel mm3) where the header is not 0.512


def etiv_synthstrip(seg_root, sub, voxel_mode="header"):
    """Non-zero voxels of the SynthStrip ICV mask x voxel volume (mm3).
    voxel_mode 'header' (default) uses each file's own voxel size; 'fixed'
    uses 0.512 mm3 for everyone, as the old creation script did."""
    import nibabel as nib
    p = os.path.join(seg_root, sub, "Segmentation", f"{sub}_ses-01_ICV_synthstrip.nii.gz")
    if not os.path.isfile(p):
        return np.nan, p
    img = nib.load(p)
    zooms = tuple(float(z) for z in img.header.get_zooms()[:3])
    vox = float(np.prod(zooms))
    if not np.isclose(vox, VOXEL_MM3, rtol=0.01):
        VOXEL_REPORT.append((sub, zooms, vox))
    used = vox if voxel_mode == "header" else VOXEL_MM3
    return float(np.count_nonzero(img.get_fdata())) * used, p


def etiv_freesurfer(fs_root, sub):
    """First row, second column of <fs_root>/<sub>/FreeSurfer/synthseg/volumetry_scores.csv (mm3)."""
    p = os.path.join(fs_root, sub, "FreeSurfer", "synthseg", "volumetry_scores.csv")
    if not os.path.isfile(p):
        return np.nan, p
    d = pd.read_csv(p)
    if d.shape[1] < 2 or d.shape[0] < 1:
        raise SystemExit(f"{p}: expected at least one row and two columns, got {d.shape}")
    v = pd.to_numeric(pd.Series([d.iloc[0, 1]]), errors="coerce").iloc[0]
    if not np.isfinite(v):
        raise SystemExit(f"{p}: first row, second column is not a number: {d.iloc[0, 1]!r}")
    return float(v), p


def collect(inc, spm_xlsx, seg_root, fs_root, plausible, allow_missing, voxel_mode="header"):
    """Long table: subject, group, method, value_mm3 (dropping nothing silently)."""
    spm = etiv_spm12(spm_xlsx)
    rows, problems = [], []
    for grp, ids in inc.items():
        for sub in ids:
            vals = {"spm12": (spm.get(sub, np.nan), spm_xlsx)}
            vals["synthstrip"] = etiv_synthstrip(seg_root, sub, voxel_mode)
            vals["freesurfer"] = etiv_freesurfer(fs_root, sub)
            for m, (v, src) in vals.items():
                if not np.isfinite(v):
                    why = "file not found" if not os.path.isfile(src) else "empty value"
                    problems.append(f"{sub} ({grp}) {m}: {why}: {src}")
                    continue
                if not (plausible[0] <= v <= plausible[1]):
                    raise SystemExit(f"{sub} {m}: eTIV {v:.4g} mm3 outside the plausible range "
                                     f"[{plausible[0]:.3g}, {plausible[1]:.3g}] -- units? ({src})")
                rows.append({"subject": sub, "group": grp, "method": m, "value_mm3": v})
    if problems:
        by_method = {m: sum(p.split(" ")[2] == m + ":" for p in problems) for m, _ in METHODS}
        n_subj = sum(len(v) for v in inc.values())
        msg = "missing values:\n  " + "\n  ".join(problems)
        for m, n in by_method.items():
            if n == n_subj:
                msg += (f"\n  -> EVERY subject lacks a {m} value: the root for that source is "
                        f"probably wrong (--seg-root / --freesurfer-root)")
        if not allow_missing:
            raise SystemExit(msg + "\n(--allow-missing-subjects to continue without them)")
        print("[!] " + msg)
    return pd.DataFrame(rows)


# ------------------------------------------------------------- statistics ---
def _stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def group_stats(df_method, method, correction):
    v = {g: df_method.loc[df_method.group == g, "value_mm3"].to_numpy() for g in GROUPS}
    for g in GROUPS:
        if len(v[g]) < 2:
            raise SystemExit(f"{method}: group {g} has {len(v[g])} value(s)")
    H, p_kw = kruskal(*(v[g] for g in GROUPS))
    rows = [{"method": method, "test": "Kruskal-Wallis", "groups": "HC, PPCS, RHI",
             "n": ", ".join(str(len(v[g])) for g in GROUPS), "statistic": H, "p": p_kw,
             "q": np.nan, "cliffs_delta": np.nan,
             "median_mm3": ", ".join(f"{np.median(v[g]):.4g}" for g in GROUPS)}]
    for a, b in PAIRS:
        U, p = mannwhitneyu(v[a], v[b], alternative="two-sided")
        rows.append({"method": method, "test": "Mann-Whitney U", "groups": f"{a} vs {b}",
                     "n": f"{len(v[a])}, {len(v[b])}", "statistic": U, "p": p, "q": np.nan,
                     "cliffs_delta": 2 * U / (len(v[a]) * len(v[b])) - 1,   # >0: a larger
                     "median_mm3": ""})
    st = pd.DataFrame(rows)
    pw = st.test == "Mann-Whitney U"
    st.loc[pw, "q"] = multipletests(st.loc[pw, "p"], method="fdr_bh")[1]
    st["p_for_stars"] = np.where(pw, st["q"] if correction == "fdr" else st["p"], st["p"])
    return st


# ------------------------------------------------------------------ plot ----
def violin_panel(ax, df_method, st, title, ylim, show_ylabel):
    sns.violinplot(data=df_method, x="group", y="value_mm3", hue="group", order=GROUPS,
                   hue_order=GROUPS, palette=COLORS, inner=None, linewidth=0.8,
                   legend=False, ax=ax)
    sns.swarmplot(data=df_method, x="group", y="value_mm3", order=GROUPS, color="yellow",
                  size=2.5, alpha=0.8, ax=ax)
    lo, hi = ax.dataLim.y0, ax.dataLim.y1
    if ylim is None:
        ylim = (lo - 0.03 * (hi - lo), hi + 0.03 * (hi - lo))
    step = 0.09 * (ylim[1] - ylim[0])
    y = hi + 0.4 * step
    drawn = []
    for a, b in PAIRS:
        s = _stars(float(st.loc[st.groups == f"{a} vs {b}", "p_for_stars"].iloc[0]))
        if not s:
            continue
        x1, x2 = GROUPS.index(a), GROUPS.index(b)
        ax.plot([x1, x2], [y, y], color="k", lw=0.9)
        ax.text((x1 + x2) / 2, y, s, ha="center", va="bottom", fontsize=13, fontweight="bold")
        y += step
        drawn.append(f"{a} vs {b}")
    ax.set_ylim(ylim[0], max(ylim[1], y))
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((0, 0))
    ax.yaxis.set_major_formatter(fmt)
    ax.yaxis.get_offset_text().set_fontsize(10)
    ax.yaxis.get_offset_text().set_fontweight("bold")
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Volume (mm$^3$)" if show_ylabel else "", fontsize=12, fontweight="bold")
    ax.set_xticks(range(len(GROUPS)))
    ax.set_xticklabels(GROUPS, fontsize=12, fontweight="bold")
    ax.tick_params(axis="y", labelsize=10)
    for t in ax.get_yticklabels():
        t.set_fontweight("bold")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(1.5)
    return drawn, ylim


def plot_three(df, stats, out_png, ylim):
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.8))
    drawn, lims = {}, []
    for ax, (m, title) in zip(axes, METHODS):
        drawn[m], lim = violin_panel(ax, df[df.method == m], stats[stats.method == m],
                                     title, ylim, show_ylabel=ax is axes[0])
        lims.append(ax.get_ylim())
    if ylim is None:                           # auto: one common range for the three
        common = (min(l[0] for l in lims), max(l[1] for l in lims))
        for ax in axes:
            ax.set_ylim(common)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(os.path.splitext(out_png)[0] + ".pdf", bbox_inches="tight")
    plt.close(fig)
    return drawn


# ------------------------------------------------------------------ main ----
def _ylim(v):
    if v.strip().lower() == "auto":
        return None
    try:
        a, b = (float(x) for x in v.split(","))
    except ValueError:
        raise argparse.ArgumentTypeError("use 'auto' or 'min,max' in mm3, e.g. 0.5e6,2.5e6")
    if not a < b:
        raise argparse.ArgumentTypeError("min must be smaller than max")
    return (a, b)


def main():
    ap = argparse.ArgumentParser(description="eTIV by SPM12, SynthStrip and FreeSurfer, by group.")
    ap.add_argument("--study-xlsx", default=DEFAULT_STUDY, help="inclusion list ('Yes' rows)")
    ap.add_argument("--spm-xlsx", default=DEFAULT_SPM, help="TIV_spm12_results.xlsx")
    ap.add_argument("--seg-root", default=DEFAULT_SEG_ROOT,
                    help="folder holding sub-*/Segmentation/<sub>_ses-01_ICV_synthstrip.nii.gz")
    ap.add_argument("--freesurfer-root", default=None,
                    help="folder holding sub-*/FreeSurfer/synthseg/volumetry_scores.csv "
                         "(default: same as --seg-root)")
    ap.add_argument("--out", required=True, help="output PNG (PDF and CSVs written next to it)")
    ap.add_argument("--correction", choices=["fdr", "none"], default="fdr",
                    help="correction of the 3 pairwise tests within each method")
    ap.add_argument("--ylim", type=_ylim, default=(0.5e6, 2.5e6),
                    help="'min,max' in mm3 (default 0.5e6,2.5e6) or 'auto'")
    ap.add_argument("--plausible", type=_ylim, default=(0.8e6, 2.5e6),
                    help="eTIV range in mm3 outside which the run stops (units check)")
    ap.add_argument("--allow-missing-subjects", action="store_true")
    ap.add_argument("--synthstrip-voxel", choices=["header", "fixed"], default="header",
                    help="voxel volume for the SynthStrip eTIV: each mask's header (default) "
                         "or 0.512 mm3 for everyone, as the old creation script")
    a = ap.parse_args()
    if a.plausible is None:
        raise SystemExit("--plausible needs 'min,max'")

    inc = load_inclusion(a.study_xlsx)
    print("included: " + ", ".join(f"{g}={len(v)}" for g, v in inc.items()))
    fs_root = a.freesurfer_root or a.seg_root
    print(f"SynthStrip masks under: {a.seg_root}\nFreeSurfer csv under:   {fs_root}")
    df = collect(inc, a.spm_xlsx, a.seg_root, fs_root, a.plausible, a.allow_missing_subjects,
                 a.synthstrip_voxel)
    if VOXEL_REPORT:
        print(f"\n[!] SynthStrip masks whose voxel volume is not {VOXEL_MM3} mm3 "
              f"(the old creation script assumed {VOXEL_MM3} for everyone):")
        for sub, z, v in VOXEL_REPORT:
            print(f"    {sub}: voxel {z[0]:.3f} x {z[1]:.3f} x {z[2]:.3f} mm = {v:.4f} mm3 "
                  f"-> old eTIV off by {100 * (VOXEL_MM3 / v - 1):+.1f}%")
        print("    using " + ("each header's voxel volume" if a.synthstrip_voxel == "header"
                              else f"{VOXEL_MM3} mm3 (--synthstrip-voxel fixed)") + "\n")

    stats = pd.concat([group_stats(df[df.method == m], m, a.correction) for m, _ in METHODS],
                      ignore_index=True)
    with pd.option_context("display.float_format", "{:.4g}".format, "display.width", 160):
        print(stats.drop(columns=["p_for_stars"]).to_string(index=False))

    wide = df.pivot(index="subject", columns="method", values="value_mm3")
    wide.insert(0, "group", df.drop_duplicates("subject").set_index("subject").loc[wide.index, "group"])
    print("\nagreement between methods (subjects with all three):")
    full = wide.dropna()
    for i, (m1, _) in enumerate(METHODS):
        for m2, _ in METHODS[i + 1:]:
            r = spearmanr(full[m1], full[m2]).correlation
            print(f"  {m1:>10s} vs {m2:<10s}: Spearman r = {r:.2f}, "
                  f"median ratio {m2}/{m1} = {np.median(full[m2] / full[m1]):.3f}")

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    drawn = plot_three(df, stats, a.out, a.ylim)
    for m, _ in METHODS:
        print(f"brackets {m}: " + (", ".join(drawn[m]) if drawn[m] else "none"))
    stem = os.path.splitext(a.out)[0]
    stats.assign(correction=a.correction, study_xlsx=a.study_xlsx, spm_xlsx=a.spm_xlsx,
                 seg_root=a.seg_root, freesurfer_root=fs_root).to_csv(stem + "_stats.csv", index=False)
    wide.to_csv(stem + "_values.csv")
    print(f"-> {a.out}, {stem}.pdf, {stem}_stats.csv, {stem}_values.csv")


if __name__ == "__main__":
    main()
