#!/usr/bin/env python3
"""
PSD volumes: group statistics and violin plots (HC vs PPCS vs RHI), restricted to the participants
listed ("Yes") in the study-info spreadsheet.

What it does
------------
1. Reads the gathered PSD table (one row per FLAIR label map; columns: Scan, Stratum, total / frontal /
   parietal / occipital PSD volume in cm3). The prefrontal column, if still present, is ignored; the total
   is expected to be frontal + parietal + occipital (checked).
2. Keeps ONLY the participants marked "Yes" in column "Will we used it?" of the sheets HC / PPCS / Boxers
   of the study-info spreadsheet. The group of a participant is the study-file sheet it is listed in
   (the CSV 'Stratum' column must agree, otherwise the script stops).
3. Per measure (total, frontal, parietal, occipital):
      - descriptives per group (n, median, IQR, mean, SD, range)
      - Kruskal-Wallis across the three groups
      - pairwise two-sided Mann-Whitney U tests: PPCS vs HC, RHI vs HC, RHI vs PPCS
      - Cliff's delta with a percentile bootstrap CI (positive = FIRST-named group larger),
        magnitude labelled with the Romano thresholds (0.147 / 0.33 / 0.474)
4. Benjamini-Hochberg FDR.  Default family = all pairwise Mann-Whitney tests together
   (4 measures x 3 contrasts = 12), which is the family behind q = 0.78 for the parietal contrasts in
   the manuscript; Kruskal-Wallis P values are corrected across the 4 measures.
   --fdr-family per-contrast corrects each contrast separately across the 4 measures instead.
5. Violin + swarm plots (one PNG per measure + one combined figure).  Following the figure conventions:
   every contrast carries a bracket with "RHI > HC: delta = +0.15 [-0.23, +0.51]" (Cliff's delta and its
   bootstrap 95% CI, bold if |delta| >= 0.474); asterisks only when q < 0.05; y-axis starting near the data
   minimum; bold axis titles.  --significant-only annotates only contrasts with q < 0.05.

Dry-run by default: the statistics, the participant audit and the n per group are printed and nothing is
written.  Add --apply to write the workbook, the filtered data table, the figures and a provenance JSON.

Example
-------
python3 step24_psd_group_stats_and_plots.py \
    --psd-csv   "$(python3 pipeline_config.py psd_csv)" \
    --study-xlsx "$(python3 pipeline_config.py study_xlsx)" \
    --out-dir   ./PSD_study-filtered \
    --apply

Importable: load_inclusion(), load_psd_table(), select_participants(), run_statistics(), make_figures().
"""
import argparse
import hashlib
import json
import os
import re
import sys
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

GROUPS = ("hc", "ppcs", "boxers")                       # internal keys (order of the x-axis)
PRETTY = {"hc": "HC", "ppcs": "PPCS", "boxers": "RHI"}  # labels shown in tables and figures
STUDY_SHEETS = {"HC": "hc", "PPCS": "ppcs", "Boxers": "boxers"}   # study-file sheet -> group key
STRATUM_TO_GROUP = {"HC": "hc", "PPCS": "ppcs", "Boxer": "boxers", "Boxers": "boxers", "RHI": "boxers"}
INCLUSION_COL = "Will we used it?"

# measure key -> (CSV column, panel title, y-axis label)
MEASURES = {
    "total":     ("PSD vol cm3",          "Total",     "Total PSD volume (cm$^3$)"),
    "frontal":   ("PSD frontal vol cm3",  "Frontal",   "Frontal PSD volume (cm$^3$)"),
    "parietal":  ("PSD parietal vol cm3", "Parietal",  "Parietal PSD volume (cm$^3$)"),
    "occipital": ("PSD occipital vol cm3", "Occipital", "Occipital PSD volume (cm$^3$)"),
}
# contrasts as (first, second): delta > 0 means FIRST larger than SECOND
CONTRASTS = (("ppcs", "hc"), ("boxers", "hc"), ("boxers", "ppcs"))
GROUP_COLORS = {"hc": "darkkhaki", "ppcs": "darkgoldenrod", "boxers": "maroon"}


# ----------------------------------------------------------------------------------------------
# Participant selection
# ----------------------------------------------------------------------------------------------
def norm_id(x):
    """'sub-105', 'Sub-105', 'sub-105_FLAIR_label.nii.gz', 105, 105.0, '105' -> 'sub-105'.
    Blank cells, 'Total:' and other non-ID text -> None."""
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return None
    s = str(x).strip()
    m = re.fullmatch(r"(\d+)(?:\.0+)?", s) or re.search(r"(?i)sub[-_]?(\d+)", s)
    return f"sub-{int(m.group(1))}" if m else None


@lru_cache(maxsize=None)
def load_inclusion(study_xlsx, col=INCLUSION_COL):
    """{'hc': frozenset, 'ppcs': frozenset, 'boxers': frozenset} of participants marked 'Yes'."""
    if not os.path.isfile(study_xlsx):
        raise FileNotFoundError(f"Study-info spreadsheet not found: {study_xlsx}")
    raw = pd.read_excel(study_xlsx, sheet_name=list(STUDY_SHEETS))
    inc = {}
    for sheet, grp in STUDY_SHEETS.items():
        df = raw[sheet]
        if col not in df.columns:
            raise KeyError(f"Column '{col}' not found in sheet '{sheet}'; columns: {list(df.columns)}")
        ids, not_yes = set(), set()
        for sid, flag in zip(df.iloc[:, 0], df[col]):
            sid = norm_id(sid)
            if sid is None:                                  # blank row / 'Total:' footer
                continue
            (ids if str(flag).strip().lower() == "yes" else not_yes).add(sid)
        if not_yes:
            print(f"[i] Study file, {sheet}: listed but not 'Yes' -> excluded: {sorted(not_yes)}")
        tot = df[df.iloc[:, 0].astype(str).str.strip().str.lower() == "total:"]
        if len(tot) and pd.notna(tot.iloc[0, 1]):
            try:
                if int(tot.iloc[0, 1]) != len(ids):
                    print(f"[!] Sheet '{sheet}': 'Total:' says {int(tot.iloc[0, 1])} but {len(ids)} rows are 'Yes'")
            except (TypeError, ValueError):
                pass
        inc[grp] = ids
    for a in GROUPS:
        for b in GROUPS:
            if a < b and inc[a] & inc[b]:
                raise ValueError(f"Participants listed in two groups of the study file: {sorted(inc[a] & inc[b])}")
    return {g: frozenset(v) for g, v in inc.items()}


def load_psd_table(psd_csv):
    """Read the gathered PSD CSV; strip header whitespace; add 'sub' (normalised ID)."""
    if not os.path.isfile(psd_csv):
        raise FileNotFoundError(f"PSD table not found: {psd_csv}")
    df = pd.read_csv(psd_csv)
    df.columns = [str(c).strip() for c in df.columns]
    need = ["Scan", "Stratum"] + [v[0] for v in MEASURES.values()]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise KeyError(f"Columns missing from {psd_csv}: {missing}; found {list(df.columns)}")
    df["sub"] = df["Scan"].map(norm_id)
    if df["sub"].isna().any():
        raise ValueError(f"Rows without a recognisable subject ID in 'Scan': {df.loc[df['sub'].isna(), 'Scan'].tolist()}")
    dup = sorted(df.loc[df["sub"].duplicated(keep=False), "sub"].unique())
    if dup:
        raise ValueError(f"Participants appearing more than once in {psd_csv}: {dup}")
    for _, (col, _, _) in MEASURES.items():
        if not pd.api.types.is_numeric_dtype(df[col]):
            raise TypeError(f"Column '{col}' is not numeric")
    # the total must be frontal + parietal + occipital (prefrontal NOT included)
    comp = df["PSD frontal vol cm3"] + df["PSD parietal vol cm3"] + df["PSD occipital vol cm3"]
    bad = df.loc[(df["PSD vol cm3"] - comp).abs() > 0.02, "sub"].tolist()
    if bad:
        print(f"[!] 'PSD vol cm3' differs from frontal+parietal+occipital (by >0.02 cm3) for: {bad}. "
              f"Is the prefrontal volume still part of the total?")
    return df


def select_participants(df, inclusion, strict_missing=False):
    """Keep only study-listed participants, assign group from the study file, return (data, audit)."""
    owner = {s: g for g in GROUPS for s in inclusion[g]}
    in_csv = set(df["sub"])
    dropped = sorted(in_csv - set(owner))
    missing = {g: sorted(inclusion[g] - in_csv) for g in GROUPS}
    out = df[df["sub"].isin(owner)].copy()
    out["group"] = out["sub"].map(owner)
    # the CSV 'Stratum' must agree with the study file
    out["stratum_group"] = out["Stratum"].astype(str).str.strip().map(STRATUM_TO_GROUP)
    clash = out[out["stratum_group"] != out["group"]]
    if len(clash):
        raise ValueError("CSV 'Stratum' disagrees with the study file for: "
                         + ", ".join(f"{r['sub']} (CSV {r['Stratum']}, study {PRETTY[r['group']]})" for _, r in clash.iterrows()))
    if strict_missing and any(missing.values()):
        raise SystemExit(f"--strict-missing: study participants without a PSD row: {missing}")
    out = out.drop(columns=["stratum_group"]).sort_values(["group", "sub"]).reset_index(drop=True)
    audit = {
        "rows_in_csv": int(len(df)),
        "rows_kept": int(len(out)),
        "n_study_list": {g: len(inclusion[g]) for g in GROUPS},
        "n_kept": {g: int((out["group"] == g).sum()) for g in GROUPS},
        "dropped_not_in_study_list": dropped,
        "study_participants_without_psd_row": missing,
    }
    return out, audit


# ----------------------------------------------------------------------------------------------
# Statistics
# ----------------------------------------------------------------------------------------------
def cliffs_delta(x, y):
    """P(x>y) - P(x<y); positive = x larger."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    return float(np.sign(x[:, None] - y[None, :]).mean())


def cliffs_delta_ci(x, y, n_boot=10000, seed=0, alpha=0.05):
    rng = np.random.default_rng(seed)
    x, y = np.asarray(x, float), np.asarray(y, float)
    d = np.empty(n_boot)
    for b in range(n_boot):
        xb = x[rng.integers(0, len(x), len(x))]
        yb = y[rng.integers(0, len(y), len(y))]
        d[b] = np.sign(xb[:, None] - yb[None, :]).mean()
    return tuple(np.percentile(d, [100 * alpha / 2, 100 * (1 - alpha / 2)]))


def magnitude(delta):
    a = abs(delta)
    return "negligible" if a < 0.147 else "small" if a < 0.33 else "medium" if a < 0.474 else "large"


def run_statistics(data, n_boot=10000, seed=0, fdr_family="all"):
    """Return dict of DataFrames: descriptives, kruskal_wallis, pairwise."""
    desc, kw_rows, pw_rows = [], [], []
    for mkey, (col, title, _) in MEASURES.items():
        vals = {g: data.loc[data["group"] == g, col].dropna().to_numpy(float) for g in GROUPS}
        for g in GROUPS:
            v = vals[g]
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            desc.append({"measure": title, "group": PRETTY[g], "n": len(v), "median": med, "Q1": q1, "Q3": q3,
                         "mean": v.mean(), "SD": v.std(ddof=1), "min": v.min(), "max": v.max()})
        H, p = kruskal(*[vals[g] for g in GROUPS])
        kw_rows.append({"measure": title, "H": H, "df": 2, "P": p, **{f"n_{PRETTY[g]}": len(vals[g]) for g in GROUPS}})
        for a, b in CONTRASTS:
            U, pu = mannwhitneyu(vals[a], vals[b], alternative="two-sided")
            d = cliffs_delta(vals[a], vals[b])
            assert abs(d - (2 * U / (len(vals[a]) * len(vals[b])) - 1)) < 1e-9      # U-based and direct delta agree
            lo, hi = cliffs_delta_ci(vals[a], vals[b], n_boot=n_boot, seed=seed)
            pw_rows.append({"measure": title, "contrast": f"{PRETTY[a]} vs {PRETTY[b]}", "first": a, "second": b,
                            "n_first": len(vals[a]), "n_second": len(vals[b]),
                            "median_first": float(np.median(vals[a])), "median_second": float(np.median(vals[b])),
                            "U": U, "P": pu, "cliffs_delta": d, "delta_CI_low": lo, "delta_CI_high": hi,
                            "magnitude": magnitude(d)})
    desc_df, kw_df, pw_df = pd.DataFrame(desc), pd.DataFrame(kw_rows), pd.DataFrame(pw_rows)

    kw_df["q_BH"] = multipletests(kw_df["P"], method="fdr_bh")[1]
    if fdr_family == "all":
        pw_df["q_BH"] = multipletests(pw_df["P"], method="fdr_bh")[1]
        fam = f"all {len(pw_df)} pairwise tests together"
    else:
        pw_df["q_BH"] = np.nan
        for c, idx in pw_df.groupby("contrast").groups.items():
            pw_df.loc[idx, "q_BH"] = multipletests(pw_df.loc[idx, "P"], method="fdr_bh")[1]
        fam = "each contrast separately across measures"
    pw_df["significant_q<0.05"] = pw_df["q_BH"] < 0.05
    kw_df["significant_q<0.05"] = kw_df["q_BH"] < 0.05
    return {"descriptives": desc_df, "kruskal_wallis": kw_df, "pairwise": pw_df, "fdr_family_text": fam}


# ----------------------------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------------------------
def _stars(q):
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""


def _draw_panel(ax, data, mkey, pairwise, ylabel, fs_tick, fs_axis, fs_ann, show_all_deltas=True):
    import seaborn as sns
    col = MEASURES[mkey][0]
    title = MEASURES[mkey][1]
    plot_df = data[["group", col]].rename(columns={col: "value"}).dropna()
    plot_df["label"] = plot_df["group"].map(PRETTY)
    order = [PRETTY[g] for g in GROUPS]
    pal = {PRETTY[g]: GROUP_COLORS[g] for g in GROUPS}
    try:
        sns.violinplot(data=plot_df, x="label", y="value", hue="label", order=order, hue_order=order,
                       palette=pal, inner=None, cut=0, legend=False, ax=ax)
    except TypeError:                                   # older seaborn
        sns.violinplot(data=plot_df, x="label", y="value", order=order, palette=pal, inner=None, cut=0, ax=ax)
    sns.swarmplot(data=plot_df, x="label", y="value", order=order, color="yellow", size=6, alpha=0.8,
                  edgecolor="dimgray", linewidth=0.5, ax=ax)

    lo, hi = plot_df["value"].min(), plot_df["value"].max()
    span = hi - lo
    rows = pairwise[pairwise["measure"] == title]
    draw = rows[(rows["q_BH"] < 0.05) | bool(show_all_deltas)]
    y_bottom = lo - 0.06 * span                         # start near the data minimum
    bracket_step = 0.24 * span                          # room for bracket + label (+ stars)
    y_top = hi + 0.06 * span + bracket_step * len(draw)
    ax.set_ylim(y_bottom, y_top)

    for k, (_, r) in enumerate(draw.iterrows()):
        sig = r["q_BH"] < 0.05
        x1, x2 = GROUPS.index(r["first"]), GROUPS.index(r["second"])
        y = hi + 0.06 * span + bracket_step * k + 0.02 * span
        ax.plot([x1, x2], [y, y], lw=2, c="k")
        strong = abs(r["cliffs_delta"]) >= 0.474
        d = r["cliffs_delta"]
        rel = ">" if d > 0.005 else "<" if d < -0.005 else "="
        label = (f"{PRETTY[r['first']]} {rel} {PRETTY[r['second']]}: \u0394 = {d:+.2f}\n"
                 f"95% CI [{r['delta_CI_low']:+.2f}, {r['delta_CI_high']:+.2f}]")
        ax.text((x1 + x2) / 2, y + 0.01 * span, label, ha="center", va="bottom",
                fontsize=fs_ann, fontweight="bold" if strong else "normal", color="k", linespacing=1.1)
        if sig:
            ax.text((x1 + x2) / 2, y + 0.14 * span, _stars(r["q_BH"]), ha="center", va="bottom",
                    fontsize=fs_ann + 6, fontweight="bold", color="k")
    n_drawn = len(draw)

    # y ticks only over the data range (not up into the bracket headroom)
    ticks = [t for t in ax.get_yticks() if y_bottom <= t <= hi + 0.02 * span]
    ax.set_yticks(ticks)
    ax.spines["left"].set_bounds(y_bottom, hi + 0.02 * span)

    ax.set_xlabel("")
    ax.set_ylabel(ylabel, fontsize=fs_axis, fontweight="bold")
    ax.tick_params(axis="both", labelsize=fs_tick, length=7, width=3)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight("bold")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("bottom", "left"):
        ax.spines[s].set_linewidth(3)
    return n_drawn


def make_figures(data, pairwise, out_dir, show_all_deltas=True, dpi=150):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    written = []
    for mkey, (col, title, ylabel) in MEASURES.items():
        fig, ax = plt.subplots(figsize=(8, 7))
        _draw_panel(ax, data, mkey, pairwise, ylabel, fs_tick=22, fs_axis=26, fs_ann=17,
                    show_all_deltas=show_all_deltas)
        fig.tight_layout()
        p = os.path.join(out_dir, f"PSD_{mkey}.png")
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        written.append(p)
    # combined 1 x 4 panel figure
    fig, axes = plt.subplots(1, 4, figsize=(26, 7.5))
    for ax, (mkey, (col, title, ylabel)) in zip(axes, MEASURES.items()):
        _draw_panel(ax, data, mkey, pairwise, "PSD volume (cm$^3$)" if mkey == "total" else "",
                    fs_tick=20, fs_axis=26, fs_ann=15, show_all_deltas=show_all_deltas)
        ax.set_title(title, fontsize=26, fontweight="bold", pad=14)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        p = os.path.join(out_dir, f"PSD_all_regions.{ext}")
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        written.append(p)
    plt.close(fig)
    return written


# ----------------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------------
def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--psd-csv", required=True, help="gathered PSD table (FLAIR_output_gathered_no-prefrontal.csv)")
    ap.add_argument("--study-xlsx", required=True, help="study-info spreadsheet (sheets HC / PPCS / Boxers)")
    ap.add_argument("--out-dir", required=True, help="output folder (created with --apply)")
    ap.add_argument("--inclusion-col", default=INCLUSION_COL, help="column holding 'Yes' (default: '%(default)s')")
    ap.add_argument("--fdr-family", choices=("all", "per-contrast"), default="all",
                    help="BH family for the pairwise tests (default: all 12 together)")
    ap.add_argument("--n-boot", type=int, default=10000, help="bootstrap resamples for the delta CI (default 10000)")
    ap.add_argument("--seed", type=int, default=0, help="random seed for the bootstrap (default 0)")
    ap.add_argument("--significant-only", action="store_true",
                    help="annotate only contrasts with q < 0.05 (default: every contrast carries "
                         "'group > group: delta [95%% CI]'; asterisks only when q < 0.05)")
    ap.add_argument("--strict-missing", action="store_true",
                    help="stop if a participant marked Yes has no row in the PSD table (e.g. sub-145)")
    ap.add_argument("--apply", action="store_true", help="write outputs (default is a dry run: print only)")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    inclusion = load_inclusion(args.study_xlsx, args.inclusion_col)
    print("[i] Study file: " + ", ".join(f"{PRETTY[g]}={len(inclusion[g])}" for g in GROUPS) + " participants marked Yes")

    table = load_psd_table(args.psd_csv)
    data, audit = select_participants(table, inclusion, strict_missing=args.strict_missing)

    print("\n--- PARTICIPANT AUDIT ---")
    print(f"PSD table rows: {audit['rows_in_csv']}; kept: {audit['rows_kept']}")
    print("n analysed:      " + ", ".join(f"{PRETTY[g]}={audit['n_kept'][g]}" for g in GROUPS))
    print(f"Dropped (in PSD table, not in the study list): {audit['dropped_not_in_study_list']}")
    for g in GROUPS:
        if audit["study_participants_without_psd_row"][g]:
            print(f"[!] {PRETTY[g]}: marked Yes but NO row in the PSD table: {audit['study_participants_without_psd_row'][g]}")

    res = run_statistics(data, n_boot=args.n_boot, seed=args.seed, fdr_family=args.fdr_family)
    pd.set_option("display.width", 200, "display.max_columns", 30)
    print("\n--- DESCRIPTIVES (cm3) ---")
    print(res["descriptives"].round(2).to_string(index=False))
    print("\n--- KRUSKAL-WALLIS (BH across the 4 measures) ---")
    print(res["kruskal_wallis"].round(3).to_string(index=False))
    print(f"\n--- PAIRWISE MANN-WHITNEY / CLIFF'S DELTA (BH: {res['fdr_family_text']}) ---")
    show = res["pairwise"][["measure", "contrast", "n_first", "n_second", "U", "P", "q_BH", "cliffs_delta",
                            "delta_CI_low", "delta_CI_high", "magnitude"]]
    print(show.round(3).to_string(index=False))
    n_sig = int(res["pairwise"]["significant_q<0.05"].sum())
    print(f"\nContrasts with q < 0.05: {n_sig} of {len(res['pairwise'])}"
          + ("  (figures will carry no significance brackets)" if n_sig == 0 and args.significant_only else ""))

    if not args.apply:
        print("\n[dry run] nothing written. Re-run with --apply to write the workbook, data table and figures.")
        return res, data, audit

    os.makedirs(args.out_dir, exist_ok=True)
    data_cols = ["sub", "group"] + [v[0] for v in MEASURES.values()]
    used = data[data_cols].copy()
    used["group"] = used["group"].map(PRETTY)
    used.to_csv(os.path.join(args.out_dir, "PSD_data_used.csv"), index=False)

    xlsx = os.path.join(args.out_dir, "PSD_statistics.xlsx")
    run_info = pd.DataFrame({
        "item": ["psd_csv", "psd_csv_sha256", "study_xlsx", "study_xlsx_sha256", "inclusion_col", "fdr_family",
                 "n_boot", "seed", "n_HC", "n_PPCS", "n_RHI", "dropped_not_in_study_list",
                 "study_participants_without_psd_row", "python", "numpy", "pandas", "scipy"],
        "value": [os.path.abspath(args.psd_csv), _sha256(args.psd_csv), os.path.abspath(args.study_xlsx),
                  _sha256(args.study_xlsx), args.inclusion_col, res["fdr_family_text"], args.n_boot, args.seed,
                  audit["n_kept"]["hc"], audit["n_kept"]["ppcs"], audit["n_kept"]["boxers"],
                  ", ".join(audit["dropped_not_in_study_list"]),
                  json.dumps(audit["study_participants_without_psd_row"]),
                  sys.version.split()[0], np.__version__, pd.__version__, __import__("scipy").__version__]})
    with pd.ExcelWriter(xlsx) as w:
        res["descriptives"].to_excel(w, sheet_name="Descriptives", index=False)
        res["kruskal_wallis"].to_excel(w, sheet_name="Kruskal-Wallis", index=False)
        res["pairwise"].drop(columns=["first", "second"]).to_excel(w, sheet_name="Pairwise", index=False)
        used.to_excel(w, sheet_name="Participants_used", index=False)
        run_info.to_excel(w, sheet_name="run_info", index=False)
    figs = make_figures(data, res["pairwise"], args.out_dir, show_all_deltas=not args.significant_only)

    with open(os.path.join(args.out_dir, "PSD_run_provenance.json"), "w") as f:
        json.dump({"args": vars(args), "audit": audit, "fdr_family": res["fdr_family_text"]}, f, indent=2, default=str)
    print(f"\nWrote: {xlsx}\n       {os.path.join(args.out_dir, 'PSD_data_used.csv')}")
    for p in figs:
        print(f"       {p}")
    return res, data, audit


if __name__ == "__main__":
    main()
