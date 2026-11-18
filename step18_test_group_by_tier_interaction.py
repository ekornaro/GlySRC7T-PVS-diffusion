#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step18_test_group_by_tier_interaction.py
========================================
CHANGELOG v11 (2026-10-04)
--------------------------
 * ONE FDR FAMILY (--fdr-families single, now the default): the low-b and
   full-range fits are two estimates of the same quantities, so they are
   corrected together (paper-wide rule, also used by the regional statistics).
   --fdr-families split reproduces v10.
 * --exclude-metrics (default: lowb_ADC). Low-b ADC duplicates low-b MD
   (r = 0.998 across subjects); it is not fitted, tested or reported. The
   family is therefore 14 metrics: MD AD RD FA MK AK RK Dpar Dperp
   Dpar_minus_Dperp lowb_MD lowb_AD lowb_RD lowb_FA.
 * Default --out MixedEffects_group_x_tier_v11.xlsx; run_info records the
   excluded metrics.

(file history: v6 statistical core -> v7 --subject-list and run_info ->
 v8 unchanged core -> v9 within-subject reference tiers -> v10 FDR families)

CHANGELOG v10 (2026-09-29)
--------------------------
Models, permutations and raw p-values are UNCHANGED; only the FDR families.
  --fdr-families split (default): two families, decided by the question each
      metric answers, not by its result --
        main : MD, AD, RD, FA, MK, AK, RK, Dpar, Dperp, Dpar_minus_Dperp
               (the metrics reported in the paper)
        lowb : lowb_ADC, lowb_MD, lowb_AD, lowb_RD, lowb_FA
               (low-b fit, reported separately as the answer to Reviewer 3,
               major concern 2: sensitivity to fast-diffusing water)
      Every q-value (interaction tests, and the within-subject q-values per
      reference tier) is computed within its family. Each row carries an
      'fdr_family' column; run_info records the setting and the members.
  --fdr-families single: one family across all metrics (= v9).

CHANGELOG v9 (2026-09-28)
-------------------------
The group x tier models, permutation tests and per-tier table are UNCHANGED.
Only the within-subject sheet is extended, so its rows for the default
reference tier reproduce v8 exactly.

  (a) --within-refs TIER [TIER ...]  reference tiers for the within-subject
      contrast (PVS minus reference, per subject). Default = --primary-fake
      only, i.e. v8 behaviour. Use '--within-refs Fake_m6 NearPVS' to add the
      PVS-vs-NearPVS contrast. Why it matters: for the orientation metrics
      (Dpar, Dperp, Dpar_minus_Dperp) the PVS axis transplanted into a distant
      control lands in DIFFERENT fibre geometry, whereas NearPVS shares the
      PVS's fibre environment, so NearPVS is the fair comparator there.
      An unknown reference tier is a hard error, not a silent skip.
  (b) PPCS columns: n_*, delta_PPCS, p_PPCS_vs_0 alongside HC and RHI.
  (c) FDR on the within-group tests against zero (q_HC_vs_0, q_PPCS_vs_0,
      q_RHI_vs_0), and q_RHI_vs_HC, each computed ACROSS METRICS WITHIN ONE
      REFERENCE TIER. Pooling tiers would change the v8 Fake_m6 q-values.
  (d) Column 'Fake_tier' renamed 'Ref_tier': a column called Fake_tier that
      holds NearPVS would be a mislabel. Downstream readers of the
      within_subject sheet must use 'Ref_tier'.
  (e) Metric/reference combinations that cannot be formed (e.g. orientation
      metrics have no NearPVSfull column) are listed in run_info and printed,
      never dropped silently.
  (f) default --out is now MixedEffects_group_x_tier_v9.xlsx (v8 still wrote
      a v6-named file).
  (g) the printed verdict names the tiers actually tested.

Mixed-effects group x tier model: the direct test of PVS spatial specificity.

WHY THIS REPLACES THE PER-TIER TESTS
------------------------------------
Every earlier analysis ran a SEPARATE RHI-vs-HC test in each tier and invited the
reader to notice that the effect looked bigger at PVS. That is not a test. Four
significance markers arranged in decreasing order do not establish that the tier
effects DIFFER, and in this dataset they demonstrably do not:

    MD, raw RHI-HC difference   PVS 6.20e-05   NearPVS 3.88e-05
                                Fake m6 5.41e-05   Fake m10 6.71e-05

Cliff's delta falls across tiers (0.66, 0.69, 0.45, 0.41) only because the
control-tier SD grows from 4.2e-05 to 1.8e-04. The quantity with physical units
does not decay at all.

The question "is the group difference LARGER at PVS than in matched distant white
matter" is a GROUP x TIER INTERACTION, and this script tests it directly, with
subject as a random intercept so the tiers are correctly treated as repeated
measures within subject.

WHAT IS FITTED
--------------
    value ~ Group * Tier + covariates,   random intercept per Subject

    H0 (no spatial specificity):  all Group x Tier terms are zero -- the group
                                  difference is the same everywhere
    H1 (PVS-localised):           the interaction is non-zero, with the group
                                  difference attenuating away from PVS

Tested by likelihood ratio against the same model without the interaction, fitted
by ML (not REML -- LR tests on fixed effects require ML).

COVARIATES
----------
Carried from registration_confounds_subject_covariates.csv, all subject-level:

    icv_ml           head size. RHI have smaller heads (1517 vs 1669 mL,
                     p = 0.033) with no difference in brain:ICV ratio (p = 0.48),
                     so this is skull size, not atrophy, and adjusting for it is
                     legitimate rather than conditioning on a mediator.
    wm_absldj_med    |log det J| of the warp: local volume change, hence how much
                     interpolation smoothing was applied. This is the mechanism
                     by which a small high-MD structure gets averaged toward its
                     surroundings. Rotation is NOT used: MD is the trace and is
                     rotation-invariant, so conditioning on rotation is
                     over-adjustment on a non-mechanistic proxy.
    scale_x          affine scaling, which differs by group (p = 0.007).

pool_frac is deliberately NOT included as a covariate. It is a property of the
control tiers only, so it does not exist for PVS or NearPVS rows; putting it in a
model spanning all tiers would be incoherent. It is reported separately.

THE WITHIN-SUBJECT CROSS-CHECK
------------------------------
The same question can be asked without any random-effects machinery:

    delta_s = median(PVS)_s - median(FakeControl)_s      (one number per subject)
    delta ~ Group + covariates

If PVS were specifically affected, delta would be larger in RHI. This is
algebraically close to the interaction contrast and is reported alongside it as a
sanity check; the two should agree.

USAGE
-----
  python step18_test_group_by_tier_interaction.py \\
      --results pvs_tier_metrics_v6.xlsx \\
      --covariates registration_confounds_subject_covariates.csv
"""

import os
import re
import argparse
import warnings
import logging
import numpy as np

warnings.filterwarnings("ignore")
logging.getLogger("statsmodels").setLevel(logging.ERROR)

try:
    import pandas as pd
    from scipy import stats
    import statsmodels.formula.api as smf
    import statsmodels.api as sm
except ImportError as e:
    raise SystemExit(f"missing dependency: {e}")

SUFFIX_RE = re.compile(r"_(Original|NearPVS|NearPVSfull|FakeControl_[dm]\d+|WholeWM)$")
ORIENTATION = {"Dpar", "Dperp", "Dpar_minus_Dperp"}
COVARS = ["icv_ml", "wm_absldj_med", "scale_x"]
GROUPS = ["HC", "PPCS", "RHI"]


def tier_label(suffix):
    if suffix == "Original":
        return "PVS"
    if suffix == "NearPVS":
        return "NearPVS"
    if suffix == "NearPVSfull":
        return "NearPVSfull"
    if suffix == "WholeWM":                 # step24 benchmark: not a matched control
        return "WholeWM"
    return "Fake_" + suffix.split("_")[-1]


def tier_margin(t):
    """Sort key in millimetres. Plain string sort puts m10 before m6."""
    if t == "PVS":
        return 0.0
    if t in ("NearPVS", "NearPVSfull"):
        return 4.5
    if t == "WholeWM":                      # step24 benchmark tier, sorts last
        return 99.0
    m = re.search(r"(\d+)", t)
    return float(m.group(1)) if m else 99.0


def build_long(results, covariates, drop_orientation=True):
    near = pd.read_excel(results, sheet_name="Case_A_Near_PVS")
    fake = pd.read_excel(results, sheet_name="Case_B_Fake_PVS")
    cov = pd.read_csv(covariates)
    keep = ["Subject", "Group"] + [c for c in COVARS if c in cov.columns]
    cov = cov[keep]

    recs = []
    for df in (near, fake):
        for c in df.columns:
            mt = SUFFIX_RE.search(c)
            if not mt:
                continue
            base = c[:mt.start()]
            if drop_orientation and base in ORIENTATION:
                continue
            t = tier_label(mt.group(1))
            if t == "PVS" and df is fake:
                continue                      # identical column in both sheets
            s = df.groupby("Subject")[c].median().reset_index()
            s.columns = ["Subject", "value"]
            s["Metric"], s["Tier"] = base, t
            recs.append(s)
    long = pd.concat(recs, ignore_index=True).merge(cov, on="Subject", how="inner")
    long = long[long.Group.isin(GROUPS)]
    return long


def fit_metric(d, metric, ref_tier="PVS", covars=(), min_n=8):
    """MixedLM with and without the Group x Tier interaction; LR test between."""
    d = d[d.Metric == metric].dropna(subset=["value"]).copy()
    if d.Subject.nunique() < min_n or d.Tier.nunique() < 2:
        return None
    sc = 1e3 if d.value.abs().median() < 1e-2 else 1.0
    d["y"] = d.value * sc
    tiers = sorted(d.Tier.unique(), key=tier_margin)
    d["Tier"] = pd.Categorical(d.Tier, categories=tiers, ordered=False)
    present = [g for g in GROUPS if (d.Group == g).any()]
    # drop unused levels: an empty category becomes an all-zero dummy column and
    # makes the design matrix singular
    d["Group"] = pd.Categorical(d.Group, categories=present, ordered=False)
    use = [c for c in covars if c in d.columns and d[c].notna().all()]
    extra = (" + " + " + ".join(use)) if use else ""

    # ML, not REML: likelihood ratio tests on FIXED effects require ML.
    #
    # MixedLM can fail here for a structural reason: the covariates are
    # SUBJECT-LEVEL and the random intercept already absorbs between-subject
    # variance, so when the intercept variance runs to the boundary the design
    # becomes near-singular. When that happens we fall back to OLS with
    # cluster-robust (by subject) standard errors, which handles the repeated
    # measures without a variance component and gives a valid Wald test on the
    # interaction. The method used is recorded per metric.
    method, p_int, lr, df_diff = None, np.nan, np.nan, np.nan
    full = None
    try:
        full = smf.mixedlm(f"y ~ C(Group) * C(Tier){extra}", d,
                           groups=d["Subject"]).fit(reml=False, method="lbfgs")
        red = smf.mixedlm(f"y ~ C(Group) + C(Tier){extra}", d,
                          groups=d["Subject"]).fit(reml=False, method="lbfgs")
        df_diff = len(full.params) - len(red.params)
        lr = 2 * (full.llf - red.llf)
        if not np.isfinite(lr) or lr < 0 or df_diff <= 0:
            raise ValueError("degenerate likelihood ratio")
        p_int = float(stats.chi2.sf(lr, df_diff))
        method = "MixedLM(ML)"
    except Exception:
        try:
            full = smf.ols(f"y ~ C(Group) * C(Tier){extra}", data=d).fit(
                cov_type="cluster", cov_kwds={"groups": d["Subject"]})
            terms = [t for t in full.params.index if ":" in t]
            if not terms:
                return None
            w = full.f_test(" = 0, ".join(terms) + " = 0")
            p_int = float(np.squeeze(w.pvalue))
            lr = float(np.squeeze(w.fvalue))
            df_diff = len(terms)
            method = "OLS+cluster-robust"
        except Exception as exc:
            print(f"  [{metric}] model failed: {exc}")
            return None

    rows = []
    for t in tiers:
        sub = d[d.Tier == t]
        hc = sub[sub.Group == "HC"].value
        rhi = sub[sub.Group == "RHI"].value
        if len(hc) < 3 or len(rhi) < 3:
            continue
        raw = float(rhi.median() - hc.median())
        sd = float(np.sqrt(((len(hc) - 1) * hc.var(ddof=1)
                            + (len(rhi) - 1) * rhi.var(ddof=1))
                           / (len(hc) + len(rhi) - 2)))
        rows.append({"Metric": metric, "Tier": t, "n_HC": len(hc), "n_RHI": len(rhi),
                     "raw_RHI_minus_HC": raw, "pooled_SD": sd,
                     "cohens_d": raw / sd if sd > 0 else np.nan,
                     "raw_rel_to_PVS": np.nan})
    if rows:
        base = rows[0]["raw_RHI_minus_HC"]
        for r in rows:
            r["raw_rel_to_PVS"] = (r["raw_RHI_minus_HC"] / base
                                   if base not in (0, np.nan) else np.nan)
    return {"metric": metric, "p_interaction": p_int, "lr": lr, "df": df_diff,
            "n_subjects": d.Subject.nunique(), "tiers": tiers, "method": method,
            "per_tier": rows, "full": full, "covars_used": use}


def within_subject(d, metric, ref_tier, covars=()):
    """delta = PVS - reference tier per subject, then delta ~ Group + covariates.

    v9: reference may be any tier (Fake_m6, Fake_m10, NearPVS, ...); returns
    per-group medians and Wilcoxon signed-rank tests of delta against zero for
    HC, PPCS and RHI, plus the HC-vs-RHI comparison of delta (as in v8).
    """
    w = d[(d.Metric == metric) & (d.Tier.isin(["PVS", ref_tier]))]
    p = w.pivot_table(index=["Subject", "Group"] + [c for c in covars if c in w],
                      columns="Tier", values="value").reset_index()
    if "PVS" not in p or ref_tier not in p:
        return None
    p["delta"] = p["PVS"] - p[ref_tier]
    p = p.dropna(subset=["delta"])
    by = {g: p[p.Group == g].delta for g in GROUPS}
    hc, rhi = by["HC"], by["RHI"]
    if len(hc) < 3 or len(rhi) < 3:
        return None

    def vs0(x):
        if len(x) < 3 or not x.abs().sum():
            return np.nan
        return float(stats.wilcoxon(x).pvalue)

    out = {"Metric": metric, "Ref_tier": ref_tier}
    for g in GROUPS:
        out[f"n_{g}"] = int(len(by[g]))
    for g in GROUPS:
        out[f"delta_{g}"] = float(by[g].median()) if len(by[g]) else np.nan
    for g in GROUPS:
        out[f"p_{g}_vs_0"] = vs0(by[g])
    out["p_RHI_vs_HC"] = float(stats.mannwhitneyu(hc, rhi).pvalue)
    use = [c for c in covars if c in p.columns and p[c].notna().all()]
    if use:
        sc = 1e3 if p.delta.abs().median() < 1e-2 else 1.0
        p["_y"] = p.delta * sc
        p["RHI"] = (p.Group == "RHI").astype(int)
        try:
            r = smf.ols("_y ~ RHI + " + " + ".join(use),
                        data=p[p.Group.isin(["HC", "RHI"])]).fit()
            out["p_RHI_adjusted"] = float(r.pvalues["RHI"])
            out["beta_RHI_adjusted"] = float(r.params["RHI"])
        except Exception:
            pass
    return out


def interaction_stat(d, extra, groups_keep=None):
    """Cluster-robust Wald F on all Group x Tier terms. Returns (F, p, df)."""
    dd = d if groups_keep is None else d[d.Group.isin(groups_keep)]
    if dd.Group.nunique() < 2:
        return np.nan, np.nan, 0
    dd = dd.copy()
    dd["Group"] = pd.Categorical(dd.Group,
                                 categories=[g for g in GROUPS
                                             if (dd.Group == g).any()])
    try:
        m = smf.ols(f"y ~ C(Group) * C(Tier){extra}", data=dd).fit(
            cov_type="cluster", cov_kwds={"groups": dd["Subject"]})
        terms = [t for t in m.params.index if ":" in t]
        if not terms:
            return np.nan, np.nan, 0
        w = m.f_test(" = 0, ".join(terms) + " = 0")
        return float(np.squeeze(w.fvalue)), float(np.squeeze(w.pvalue)), len(terms)
    except Exception:
        return np.nan, np.nan, 0


def permutation_interaction(d, extra, n_perm, rng, groups_keep=None):
    """Exact-under-null p for the interaction.

    Cluster-robust Wald tests are known to OVER-REJECT when the number of
    clusters is modest (here 55 subjects) and several restrictions are tested at
    once (6 df with three groups). Permuting the GROUP LABEL BETWEEN SUBJECTS --
    keeping all of a subject's tiers together -- preserves the within-subject
    covariance and the tier-specific variances exactly, and is valid regardless
    of heteroscedasticity or cluster count.
    """
    obs, _, df = interaction_stat(d, extra, groups_keep)
    if not np.isfinite(obs) or n_perm <= 0:
        return obs, np.nan, df
    lab = d[["Subject", "Group"]].drop_duplicates().reset_index(drop=True)
    hits, done = 0, 0
    for _ in range(n_perm):
        lab2 = lab.copy()
        lab2["Group"] = rng.permutation(lab2.Group.values)
        dp = d.drop(columns=["Group"]).merge(lab2, on="Subject", how="left")
        f, _, _ = interaction_stat(dp, extra, groups_keep)
        if np.isfinite(f):
            done += 1
            hits += (f >= obs)
    return obs, ((hits + 1) / (done + 1) if done else np.nan), df


def fdr_family_of(metric, mode):
    """v10: FDR family of a metric ('main' / 'lowb', or 'all' when single)."""
    if mode == "single":
        return "all"
    return "lowb" if str(metric).startswith("lowb_") else "main"


def fdr_within(df, pcol, keys):
    """BH q-values of df[pcol] computed separately within each group of keys."""
    q = pd.Series(np.nan, index=df.index)
    for _, idx in df.groupby(keys).groups.items():
        q.loc[idx] = fdr(df.loc[idx, pcol].values)
    return q.values


def fdr(p):
    p = np.asarray(p, float)
    ok = np.isfinite(p)
    q = np.full_like(p, np.nan)
    if ok.sum() == 0:
        return q
    v = p[ok]
    o = np.argsort(v)
    n = len(v)
    adj = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        prev = min(prev, v[o[i]] * n / (i + 1))
        adj[o[i]] = prev
    q[ok] = adj
    return q


def main():
    ap = argparse.ArgumentParser(
        description="Group x tier interaction tests on the tier workbook. The "
                    "permutation seed is fixed at 0 (np.random.default_rng(0)); "
                    "the reported runs used --n-perm 10000.")
    ap.add_argument("--results", required=True)
    ap.add_argument("--covariates", required=True)
    ap.add_argument("--out", default="MixedEffects_group_x_tier_v11.xlsx")
    ap.add_argument("--primary-fake", default="Fake_m6",
                    help="control tier for the within-subject check; m6 is "
                         "primary because m10 relaxes matching heavily in the "
                         "high-burden subjects")
    ap.add_argument("--tiers", nargs="*", default=None,
                    help="restrict to these tiers, e.g. --tiers PVS Fake_m6. "
                         "With exactly two tiers and HC vs RHI the interaction "
                         "becomes a SINGLE-DF contrast: is the PVS-minus-control "
                         "difference larger in RHI? That is the focused question, "
                         "free of the NearPVS tier that drives the 4-tier omnibus.")
    ap.add_argument("--n-perm", type=int, default=1000,
                    help="subject-level label permutations for the exact "
                         "interaction test (0 to skip); the reported runs "
                         "used 10000; seed fixed at 0")
    ap.add_argument("--subject-list", default=None,
                    help="v7: restrict the analysis to the subjects named in "
                         "this file (one sub-XXX per line; '#' comments and "
                         "blank lines ignored). Unknown IDs are an ERROR, not "
                         "a silent drop. Use for sensitivity analyses: the "
                         "u_rel_det reflection split, the pre-revision cohort, "
                         "single-block runs.")
    ap.add_argument("--within-refs", nargs="+", default=None, metavar="TIER",
                    help="v9: reference tier(s) for the within-subject contrast "
                         "(PVS minus reference). Default: --primary-fake only "
                         "(v8 behaviour). E.g. --within-refs Fake_m6 NearPVS. "
                         "Unknown tiers are an ERROR.")
    ap.add_argument("--exclude-metrics", nargs="*", default=["lowb_ADC"],
                    help="v11: metrics not fitted or reported (default: lowb_ADC, a "
                         "duplicate of lowb_MD); pass with no names to keep all")
    ap.add_argument("--fdr-families", choices=["split", "single"], default="single",
                    help="v10: 'split' = main metrics and lowb_* metrics are separate "
                         "FDR families (default); 'single' = one family (v9)")
    ap.add_argument("--no-covariates", action="store_true")
    ap.add_argument("--include-orientation", action="store_true")
    a = ap.parse_args()

    rng = np.random.default_rng(0)
    covars = [] if a.no_covariates else COVARS
    d = build_long(a.results, a.covariates,
                   drop_orientation=not a.include_orientation)
    if a.subject_list:                                    # v7
        with open(a.subject_list) as fh:
            wanted = [ln.split("#")[0].strip() for ln in fh]
        wanted = [w for w in wanted if w]
        if not wanted:
            raise SystemExit(f"--subject-list {a.subject_list} is empty")
        have = set(d.Subject.unique())
        unknown = [w for w in wanted if w not in have]
        if unknown:
            raise SystemExit(
                f"--subject-list contains {len(unknown)} ID(s) absent from the "
                f"results/covariates join: {unknown[:8]}"
                f"{' ...' if len(unknown) > 8 else ''}. Refusing to run on a "
                f"silently reduced cohort -- fix the list or the inputs.")
        before = d.Subject.nunique()
        d = d[d.Subject.isin(wanted)]
        print(f"SUBJECT LIST: {os.path.basename(a.subject_list)} -> "
              f"{d.Subject.nunique()}/{before} subjects "
              f"(SENSITIVITY RUN, not the primary analysis)")

    if a.tiers:
        avail = sorted(d.Tier.unique(), key=tier_margin)
        missing = [t for t in a.tiers if t not in avail]
        if missing:
            raise SystemExit(f"tier(s) {missing} not found. Available: {avail}")
        d = d[d.Tier.isin(a.tiers)]
        print(f"RESTRICTED to tiers: {a.tiers}")
        if len(a.tiers) == 2:
            print("  -> single-df focused contrast (HC vs RHI)")
    # v9: resolve the within-subject reference tiers, hard error if absent
    avail_now = sorted(d.Tier.unique(), key=tier_margin)
    if a.within_refs:
        refs = list(dict.fromkeys(a.within_refs))
    else:
        pf0 = a.primary_fake
        if a.tiers and pf0 not in a.tiers:            # v8 fallback, kept
            alt = [t for t in a.tiers if t != "PVS"]
            pf0 = alt[0] if alt else pf0
        refs = [pf0]
    bad = [t for t in refs if t == "PVS" or t not in avail_now]
    if bad:
        raise SystemExit(f"--within-refs: {bad} not usable. PVS is the "
                         f"minuend; available references: "
                         f"{[t for t in avail_now if t != 'PVS']}")
    # v11: drop excluded metrics before anything is fitted or corrected
    excl = sorted(set(a.exclude_metrics or []))
    unknown = [m for m in excl if m not in set(d.Metric)]
    if unknown:
        raise SystemExit(f"--exclude-metrics {unknown}: not in the results "
                         f"({sorted(d.Metric.unique())})")
    if excl:
        d = d[~d.Metric.isin(excl)]
        print(f"excluded metrics: {excl}")
    print(f"subjects: {d.Subject.nunique()}   metrics: {d.Metric.nunique()}   "
          f"tiers: {sorted(d.Tier.unique(), key=tier_margin)}")
    print(f"groups: " + ", ".join(f"{g}={d[d.Group==g].Subject.nunique()}"
                                  for g in GROUPS))
    print(f"covariates: {covars if covars else 'none'}")
    run_info = {"results": os.path.basename(a.results),
                "covariates_file": os.path.basename(a.covariates),
                "covariates_used": ",".join(covars) if covars else "none",
                "subject_list": (os.path.basename(a.subject_list)
                                 if a.subject_list else "ALL"),
                "n_subjects": int(d.Subject.nunique()),
                "tiers": ",".join(sorted(d.Tier.unique(), key=tier_margin)),
                "n_perm": a.n_perm,
                "include_orientation": bool(a.include_orientation),
                "within_refs": ",".join(refs),
                "excluded_metrics": ",".join(excl) if excl else "none",
                "within_q_scope": "per reference tier, across metrics",
                "script": os.path.basename(__file__)}
    print("-" * 78)

    fits, per_tier, within, ws_skipped = [], [], [], []
    skipped = []
    for m in sorted(d.Metric.unique()):
        r = fit_metric(d, m, covars=covars)
        if r is None:                                   # v7: say why
            dm_ = d[d.Metric == m].dropna(subset=["value"])
            skipped.append(f"{m} (n_subj={dm_.Subject.nunique()}, "
                           f"n_tiers={dm_.Tier.nunique()})")
            continue
        fits.append({"Metric": m, "n_subjects": r["n_subjects"],
                     "method": r["method"], "stat": r["lr"], "df": r["df"],
                     "p_group_x_tier": r["p_interaction"]})
        per_tier.extend(r["per_tier"])
        for ref in refs:                                # v9
            w = within_subject(d, m, ref, covars)
            if w:
                within.append(w)
            else:
                ws_skipped.append(f"{m} vs {ref}")
        # robust refits: permutation, and HC-vs-RHI only (drops PPCS contrasts,
        # which are in the omnibus 6 df but not in the question being asked)
        dm = d[d.Metric == m].dropna(subset=["value"]).copy()
        sc = 1e3 if dm.value.abs().median() < 1e-2 else 1.0
        dm["y"] = dm.value * sc
        tiers_m = sorted(dm.Tier.unique(), key=tier_margin)
        dm["Tier"] = pd.Categorical(dm.Tier, categories=tiers_m)
        use_m = [c for c in covars if c in dm.columns and dm[c].notna().all()]
        ex = (" + " + " + ".join(use_m)) if use_m else ""
        f_all, p_perm, df_all = permutation_interaction(dm, ex, a.n_perm, rng)
        f_hr, p_hr, df_hr = interaction_stat(dm, ex, groups_keep=["HC", "RHI"])
        _, p_perm_hr, _ = permutation_interaction(dm, ex, a.n_perm, rng,
                                                  groups_keep=["HC", "RHI"])
        fits[-1].update({"p_perm_allgroups": p_perm,
                         "df_HCvRHI": df_hr, "p_HCvRHI_wald": p_hr,
                         "p_HCvRHI_perm": p_perm_hr})
        print(f"  {m:12s} omnibus p={r['p_interaction']:.4f} [{r['method']}]"
              f"   perm p={p_perm:.4f}   HCvRHI wald p={p_hr:.4f} "
              f"perm p={p_perm_hr:.4f}")

    if skipped:
        print(f"\n  [!] {len(skipped)} metric(s) could not be fitted: "
              + ", ".join(skipped[:6]) + (" ..." if len(skipped) > 6 else ""))
    F = pd.DataFrame(fits)
    if len(F) == 0:
        raise SystemExit(
            "No metric could be fitted on this cohort, so there is nothing to "
            "report.\n"
            f"  subjects={d.Subject.nunique()}  metrics={d.Metric.nunique()}  "
            f"tiers={sorted(d.Tier.unique(), key=tier_margin)}  "
            f"covariates={covars if covars else 'none'}\n"
            "  Most likely causes, in order:\n"
            "   1. a covariate is collinear with Group in this subset (the "
            "reflection\n"
            "      split is heavily block-structured, and scale_x differs by "
            "block) --\n"
            "      re-run with --no-covariates to see whether the fits then "
            "succeed;\n"
            "   2. too few subjects per group for the group x tier design;\n"
            "   3. the subject list and the covariates file do not overlap as "
            "expected.\n"
            "  This is an underpowered/ill-conditioned subset, not a null "
            "result -- report it as such.")
    if F.empty:
        raise SystemExit("no metric could be fitted -- check groups and tiers")
    # v10: q-values within FDR families (main / lowb, or one family)
    F["fdr_family"] = F.Metric.map(lambda m: fdr_family_of(m, a.fdr_families))
    F["q_group_x_tier"] = fdr_within(F, "p_group_x_tier", "fdr_family")
    for c in ("p_perm_allgroups", "p_HCvRHI_wald", "p_HCvRHI_perm"):
        if c in F:
            F["q_" + c.replace("p_", "", 1)] = fdr_within(F, c, "fdr_family")
    fams = {f: sorted(F.loc[F.fdr_family == f, "Metric"]) for f in F.fdr_family.unique()}
    run_info["fdr_families"] = a.fdr_families
    for f, mem in fams.items():
        run_info[f"fdr_family_{f}"] = ", ".join(mem)
    print("\nFDR families (" + a.fdr_families + "): " +
          "; ".join(f"{f} = {len(m)} metrics" for f, m in fams.items()))
    T = pd.DataFrame(per_tier)
    W = pd.DataFrame(within)
    if len(W):                       # v10: FDR per reference tier AND family
        W["fdr_family"] = W.Metric.map(lambda m: fdr_family_of(m, a.fdr_families))
        for col in ["p_HC_vs_0", "p_PPCS_vs_0", "p_RHI_vs_0", "p_RHI_vs_HC"]:
            W["q_" + col[2:]] = fdr_within(W, col, ["Ref_tier", "fdr_family"])
    run_info["within_skipped"] = "; ".join(ws_skipped) if ws_skipped else "none"
    if ws_skipped:
        print(f"\n  [i] within-subject contrasts not formed (no data for that "
              f"metric/tier pair): {', '.join(ws_skipped)}")

    print("\n" + "=" * 78)
    print("IS THE GROUP DIFFERENCE LARGER AT PVS THAN IN MATCHED DISTANT WM?")
    print("=" * 78)
    pcol = "q_HCvRHI_perm" if "q_HCvRHI_perm" in F else "q_group_x_tier"
    print(f"  primary test: {pcol} (permutation, HC vs RHI only)")
    sig = F[F[pcol] < 0.05] if pcol in F else F[F.q_group_x_tier < 0.05]
    if len(sig) == 0:
        print("  NO. Not one metric shows a group x tier interaction surviving FDR.")
        print("  The RHI-HC difference is statistically indistinguishable across")
        print("  the tiers tested: " + ", ".join(sorted(d.Tier.unique(), key=tier_margin)) + ".")
        print("  Nothing here says the group effect is PVS-specific; whether the")
        print("  group effect itself is significant is a separate test.")
        col = pcol.replace("q_", "p_", 1)
        if col in F and F[col].notna().any():
            print(f"  smallest p = {F[col].min():.4f} "
                  f"({F.loc[F[col].idxmin(), 'Metric']})")
    else:
        print(f"  YES for {len(sig)} metric(s): "
              + ", ".join(f"{r.Metric} (q={getattr(r, pcol):.4f})"
                          for r in sig.itertuples()))
        print("  Check the per-tier raw differences: the interaction must reflect")
        print("  a DECAYING raw difference, not merely a growing control-tier SD.")

    print("\nPER-TIER RAW DIFFERENCES (RHI - HC, metric units)")
    if len(T) == 0 or "Metric" not in T.columns:
        # v7: an empty per-tier table means no tier had enough HC AND RHI
        # subjects for a tier-wise summary, even though the models converged.
        print("  (none: no tier had enough HC and RHI subjects in this cohort "
              "for a\n   tier-wise summary -- expected in the reflection-split "
              "sensitivity arms.\n   The interaction test above is still valid; "
              "only this descriptive table\n   is unavailable.)")
        T = pd.DataFrame(columns=["Metric", "Tier", "raw_RHI_minus_HC",
                                  "raw_rel_to_PVS", "cohens_d"])
    for m in ["MD", "AD", "RD"]:
        s = T[T.Metric == m]
        if len(s):
            print(f"  {m:4s} " + "   ".join(
                f"{r.Tier}: {r.raw_RHI_minus_HC:.2e} ({r.raw_rel_to_PVS:.2f}x, "
                f"d={r.cohens_d:+.2f})" for r in s.itertuples()))

    if len(W) and "Metric" in W.columns:
        for ref in refs:                                # v9: one block per ref
            Wr = W[W.Ref_tier == ref]
            if not len(Wr):
                continue
            print(f"\nWITHIN-SUBJECT  (PVS minus own {ref}; medians of per-subject deltas)")
            print(f"  {'metric':17s} {'HC':>10s} {'PPCS':>10s} {'RHI':>10s} "
                  f"{'q_HC0':>8s} {'q_PPCS0':>8s} {'q_RHI0':>8s} {'q_RvH':>7s}")
            for fam in [f for f in ("main", "all", "lowb") if (Wr.fdr_family == f).any()]:
                if a.fdr_families == "split":
                    print(f"  -- FDR family: {fam} --")
                for r in Wr[Wr.fdr_family == fam].sort_values("p_RHI_vs_HC").itertuples():
                    print(f"  {r.Metric:17s} {r.delta_HC:10.2e} {r.delta_PPCS:10.2e} "
                          f"{r.delta_RHI:10.2e} {r.q_HC_vs_0:8.1e} {r.q_PPCS_vs_0:8.1e} "
                          f"{r.q_RHI_vs_0:8.1e} {r.q_RHI_vs_HC:7.3f}")
            if (Wr.q_RHI_vs_HC >= 0.05).all():
                print(f"  The PVS-minus-{ref} excess does not differ between HC "
                      f"and RHI for any metric.")
            print("  q_*0: is PVS different from the reference WITHIN that group "
                  "(FDR across the metrics of one family).")

    with pd.ExcelWriter(a.out) as w:
        F.to_excel(w, sheet_name="interaction_tests", index=False)
        T.to_excel(w, sheet_name="per_tier_raw", index=False)
        if len(W) and "Metric" in W.columns:
            W.to_excel(w, sheet_name="within_subject", index=False)
        # v7: every workbook records the exact run that produced it, so two
        # files with different cohorts/covariates can never be confused.
        pd.DataFrame([run_info]).T.reset_index().rename(
            columns={"index": "setting", 0: "value"}).to_excel(
            w, sheet_name="run_info", index=False)
    print(f"\n-> {a.out}")
    print("\nReport the interaction p-value, not four per-tier p-values. The")
    print("per-tier table is context; the interaction is the test.")


if __name__ == "__main__":
    main()
