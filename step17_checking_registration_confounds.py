#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step17_checking_registration_confounds.py
=========================================
Follow-up questions from the b0 -> T1w rotation audit (v2).

WHAT v2 ADDS
------------
Q1 now also computes |log det J|, the local volume change. Rotation affects only
u^T D u; it cannot change MD, which is the rotation-invariant trace. The route by
which the warp CAN change MD is interpolation smoothing: where the deformation
compresses, several source voxels collapse into one target voxel and the map is
locally averaged, pulling a small high-MD PVS voxel toward its parenchymal
surroundings. |log det J| measures that directly and is the mechanistically
correct covariate.

Q3 tests head size as the upstream variable. The RHI group is age-matched but has
smaller heads, which gives a causal ordering:

    head size -> warp deformation -> interpolation smoothing -> measured MD

Nothing on that path is biology, so adjusting for it is legitimate rather than
over-adjustment. Intracranial volume is also a far more defensible covariate in a
manuscript than a registration diagnostic.

v2_2: STUDY_XLSX -> 2026-08-28 authority (sub-127/130 excluded at source).

Q4 (optional, --metrics-xlsx) refits the group comparison stepwise against that
chain, so you can see exactly where each metric's effect dies -- or does not.

Q1. SPATIAL COINCIDENCE. The SyN warp's local rotation has a median of ~0.9 deg
    but a p90 of ~9.4 deg and a p99 of ~19 deg, consistently across subjects. A
    19 deg local rotation induces |dDpar| ~ 2e-4, which is half the cluster-level
    IQR of Dpar-Dperp and several times the raw RHI-HC difference. That only
    matters if the high-rotation voxels sit where the PVS are. If they sit at the
    frontal/temporal EPI dropout margins, it is irrelevant.

    This script computes the local rotation AT the PVS voxels and compares it to
    a matched random sample of cerebral WM in the same subject, and breaks it
    down by region (CSO / basal ganglia / midbrain).

    Read the enrichment number: the fraction of PVS voxels above that subject's
    WM p90. Under no spatial coincidence it is 0.10 by construction. Values near
    0.10 mean the tail is somewhere else and you can stop worrying; values of
    0.2-0.3 mean the warp is roughest exactly where you are measuring.

Q2. AFFINE SCALING BY GROUP. Residual EPI geometric scaling absorbed by the
    affine is harmless in itself -- the affine corrects it -- but if its
    magnitude differs by group, then the PVS mask aligns to the DKI data
    differently in different groups, which is a genuine confound.

    IMPORTANT METHOD NOTE. An earlier version of this check read the SVD
    singular values of the affine matrix and labelled them sx/sy/sz. Singular
    values come back sorted descending, so s1 > s3 always, and the apparent
    "every subject stretches in x and compresses in z" was an artifact of the
    sort with no axis meaning. This script instead takes the right polar
    decomposition A = R P and reads the DIAGONAL OF P, which is scaling along
    the actual x, y and z axes. Off-diagonal entries of P are reported as shear.

CAVEAT (both sections)
----------------------
The displacement field is in ITK physical (LPS) coordinates while the array axes
are i,j,k. Differentiation here uses voxel spacing without applying the direction
cosine matrix, so angles are accurate to roughly a degree for near-axial
acquisitions. That is ample for a coincidence test; it is not a correction.

USAGE
-----
  python check_registration_confounds.py --root <ROOT> --q1
  python check_registration_confounds.py --root <ROOT> --q2
  python check_registration_confounds.py --root <ROOT>          # both
  python check_registration_confounds.py --root <ROOT> --q1 --subjects sub-109
"""

import os
import glob
import argparse
import numpy as np
import pandas as pd
from pipeline_config import cfg

try:
    import nibabel as nib
    from scipy.io import loadmat
    from scipy import stats
except ImportError as e:
    raise SystemExit(f"missing dependency: {e}")

try:
    import statsmodels.formula.api as smf
    _HAVE_SM = True
except ImportError:
    _HAVE_SM = False

# ==========================================
# CONFIGURATION
# ==========================================
ROOT_DIR = cfg("pvs_root")
STUDY_XLSX = cfg("study_xlsx")
SEG_SUBDIR = "Segmentation"

WM_LABELS = [2, 41]
N_WM_SAMPLE = 200000        # random WM voxels for the reference distribution
N_PVS_CAP = 200000          # cap PVS voxels (they are far fewer in practice)
SHEETS = {"HC": "HC", "PPCS": "PPCS", "Boxers": "RHI"}


# ==========================================
# SHARED
# ==========================================
def jacobian_at(warp_img, idx, zooms):
    """Local Jacobian statistics of F = I + grad(u) at voxel indices idx (N,3).

    Central differences at scattered points; far cheaper than a whole-volume
    Jacobian when PVS voxels are sparse.

    Returns (angle_deg, log_det, n_folded):
      angle_deg  rotation angle of the polar factor -- affects u^T D u
      log_det    log|det F| -- local volume change. Negative = compression, so
                 several source voxels collapse into one target voxel and the
                 map is locally averaged; positive = expansion, interpolated up.
                 |log_det| is therefore the natural proxy for how much
                 interpolation smoothing the warp imposed at that location,
                 which is the mechanism by which a high-MD PVS voxel gets pulled
                 toward its parenchymal surroundings.
      n_folded   voxels with det F <= 0, i.e. the deformation folds. These are a
                 topology violation and are reported rather than silently
                 dropped from the log.
    """
    d = np.asarray(warp_img.dataobj, np.float32)
    d = np.squeeze(d)
    if d.ndim != 4 or d.shape[-1] != 3:
        raise ValueError(f"warp shape {d.shape} unexpected; want (X,Y,Z,3)")
    shape = d.shape[:3]

    keep = np.all((idx >= 1) & (idx < np.array(shape) - 1), axis=1)
    idx = idx[keep]
    if len(idx) == 0:
        return np.array([]), np.array([]), 0

    J = np.empty((len(idx), 3, 3), np.float32)
    for ax in range(3):
        step = np.zeros(3, int)
        step[ax] = 1
        p = idx + step
        m = idx - step
        up = d[p[:, 0], p[:, 1], p[:, 2], :]
        dn = d[m[:, 0], m[:, 1], m[:, 2], :]
        J[:, :, ax] = (up - dn) / (2.0 * zooms[ax])
    J += np.eye(3, dtype=np.float32)

    det = np.linalg.det(J)
    folded = int((det <= 0).sum())
    ldj = np.full(det.shape, np.nan, np.float64)
    ok = det > 0
    ldj[ok] = np.log(det[ok])

    W, s, Vt = np.linalg.svd(J)
    R = W @ Vt
    bad = np.linalg.det(R) < 0
    if bad.any():
        Wb = W[bad].copy()
        Wb[:, :, -1] *= -1
        R[bad] = Wb @ Vt[bad]
    tr = np.trace(R, axis1=1, axis2=2)
    ang = np.degrees(np.arccos(np.clip((tr - 1.0) / 2.0, -1.0, 1.0)))
    return ang, ldj, folded


def map_indices(src_img, dst_img, idx):
    """Voxel indices in src grid -> nearest voxel indices in dst grid."""
    if src_img.shape[:3] == dst_img.shape[:3] and \
            np.allclose(src_img.affine, dst_img.affine, atol=1e-4):
        return idx
    M = np.linalg.inv(dst_img.affine) @ src_img.affine
    h = np.c_[idx, np.ones(len(idx))]
    out = (h @ M.T)[:, :3]
    return np.rint(out).astype(int)


def polar_stretch(A):
    """A = R P. Returns R, per-axis scaling diag(P), max off-diagonal shear."""
    W, s, Vt = np.linalg.svd(A)
    R = W @ Vt
    if np.linalg.det(R) < 0:
        W[:, -1] *= -1
        R = W @ Vt
    P = R.T @ A
    scale = np.diag(P).copy()
    shear = float(np.abs(P - np.diag(np.diag(P))).max())
    return R, scale, shear


def find(subj_dir, pattern, recursive=True):
    hits = sorted(glob.glob(os.path.join(subj_dir, "**", pattern),
                            recursive=recursive))
    hits = [h for h in hits if "InverseWarp" not in os.path.basename(h)]
    return hits[0] if hits else None


def group_map(xlsx=STUDY_XLSX):
    x = pd.ExcelFile(xlsx)
    out = {}
    for sheet, label in SHEETS.items():
        t = x.parse(sheet)
        u = t["Will we used it?"].astype(str).str.strip().str.lower()
        for n in t.loc[u == "yes", "Subj"]:
            out[f"sub-{int(n)}"] = label
    return out


def kw_and_pairs(df, col, gcol="Group"):
    """Kruskal-Wallis across groups plus pairwise Mann-Whitney."""
    groups = [g for g in ("HC", "PPCS", "RHI") if (df[gcol] == g).any()]
    arrs = [df.loc[df[gcol] == g, col].dropna().values for g in groups]
    arrs = [a for a in arrs if len(a) > 2]
    if len(arrs) < 2:
        return None
    if np.ptp(np.concatenate(arrs)) < 1e-12:
        return f"  {col:<12s} constant across all subjects -- no test"
    H, p = stats.kruskal(*arrs)
    line = "  %-12s " % col + "  ".join(
        f"{g}={np.median(df.loc[df[gcol]==g, col].dropna()):.4f}" for g in groups)
    line += f"   KW p={p:.3f}"
    if p < 0.05:
        line += "  <-- differs by group"
    pairs = []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            a = df.loc[df[gcol] == groups[i], col].dropna()
            b = df.loc[df[gcol] == groups[j], col].dropna()
            if len(a) > 2 and len(b) > 2:
                _, pp = stats.mannwhitneyu(a, b)
                pairs.append(f"{groups[i]}v{groups[j]} p={pp:.3f}")
    return line + "\n      " + "  ".join(pairs)


# ==========================================
# Q1 -- ROTATION AT PVS vs WM
# ==========================================
def q1_subject(subj_dir, subj, rng):
    warp = find(subj_dir, "*from_b0_1Warp.nii.gz")
    pvs_p = find(subj_dir, "*_ses-01_PVS_postprocessed.nii*")
    seg_p = find(subj_dir, "*_ses-01_synthseg.nii.gz")
    if not (warp and pvs_p and seg_p):
        return None

    wimg = nib.load(warp)
    pimg = nib.load(pvs_p)
    simg = nib.load(seg_p)
    zooms = np.asarray(wimg.header.get_zooms()[:3], float)

    pvs = np.asarray(pimg.dataobj) > 0
    seg = np.asarray(simg.dataobj)
    wm = np.isin(seg, WM_LABELS)
    if not pvs.any() or not wm.any():
        return None

    pidx = np.argwhere(pvs)
    if len(pidx) > N_PVS_CAP:
        pidx = pidx[rng.choice(len(pidx), N_PVS_CAP, replace=False)]
    widx = np.argwhere(wm)
    if len(widx) > N_WM_SAMPLE:
        widx = widx[rng.choice(len(widx), N_WM_SAMPLE, replace=False)]

    pidx_w = map_indices(pimg, wimg, pidx)
    widx_w = map_indices(simg, wimg, widx)

    ang_p, ldj_p, fold_p = jacobian_at(wimg, pidx_w, zooms)
    ang_w, ldj_w, fold_w = jacobian_at(wimg, widx_w, zooms)
    if ang_p.size == 0 or ang_w.size == 0:
        return None

    abs_p = np.abs(ldj_p[np.isfinite(ldj_p)])
    abs_w = np.abs(ldj_w[np.isfinite(ldj_w)])
    if abs_p.size == 0 or abs_w.size == 0:
        return None

    wm_p90 = float(np.percentile(ang_w, 90))
    wm_ldj_p90 = float(np.percentile(abs_w, 90))
    vox_ml = float(np.prod(np.asarray(simg.header.get_zooms()[:3], float))) / 1000.0

    row = {"Subject": subj,
           "n_pvs": int(len(ang_p)), "n_wm": int(len(ang_w)),
           # --- rotation (affects u^T D u) ---
           "pvs_med": float(np.median(ang_p)),
           "pvs_p90": float(np.percentile(ang_p, 90)),
           "pvs_p99": float(np.percentile(ang_p, 99)),
           "wm_med": float(np.median(ang_w)),
           "wm_p90": wm_p90,
           "wm_p99": float(np.percentile(ang_w, 99)),
           "delta_med": float(np.median(ang_p) - np.median(ang_w)),
           # 0.10 under no spatial coincidence, by construction
           "enrichment_p90": float((ang_p > wm_p90).mean()),
           "frac_pvs_over_10deg": float((ang_p > 10).mean()),
           # --- |log det J| (proxy for interpolation smoothing) ---
           "pvs_absldj_med": float(np.median(abs_p)),
           "pvs_absldj_p90": float(np.percentile(abs_p, 90)),
           "wm_absldj_med": float(np.median(abs_w)),
           "wm_absldj_p90": wm_ldj_p90,
           "delta_absldj": float(np.median(abs_p) - np.median(abs_w)),
           "enrichment_absldj_p90": float((abs_p > wm_ldj_p90).mean()),
           # signed: negative = net local compression
           "pvs_ldj_med": float(np.nanmedian(ldj_p)),
           "wm_ldj_med": float(np.nanmedian(ldj_w)),
           "pvs_folded_frac": fold_p / max(len(ang_p), 1),
           "wm_folded_frac": fold_w / max(len(ang_w), 1),
           # --- head size ---
           "brain_ml": float(np.count_nonzero(seg) * vox_ml),
           "wm_ml": float(int(wm.sum()) * vox_ml)}

    icv_p = find(subj_dir, "*_ses-01_ICV_synthstrip.nii.gz")
    if icv_p:
        icv = np.asarray(nib.load(icv_p).dataobj) > 0
        row["icv_ml"] = float(int(icv.sum()) * vox_ml)
    else:
        row["icv_ml"] = np.nan

    # regional breakdown, if the template ROIs are present
    for key, tag in [("template_LBG", "LBG"), ("template_RBG", "RBG"),
                     ("template_MB", "MB")]:
        p = find(subj_dir, f"*_ses-01_{key}.nii.gz")
        if not p:
            continue
        m = np.asarray(nib.load(p).dataobj) > 0
        if m.shape != pvs.shape:
            continue
        sub = pvs & m
        if sub.sum() < 20:
            continue
        si = np.argwhere(sub)
        a, l, _ = jacobian_at(wimg, map_indices(pimg, wimg, si), zooms)
        if a.size:
            row[f"pvs_med_{tag}"] = float(np.median(a))
            la = np.abs(l[np.isfinite(l)])
            if la.size:
                row[f"pvs_absldj_{tag}"] = float(np.median(la))
    return row


def run_q1(root, subjects, gmap, rng):
    rows = []
    for s in subjects:
        try:
            r = q1_subject(os.path.join(root, s), s, rng)
        except Exception as exc:
            print(f"[{s}] Q1 failed: {exc}")
            continue
        if r is None:
            print(f"[{s}] Q1 skipped (missing warp/PVS/SynthSeg)")
            continue
        r["Group"] = gmap.get(s, "UNKNOWN")
        rows.append(r)
        print(f"[{s}] PVS med {r['pvs_med']:.2f} vs WM med {r['wm_med']:.2f} deg "
              f"| enrichment {r['enrichment_p90']:.3f} "
              f"| PVS>10deg {100*r['frac_pvs_over_10deg']:.1f}%", flush=True)
    if not rows:
        return None
    d = pd.DataFrame(rows)

    print("\n" + "=" * 72)
    print("Q1  LOCAL WARP ROTATION AT PVS vs CEREBRAL WM")
    print("=" * 72)
    print(f"  subjects: {len(d)}")
    print(f"  PVS  median {d.pvs_med.median():.2f}  p90 {d.pvs_p90.median():.2f}"
          f"  p99 {d.pvs_p99.median():.2f} deg")
    print(f"  WM   median {d.wm_med.median():.2f}  p90 {d.wm_p90.median():.2f}"
          f"  p99 {d.wm_p99.median():.2f} deg")
    w = None
    if len(d) > 5 and np.ptp(d.pvs_med - d.wm_med) > 1e-12:
        w = stats.wilcoxon(d.pvs_med, d.wm_med)
    print(f"  paired difference in medians: {d.delta_med.median():+.3f} deg"
          + (f"   Wilcoxon p={w.pvalue:.4f}" if w else ""))
    e = d.enrichment_p90.median()
    print(f"\n  ENRICHMENT (fraction of PVS voxels above that subject's WM p90)")
    print(f"    median {e:.3f}   range {d.enrichment_p90.min():.3f}"
          f"-{d.enrichment_p90.max():.3f}   (0.100 = no coincidence)")
    if e < 0.15:
        print("    -> PVS are NOT in the high-rotation tail. The p99 warp tail")
        print("       sits elsewhere; no action needed beyond a methods sentence.")
    elif e < 0.25:
        print("    -> mild enrichment. Worth a limitations sentence; unlikely to")
        print("       drive a group effect unless it also differs by group (below).")
    else:
        print("    -> substantial enrichment. The warp is roughest where you are")
        print("       measuring; reorientation or a rotation covariate is needed.")
    print(f"\n  PVS voxels above 10 deg: median {100*d.frac_pvs_over_10deg.median():.1f}%")

    print("\n  |log det J| -- local volume change, the interpolation-smoothing proxy")
    print(f"    PVS median {d.pvs_absldj_med.median():.4f}  p90 "
          f"{d.pvs_absldj_p90.median():.4f}")
    print(f"    WM  median {d.wm_absldj_med.median():.4f}  p90 "
          f"{d.wm_absldj_p90.median():.4f}")
    wl = None
    if len(d) > 5 and np.ptp(d.pvs_absldj_med - d.wm_absldj_med) > 1e-12:
        wl = stats.wilcoxon(d.pvs_absldj_med, d.wm_absldj_med)
    print(f"    paired difference: {d.delta_absldj.median():+.4f}"
          + (f"   Wilcoxon p={wl.pvalue:.4f}" if wl else ""))
    el = d.enrichment_absldj_p90.median()
    print(f"    enrichment: {el:.3f}  range {d.enrichment_absldj_p90.min():.3f}"
          f"-{d.enrichment_absldj_p90.max():.3f}   (0.100 = no coincidence)")
    print(f"    signed log det J (PVS): {d.pvs_ldj_med.median():+.4f}  "
          f"(negative = net compression -> averaging)")
    fold = d.wm_folded_frac.median()
    print(f"    folded voxels (det J <= 0) in WM: {100*fold:.3f}%"
          + ("   [!] non-trivial folding" if fold > 0.001 else ""))

    if d.Group.nunique() > 1:
        print("\n  does |log det J| differ by group? (this is the covariate that")
        print("  matters -- it is the mechanism by which a high-MD PVS voxel gets")
        print("  averaged toward its surroundings)")
        for c in ["pvs_absldj_med", "wm_absldj_med", "enrichment_absldj_p90"]:
            r = kw_and_pairs(d, c)
            if r:
                print(r)

    regs = [c for c in d.columns if c.startswith("pvs_med_")]
    if regs:
        print("\n  by region (median across subjects):")
        for c in regs:
            print(f"    {c.replace('pvs_med_',''):5s} {d[c].median():.2f} deg")

    if d.Group.nunique() > 1:
        print("\n  does the PVS-region rotation differ by group?")
        for c in ["pvs_med", "enrichment_p90", "frac_pvs_over_10deg"]:
            r = kw_and_pairs(d, c)
            if r:
                print(r)
    return d


# ==========================================
# Q2 -- AFFINE SCALING BY GROUP
# ==========================================
def q2_subject(subj_dir, subj):
    aff = find(subj_dir, "*from_b0_0GenericAffine.mat")
    if not aff:
        return None
    m = loadmat(aff)
    key = next((k for k in m if k.startswith("AffineTransform")), None)
    if key is None:
        return None
    p = np.asarray(m[key]).ravel()
    A = p[:9].reshape(3, 3)
    R, scale, shear = polar_stretch(A)
    ang = float(np.degrees(np.arccos(
        np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0))))
    return {"Subject": subj, "scale_x": float(scale[0]),
            "scale_y": float(scale[1]), "scale_z": float(scale[2]),
            "shear": shear, "affine_deg": ang,
            "aniso": float(scale.max() / scale.min()),
            "volume_ratio": float(np.abs(np.linalg.det(A)))}


def run_q2(root, subjects, gmap):
    rows = []
    for s in subjects:
        r = q2_subject(os.path.join(root, s), s)
        if r is None:
            print(f"[{s}] Q2 skipped (no affine)")
            continue
        r["Group"] = gmap.get(s, "UNKNOWN")
        rows.append(r)
    if not rows:
        return None
    d = pd.DataFrame(rows)

    print("\n" + "=" * 72)
    print("Q2  PER-AXIS AFFINE SCALING (diagonal of the stretch factor P)")
    print("=" * 72)
    print("  NOTE: this is axis-associated. SVD singular values are sorted and")
    print("        cannot be read as x/y/z -- an earlier check made that error.")
    print(f"\n  subjects: {len(d)}  |  groups: "
          + ", ".join(f"{g}={int((d.Group==g).sum())}" for g in d.Group.unique()))
    print("\n  cohort medians:")
    for c in ["scale_x", "scale_y", "scale_z", "aniso", "shear",
              "volume_ratio", "affine_deg"]:
        print(f"    {c:13s} {d[c].median():.4f}   "
              f"[{d[c].min():.4f}, {d[c].max():.4f}]")

    if d.Group.nunique() > 1:
        print("\n  by group (median per group, Kruskal-Wallis, then pairwise):")
        for c in ["scale_x", "scale_y", "scale_z", "aniso", "shear",
                  "volume_ratio", "affine_deg"]:
            r = kw_and_pairs(d, c)
            if r:
                print(r)
        sig = []
        for c in ["scale_x", "scale_y", "scale_z", "aniso", "volume_ratio"]:
            arrs = [d.loc[d.Group == g, c].dropna().values
                    for g in ("HC", "PPCS", "RHI") if (d.Group == g).any()]
            arrs = [a for a in arrs if len(a) > 2]
            if len(arrs) > 1 and np.ptp(np.concatenate(arrs)) > 1e-12 \
                    and stats.kruskal(*arrs)[1] < 0.05:
                sig.append(c)
        print()
        if sig:
            print(f"  -> {sig} differ by group. The PVS mask therefore aligns to")
            print("     the DKI data differently across groups. Report this and")
            print("     covary it, or show the tier results are unchanged when it")
            print("     is included.")
        else:
            print("  -> no group difference in affine geometry. This is the null")
            print("     you want: state it in one sentence in the rebuttal, since")
            print("     it forecloses a PVS-mask-to-DKI alignment confound.")
    else:
        print("\n  [!] no group labels resolved -- check the study spreadsheet.")
    return d


# ==========================================
# Q3 -- HEAD SIZE AS THE UPSTREAM VARIABLE
# ==========================================
def run_q3(d1, d2):
    """Is head size the common cause of the group difference in registration?

    The RHI group is age-matched but has smaller heads. That ordering matters
    causally. If head size drives how much the warp deforms, then registration
    geometry is a MEDIATOR on a purely technical path

        head size -> warp deformation -> interpolation smoothing -> measured MD

    and adjusting for it is the right thing to do, because nothing on that path
    is biology. If instead head size and deformation are independent, they are
    two separate covariates and both belong in the model.

    Adjusting for head size directly is also more defensible to a reviewer than
    adjusting for a registration diagnostic: intracranial volume is a standard
    covariate, whereas "median local Jacobian rotation" needs three paragraphs
    of explanation.
    """
    if d1 is None or d2 is None:
        print("\n[Q3] needs both Q1 and Q2 -- run without --q1/--q2 to get both")
        return None
    d = d1.merge(d2.drop(columns=[c for c in ("Group",) if c in d2]),
                 on="Subject", how="inner")
    print("\n" + "=" * 72)
    print("Q3  HEAD SIZE, AND WHETHER IT EXPLAINS THE REGISTRATION DIFFERENCES")
    print("=" * 72)
    if "icv_ml" not in d or d.icv_ml.isna().all():
        print("  [!] no ICV masks found; falling back to brain_ml from SynthSeg")
        d["icv_ml"] = d["brain_ml"]

    print("\n  head size by group:")
    for c in ["icv_ml", "brain_ml", "wm_ml"]:
        r = kw_and_pairs(d, c)
        if r:
            print(r)

    print("\n  does head size predict the registration geometry?")
    targets = ["wm_med", "pvs_med", "wm_absldj_med", "pvs_absldj_med",
               "scale_x", "scale_y", "scale_z", "volume_ratio", "affine_deg"]
    for c in targets:
        if c not in d:
            continue
        ok = d[["icv_ml", c]].dropna()
        if len(ok) < 8 or np.ptp(ok[c]) < 1e-12 or np.ptp(ok.icv_ml) < 1e-12:
            print(f"    spearman(icv_ml, {c:16s}) = constant / too few -- skipped")
            continue
        r, p = stats.spearmanr(ok.icv_ml, ok[c])
        flag = "  <-- head size tracks this" if p < 0.05 else ""
        print(f"    spearman(icv_ml, {c:16s}) = {r:+.3f}  p={p:.4f}{flag}")

    if _HAVE_SM and d.Group.nunique() > 1:
        print("\n  does the GROUP difference in registration survive adjusting")
        print("  for head size? (if it does not, head size is the whole story)")
        dd = d[d.Group.isin(["HC", "RHI"])].copy()
        dd["RHI"] = (dd.Group == "RHI").astype(int)
        for c in ["wm_med", "wm_absldj_med", "scale_x"]:
            if c not in dd:
                continue
            sub = dd[[c, "RHI", "icv_ml"]].dropna()
            if len(sub) < 10:
                continue
            sub = sub.rename(columns={c: "_y"})
            r1 = smf.ols("_y ~ RHI", data=sub).fit()
            r2 = smf.ols("_y ~ RHI + icv_ml", data=sub).fit()
            shrink = (100 * (1 - r2.params.RHI / r1.params.RHI)
                      if abs(r1.params.RHI) > 1e-12 else np.nan)
            print(f"    {c:16s} beta_RHI {r1.params.RHI:+.4f} (p={r1.pvalues.RHI:.4f})"
                  f"  ->  {r2.params.RHI:+.4f} (p={r2.pvalues.RHI:.4f})"
                  f"   shrink {shrink:.0f}%")
    return d


# ==========================================
# Q4 -- OUTCOMES ADJUSTED FOR THE FULL CHAIN
# ==========================================
def run_q4(cov, xlsx, sheet="Case_A_Near_PVS"):
    """Group effect on each metric, adjusted stepwise for the technical chain.

    Steps are nested deliberately, so you can see WHERE an effect dies:
        y ~ Group
        y ~ Group + icv_ml                     head size
        y ~ Group + icv_ml + wm_absldj_med     + interpolation smoothing
        y ~ Group + icv_ml + wm_absldj_med + scale_x   + affine geometry

    A metric whose beta is stable across all four is orthogonal to the
    acquisition-and-registration chain. A metric whose beta collapses is not
    necessarily an artifact -- with n ~ 40 and correlated covariates the model
    often cannot separate them -- but it can no longer be reported as a clean
    group difference either.
    """
    if not _HAVE_SM:
        print("\n[Q4] statsmodels not available; skipping")
        return None
    if not os.path.exists(xlsx):
        print(f"\n[Q4] {xlsx} not found; skipping")
        return None
    t = pd.read_excel(xlsx, sheet_name=sheet)
    mets = [c for c in t.columns if c.endswith("_Original")]
    if not mets or "Subject" not in t:
        print(f"\n[Q4] no *_Original columns in {sheet}; skipping")
        return None
    sub = t.groupby("Subject")[mets].median()
    d = cov.set_index("Subject").join(sub, how="inner").reset_index()
    d = d[d.Group.isin(["HC", "RHI"])].copy()
    d["RHI"] = (d.Group == "RHI").astype(int)

    print("\n" + "=" * 72)
    print("Q4  RHI vs HC, ADJUSTED STEPWISE FOR THE TECHNICAL CHAIN")
    print("=" * 72)
    print(f"  HC n={(d.RHI==0).sum()}  RHI n={(d.RHI==1).sum()}  "
          f"metrics={len(mets)}")
    models = [("unadjusted", "_y ~ RHI"),
              ("+ head size", "_y ~ RHI + icv_ml"),
              ("+ |log det J|", "_y ~ RHI + icv_ml + wm_absldj_med"),
              ("+ scale_x", "_y ~ RHI + icv_ml + wm_absldj_med + scale_x")]
    out = []
    print(f"\n  {'metric':<24s}" + "".join(f"{n:>16s}" for n, _ in models))
    for m in sorted(mets):
        cols = ["RHI", "icv_ml", "wm_absldj_med", "scale_x", m]
        s = d[[c for c in cols if c in d]].dropna()
        if len(s) < 12:
            continue
        # rescale for numerical conditioning only; p-values are unaffected
        sc = 1e3 if s[m].abs().median() < 1e-2 else 1.0
        s = s.assign(_y=s[m] * sc)
        row, ps = {"metric": m}, []
        for name, f in models:
            try:
                r = smf.ols(f, data=s).fit()
                ps.append(f"{r.params.RHI:+.3f}/{r.pvalues.RHI:.3f}")
                row[name] = float(r.pvalues.RHI)
                row[name + "_beta"] = float(r.params.RHI)
            except Exception:
                ps.append("     --     ")
                row[name] = np.nan
        mark = ""
        if np.isfinite(row.get(models[0][0], np.nan)) and \
                np.isfinite(row.get(models[-1][0], np.nan)):
            if row[models[0][0]] < 0.05 and row[models[-1][0]] < 0.05:
                mark = "  SURVIVES"
            elif row[models[0][0]] < 0.05:
                mark = "  dies"
        print(f"  {m:<24s}" + "".join(f"{v:>16s}" for v in ps) + mark)
        out.append(row)
    print("\n  cells are beta_RHI/p. 'SURVIVES' = significant before and after.")
    print("  Caveat: if a covariate is itself non-significant while beta shrinks,")
    print("  group and covariate are collinear and the data cannot separate them.")
    print("  Report that as 'cannot be distinguished from', not 'is an artifact'.")
    return pd.DataFrame(out)


# ==========================================
# MAIN
# ==========================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=ROOT_DIR)
    ap.add_argument("--out-prefix", default="registration_confounds")
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--q1", action="store_true", help="rotation at PVS vs WM")
    ap.add_argument("--q2", action="store_true", help="affine scaling by group")
    ap.add_argument("--q3", action="store_true", help="head size (needs q1+q2)")
    ap.add_argument("--metrics-xlsx", default="",
                    help="v4/v5 spreadsheet; enables Q4 adjusted outcome models")
    ap.add_argument("--metrics-sheet", default="Case_A_Near_PVS")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    any_flag = a.q1 or a.q2 or a.q3
    do_q1 = a.q1 or a.q3 or not any_flag
    do_q2 = a.q2 or a.q3 or not any_flag
    do_q3 = a.q3 or not any_flag

    try:
        gmap = group_map()
        print(f"group labels: {len(gmap)} subjects from the study spreadsheet")
    except Exception as exc:
        print(f"[!] group labels unavailable ({exc}); group tests will be skipped")
        gmap = {}

    if a.subjects:
        subjects = a.subjects
    else:
        subjects = sorted(s for s in os.listdir(a.root)
                          if s.startswith("sub-")
                          and os.path.isdir(os.path.join(a.root, s)))
        if gmap:
            subjects = [s for s in subjects if s in gmap]
    print(f"processing {len(subjects)} subjects\n")

    rng = np.random.default_rng(a.seed)
    d1 = d2 = None
    if do_q1:
        d1 = run_q1(a.root, subjects, gmap, rng)
        if d1 is not None:
            p = f"{a.out_prefix}_q1_rotation_at_pvs.csv"
            d1.to_csv(p, index=False)
            print(f"\n-> {p}")
    if do_q2:
        d2 = run_q2(a.root, subjects, gmap)
        if d2 is not None:
            p = f"{a.out_prefix}_q2_affine_scaling.csv"
            d2.to_csv(p, index=False)
            print(f"\n-> {p}")
    if do_q3:
        cov = run_q3(d1, d2)
        if cov is not None:
            keep = [c for c in ("Subject", "Group", "icv_ml", "brain_ml",
                                "wm_ml", "wm_med", "pvs_med", "wm_absldj_med",
                                "pvs_absldj_med", "enrichment_p90",
                                "enrichment_absldj_p90", "scale_x", "scale_y",
                                "scale_z", "aniso", "shear", "volume_ratio",
                                "affine_deg") if c in cov]
            p = f"{a.out_prefix}_subject_covariates.csv"
            cov[keep].to_csv(p, index=False)
            print(f"\n-> {p}   (join this on Subject in the stats script)")
            if a.metrics_xlsx:
                q4 = run_q4(cov[keep], a.metrics_xlsx, a.metrics_sheet)
                if q4 is not None:
                    p = f"{a.out_prefix}_q4_adjusted_models.csv"
                    q4.to_csv(p, index=False)
                    print(f"\n-> {p}")


if __name__ == "__main__":
    main()
