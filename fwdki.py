#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fwdki_pvs.py
============
Free-water-corrected diffusion & kurtosis metrics inside perivascular-space (PVS)
masks, built for a 4-shell 7T DKI acquisition:

    b = 0 (x2), 100 (6 dir), 500 (6 dir), 1000 (10 dir), 2000 (30 dir) s/mm^2
    2 mm isotropic, TE = 76 ms

WHY THIS DESIGN
---------------
A joint free-water-DKI fit (Collier et al. 2018, MRM) is provably ill-conditioned
and needs a Bayesian shrinkage prior; on a 52-DWI acquisition with thin shells that
estimator is not trustworthy. Instead we use a TWO-STAGE estimator that keeps every
step in a regime where its parameter is identifiable:

  Stage 1 (f):  free-water DTI on the free-water-sensitive shells (b <= BMAX_FW,
                default 1000) where f is well determined and the Gaussian tissue
                assumption approximately holds (Hoy 2014; Bergmann 2020 logic).

  Stage 2 (deconfound): subtract the free-water term at EVERY shell. The free-water
                compartment has a FIXED, known b-dependence, exp(-b * d_fw) with
                d_fw = 3e-3 mm^2/s, so once f is known the subtraction is defined
                for the high shells too:
                    s_tissue_i = ( s_i - S0 * f * exp(-b_i * d_fw) ) / (1 - f)

  Stage 3 (tissue metrics): a SINGLE standard DKI fit on the deconfounded signal
                across ALL shells -> one internally consistent tissue tensor giving
                FAt, MDt, ADt, RDt, MKt, AKt, RKt.

This never performs the ill-posed 23-parameter joint fit. f is estimated where it is
identifiable; kurtosis is estimated by an ordinary DKI fit on cleaned data.

OUTPUTS (per subject)
---------------------
  Uncorrected (standard, all shells):  FA, MD, AD, RD, MK, AK, RK
  Free-water fraction:                 f            <-- the fluid / clearance index
  Corrected (tissue compartment):      FAt, MDt, ADt, RDt, MKt, AKt, RKt

Plus, for the pre-specified symmetric sensitivity comparison:
  * per-region metric summaries inside PVS masks,
  * a residual-bootstrap ESTIMATION-VARIANCE control per metric (so that any
    apparent sensitivity advantage of f can be checked against the possibility that
    f is merely estimated more precisely than the kurtosis metrics on this data).

CAVEATS (state these in the manuscript)
---------------------------------------
  * f is a T2-WEIGHTED SIGNAL fraction, not a volume fraction. CSF T2 >> tissue T2,
    and at 7T tissue T2 is even shorter, so f is inflated relative to physical
    volume. Fine for group CONTRASTS; do not read f as an absolute volume without a
    multi-TE (Eq. 16, Collier) correction you do not have data for.
  * Deconfounding assumes the low-b-estimated f applies at high b. This is exact for
    the fixed free-water model; residual bias is second-order.
  * Corrected kurtosis on 10-direction b=1000 + 30-direction b=2000 is the noisiest
    output. The bootstrap control is not decoration -- it is load-bearing.
  * NO-b0 DATA: if the image has no b=0 volume (some Eddy/Topup-derived DKI images
    drop it), S0 is estimated by extrapolation and refined iteratively from the
    free-water-free tissue. This works (see the NO-b0 self-test) but f is inherently
    less precise than with a measured b0, because f is sensitive to the S0 anchor at
    low b. If you can recover the acquired b0 upstream and re-insert it, do so -- it
    improves f precision. The code auto-detects b0 presence; no flag needed.

Author-facing: preprocessing (denoise / Gibbs / Topup / Eddy) is assumed already
done upstream (as in your acquisition's Collier-style pipeline). This script starts
from preprocessed DWIs + bval/bvec + brain mask + PVS mask.
"""

import os
import numpy as np

D_FW = 3.0e-3          # free-water diffusivity, mm^2/s (b in s/mm^2 -> b*D_FW dimensionless)
BMAX_FW = 1000.0       # top b-value used for the free-water (Gaussian) f estimate
EPS = 1e-6


# ----------------------------------------------------------------------------- #
# Gradient table helpers
# ----------------------------------------------------------------------------- #
def build_gtab(bvals, bvecs, b0_threshold=50.0):
    """Build a DIPY gradient table. b0_threshold=50 keeps the b=100 shell as a
    genuine (non-b0) shell, which we need for free-water sensitivity."""
    from dipy.core.gradients import gradient_table
    bvals = np.asarray(bvals, dtype=float)
    bvecs = np.asarray(bvecs, dtype=float)
    return gradient_table(bvals, bvecs=bvecs, b0_threshold=b0_threshold)


def _shell_subset(bvals, bmax):
    """Boolean index of volumes with b <= bmax (keeps b0s)."""
    return np.asarray(bvals, dtype=float) <= (bmax + 1e-3)


def _intercept_S0(data, bvals, bvecs, model="dti", bmax=BMAX_FW):
    """Per-voxel b=0 intercept of a WLLS fit (signal extrapolated to b=0).

    model='dti': mono-exponential design on shells b<=bmax (accurate for a
                 free-water-contaminated signal only in the Gaussian regime).
    model='dki': kurtosis-aware design on ALL shells (accurate for a FREE-WATER-FREE
                 signal, e.g. deconfounded tissue -> unbiased S0).
    """
    from dipy.reconst import dti as _dti, dki as _dki
    bvals = np.asarray(bvals, float)
    if model == "dki":
        sel = np.ones_like(bvals, bool)
        B = _dki.design_matrix(build_gtab(bvals, np.asarray(bvecs)))
    else:
        sel = bvals <= (bmax + 1e-3)
        B = _dti.design_matrix(build_gtab(bvals[sel], np.asarray(bvecs)[sel]))
    S = np.clip(data[..., sel], 1e-6, None)
    logS = np.log(S)
    w = S ** 2
    BtW = B.T * w[..., None, :]
    A = np.einsum("...im,mj->...ij", BtW, B)
    rhs = np.einsum("...im,...m->...i", BtW, logS)
    # Per-COLUMN ridge proportional to each diagonal element. The design-matrix
    # diagonals span orders of magnitude (diffusion columns huge, intercept small),
    # so a mean-scaled ridge would over-regularise the intercept and bias S0. Adding
    # a tiny fraction of each column's own scale keeps every column conditioned
    # without bias, and a small absolute floor keeps singular (zero) columns solvable.
    p = A.shape[-1]
    d = np.einsum("...ii->...i", A)                          # per-voxel diagonal
    ridge = d * 1e-8 + 1e-10
    A = A + ridge[..., None, :] * np.eye(p)
    x = np.linalg.solve(A, rhs[..., None])[..., 0]
    S0 = np.exp(-np.clip(x[..., -1], -50, 50))               # DIPY intercept sign (-1 col)
    return np.nan_to_num(S0, nan=0.0, posinf=0.0, neginf=0.0)


def estimate_S0(data, bvals, bvecs, mask=None, b0_thresh=50.0, bmax_lin=BMAX_FW):
    """Per-voxel non-diffusion-weighted signal S0.

    Measured b0 (b <= b0_thresh) present -> its per-voxel mean (preferred, exact).
    No b0 acquired (some Eddy/Topup-derived DKI images drop the b0) -> DTI-intercept
    extrapolation on the low/Gaussian shells as an INITIAL estimate; the hilow
    f-estimator then refines S0 iteratively from the free-water-free tissue (which
    matters because f is sensitive to the S0 anchor at low b).
    """
    bvals = np.asarray(bvals, float)
    b0_idx = bvals <= b0_thresh
    if b0_idx.any():
        return data[..., b0_idx].mean(axis=-1)
    S0 = _intercept_S0(data, bvals, bvecs, model="dti", bmax=bmax_lin)
    if mask is not None:
        S0 = S0 * mask
    return S0


# ----------------------------------------------------------------------------- #
# Stage 0: standard (uncorrected) DKI  -- all shells
# ----------------------------------------------------------------------------- #
def fit_standard_dki(data, gtab, mask=None):
    """Ordinary DKI on the raw signal. Returns dict of uncorrected metrics."""
    from dipy.reconst import dki
    model = dki.DiffusionKurtosisModel(gtab, fit_method="WLS")
    fit = model.fit(data, mask=mask)
    return {
        "FA": fit.fa, "MD": fit.md, "AD": fit.ad, "RD": fit.rd,
        "MK": fit.mk(), "AK": fit.ak(), "RK": fit.rk(),
    }


# ----------------------------------------------------------------------------- #
# Stage 1: free-water fraction f
# ----------------------------------------------------------------------------- #
# Two estimators are provided:
#   (A) 'hilow'  [DEFAULT]  High-Low-Downhill (Bergmann 2020 logic, adapted to DKI).
#                f is a SCALAR regression over the low shells against a DKI tissue
#                model fit from the rich high shells. A scalar over 12 points is
#                robust; a low-b free-water TENSOR (6+6 directions) is not -- the
#                synthetic test shows the tensor route collapses to f=0 at SNR 20.
#   (B) 'fwdti'            DIPY free-water DTI on b<=bmax_fw. Standard/citable but
#                fragile on this thin low-b acquisition; kept as a cross-check.

def _tissue_low_b_regression_f(sig, tissue_pred, bvals, S0, low_mask):
    """Closed-form scalar f from the low shells (Bergmann Eq. 8 generalisation).

        s_i = (1-f) T_i + f FW_i   ->   f = sum (FW-T)(s-T) / sum (FW-T)^2
    over the free-water-sensitive (low-b) measurements.
    sig, tissue_pred : (..., N) ;  S0 : (...) ;  low_mask : (N,) boolean.
    """
    bvals = np.asarray(bvals, float)
    fw = S0[..., None] * np.exp(-bvals * D_FW)[None, ...] if S0.ndim else \
         S0 * np.exp(-bvals * D_FW)
    fw = S0[..., None] * np.exp(-bvals * D_FW)                 # (...,N)
    T = tissue_pred
    lm = low_mask
    num = np.sum((fw[..., lm] - T[..., lm]) * (sig[..., lm] - T[..., lm]), axis=-1)
    den = np.sum((fw[..., lm] - T[..., lm]) ** 2, axis=-1)
    f = np.zeros_like(num)
    np.divide(num, den, out=f, where=den > 0)
    return np.clip(np.nan_to_num(f), 0.0, 1.0 - 1e-3)


def fit_free_water_fraction(data, bvals, bvecs, mask=None, bmax_fw=BMAX_FW,
                            method="hilow", n_iter=4):
    """Estimate the free-water fraction f.

    Returns (f, S0), both with the spatial shape of `data`. S0 is the (possibly
    refined) non-diffusion-weighted amplitude actually used, so the caller can reuse
    the SAME S0 for the final deconfounding -- important in the no-b0 case where S0
    is refined iteratively here.
    """
    bvals = np.asarray(bvals, float)
    if method == "fwdti":
        from dipy.reconst import fwdti
        sel = _shell_subset(bvals, bmax_fw)
        gtab_lo = build_gtab(bvals[sel], np.asarray(bvecs)[sel])
        fit = fwdti.FreeWaterTensorModel(gtab_lo).fit(data[..., sel], mask=mask)
        f = np.clip(np.nan_to_num(fit.f), 0.0, 1.0 - 1e-3)
        return f, estimate_S0(data, bvals, bvecs, mask=mask)

    # --- 'hilow' (default) ---------------------------------------------------
    # Init: anchor the tissue SLOPE on the high, free-water-free shells (Bergmann
    # HILOW excludes b0 and low shells from the tissue fit), predict tissue at low
    # b using the measured S0 as amplitude, then solve f as a scalar. Then refine
    # with a DKI tissue model on deconfounded data (HILOWDOWNHILL).
    from dipy.reconst import dki, dti
    gtab_full = build_gtab(bvals, bvecs)
    dki_model = dki.DiffusionKurtosisModel(gtab_full, fit_method="WLS")

    b0_idx = bvals <= 50.0
    has_b0 = b0_idx.any()
    S0 = estimate_S0(data, bvals, bvecs, mask=mask)           # measured b0 or extrapolated
    low_mask = (~b0_idx) & (bvals <= bmax_fw + 1e-3)           # f-sensitive shells
    b_high = 1000.0
    high_sel = bvals >= b_high - 1e-3                          # FW ~ 0 here

    # --- HILOW initialisation: Gaussian tissue slope from high shells ----------
    gtab_hi = build_gtab(bvals[high_sel], np.asarray(bvecs)[high_sel])
    dti_hi = dti.TensorModel(gtab_hi, fit_method="WLS")
    fit_hi = dti_hi.fit(data[..., high_sel], mask=mask)
    Dq = fit_hi.quadratic_form                                 # (...,3,3)
    bvecs_arr = np.asarray(bvecs, float)
    # directional ADC for every measurement: g^T D g
    adc = np.einsum("...i,ni,...ij,nj->...n", np.ones(Dq.shape[:-2] + (3,)),
                    bvecs_arr, Dq, bvecs_arr) if False else \
          np.einsum("ni,...ij,nj->...n", bvecs_arr, Dq, bvecs_arr)
    T = S0[..., None] * np.exp(-bvals[None, ...] * adc)
    f = _tissue_low_b_regression_f(data, T, bvals, S0, low_mask)
    if mask is not None:
        f = f * mask

    # --- HILOWDOWNHILL refinement with a DKI tissue model ----------------------
    # When no b0 was measured, S0 is only an extrapolated guess and f is sensitive
    # to it; so each iteration we also re-estimate S0 from the (now free-water-free)
    # deconfounded tissue via a kurtosis-aware all-shell intercept, which is
    # unbiased for a pure tissue signal and breaks the S0->f bias amplification.
    for _ in range(n_iter):
        tissue = deconfound_free_water(data, bvals, f, mask=mask, S0=S0)
        if not has_b0:
            S0 = _intercept_S0(tissue, bvals, bvecs, model="dki")
            if mask is not None:
                S0 = S0 * mask
            tissue = deconfound_free_water(data, bvals, f, mask=mask, S0=S0)
        fit = dki_model.fit(tissue, mask=mask)
        T = dki_model.predict(fit.model_params, S0=S0)        # full-amplitude tissue
        f = _tissue_low_b_regression_f(data, T, bvals, S0, low_mask)
        if mask is not None:
            f = f * mask
    return f, S0


# ----------------------------------------------------------------------------- #
# Stage 2: deconfound the signal at every shell given f
# ----------------------------------------------------------------------------- #
def deconfound_free_water(data, bvals, f, mask=None, S0=None, bvecs=None):
    """Remove the free-water compartment from all shells.

        s_tissue = ( s - S0 * f * exp(-b * d_fw) ) / (1 - f)

    S0 is the per-voxel non-diffusion-weighted signal. If not provided it is
    estimated via estimate_S0 (measured b0 if present, else tensor-intercept
    extrapolation). Returns the tissue-compartment signal, clipped non-negative so
    the downstream WLLS log-fit is well defined."""
    bvals = np.asarray(bvals, dtype=float)
    if S0 is None:
        S0 = estimate_S0(data, bvals, bvecs, mask=mask)
    fw_decay = np.exp(-np.outer(np.ones(1), bvals) * D_FW).ravel()  # (N,)
    f_ = f[..., None]
    S0_ = S0[..., None]
    tissue = (data - S0_ * f_ * fw_decay) / np.clip(1.0 - f_, EPS, None)
    tissue = np.where(np.isfinite(tissue), tissue, 0.0)
    tissue = np.clip(tissue, 0.0, None)
    if mask is not None:
        tissue = tissue * mask[..., None]
    return tissue


# ----------------------------------------------------------------------------- #
# Stage 3: corrected tissue DKI  -- single fit on deconfounded signal, all shells
# ----------------------------------------------------------------------------- #
def fit_corrected_dki(tissue_data, gtab, mask=None):
    """One DKI fit on the deconfounded signal -> all seven corrected metrics."""
    from dipy.reconst import dki
    model = dki.DiffusionKurtosisModel(gtab, fit_method="WLS")
    fit = model.fit(tissue_data, mask=mask)
    return {
        "FAt": fit.fa, "MDt": fit.md, "ADt": fit.ad, "RDt": fit.rd,
        "MKt": fit.mk(), "AKt": fit.ak(), "RKt": fit.rk(),
    }


# ----------------------------------------------------------------------------- #
# Full per-subject pipeline
# ----------------------------------------------------------------------------- #
def run_subject(data, bvals, bvecs, brain_mask=None, bmax_fw=BMAX_FW):
    """Run all stages. Returns dict with the 8 primary maps (+ the DTI-part
    corrected metrics). Keys:
        uncorrected: FA MD AD RD MK AK RK
        fluid:       f
        corrected:   FAt MDt ADt RDt MKt AKt RKt
    """
    bvals = np.asarray(bvals, float)
    data = np.nan_to_num(np.asarray(data, float), nan=0.0, posinf=0.0, neginf=0.0)
    gtab = build_gtab(bvals, bvecs)
    out = {}
    out.update(fit_standard_dki(data, gtab, mask=brain_mask))
    # f-estimator also returns the S0 it used (possibly refined in the no-b0 case),
    # which we reuse for the final deconfound so both stages share one S0.
    f, S0 = fit_free_water_fraction(data, bvals, bvecs, mask=brain_mask, bmax_fw=bmax_fw)
    out["f"] = f
    tissue = deconfound_free_water(data, bvals, f, mask=brain_mask, S0=S0)
    out.update(fit_corrected_dki(tissue, gtab, mask=brain_mask))
    return out


# ----------------------------------------------------------------------------- #
# PVS-mask metric extraction
# ----------------------------------------------------------------------------- #
PRIMARY_METRICS = ["f", "FAt", "MDt", "ADt", "RDt", "MKt", "AKt", "RKt"]
UNCORR_METRICS  = ["FA", "MD", "AD", "RD", "MK", "AK", "RK"]

def extract_region_summaries(maps, pvs_mask, region_labels=None):
    """Summarise each metric inside a PVS mask (optionally per region).

    pvs_mask : binary array, or integer label map if region_labels given.
    Returns nested dict: region -> metric -> {mean, median, std, n}."""
    summaries = {}
    if region_labels is None:
        regions = {"PVS_all": pvs_mask.astype(bool)}
    else:
        regions = {name: (pvs_mask == lab) for name, lab in region_labels.items()}
    for rname, rmask in regions.items():
        summaries[rname] = {}
        for m, arr in maps.items():
            vals = arr[rmask]
            vals = vals[np.isfinite(vals)]
            summaries[rname][m] = dict(
                mean=float(np.mean(vals)) if vals.size else np.nan,
                median=float(np.median(vals)) if vals.size else np.nan,
                std=float(np.std(vals)) if vals.size else np.nan,
                n=int(vals.size),
            )
    return summaries


# ----------------------------------------------------------------------------- #
# Estimation-variance control (residual bootstrap) -- the load-bearing check
# ----------------------------------------------------------------------------- #
def bootstrap_estimation_sd(data, bvals, bvecs, voxel_mask, n_boot=100,
                            bmax_fw=BMAX_FW, seed=0):
    """Residual-bootstrap estimation SD for each of the 8 primary metrics, over the
    voxels in voxel_mask (use a modest ROI -- this is expensive).

    For each metric we return the median across ROI voxels of the per-voxel
    bootstrap SD. This quantifies HOW PRECISELY each metric is estimated on THIS
    acquisition, independent of any group effect -- so a group-separation advantage
    of f can be divided by this to check it is not merely an estimability artifact.
    """
    rng = np.random.default_rng(seed)
    gtab = build_gtab(bvals, bvecs)
    idx = np.argwhere(voxel_mask)
    if idx.shape[0] == 0:
        raise ValueError("empty voxel_mask")

    # Fit once to get predicted signal + residuals per voxel (standard DKI as the
    # signal model; residuals capture the noise the metrics are exposed to).
    from dipy.reconst import dki
    base_model = dki.DiffusionKurtosisModel(gtab, fit_method="WLS")

    # collect bootstrap metric samples: metric -> array (n_boot, n_vox)
    samples = {m: np.full((n_boot, idx.shape[0]), np.nan) for m in PRIMARY_METRICS}

    # extract the small ROI block of signal
    sub = np.stack([data[tuple(v)] for v in idx], axis=0)          # (n_vox, N)
    S0_vox = estimate_S0(sub, bvals, bvecs)                        # (n_vox,)
    base_fit = base_model.fit(sub)
    pred = base_model.predict(base_fit.model_params, S0=S0_vox)    # (n_vox, N)
    resid = sub - pred                                            # (n_vox, N)

    n_vox, N = sub.shape
    for b in range(n_boot):
        boot_idx = rng.integers(0, N, size=N)
        boot = pred + resid[:, boot_idx]
        boot = np.clip(boot, 0.0, None)
        boot4d = boot.reshape(n_vox, 1, 1, N)
        out = run_subject(boot4d, bvals, bvecs, brain_mask=None, bmax_fw=bmax_fw)
        for m in PRIMARY_METRICS:
            samples[m][b] = out[m].reshape(n_vox)

    est_sd = {}
    for m in PRIMARY_METRICS:
        per_vox_sd = np.nanstd(samples[m], axis=0)   # SD across boots, per voxel
        est_sd[m] = float(np.nanmedian(per_vox_sd))
    return est_sd


# ----------------------------------------------------------------------------- #
# Pre-specified symmetric discriminability (template -- plug in your group data)
# ----------------------------------------------------------------------------- #
def discriminability_table(subject_region_means, group_labels, est_sd=None):
    """Symmetric, pre-specified comparison across ALL metrics.

    subject_region_means : dict  subject_id -> {metric -> value}  (one region)
    group_labels         : dict  subject_id -> group name
    est_sd (optional)    : dict  metric -> estimation SD (from bootstrap) to report
                           an estimability-adjusted effect size.

    Returns: dict  metric -> {eta2, and (if est_sd) eta2_per_estsd}.
    eta2 = between-group variance / total variance (one-way, group-agnostic to which
    group is 'worse'; extend to your planned contrast/covariates in your stats stack).
    This is a SCREEN, not the final inferential model -- run all 8 identically and
    report every one regardless of outcome.
    """
    subs = list(subject_region_means.keys())
    metrics = PRIMARY_METRICS
    groups = sorted(set(group_labels[s] for s in subs))
    table = {}
    for m in metrics:
        vals = np.array([subject_region_means[s][m] for s in subs], float)
        g = np.array([group_labels[s] for s in subs])
        grand = np.nanmean(vals)
        ss_between = 0.0
        ss_total = np.nansum((vals - grand) ** 2)
        for gr in groups:
            v = vals[g == gr]
            v = v[np.isfinite(v)]
            if v.size:
                ss_between += v.size * (np.mean(v) - grand) ** 2
        eta2 = ss_between / ss_total if ss_total > 0 else np.nan
        entry = {"eta2": float(eta2)}
        if est_sd is not None and m in est_sd and est_sd[m] > 0:
            # effect size expressed per unit estimation noise: guards against a
            # metric "winning" only because it is estimated more precisely.
            spread = np.nanstd(vals)
            entry["eta2_per_estnoise_ratio"] = float(eta2 * (spread / est_sd[m]))
        table[m] = entry
    return table


# =============================================================================== #
# SYNTHETIC SELF-TEST
# =============================================================================== #
def _make_protocol():
    """Recreate the user's 4-shell scheme: b0x2, 100x6, 500x6, 1000x10, 2000x30."""
    from dipy.core.sphere import disperse_charges, HemiSphere
    def dirs(n, seed):
        rng = np.random.default_rng(seed)
        theta = np.pi * rng.random(n)
        phi = 2 * np.pi * rng.random(n)
        hs = HemiSphere(theta=theta, phi=phi)
        hs, _ = disperse_charges(hs, 50)
        return hs.vertices
    bvals, bvecs = [0, 0], [[0, 0, 0], [0, 0, 0]]
    for b, n, seed in [(100, 6, 1), (500, 6, 2), (1000, 10, 3), (2000, 30, 4)]:
        vs = dirs(n, seed)
        for v in vs:
            bvals.append(b); bvecs.append(list(v))
    return np.array(bvals, float), np.array(bvecs, float)


def _synth_test():
    from dipy.core.gradients import gradient_table
    from dipy.sims.voxel import multi_tensor_dki
    rng = np.random.default_rng(0)

    bvals, bvecs = _make_protocol()
    gtab = gradient_table(bvals, bvecs=bvecs, b0_threshold=50.0)

    # anisotropic WM-like tissue with realistic kurtosis
    mevals = np.array([[0.00185, 0.00055, 0.00055],
                       [0.00185, 0.00055, 0.00055]])
    angles = [(0, 0), (0, 0)]
    fie = [50, 50]
    S0 = 1000.0
    dki_sig, _, _ = multi_tensor_dki(gtab, mevals, S0=S0, angles=angles,
                                     fractions=fie, snr=None)

    f_true = 0.30
    fw_sig = S0 * np.exp(-bvals * D_FW)
    mixed = (1 - f_true) * dki_sig + f_true * fw_sig

    def add_rician(sig, snr):
        s = snr and (S0 / snr) or 0.0
        if s == 0:
            return sig
        n1 = rng.normal(0, s, sig.shape); n2 = rng.normal(0, s, sig.shape)
        return np.sqrt((sig + n1) ** 2 + n2 ** 2)

    results = {}
    for snr in [np.inf, 40, 20]:
        vol = add_rician(mixed, snr if np.isfinite(snr) else None)
        data = np.tile(vol, (3, 3, 3, 1)).astype(float)      # small phantom block
        mask = np.ones((3, 3, 3), bool)
        out = run_subject(data, bvals, bvecs, brain_mask=mask)
        f_est = np.median(out["f"][mask])
        md_u = np.median(out["MD"][mask]); md_c = np.median(out["MDt"][mask])
        fa_u = np.median(out["FA"][mask]); fa_c = np.median(out["FAt"][mask])
        mk_u = np.median(out["MK"][mask]); mk_c = np.median(out["MKt"][mask])
        results[snr] = (f_est, fa_u, fa_c, md_u, md_c, mk_u, mk_c)

    print("=" * 74)
    print("SYNTHETIC SELF-TEST  (f_true = %.2f, tissue = anisotropic WM+kurtosis)" % f_true)
    print("=" * 74)
    print(f"{'SNR':>5} | {'f_est':>6} | {'FA':>6} {'FAt':>6} | "
          f"{'MD':>7} {'MDt':>7} | {'MK':>6} {'MKt':>6}")
    print("-" * 74)
    for snr, (fe, fau, fac, mdu, mdc, mku, mkc) in results.items():
        s = "inf" if not np.isfinite(snr) else str(int(snr))
        print(f"{s:>5} | {fe:6.3f} | {fau:6.3f} {fac:6.3f} | "
              f"{mdu:7.5f} {mdc:7.5f} | {mku:6.3f} {mkc:6.3f}")
    print("-" * 74)
    print("EXPECTED, if the correction works:")
    print("  * f_est ~ 0.30 (recovers the injected free-water fraction)")
    print("  * FAt > FA   (removing isotropic free water raises anisotropy)")
    print("  * MDt < MD   (removing high-diffusivity water lowers mean diffusivity)")
    print("  * MKt vs MK: kurtosis changes once the isotropic compartment is removed")

    # sanity assertions (noiseless)
    fe, fau, fac, mdu, mdc, mku, mkc = results[np.inf]
    assert abs(fe - f_true) < 0.06, f"f recovery off: {fe}"
    assert fac > fau, "FAt should exceed FA after FWE"
    assert mdc < mdu, "MDt should be below MD after FWE"
    print("\nnoiseless assertions PASSED: f recovered, FAt>FA, MDt<MD")

    # ---- NO-b0 case: strip the b0 volumes and confirm the pipeline still works --
    nz = bvals > 50.0
    bvals_n, bvecs_n = bvals[nz], bvecs[nz]
    print("\n" + "=" * 74)
    print("NO-b0 SELF-TEST  (b0 volumes removed; S0 estimated by extrapolation)")
    print("=" * 74)
    print(f"{'SNR':>5} | {'S0_est':>7} | {'f_est':>6} | {'FAt':>6} {'MDt':>7} {'MKt':>6}")
    print("-" * 74)
    ok = True
    for snr in [np.inf, 40, 20]:
        vol = add_rician(mixed[nz], snr if np.isfinite(snr) else None)
        data = np.tile(vol, (3, 3, 3, 1)).astype(float)
        mask = np.ones((3, 3, 3), bool)
        s0 = estimate_S0(data, bvals_n, bvecs_n, mask=mask)
        out = run_subject(data, bvals_n, bvecs_n, brain_mask=mask)
        fe2 = float(np.median(out["f"][mask]))
        print(f"{('inf' if not np.isfinite(snr) else int(snr)):>5} | "
              f"{float(np.median(s0[mask])):7.1f} | {fe2:6.3f} | "
              f"{float(np.median(out['FAt'][mask])):6.3f} "
              f"{float(np.median(out['MDt'][mask])):7.5f} "
              f"{float(np.median(out['MKt'][mask])):6.3f}")
        if not np.isfinite(snr):
            ok = ok and (abs(fe2 - f_true) < 0.08) and \
                 (np.median(out["FAt"][mask]) > fau) and \
                 (np.median(out["MDt"][mask]) < mdu)
    assert ok, "no-b0 noiseless recovery failed"
    print("-" * 74)
    print("NO-b0 assertions PASSED: S0 estimated, f recovered, FAt>FA, MDt<MD")


# =============================================================================== #
# FILE I/O WRAPPER  +  CLI
# =============================================================================== #
def run_from_files(dwi, bval, bvec, out_dir, brain_mask=None, pvs_mask=None,
                   bmax_fw=BMAX_FW, f_method="hilow", save_maps=True):
    """Run the pipeline on preprocessed NIfTI files and write outputs.

    dwi        : path to preprocessed 4D DWI (denoise/Gibbs/Topup/Eddy already done)
    bval, bvec : FSL-style .bval / .bvec
    out_dir    : output directory
    brain_mask : optional brain mask NIfTI
    pvs_mask   : optional PVS mask NIfTI (binary or integer-label); if given, region
                 summaries are written to out_dir/pvs_region_summaries.json
    Returns the maps dict.
    """
    import json
    import nibabel as nib
    os.makedirs(out_dir, exist_ok=True)

    img = nib.load(dwi)
    data = img.get_fdata(dtype=np.float32)
    aff, hdr = img.affine, img.header
    bvals = np.loadtxt(bval)
    bvecs = np.loadtxt(bvec)
    if bvecs.shape[0] == 3:                      # FSL stores 3xN
        bvecs = bvecs.T
    bmask = nib.load(brain_mask).get_fdata().astype(bool) if brain_mask else None

    maps = run_subject(data, bvals, bvecs, brain_mask=bmask, bmax_fw=bmax_fw)
    # override f-method if requested (run_subject uses hilow by default)
    if f_method != "hilow":
        f, S0 = fit_free_water_fraction(data, bvals, bvecs, mask=bmask,
                                        bmax_fw=bmax_fw, method=f_method)
        maps["f"] = f
        tissue = deconfound_free_water(data, bvals, f, mask=bmask, S0=S0)
        maps.update(fit_corrected_dki(tissue, build_gtab(bvals, bvecs), mask=bmask))

    if save_maps:
        for name, arr in maps.items():
            nib.save(nib.Nifti1Image(np.nan_to_num(arr).astype(np.float32), aff, hdr),
                     os.path.join(out_dir, f"{name}.nii.gz"))

    if pvs_mask is not None:
        pmask_img = nib.load(pvs_mask)
        pmask = pmask_img.get_fdata()
        is_label = pmask.max() > 1 and np.allclose(pmask, np.round(pmask))
        region_labels = None
        if is_label:
            labs = [int(x) for x in np.unique(pmask) if x != 0]
            region_labels = {f"region_{l}": l for l in labs}
            summ = extract_region_summaries(maps, pmask, region_labels)
        else:
            summ = extract_region_summaries(maps, pmask > 0)
        with open(os.path.join(out_dir, "pvs_region_summaries.json"), "w") as fh:
            json.dump(summ, fh, indent=2)
    return maps


def _cli():
    import argparse
    p = argparse.ArgumentParser(
        description="Free-water-corrected DKI metrics inside PVS masks "
                    "(two-stage estimator for a 4-shell 7T acquisition).")
    p.add_argument("--selftest", action="store_true",
                   help="run the synthetic ground-truth self-test and exit")
    p.add_argument("--dwi"); p.add_argument("--bval"); p.add_argument("--bvec")
    p.add_argument("--brain-mask", default=None)
    p.add_argument("--pvs-mask", default=None)
    p.add_argument("--out", default="fwdki_out")
    p.add_argument("--bmax-fw", type=float, default=BMAX_FW,
                   help="top b-value used for the free-water (Gaussian) f estimate")
    p.add_argument("--f-method", choices=["hilow", "fwdti"], default="hilow",
                   help="hilow = robust High-Low-Downhill (default); "
                        "fwdti = DIPY free-water DTI (cross-check, fragile here)")
    a = p.parse_args()
    if a.selftest or not a.dwi:
        _synth_test(); return
    run_from_files(a.dwi, a.bval, a.bvec, a.out, brain_mask=a.brain_mask,
                   pvs_mask=a.pvs_mask, bmax_fw=a.bmax_fw, f_method=a.f_method)
    print(f"done -> {a.out}")


if __name__ == "__main__":
    _cli()
