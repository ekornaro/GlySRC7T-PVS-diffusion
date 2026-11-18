#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
step16_pvs_reference_tissue_extraction.py
=========================================
Three-tier spatial-specificity extraction with ANATOMICALLY MATCHED controls.

WHAT CHANGED FROM v8  (v9, 2026-10-04)
--------------------
 * FIBRE ANGLE IN EVERY TIER. v8 measured the angle between the PVS axis u and
   the principal fibre direction (V1 of the median tensor) at the PVS only
   (Theta_deg, cos2_theta); the matched-control sheet carried copies of that
   PVS value. v9 measures it, with the SAME u, in the reference tissue too:
       Case_A_Near_PVS : NearPVS_Theta_deg, NearPVS_cos2_theta
       Case_B_Fake_PVS : FakeControl_m<mm>_Theta_deg, FakeControl_m<mm>_cos2_theta
                         (mean over the control draws)
   Named tier-first (like NearPVS_dist_min), so the statistics and plotting
   scripts, which read '<metric>_<tier>' columns, never mistake them for
   diffusion metrics.
   The adjacent-tissue angle uses the same ROI as Dpar/Dperp_NearPVS (see the
   next point), so the angle and the projected diffusivities always describe
   the same voxels. A summary 'ORIENTATION BY TIER' is printed.
 * --near-orientation-roi full|matched. In v8 the orientation-referenced
   adjacent values (Dpar/Dperp/Dpar_minus_Dperp_NearPVS) came from the FULL
   3-6 mm shell (near_roi = near_full), while every scalar *_NearPVS metric used
   the SIZE-MATCHED draws. 'full' (default) keeps v8 exactly; 'matched' uses
   the size-matched draws for the orientation-referenced values and the angle
   as well (mean over draws), consistent with the scalar metrics.
 * No random numbers are drawn by the new code, so with the default every
   existing column is identical to v8 (same seed, same draws).
 * Default output: pvs_tier_metrics_v9.xlsx.

WHAT CHANGED FROM v7  (v8, post frame-audit)
--------------------
(I) HEADER-RELATION CORRECTION OF THE CLUSTER AXIS u (--rotate-u header,
    default). The tensor components are expressed in the NATIVE index frame
    (warped as scalars); u comes from T1w index coordinates. For the 27
    subjects whose two image headers differ by a REFLECTION (step20b flag ==
    header det=-1, verified identical lists), projecting the native-frame
    tensor onto an uncorrected u scrambles every orientation column to the
    uniform null -- exactly the artifact diagnosed and fixed in step16 v3
    (convergence test passed: 44.1 vs 42.6 deg, p=0.905). v8 ports that fix:
        R_rel = orth( inv(C_native) @ C_t1w )     [det preserved, NO repair]
        u     = R_rel @ u_pca
    u_rel_det is written to Subject_QC. '--rotate-u off' reproduces the
    (known-scrambled) v7 behaviour for audit only. Scalar metric columns are
    frame-invariant and unaffected by this choice.
(J) PATHS: DKI_SUBDIR is now 'DKI_in_T1w_space' (the 'revised_' tree no
    longer exists); STUDY_XLSX updated to the 2026-08-28 authority; the
    minimum mtime of the six tensor components is written to Subject_QC
    (tensor_mtime_min) so a stale-tensor run is visible in the spreadsheet.

WHAT CHANGED FROM v6
--------------------
(G) MEDIAL-AXIS THICKNESS (3D skeleton x 3D EDT), replacing 2*max(EDT).

    v6 defined Thickness_mm = 2 * max(EDT) over the cluster. Two problems:
      * max() is set by the single widest point -- one junction blob inflates
        the whole cluster's thickness;
      * the EDT was taken from the UNION PVS mask, so where two clusters touch
        diagonally (3x3x3 labelling keeps them separate, the EDT does not) the
        inscribed distance bleeds across cluster boundaries.

    v7 computes, per cluster, on the cluster's OWN local mask:
        1. 3D EDT with anisotropic voxel sampling  (radius of the inscribed
           sphere at every voxel);
        2. 3D skeletonization (skimage) -> the medial axis / centerline;
        3. EDT sampled ALONG the skeleton -> a thickness profile.
    Written columns:
        ThicknessSkel_med_mm   2 * median(EDT on skeleton)  <- primary
        ThicknessSkel_max_mm   2 * max(EDT on skeleton)
        Skel_len_mm            curved centerline length (voxel-count approx.)
        Elong_skel             Skel_len_mm / ThicknessSkel_med_mm
        Skel_voxels            skeleton size (0 -> fallback used)
    Fallback for clusters too small to skeletonize (skeleton empty): the
    single voxel of maximal local EDT is used as the skeleton.
    Thickness_mm / Thickness_Bin (the v6 definitions) are RETAINED unchanged
    for cross-version continuity.

    HONESTY NOTE: this does not un-quantize thickness. At ~0.8 mm T1w
    resolution most clusters are 1-2 voxels across, so the median-EDT
    thickness still sits on a discrete floor for the majority of clusters
    (the v4/v5 degeneracy). What the skeleton genuinely adds is robustness of
    the thickness estimate to junction blobs, a true CURVED length (Skel_len,
    vs the straight-line PCA extent), and a length-based elongation. Do not
    re-frame thickness as a dose-response axis on the strength of this change.

(H) EXPLICIT --min-cluster-size FLAG (sensitivity analysis support).

    MIN_CLUSTER_SIZE was a hidden module constant (3). Per the no-silent-
    defaults rule it is now a CLI flag; the value used is echoed in the run
    header and written to Subject_QC (min_cluster_size, n_below_min).
    Intended use: run once at the default 3 (primary analysis), once at e.g.
    20 (supplementary sensitivity check restricted to clusters large enough
    for stable per-tier medians -- note this selects the largest PVS, a
    different subpopulation, so it supplements rather than replaces the
    primary run):
        python step16_pvs_reference_tissue_extraction.py \
            --min-cluster-size 20 --out pvs_tier_metrics_v7_minsize20.xlsx

WHY v6 EXISTS
-------------
Diagnostics on the v4 spreadsheet (33,333 clusters, 57 subjects) showed that the
apparent PVS -> NearPVS -> FakeControl gradient was confounded on BOTH sides of
the effect-size ratio:

  NUMERATOR.  The control pool was "cerebral WM minus a dilation of the PVS
              mask". PVS burden therefore determined how much WM survived:

                  corr(PVS_voxels, validWM_frac_d10)         = -0.91
                  corr(validWM_frac_d10, MD_FakeControl_d10) = -0.48

              Burden is elevated in RHI, so the control tier acquired an upward
              RHI bias with no biological content.

  DENOMINATOR. Subject-level SD grew with the margin (d10 ~ 1.6x, d15 ~ 2.6x the
              PVS-tier SD), because larger margins leave a smaller and more
              idiosyncratic remnant of WM. Cohen's d fell across tiers largely
              because its denominator grew. Back-converting to raw units, the
              RHI-HC difference in AD was FLAT from PVS out to 12 mm.

  NearPVS.    near_roi excluded only the PARENT cluster's dilation, so with a
              median of 459 clusters per subject the shell landed next to
              neighbouring PVS: 88.1% of clusters had a shell voxel within 2 mm
              of some PVS, median distance 3.19 mm -- inside the ~4 mm
              interpolation point-spread of a 2 mm native grid.

WHAT CHANGED FROM v4
--------------------
(1) TIERS ARE DISTANCE BANDS IN MILLIMETRES, FROM THE UNION PVS MASK.

    A single EDT of the whole PVS mask defines every tier:

        PVS         inside the mask
        NearPVS     NEAR_BAND_MM[0] <= d < NEAR_BAND_MM[1]   (3-6 mm)
        FakeControl d >= margin                              (6 mm, 10 mm)

    A NearPVS voxel adjacent to a NEIGHBOURING cluster now fails its own
    distance test, which is the leak v4 could not catch. Tier ordering is
    guaranteed by construction rather than by dilation bookkeeping, and the
    margins are declared in mm so the PSF argument is auditable.

    Margins fell from 8/12 mm to 6/10 mm. The point-spread of linear
    interpolation from a 2 mm grid is ~4 mm; 6 mm clears it by 1.5x. The extra
    6 mm in v4 bought no additional PSF safety and cost up to 90% of the
    eligible WM in high-burden subjects.

(2) CONTROLS ARE MATCHED ON REGION x DEPTH, WITHIN SUBJECT, NO REGISTRATION.

    What varies with PVS burden is WHERE the surviving WM sits: burden is
    highest in superficial centrum semiovale, so exclusion strips superficial WM
    preferentially, and more so in high-burden subjects. Both coordinates needed
    to correct that are computable natively:

        depth  = d_cortex / (d_cortex + d_ventricle)     two EDTs, 0=cortex,
                                                          1=ventricle
        region = template_LBG / template_RBG / template_MB, else CSO by
                 hemisphere from SynthSeg

    Each cluster's control is drawn from the SAME region and the SAME depth bin
    in the SAME brain. This is propensity matching on subject-intrinsic
    geometry; no template warp is involved, so no group-dependent normalisation
    step is introduced into an analysis whose purpose is to remove a
    group-dependent artifact.

    v4 drew every control from SynthSeg labels 2/41 (cerebral WM) regardless of
    where the cluster lived, so basal-ganglia and midbrain clusters were
    compared against centrum semiovale. v6 keeps them in their own tissue.

(3) WMH AND ARTEFACT VOXELS ARE EXCLUDED FROM CONTROL POOLS. WMH carry elevated
    MD and were eligible controls in v4.

(4) N_FAKE_DRAWS DEFAULTS TO 5. Per-cluster draw noise is a minor variance term
    (~6% of between-subject variance) but it is free to remove.

(5) TRUE ORIENTATION ANGLE. v4 could only support an inversion of
    Dpar ~ RD + (AD-RD)cos^2(theta), which returned cos^2 = 0.354 against 0.333
    for a random direction, with 10.8% of clusters physically impossible. v6
    eigendecomposes each cluster's median tensor and reports theta, cos^2(theta)
    and tensor FA directly, so the question "is the PVS axis related to the
    fibre axis, or is the tensor frame simply wrong?" becomes answerable.

    NOTE: this does NOT settle tensor reorientation. Component maps warped as
    scalars are correct for MD/FA/AD/RD (rotation-invariant) but NOT for
    u^T D u. Check the rotation in your DKI->T1w transform independently; if it
    is non-trivial, reorient (ANTs ReorientTensorImage) and rerun.

(6) VENTRICULAR CSF REFERENCE TIER, per subject, from an eroded ventricle mask.
    Gives the dynamic range of every metric between parenchyma and free fluid.

(7) MATCH BALANCE IS REPORTED. Match_Balance sheet: achieved depth and region
    distribution for PVS vs controls, and how often matching had to be relaxed.
    Relaxation will concentrate in high-burden subjects; that must be visible
    rather than silent.

WHAT v6 DELIBERATELY DOES NOT DO
--------------------------------
Thickness binning is retained as a COLUMN but no longer framed as a
dose-response test. In v4 the bins were T1_thin=0, T2_mid=33,300, T3_thick=33.
There is no dynamic range: 93.2% of clusters are smaller than one native DKI
voxel (median 2.55 mm^3 vs 8 mm^3). Report that as a limitation; do not plot it
as a partial-volume result.

OUTPUTS
-------
  <out>.xlsx
      Case_A_Near_PVS   per-cluster PVS and NearPVS values
      Case_B_Fake_PVS   per-cluster matched-control values, per margin
      Subject_QC        burden, pool fractions, ventricular CSF reference
      Tier_Distances    measured distance-to-nearest-PVS per tier
      Match_Balance     depth/region balance and match relaxation counts
  <out>_provenance.csv

USAGE
-----
  python step16_pvs_reference_tissue_extraction.py --check
  python step16_pvs_reference_tissue_extraction.py
  python step16_pvs_reference_tissue_extraction.py --margins-mm 6 10 --fake-draws 5
"""

import os
import glob
import csv
import argparse
import datetime
import numpy as np
import nibabel as nib
from scipy import ndimage
import pandas as pd
from pipeline_config import cfg

try:
    from skimage.morphology import skeletonize as _skeletonize_3d
    HAVE_SKIMAGE = True
except ImportError:
    HAVE_SKIMAGE = False
    print("[!] scikit-image not available: ThicknessSkel_* columns will be NaN. "
          "pip install scikit-image to enable medial-axis thickness.")

# ==========================================
# CONFIGURATION
# ==========================================
ROOT_DIR = cfg("pvs_root")
SPRD_DIR = cfg("workbooks_dir")   # the folder this script was run from during the revision
OUTPUT_XLSX = os.path.join(SPRD_DIR, "pvs_tier_metrics_v9.xlsx")
NEAR_ORIENTATION_ROI = "full"     # v9: 'full' (= v8) or 'matched'; set in main

# Read the ORIGINAL warped maps. The reorientation attempt is not used: the
# native-space check showed the gradient sign convention differs between the two
# acquisition blocks (block A needs a y/z flip, block B does not; each reaches
# coherence ~0.73 under its own convention, p=0.47). That is a source-level
# difference, so no downstream rotation fixes it. It affects only
# orientation-dependent quantities (V1, theta, Dpar, Dperp); MD, FA, AD, RD and
# the kurtosis metrics are rotation-invariant and unaffected.
DKI_SUBDIR = "DKI_in_T1w_space"
SEG_SUBDIR = "Segmentation"
MAP_SUFFIX = "_toT1w_Warped.nii.gz"

METRICS_OF_INTEREST = ['MD', 'AD', 'RD', 'FA', 'MK_mf', 'AK_mf', 'RK_mf',
                       'lowb_ADC', 'lowb_MD', 'lowb_AD', 'lowb_RD', 'lowb_FA']

TENSOR_COMPONENTS = ['Dxx', 'Dyy', 'Dzz', 'Dxy', 'Dxz', 'Dyz']
MIN_LINEARITY = 2.0

STUDY_XLSX = cfg("study_xlsx")

# --- SynthSeg labels ---
WM_LABELS = [2, 41]                 # cerebral WM  (L, R)
# SynthSeg run with --parc emits Desikan-Killiany cortical labels (1001-1035
# left, 2001-2035 right) INSTEAD of the aggregate 3/42. Accepting only 3/42 left
# the cortex mask empty, which silently zeroed the depth coordinate and reduced
# control matching to region alone.
CORTEX_LABELS = [3, 42]
CORTEX_RANGES = [(1000, 1099), (2000, 2099)]

# DK parcel -> lobe. Propagated into white matter by nearest-cortex assignment,
# which the depth EDT already computes, so lobar matching is nearly free.
DK_LOBE = {}
for _lo, _ids in {
    "FRO": [3, 12, 14, 17, 18, 19, 20, 24, 27, 28, 32],
    "PAR": [8, 22, 25, 29, 31],
    "TMP": [1, 6, 7, 9, 15, 16, 30, 33, 34],
    "OCC": [5, 11, 13, 21],
    "CIN": [2, 10, 23, 26],
    "INS": [35],
}.items():
    for _i in _ids:
        DK_LOBE[1000 + _i] = _lo
        DK_LOBE[2000 + _i] = _lo
LOBES = ["FRO", "PAR", "TMP", "OCC", "CIN", "INS"]
VENTRICLE_LABELS = [4, 43, 5, 44, 14, 15]

# --- TIER GEOMETRY, IN MILLIMETRES (from the EDT of the union PVS mask) ---
# PSF of linear interpolation from a 2 mm native grid is ~4 mm.
NEAR_BAND_MM = (3.0, 6.0)           # immediately adjacent; declared PSF-adjacent
FAKE_MARGINS_MM = [6.0, 10.0]       # primary, sensitivity

# --- CONTROL MATCHING ---
# normalized depth: 0 at the cortical GM/WM interface, 1 at the ventricle
DEPTH_BIN_EDGES = [0.20, 0.40, 0.60, 0.80]
N_FAKE_DRAWS = 5
# Match the NearPVS ROI to the PVS cluster's own voxel count.
#
# WHY: the NearPVS shell spans ~683 voxels against ~5 for the cluster it
# surrounds -- roughly 137x more -- so its median is far more spatially smoothed.
# Smoothing pulls a local elevation toward the surrounding mean and suppresses
# the group difference, which is exactly what the raw differences showed:
#
#     AD   PVS 1.00x   NearPVS 0.55x   Fake_m6 0.89x   Fake_m10 0.94x
#     MD   PVS 1.00x   NearPVS 0.63x   Fake_m6 1.09x   Fake_m10 1.46x
#
# A V-shaped profile with its minimum at 3-6 mm is not what a PVS-localised
# effect looks like; it is what unequal smoothing looks like. Drawing a
# CONTIGUOUS ROI of the cluster's own size inside the shell -- the same way
# FakeControl is built -- makes all three tiers comparable in both size and
# compactness, so they differ only in distance from PVS.
#
# The unmatched full-shell value is still written as *_NearPVSfull, so the two
# can be compared directly: if the dip disappears under matching, that is
# positive evidence for the smoothing account rather than an assertion.
MATCH_NEAR_SIZE = True
N_NEAR_DRAWS = 5
MIN_POOL_FOR_STRATUM = 50           # voxels; below this, relax the stratum

THICKNESS_BINS_MM = [1.0, 2.0]      # retained as a column only (see docstring)
NATIVE_VOXEL_MM = 2.0
NATIVE_VOXEL_MM3 = NATIVE_VOXEL_MM ** 3

USE_LOBES = True          # set False (or --no-lobes) for CSO_L/CSO_R only
MIN_CLUSTER_SIZE = 3
VENTRICLE_ERODE = 2                 # voxels, to keep the CSF tier off the wall
RANDOM_SEED = 0

REGIONS = ["CSO_L", "CSO_R", "LBG", "RBG", "MB"]   # rebuilt per subject if lobes are used


# ==========================================
# FILE DISCOVERY
# ==========================================
def _unique(pattern, what, subj_id, errors, required=True):
    hits = sorted(glob.glob(pattern, recursive=True))
    if len(hits) == 1:
        return hits[0]
    if not required:
        return None
    errors.append(f"{subj_id}: " + ("MISSING " + what if not hits else
                  f"AMBIGUOUS {what} -> {[os.path.basename(h) for h in hits][:4]}"))
    return None


def find_files_for_subject(subj_dir, subj_id, errors):
    found = {"dki_paths": {}, "seg": {}}
    map_dir = os.path.join(subj_dir, DKI_SUBDIR)
    if not os.path.isdir(map_dir):
        errors.append(f"{subj_id}: no '{DKI_SUBDIR}' folder")
        return None
    for metric in METRICS_OF_INTEREST:
        p = os.path.join(map_dir, f"{metric}{MAP_SUFFIX}")
        if os.path.exists(p):
            found["dki_paths"][metric] = p
        else:
            errors.append(f"{subj_id}: MISSING map {metric}{MAP_SUFFIX}")

    found["pvs_path"] = _unique(
        os.path.join(subj_dir, "**", "*_ses-01_PVS_postprocessed.nii*"),
        "PVS mask", subj_id, errors)
    found["wm_path"] = _unique(
        os.path.join(subj_dir, "**", "*_ses-01_synthseg.nii.gz"),
        "SynthSeg mask", subj_id, errors)

    # optional segmentation inputs -- absence degrades matching, never fatal
    seg_glob = os.path.join(subj_dir, SEG_SUBDIR, "*_ses-01_%s.nii.gz")
    for key, tag in [("ventricles", "ventricles_final"), ("wmh", "WMH_samseg"),
                     ("artefact", "artefact"), ("LBG", "template_LBG"),
                     ("RBG", "template_RBG"), ("MB", "template_MB"),
                     ("gm", "GM_fast")]:
        found["seg"][key] = _unique(seg_glob % tag, tag, subj_id, errors,
                                    required=False)

    if found["pvs_path"] and found["wm_path"] and found["dki_paths"]:
        return found
    return None


def included_subjects(xlsx=STUDY_XLSX):
    x = pd.ExcelFile(xlsx)
    out = set()
    for sheet in ("HC", "PPCS", "Boxers"):
        t = x.parse(sheet)
        u = t["Will we used it?"].astype(str).str.strip().str.lower()
        out |= {f"sub-{int(n)}" for n in t.loc[u == "yes", "Subj"]}
    return out


def _load_bool(path, shape):
    """Load an optional mask, returning None if absent or grid-mismatched."""
    if not path:
        return None
    arr = np.asarray(nib.load(path).dataobj)
    if arr.shape != shape:
        return None
    return arr > 0


# ==========================================
# GEOMETRY
# ==========================================
def grow_roi(valid_mask, seed_coords, target_size, rng, max_attempts=40):
    """Grow a contiguous ROI of target_size voxels entirely inside valid_mask.

    Seeds are drawn from seed_coords (the matched stratum); growth is permitted
    anywhere in valid_mask (the region), so a tight depth stratum cannot make
    small ROIs unbuildable.
    """
    if seed_coords is None or len(seed_coords) == 0:
        return None
    shape = valid_mask.shape
    for _ in range(max_attempts):
        seed = tuple(int(v) for v in seed_coords[rng.integers(0, len(seed_coords))])
        if not valid_mask[seed]:
            continue
        roi, seen, queue = [], {seed}, [seed]
        while queue and len(roi) < target_size:
            x, y, z = queue.pop(0)
            roi.append((x, y, z))
            for nb in ((x + 1, y, z), (x - 1, y, z), (x, y + 1, z),
                       (x, y - 1, z), (x, y, z + 1), (x, y, z - 1)):
                if nb in seen:
                    continue
                nx, ny, nz = nb
                if 0 <= nx < shape[0] and 0 <= ny < shape[1] \
                        and 0 <= nz < shape[2] and valid_mask[nx, ny, nz]:
                    seen.add(nb)
                    queue.append(nb)
        if len(roi) == target_size:
            return tuple(np.array(c) for c in zip(*roi))
    return None


def axis_precision_deg(coords_mm, u, rng, max_leave=30):
    """Jackknife angular standard error of the cluster axis, in degrees.

    Leave-one-voxel-out, refit the coordinate PCA, measure the angle to the full
    axis. This is the principled variable on which to select well-conditioned
    clusters, because it is computed WITHOUT group labels: filtering on a
    group-blind statistic is exchangeable across groups under the null and so
    does not inflate type I error, whereas filtering on anything touching the
    contrast does.

    Small values mean the axis is determined by the cluster shape; large values
    mean it is determined by which voxels happened to survive segmentation.
    """
    n = len(coords_mm)
    if n < 4:
        return np.nan
    idxs = (range(n) if n <= max_leave
            else rng.choice(n, max_leave, replace=False))
    angs = []
    for i in idxs:
        c = np.delete(coords_mm, i, axis=0)
        c = c - c.mean(axis=0)
        try:
            w, V = np.linalg.eigh((c.T @ c) / max(len(c) - 1, 1))
        except np.linalg.LinAlgError:
            continue
        angs.append(np.degrees(np.arccos(
            np.clip(abs(float(np.dot(u, V[:, -1]))), 0.0, 1.0))))
    return float(np.mean(angs)) if angs else np.nan



# --- v8: native/T1w header relation (ported from step16 v3) ---------------
DKI_ROOT_DEFAULT = cfg("data_root")
NATIVE_REF_SUBPATH = "pipeline_output/DATASINK/DKISpace/tensor/Dxx.nii.gz"


def _dircos(affine):
    M = np.asarray(affine, float)[:3, :3]
    return M / np.linalg.norm(M, axis=0, keepdims=True)


def _orthogonal_part(M):
    """Orthogonal factor via SVD, determinant PRESERVED (no reflection repair)."""
    W, _, Vt = np.linalg.svd(M)
    return W @ Vt


def build_u_header_relation(sub, t1w_affine, dki_root=DKI_ROOT_DEFAULT):
    """R_rel mapping T1w-index directions into the native tensor frame.
    Reflections kept: det=-1 subjects are the step20b-reflected ones."""
    ref = os.path.join(dki_root, sub, NATIVE_REF_SUBPATH)
    if not os.path.exists(ref):
        raise RuntimeError(f"native reference missing for header relation: {ref}")
    R = _orthogonal_part(np.linalg.inv(_dircos(nib.load(ref).affine))
                         @ _dircos(t1w_affine))
    return R.astype(np.float64), float(np.sign(np.linalg.det(R)))


def skeleton_thickness(single_loc, zooms, vox_mm):
    """
    Medial-axis thickness of ONE cluster, on its own local mask (v7, item G).

        EDT_loc  = 3D EDT of the cluster mask, anisotropic sampling
                   (radius of the inscribed sphere at each voxel)
        skeleton = 3D skeletonization of the mask (centerline)
        profile  = EDT_loc sampled along the skeleton

    Returns dict:
        med_mm   2 * median(profile)   -- primary, robust to junction blobs
        max_mm   2 * max(profile)
        len_mm   centerline length: total edge length of the minimum spanning
                 tree over skeleton voxels in mm (handles curves and branches;
                 a single-voxel skeleton is assigned one voxel length)
        n_skel   number of skeleton voxels (0 => fallback voxel was used)

    Fallback: if the skeleton is empty (cluster too small/thin to
    skeletonize), the voxel of maximal EDT stands in as a 1-point skeleton.
    """
    out = {"med_mm": np.nan, "max_mm": np.nan, "len_mm": np.nan, "n_skel": 0}
    if not HAVE_SKIMAGE:
        return out
    edt_loc = ndimage.distance_transform_edt(single_loc, sampling=zooms)
    skel = _skeletonize_3d(single_loc)
    skel = skel.astype(bool) & single_loc
    n_skel = int(skel.sum())
    if n_skel == 0:                              # fallback: 1-point skeleton
        idx = np.unravel_index(np.argmax(edt_loc), edt_loc.shape)
        skel = np.zeros_like(single_loc)
        skel[idx] = True
    prof = edt_loc[skel]
    out["med_mm"] = float(2.0 * np.median(prof))
    out["max_mm"] = float(2.0 * np.max(prof))
    out["n_skel"] = n_skel
    pts = np.argwhere(skel) * np.asarray(zooms)
    if len(pts) == 1:
        out["len_mm"] = float(vox_mm)
    elif len(pts) > 2000:                        # degenerate 'skeleton'; approx
        out["len_mm"] = float(len(pts) * vox_mm)
    else:
        from scipy.spatial.distance import pdist, squareform
        from scipy.sparse.csgraph import minimum_spanning_tree
        mst = minimum_spanning_tree(squareform(pdist(pts)))
        out["len_mm"] = float(mst.sum())
    return out


def thickness_bin(mm, edges=THICKNESS_BINS_MM):
    if mm < edges[0]:
        return "T1_thin"
    if mm < edges[1]:
        return "T2_mid"
    return "T3_thick"


def size_bin(mm3, native=NATIVE_VOXEL_MM3):
    r = mm3 / native
    if r < 0.5:
        return "A_subvoxel"
    if r < 1.0:
        return "B_nearvoxel"
    return "C_supravoxel"


def depth_bin_index(d, edges=DEPTH_BIN_EDGES):
    return int(np.searchsorted(edges, d, side="right"))


def _dist_stats(dist_map, sel, prefix):
    v = dist_map[sel]
    v = v[np.isfinite(v)]
    if v.size == 0:
        return {f"{prefix}_dist_min": np.nan, f"{prefix}_dist_med": np.nan,
                f"{prefix}_dist_mean": np.nan}
    return {f"{prefix}_dist_min": float(v.min()),
            f"{prefix}_dist_med": float(np.median(v)),
            f"{prefix}_dist_mean": float(v.mean())}


def project_tensor(comp, sel, u):
    """Median Dpar and Dperp over `sel`, referenced to unit axis u."""
    xx, yy, zz = comp["Dxx"][sel], comp["Dyy"][sel], comp["Dzz"][sel]
    xy, xz, yz = comp["Dxy"][sel], comp["Dxz"][sel], comp["Dyz"][sel]
    if xx.size == 0:
        return np.nan, np.nan
    d_par = (u[0] ** 2 * xx + u[1] ** 2 * yy + u[2] ** 2 * zz
             + 2 * u[0] * u[1] * xy + 2 * u[0] * u[2] * xz
             + 2 * u[1] * u[2] * yz)
    d_perp = ((xx + yy + zz) - d_par) / 2.0
    return float(np.nanmedian(d_par)), float(np.nanmedian(d_perp))


def tensor_orientation(comp, sel, u):
    """Eigendecompose the median tensor over `sel`; angle to axis u.

    Returns (theta_deg, cos2_theta, tensor_FA, lam1, lam2, lam3).
    This is the direct measurement that v4 could only approximate by inverting
    Dpar ~ RD + (AD-RD)cos^2(theta).
    """
    try:
        m = {c: float(np.nanmedian(comp[c][sel])) for c in TENSOR_COMPONENTS}
    except (ValueError, IndexError):
        return (np.nan,) * 6
    if any(not np.isfinite(v) for v in m.values()):
        return (np.nan,) * 6
    D = np.array([[m["Dxx"], m["Dxy"], m["Dxz"]],
                  [m["Dxy"], m["Dyy"], m["Dyz"]],
                  [m["Dxz"], m["Dyz"], m["Dzz"]]], dtype=float)
    try:
        w, V = np.linalg.eigh(D)
    except np.linalg.LinAlgError:
        return (np.nan,) * 6
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]
    v1 = V[:, 0]
    c = float(np.clip(abs(float(np.dot(u, v1))), 0.0, 1.0))
    md = float(w.mean())
    num = float(np.sqrt(((w - md) ** 2).sum()))
    den = float(np.sqrt((w ** 2).sum()))
    fa = float(np.sqrt(1.5) * num / den) if den > 1e-12 else np.nan
    return (float(np.degrees(np.arccos(c))), c ** 2, fa,
            float(w[0]), float(w[1]), float(w[2]))


def orientation_over_rois(comp, rois, u):
    """v9: angle between u and V1 in each ROI of `rois`; returns the mean
    theta (deg) and mean cos^2 over the ROIs with a finite result."""
    th, c2 = [], []
    for roi in rois:
        t, c, *_ = tensor_orientation(comp, roi, u)
        if np.isfinite(t):
            th.append(t); c2.append(c)
    if not th:
        return np.nan, np.nan
    return float(np.mean(th)), float(np.mean(c2))


def projection_over_rois(comp, rois, u):
    """v9: mean over ROIs of the median Dpar / Dperp (as for the controls)."""
    vals = [project_tensor(comp, r, u) for r in rois]
    vals = [v for v in vals if np.isfinite(v[0])]
    if not vals:
        return np.nan, np.nan
    return float(np.mean([v[0] for v in vals])), float(np.mean([v[1] for v in vals]))


# ==========================================
# PER-SUBJECT PROCESSING
# ==========================================
def build_anatomy(paths, pvs_mask, zooms, errors, subj_id):
    """Native-space region labels and normalized depth. No registration."""
    shape = pvs_mask.shape
    seg = np.asarray(nib.load(paths["wm_path"]).dataobj)
    if seg.shape != shape:
        errors.append(f"{subj_id}: SynthSeg grid {seg.shape} != PVS {shape}")
        return None

    wm_l = seg == WM_LABELS[0]
    wm_r = seg == WM_LABELS[1]
    cortex = np.isin(seg, CORTEX_LABELS)
    for lo, hi in CORTEX_RANGES:
        cortex |= (seg >= lo) & (seg <= hi)
    if not cortex.any():
        gm = _load_bool(paths["seg"].get("gm"), shape)
        if gm is not None and gm.any():
            cortex = gm
            print(f"  [!] {subj_id}: no cortex labels in SynthSeg; using GM_fast")
    if not cortex.any():
        errors.append(f"{subj_id}: NO CORTEX MASK -- depth cannot be computed. "
                      f"Labels present: {sorted(np.unique(seg))[:20]}")
        return None          # fatal: silently zeroing depth disables the matching

    vent = _load_bool(paths["seg"].get("ventricles"), shape)
    if vent is None:
        vent = np.isin(seg, VENTRICLE_LABELS)
        print(f"  [i] {subj_id}: no ventricles file; derived from SynthSeg "
              f"labels {VENTRICLE_LABELS} ({int(vent.sum())} vox)")

    # normalized depth: 0 at the cortical interface, 1 at the ventricle
    d_cx, near_idx = ndimage.distance_transform_edt(
        ~cortex, sampling=zooms, return_indices=True)
    # each voxel inherits the label of its nearest cortical parcel
    parcel = seg[near_idx[0], near_idx[1], near_idx[2]]
    if vent.any():
        d_vt = ndimage.distance_transform_edt(~vent, sampling=zooms)
    else:
        d_vt = np.ones(shape, np.float32)
    depth = (d_cx / (d_cx + d_vt + 1e-6)).astype(np.float32)

    lbg = _load_bool(paths["seg"].get("LBG"), shape)
    rbg = _load_bool(paths["seg"].get("RBG"), shape)
    mb = _load_bool(paths["seg"].get("MB"), shape)
    deep = np.zeros(shape, bool)
    for m in (lbg, rbg, mb):
        if m is not None:
            deep |= m

    pools = {"LBG": lbg if lbg is not None else np.zeros(shape, bool),
             "RBG": rbg if rbg is not None else np.zeros(shape, bool),
             "MB": mb if mb is not None else np.zeros(shape, bool)}
    lobe_of = np.zeros(shape, np.int8)
    if USE_LOBES:
        for li, lo in enumerate(LOBES, start=1):
            ids = [k for k, v in DK_LOBE.items() if v == lo]
            if ids:
                lobe_of[np.isin(parcel, ids)] = li
        for li, lo in enumerate(LOBES, start=1):
            m = lobe_of == li
            pools[f"{lo}_L"] = wm_l & ~deep & m
            pools[f"{lo}_R"] = wm_r & ~deep & m
        # anything whose nearest parcel was not a DK cortical label
        pools["OTH_L"] = wm_l & ~deep & (lobe_of == 0)
        pools["OTH_R"] = wm_r & ~deep & (lobe_of == 0)
    else:
        pools["CSO_L"] = wm_l & ~deep
        pools["CSO_R"] = wm_r & ~deep
    pools = {k: v for k, v in pools.items() if v.any()}

    bad = np.zeros(shape, bool)
    for key in ("wmh", "artefact"):
        m = _load_bool(paths["seg"].get(key), shape)
        if m is not None:
            bad |= m

    return {"seg": seg, "pools": pools, "depth": depth, "bad": bad,
            "vent": vent, "any_wm": wm_l | wm_r,
            "hemi": {k: ("L" if k.endswith("_L") else
                         "R" if k.endswith("_R") else k) for k in pools}}


def process_subject(subj_id, paths, rng, margins_mm, n_draws, errors,
                    n_near_draws=N_NEAR_DRAWS,
                    min_cluster_size=MIN_CLUSTER_SIZE,
                    rotate_u="header", dki_root=DKI_ROOT_DEFAULT):
    pvs_img = nib.load(paths["pvs_path"])
    pvs_mask = np.asarray(pvs_img.dataobj) > 0
    shape = pvs_mask.shape
    zooms = tuple(float(z) for z in pvs_img.header.get_zooms()[:3])
    vox_mm3 = float(np.prod(zooms))
    vox_mm = float(np.mean(zooms))

    anat = build_anatomy(paths, pvs_mask, zooms, errors, subj_id)
    if anat is None:
        return [], [], None, []

    dki_maps = {}
    for metric, p in paths["dki_paths"].items():
        d = np.asarray(nib.load(p).dataobj, np.float32)
        if d.shape != shape:
            errors.append(f"{subj_id}: {metric} grid {d.shape} != PVS {shape}")
            continue
        dki_maps[metric] = d
    if not dki_maps:
        return [], [], None, []

    comp = {}
    map_dir = os.path.dirname(list(paths["dki_paths"].values())[0])
    for c in TENSOR_COMPONENTS:
        p = os.path.join(map_dir, f"{c}{MAP_SUFFIX}")
        if os.path.exists(p):
            arr = np.asarray(nib.load(p).dataobj, np.float32)
            if arr.shape == shape:
                comp[c] = arr
    have_tensor = len(comp) == len(TENSOR_COMPONENTS)
    if not have_tensor:
        print(f"  [!] {subj_id}: tensor components incomplete "
              f"({len(comp)}/{len(TENSOR_COMPONENTS)}); Dpar/Dperp skipped")

    # ---- v8: frame relation for u (reflections kept) ----
    R_rel, u_rel_det = np.eye(3), np.nan
    if have_tensor and rotate_u == "header":
        try:
            R_rel, u_rel_det = build_u_header_relation(subj_id, pvs_img.affine,
                                                       dki_root)
        except RuntimeError as e:
            errors.append(f"{subj_id}: {e} -- orientation columns skipped")
            print(f"  [!] {subj_id}: {e}; Dpar/Dperp/theta skipped")
            have_tensor = False
    tensor_mtime = (min(os.path.getmtime(os.path.join(map_dir, f"{c}{MAP_SUFFIX}"))
                        for c in TENSOR_COMPONENTS) if have_tensor else np.nan)

    # ---- ONE EDT defines every tier ----
    edt_out = ndimage.distance_transform_edt(~pvs_mask, sampling=zooms)
    near_band = (edt_out >= NEAR_BAND_MM[0]) & (edt_out < NEAR_BAND_MM[1])

    depth = anat["depth"]
    dbin_map = np.searchsorted(DEPTH_BIN_EDGES, depth, side="right").astype(np.int8)
    n_dbins = len(DEPTH_BIN_EDGES) + 1

    # ---- control pools, per margin, per region, per depth bin ----
    regions = list(anat["pools"].keys())
    pool_mask, seed_coords = {}, {}
    for mm in margins_mm:
        far = edt_out >= mm
        for rname, rmask in anat["pools"].items():
            vm = rmask & far & ~anat["bad"] & ~pvs_mask
            pool_mask[(mm, rname)] = vm
            for b in range(n_dbins):
                sel = vm & (dbin_map == b)
                n = int(sel.sum())
                seed_coords[(mm, rname, b)] = (np.argwhere(sel)
                                               if n >= MIN_POOL_FOR_STRATUM
                                               else None)
            seed_coords[(mm, rname, None)] = np.argwhere(vm) if vm.any() else None

    n_wm = int(anat["any_wm"].sum())
    qc = {"Subject": subj_id, "voxel_mm": round(vox_mm, 3),
          "voxel_mm3": round(vox_mm3, 4),
          "PVS_voxels": int(pvs_mask.sum()), "WM_voxels": n_wm,
          "NearBand_voxels": int((near_band & anat["any_wm"]).sum())}
    for mm in margins_mm:
        tot = sum(int(pool_mask[(mm, r)].sum()) for r in regions)
        qc[f"pool_vox_m{mm:g}"] = tot
        qc[f"pool_frac_m{mm:g}"] = round(tot / n_wm, 4) if n_wm else np.nan
        qc[f"fakeFail_m{mm:g}"] = 0
        qc[f"relax_region_m{mm:g}"] = 0
        qc[f"relax_any_m{mm:g}"] = 0

    # ---- ventricular CSF reference tier ----
    vent = anat["vent"]
    if vent.any():
        st1 = ndimage.generate_binary_structure(3, 1)
        core = ndimage.binary_erosion(vent, structure=st1,
                                      iterations=VENTRICLE_ERODE)
        if not core.any():
            core = vent
        qc["vCSF_voxels"] = int(core.sum())
        for metric, m in dki_maps.items():
            qc[f"vCSF_{metric}"] = float(np.nanmedian(m[core]))
    else:
        qc["vCSF_voxels"] = 0

    st = ndimage.generate_binary_structure(3, 1)
    labeled, n_feat = ndimage.label(pvs_mask, structure=np.ones((3, 3, 3)))
    objs = ndimage.find_objects(labeled)
    edt_in = ndimage.distance_transform_edt(pvs_mask, sampling=zooms)

    pad = int(np.ceil(NEAR_BAND_MM[1] / vox_mm)) + 2
    near_rows, fake_rows, bal_rows = [], [], []
    n_below_min = 0

    for lab in range(1, n_feat + 1):
        sl = objs[lab - 1]
        if sl is None:
            continue
        # local bounding box, padded -- avoids whole-volume dilation per cluster
        loc = tuple(slice(max(0, s.start - pad), min(shape[i], s.stop + pad))
                    for i, s in enumerate(sl))
        sub_lab = labeled[loc]
        single_loc = sub_lab == lab
        size = int(single_loc.sum())
        if size < min_cluster_size:
            n_below_min += 1
            continue
        mm3 = size * vox_mm3

        single = np.zeros(shape, bool)
        single[loc] = single_loc

        thick_mm = float(2.0 * edt_in[single].max())   # v6 definition, retained
        skel = skeleton_thickness(single_loc, zooms, vox_mm)   # v7, item (G)
        # NOTE: the old Elongation = mm3 / ((pi/6) thickness^3) was NOT a shape
        # measure. Thickness is pinned at 2*voxel for 99.1% of clusters, so the
        # denominator is constant and Elongation correlated with Size_Voxels at
        # r = 0.9997 while correlating only 0.20 with Linearity. Shape now comes
        # from the coordinate PCA eigenvalues, below.

        # cluster principal axis -- the SAME u is used for every tier
        cc = np.argwhere(single_loc) * np.array(zooms)
        cc = cc - cc.mean(axis=0)
        wpca, Vpca = np.linalg.eigh((cc.T @ cc) / max(len(cc) - 1, 1))
        u_axis = Vpca[:, -1]
        axis_se_u = u_axis                       # precision: frame-internal
        if have_tensor and rotate_u == "header":
            u_axis = R_rel @ u_axis              # v8: into the native frame
            u_axis = u_axis / np.linalg.norm(u_axis)
        linearity = float(wpca[-1] / wpca[-2]) if wpca[-2] > 1e-12 else np.inf
        # genuine shape: ratio of longest to shortest principal extent
        elong = (float(np.sqrt(wpca[-1] / wpca[0]))
                 if wpca[0] > 1e-12 else np.inf)
        planarity = (float(wpca[-2] / wpca[0]) if wpca[0] > 1e-12 else np.inf)
        axis_se = axis_precision_deg(cc, axis_se_u, rng)

        # ---- region and depth of THIS cluster ----
        reg = None
        best = 0
        for rname in regions:
            ov = int((anat["pools"][rname] & single).sum())
            if ov > best:
                best, reg = ov, rname
        if reg is None:
            reg = regions[0]
        cl_depth = float(np.nanmedian(depth[single]))
        cl_dbin = depth_bin_index(cl_depth)

        # ---- NearPVS: local neighbourhood INTERSECTED with the global band ----
        dil_loc = ndimage.binary_dilation(single_loc, structure=st,
                                          iterations=pad)
        near_loc = (dil_loc & near_band[loc] & anat["pools"][reg][loc]
                    & ~anat["bad"][loc] & ~pvs_mask[loc])
        near_full = np.zeros(shape, bool)
        near_full[loc] = near_loc

        # size- and shape-matched sub-ROIs inside the shell
        near_draws, near_matched = [], False
        if MATCH_NEAR_SIZE and near_full.sum() >= size:
            nseeds = np.argwhere(near_full)
            for _ in range(n_near_draws):
                g = grow_roi(near_full, nseeds, size, rng)
                if g is not None:
                    near_draws.append(g)
            near_matched = len(near_draws) > 0
        if not near_matched:
            near_draws = [tuple(np.argwhere(near_full).T)] if near_full.any() else []
        near_roi = near_full        # kept for distance stats and tensor projection

        common = {"Subject": subj_id, "PVS_ID": lab, "Size_Voxels": size,
                  "Size_mm3": round(mm3, 3), "Size_Bin": size_bin(mm3),
                  "Thickness_mm": round(thick_mm, 3),
                  "Thickness_Bin": thickness_bin(thick_mm),
                  "ThicknessSkel_med_mm": (round(skel["med_mm"], 3)
                                           if np.isfinite(skel["med_mm"]) else np.nan),
                  "ThicknessSkel_max_mm": (round(skel["max_mm"], 3)
                                           if np.isfinite(skel["max_mm"]) else np.nan),
                  "Skel_len_mm": (round(skel["len_mm"], 3)
                                  if np.isfinite(skel["len_mm"]) else np.nan),
                  "Elong_skel": (round(skel["len_mm"] / skel["med_mm"], 3)
                                 if np.isfinite(skel["len_mm"])
                                 and np.isfinite(skel["med_mm"])
                                 and skel["med_mm"] > 0 else np.nan),
                  "Skel_voxels": skel["n_skel"],
                  "Elongation": round(elong, 3) if np.isfinite(elong) else np.nan,
                  "Planarity": round(planarity, 3) if np.isfinite(planarity) else np.nan,
                  "Linearity": round(linearity, 3) if np.isfinite(linearity) else np.nan,
                  "Axis_SE_deg": round(axis_se, 2) if np.isfinite(axis_se) else np.nan,
                  "Region": reg, "Depth": round(cl_depth, 4), "Depth_Bin": cl_dbin,
                  "axis_x": round(float(u_axis[0]), 4),
                  "axis_y": round(float(u_axis[1]), 4),
                  "axis_z": round(float(u_axis[2]), 4)}
        r_near = dict(common,
                      NearPVS_Voxels=(size if near_matched
                                      else int(near_full.sum())),
                      NearPVS_Full_Voxels=int(near_full.sum()),
                      NearPVS_SizeMatched=bool(near_matched),
                      NearPVS_Draws=len(near_draws))
        r_fake = dict(common)
        r_near.update(_dist_stats(edt_out, near_roi, "NearPVS"))

        # ---- matched controls, with a recorded relaxation ladder ----
        fake_coords, bal = {}, {"Subject": subj_id, "PVS_ID": lab,
                                "Region": reg, "PVS_Depth": round(cl_depth, 4)}
        for mm in margins_mm:
            vm = pool_mask[(mm, reg)]
            level, sc = 0, seed_coords.get((mm, reg, cl_dbin))
            if sc is None:                       # 1: relax depth, keep region
                level, sc = 1, seed_coords.get((mm, reg, None))
            if sc is None:                       # 2: same hemisphere, any lobe
                hemi = anat["hemi"].get(reg)
                for alt in regions:
                    if anat["hemi"].get(alt) == hemi and \
                            seed_coords.get((mm, alt, None)) is not None:
                        level = 2
                        vm, sc = pool_mask[(mm, alt)], seed_coords[(mm, alt, None)]
                        break
            if sc is None:                       # 3: anywhere
                for alt in regions:
                    if seed_coords.get((mm, alt, None)) is not None:
                        level = 3
                        vm, sc = pool_mask[(mm, alt)], seed_coords[(mm, alt, None)]
                        break
            draws = [grow_roi(vm, sc, size, rng) for _ in range(n_draws)]
            draws = [d for d in draws if d is not None]
            fake_coords[mm] = draws
            tag = f"m{mm:g}"
            r_fake[f"FakeControl_Voxels_{tag}"] = size if draws else 0
            r_fake[f"Match_Level_{tag}"] = level if draws else -1
            if not draws:
                qc[f"fakeFail_m{mm:g}"] += 1
                bal[f"ctrl_depth_{tag}"] = np.nan
            else:
                if level == 1:
                    qc[f"relax_region_m{mm:g}"] += 1
                elif level >= 2:
                    qc[f"relax_any_m{mm:g}"] += 1
                r_fake.update(_dist_stats(edt_out, draws[0], f"FakeControl_{tag}"))
                bal[f"ctrl_depth_{tag}"] = float(np.mean(
                    [np.nanmedian(depth[c]) for c in draws]))
            bal[f"match_level_{tag}"] = r_fake[f"Match_Level_{tag}"]
        bal_rows.append(bal)

        # ---- metric extraction ----
        for metric, m in dki_maps.items():
            o = m[single]
            val_o = float(np.nanmedian(o)) if o.size else np.nan
            r_near[f"{metric}_Original"] = val_o
            r_fake[f"{metric}_Original"] = val_o
            nf = m[near_full]
            r_near[f"{metric}_NearPVSfull"] = (float(np.nanmedian(nf))
                                               if nf.size else np.nan)
            vals = [np.nanmedian(m[c]) for c in near_draws if len(c[0])]
            r_near[f"{metric}_NearPVS"] = float(np.mean(vals)) if vals else np.nan
            for mm in margins_mm:
                vals = [np.nanmedian(m[c]) for c in fake_coords[mm]]
                r_fake[f"{metric}_FakeControl_m{mm:g}"] = (float(np.mean(vals))
                                                           if vals else np.nan)

        # ---- orientation-referenced diffusivity + true theta ----
        if have_tensor:
            p_o, p_e = project_tensor(comp, single, u_axis)
            th, c2, tfa, l1, l2, l3 = tensor_orientation(comp, single, u_axis)
            for r in (r_near, r_fake):
                r["Dpar_Original"] = p_o
                r["Dperp_Original"] = p_e
                r["Dpar_minus_Dperp_Original"] = p_o - p_e
                r["Theta_deg"] = round(th, 2) if np.isfinite(th) else np.nan
                r["cos2_theta"] = round(c2, 4) if np.isfinite(c2) else np.nan
                r["Tensor_FA"] = round(tfa, 4) if np.isfinite(tfa) else np.nan
                r["lam1"] = l1
                r["lam2"] = l2
                r["lam3"] = l3

            if NEAR_ORIENTATION_ROI == "matched" and near_draws:
                n_o, n_e = projection_over_rois(comp, near_draws, u_axis)
                nth, nc2 = orientation_over_rois(comp, near_draws, u_axis)
            else:                                   # v8 behaviour: full shell
                n_o, n_e = project_tensor(comp, near_roi, u_axis)
                nth, nc2, *_ = tensor_orientation(comp, near_roi, u_axis)
            r_near["Dpar_NearPVS"] = n_o
            r_near["Dperp_NearPVS"] = n_e
            r_near["Dpar_minus_Dperp_NearPVS"] = n_o - n_e
            r_near["NearPVS_Theta_deg"] = round(nth, 2) if np.isfinite(nth) else np.nan
            r_near["NearPVS_cos2_theta"] = round(nc2, 4) if np.isfinite(nc2) else np.nan

            for mm in margins_mm:
                tag = f"m{mm:g}"
                vals = [project_tensor(comp, c, u_axis) for c in fake_coords[mm]]
                if vals:
                    fo = float(np.mean([v[0] for v in vals]))
                    fe = float(np.mean([v[1] for v in vals]))
                else:
                    fo = fe = np.nan
                r_fake[f"Dpar_FakeControl_{tag}"] = fo
                r_fake[f"Dperp_FakeControl_{tag}"] = fe
                r_fake[f"Dpar_minus_Dperp_FakeControl_{tag}"] = fo - fe
                fth, fc2 = orientation_over_rois(comp, fake_coords[mm], u_axis)
                r_fake[f"FakeControl_{tag}_Theta_deg"] = round(fth, 2) if np.isfinite(fth) else np.nan
                r_fake[f"FakeControl_{tag}_cos2_theta"] = round(fc2, 4) if np.isfinite(fc2) else np.nan

        near_rows.append(r_near)
        fake_rows.append(r_fake)

    qc["n_clusters"] = len(near_rows)
    qc["rotate_u"] = rotate_u
    qc["u_rel_det"] = u_rel_det
    qc["tensor_mtime_min"] = (datetime.datetime.fromtimestamp(tensor_mtime)
                              .isoformat(timespec="seconds")
                              if np.isfinite(tensor_mtime) else "")
    qc["min_cluster_size"] = min_cluster_size
    qc["n_below_min"] = n_below_min
    if near_rows:
        qc["PVS_depth_med"] = round(float(np.median([r["Depth"] for r in near_rows])), 4)
        for rname in regions:
            qc[f"n_{rname}"] = sum(1 for r in near_rows if r["Region"] == rname)
    return near_rows, fake_rows, qc, bal_rows


# ==========================================
# MAIN
# ==========================================
def main():
    ap = argparse.ArgumentParser(description="Three-tier PVS extraction v8.")
    ap.add_argument("--root", default=ROOT_DIR)
    ap.add_argument("--out", default=OUTPUT_XLSX)
    ap.add_argument("--near-orientation-roi", choices=["full", "matched"], default="full",
                    help="v9: ROI for the orientation-referenced adjacent values and the "
                         "adjacent fibre angle: 'full' 3-6 mm shell (= v8, default) or the "
                         "size-matched draws (as the scalar *_NearPVS metrics)")
    ap.add_argument("--margins-mm", type=float, nargs="*", default=FAKE_MARGINS_MM)
    ap.add_argument("--fake-draws", type=int, default=N_FAKE_DRAWS)
    ap.add_argument("--near-draws", type=int, default=N_NEAR_DRAWS,
                    help="contiguous size-matched ROIs drawn inside the NearPVS "
                         "shell; the full-shell value is still written as "
                         "*_NearPVSfull for comparison")
    ap.add_argument("--rotate-u", choices=["header", "off"], default="header",
                    help="v8: map the cluster axis u into the native tensor "
                         "frame via the image-header relation, reflections "
                         "kept (default; step16 v3 convergence-validated). "
                         "'off' reproduces the scrambled v7 behaviour -- "
                         "audit only.")
    ap.add_argument("--dki-root", default=DKI_ROOT_DEFAULT,
                    help="root holding native tensor maps (header relation)")
    ap.add_argument("--min-cluster-size", type=int, default=MIN_CLUSTER_SIZE,
                    help="minimum PVS cluster size in T1w voxels (default "
                         f"{MIN_CLUSTER_SIZE} = primary analysis; e.g. 20 for "
                         "the large-cluster sensitivity check -- use a "
                         "distinct --out so spreadsheets cannot be confused)")
    ap.add_argument("--no-inclusion-filter", action="store_true")
    ap.add_argument("--no-lobes", action="store_true",
                    help="match on CSO_L/CSO_R instead of lobe x hemisphere")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()

    bad = [m for m in a.margins_mm if m < NEAR_BAND_MM[1]]
    if bad:
        raise SystemExit(
            f"Margin(s) {bad} mm < NearPVS outer bound {NEAR_BAND_MM[1]} mm.\n"
            f"FakeControl would then overlap NearPVS in distance and the tiers\n"
            f"would not be ordered. Raise the margin or lower NEAR_BAND_MM.")
    if min(a.margins_mm) < 4.0:
        print(f"[!] margin {min(a.margins_mm)} mm is inside the ~4 mm "
              f"interpolation point-spread of a 2 mm native grid.")

    global USE_LOBES
    if a.no_lobes:
        USE_LOBES = False
    global NEAR_ORIENTATION_ROI
    NEAR_ORIENTATION_ROI = a.near_orientation_roi
    print(f"adjacent-tissue orientation ROI: {NEAR_ORIENTATION_ROI}"
          + ("  (= v8)" if NEAR_ORIENTATION_ROI == "full" else ""))
    rng = np.random.default_rng(RANDOM_SEED)
    if not os.path.isdir(a.root):
        raise SystemExit(f"root not found: {a.root}")
    subjects = sorted(d for d in os.listdir(a.root)
                      if os.path.isdir(os.path.join(a.root, d))
                      and d.startswith("sub-"))
    if not a.no_inclusion_filter:
        try:
            inc = included_subjects()
            before = len(subjects)
            subjects = [s for s in subjects if s in inc]
            print(f"study inclusion filter: {before} -> {len(subjects)} subjects")
        except Exception as exc:
            print(f"[!] inclusion list unreadable ({exc}); processing ALL subjects")

    print(f"Subjects: {len(subjects)}  |  maps from <sub>/{DKI_SUBDIR}/")
    print(f"Tiers (mm from union PVS mask): PVS=inside | "
          f"NearPVS={NEAR_BAND_MM[0]}-{NEAR_BAND_MM[1]} | "
          f"FakeControl>={a.margins_mm}")
    print(f"Controls matched on Region x Depth, {a.fake_draws} draws/cluster")
    print(f"u frame handling: --rotate-u {a.rotate_u}"
          + ("   <-- AUDIT MODE, orientation columns scrambled for reflected "
             "subjects" if a.rotate_u == "off" else " (header relation, v8)"))
    tag = "" if a.min_cluster_size == MIN_CLUSTER_SIZE else \
          "   <-- SENSITIVITY RUN, not the primary analysis"
    print(f"Min cluster size: {a.min_cluster_size} voxels{tag}")
    print("-" * 80)

    errors, provenance = [], []
    all_near, all_fake, all_qc, all_bal = [], [], [], []

    for subj in subjects:
        subj_dir = os.path.join(a.root, subj)
        paths = find_files_for_subject(subj_dir, subj, errors)
        if paths is None:
            print(f"[SKIP ] {subj}")
            continue
        for metric, p in paths["dki_paths"].items():
            provenance.append([subj, metric, os.path.relpath(p, subj_dir)])
        for key, p in [("PVS_mask", paths["pvs_path"]),
                       ("SynthSeg", paths["wm_path"])]:
            provenance.append([subj, key, os.path.relpath(p, subj_dir)])
        for key, p in paths["seg"].items():
            provenance.append([subj, key,
                               os.path.relpath(p, subj_dir) if p else "ABSENT"])
        if a.check:
            miss = [k for k, v in paths["seg"].items() if v is None]
            print(f"[OK   ] {subj}: {len(paths['dki_paths'])}/"
                  f"{len(METRICS_OF_INTEREST)} maps"
                  + (f" | missing seg: {miss}" if miss else ""))
            continue

        near, fake, qc, bal = process_subject(subj, paths, rng, a.margins_mm,
                                              a.fake_draws, errors,
                                              a.near_draws,
                                              a.min_cluster_size,
                                              a.rotate_u, a.dki_root)
        if qc is None:
            print(f"[SKIP ] {subj}: no usable maps")
            continue
        all_near.extend(near); all_fake.extend(fake)
        all_qc.append(qc); all_bal.extend(bal)
        pf = "  ".join(f"m{m:g}:{qc[f'pool_frac_m{m:g}']*100:.0f}%"
                       for m in a.margins_mm)
        rl = "  ".join(f"m{m:g}:{qc[f'relax_region_m{m:g}']}/"
                       f"{qc[f'relax_any_m{m:g}']}/{qc[f'fakeFail_m{m:g}']}"
                       for m in a.margins_mm)
        print(f"[OK   ] {subj}: {qc['n_clusters']} clusters | pool {pf} | "
              f"relax(depth/region/fail) {rl}", flush=True)

    prov_csv = os.path.splitext(a.out)[0] + "_provenance.csv"
    with open(prov_csv, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["subject", "metric", "file_used"])
        w.writerows(provenance)
    print("-" * 80)
    print(f"Provenance -> {prov_csv}")
    if errors:
        print(f"\n{len(errors)} file problems (first 20):")
        for e in errors[:20]:
            print("   ", e)
    if a.check:
        return
    if not all_near:
        raise SystemExit("No data extracted -- see provenance/errors above.")

    dn = pd.DataFrame(all_near)
    df_ = pd.DataFrame(all_fake)
    dq = pd.DataFrame(all_qc)
    db = pd.DataFrame(all_bal)

    tier_rows = [{"Tier": "NearPVS",
                  "dist_min_mm": dn["NearPVS_dist_min"].median(),
                  "dist_med_mm": dn["NearPVS_dist_med"].median(),
                  "dist_mean_mm": dn["NearPVS_dist_mean"].median()}]
    for mm in a.margins_mm:
        c = f"FakeControl_m{mm:g}_dist_min"
        if c in df_:
            tier_rows.append({"Tier": f"FakeControl_m{mm:g}",
                              "dist_min_mm": df_[c].median(),
                              "dist_med_mm": df_[f"FakeControl_m{mm:g}_dist_med"].median(),
                              "dist_mean_mm": df_[f"FakeControl_m{mm:g}_dist_mean"].median()})
    dt = pd.DataFrame(tier_rows)

    with pd.ExcelWriter(a.out) as writer:
        dn.to_excel(writer, sheet_name="Case_A_Near_PVS", index=False)
        df_.to_excel(writer, sheet_name="Case_B_Fake_PVS", index=False)
        dq.to_excel(writer, sheet_name="Subject_QC", index=False)
        dt.to_excel(writer, sheet_name="Tier_Distances", index=False)
        db.to_excel(writer, sheet_name="Match_Balance", index=False)
    print(f"Spreadsheet -> {a.out}")

    print("\nMEASURED TIER SEPARATION (distance to nearest PVS, mm, cohort medians)")
    print(dt.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    if not dt["dist_min_mm"].is_monotonic_increasing:
        print("  [!] tiers are NOT ordered by distance.")
    else:
        print("  tiers are ordered by distance.")

    if "Depth" in dn and dn.Depth.nunique() <= 1:
        print("\n[!] DEPTH IS DEGENERATE -- every cluster has the same value, so")
        print("    controls were matched on region alone. Check the cortex mask.")
    print("\nMATCH BALANCE (PVS vs control normalized depth; 0=cortex, 1=ventricle)")
    for mm in a.margins_mm:
        c = f"ctrl_depth_m{mm:g}"
        if c in db:
            ok = db[[c, "PVS_Depth"]].dropna()
            if len(ok):
                bias = float((ok[c] - ok.PVS_Depth).mean())
                print(f"  m{mm:g}: PVS {ok.PVS_Depth.mean():.3f}  "
                      f"control {ok[c].mean():.3f}  bias {bias:+.4f}")
                if abs(bias) > 0.05:
                    print("      [!] residual depth imbalance > 0.05 -- inspect "
                          "Match_Balance by group before interpreting.")
        lv = db.get(f"match_level_m{mm:g}")
        if lv is not None:
            n = len(lv)
            print(f"      match level: exact {100*(lv==0).mean():.1f}%  "
                  f"depth-relaxed {100*(lv==1).mean():.1f}%  "
                  f"region-relaxed {100*(lv==2).mean():.1f}%  "
                  f"failed {100*(lv==-1).mean():.1f}%  (n={n})")

    if "cos2_theta" in dn:
        c2 = dn["cos2_theta"].dropna()
        print(f"\nORIENTATION: cos^2(theta) between PVS axis and V1 -- "
              f"mean {c2.mean():.3f}, median {c2.median():.3f}")
        print("  a uniformly random direction gives 0.333. If this sits at 0.333")
        print("  irrespective of Linearity and cluster size, either the PVS axis")
        print("  is genuinely unrelated to the fibre axis OR the tensor frame is")
        print("  wrong. Check the rotation in your DKI->T1w transform before")
        print("  interpreting Dpar; component maps warped as scalars are valid")
        print("  for MD/FA/AD/RD but NOT for u^T D u.")
        rel = dn[dn["Linearity"] >= MIN_LINEARITY]
        if len(rel):
            print(f"  Linearity>={MIN_LINEARITY} ({len(rel)}/{len(dn)}): "
                  f"cos^2 {rel['cos2_theta'].mean():.3f}")

    # v9: fibre angle in every tier, same PVS axis u
    if "cos2_theta" in dn and "NearPVS_cos2_theta" in dn:
        print("\nORIENTATION BY TIER: angle between the PVS axis u and the fibre "
              "direction (V1)\n  measured with the same u in each tissue tier "
              "(random direction: cos^2 0.333, median angle 60 deg)")
        tiers_or = [("PVS", dn["cos2_theta"], dn["Theta_deg"]),
                    ("Adjacent WM (3-6 mm)", dn["NearPVS_cos2_theta"], dn["NearPVS_Theta_deg"])]
        for mm in a.margins_mm:
            c = f"FakeControl_m{mm:g}_cos2_theta"
            if c in df_:
                tiers_or.append((f"Matched WM (>={mm:g} mm)", df_[c],
                                 df_[f"FakeControl_m{mm:g}_Theta_deg"]))
        good = dn["Axis_SE_deg"] <= 5 if "Axis_SE_deg" in dn else None
        for name, c2s, ths in tiers_or:
            line = (f"  {name:24s} cos^2 mean {c2s.mean():.3f}   angle median "
                    f"{ths.median():5.1f} deg")
            if good is not None and len(c2s) == len(good):
                line += (f"   | axis SE<=5 deg: cos^2 {c2s[good.values].mean():.3f}, "
                         f"angle {ths[good.values].median():5.1f} deg")
            print(line)

    if "NearPVS_SizeMatched" in dn:
        ok = dn.NearPVS_SizeMatched.mean()
        print(f"\nNEARPVS SIZE MATCHING: {100*ok:.1f}% of clusters matched to "
              f"their own voxel count")
        print(f"  matched ROI median {dn.NearPVS_Voxels.median():.0f} vox  vs  "
              f"full shell median {dn.NearPVS_Full_Voxels.median():.0f} vox")
        print("  compare *_NearPVS against *_NearPVSfull: if the NearPVS dip in")
        print("  the raw group difference is smaller under matching, unequal")
        print("  smoothing was the cause.")

    if "Axis_SE_deg" in dn:
        se = dn["Axis_SE_deg"].dropna()
        if len(se):
            print(f"\nAXIS PRECISION (jackknife angular SE of the PVS axis)")
            print(f"  median {se.median():.1f} deg  p10 {se.quantile(0.1):.1f}  "
                  f"p90 {se.quantile(0.9):.1f}")
            for thr in (5, 10, 20):
                print(f"    <= {thr:2d} deg: {100*(se <= thr).mean():5.1f}% "
                      f"({int((se <= thr).sum())} clusters)")
            print("  Choose ONE threshold from this table, write it down, and use")
            print("  it for the single pre-specified orientation test. This column")
            print("  is group-blind, so filtering on it does not inflate type I")
            print("  error -- but only if the threshold is fixed before the")
            print("  contrast is run.")

    if "Thickness_mm" in dn:
        frac_thick = float((dn["Thickness_mm"] >= NATIVE_VOXEL_MM).mean())
        print(f"\nThickness: median {dn['Thickness_mm'].median():.2f} mm, "
              f"p90 {dn['Thickness_mm'].quantile(0.9):.2f} mm | "
              f">= one native voxel: {frac_thick*100:.1f}%")
        sub_native = float((dn["Size_mm3"] < NATIVE_VOXEL_MM3).mean())
        print(f"  clusters smaller than ONE native DKI voxel: {sub_native*100:.1f}%")
        if frac_thick < 0.05:
            print("  [!] no dynamic range -- report as a limitation, do not plot "
                  "thickness as a partial-volume dose-response.")

    print("\nNEXT: (1) mixed-effects group x tier interaction with subject as a "
          "random\n          effect -- not per-tier Mann-Whitney;"
          "\n      (2) report RAW tier differences alongside Cohen's d, since the"
          "\n          v4 gradient was largely denominator inflation;"
          "\n      (3) test pool_frac_m* and match-level counts across groups;"
          "\n      (4) resolve tensor reorientation before reporting Dpar.")


if __name__ == "__main__":
    main()
