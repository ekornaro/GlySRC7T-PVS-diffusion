# GlySRC7T — diffusion MRI of perivascular spaces: processing and analysis pipeline

This document lists the scripts in this repository in the order they were run for the
manuscript, with the exact command used for each. Paths are given relative to the two
data roots defined in `config.yaml`:

| placeholder | meaning |
|---|---|
| `<DATA_ROOT>` | EvgBrainPipe output tree, one folder per participant: `sub-XXX/pipeline_output/DATASINK/DKISpace/` |
| `<PVS_ROOT>` | T1w-space tree with the PVS segmentations and the maps warped into T1w space: `sub-XXX/DKI_in_T1w_space/` |

Participant inclusion is read everywhere from the study spreadsheet (sheets `HC`, `PPCS`,
`Boxers`; column *Will we used it?* = `Yes`, spelled as in the spreadsheet), giving HC n = 20, PPCS n = 18, RHI n = 21.

Vocabulary used in the code versus the manuscript:

| code | manuscript |
|---|---|
| `Original` tier | PVS |
| `NearPVS` | adjacent white matter, 3–6 mm, same bundle (size-matched draws) |
| `NearPVSfull` | the full 3–6 mm shell (not used in the primary analysis) |
| `FakeControl_m6` / `Fake_m6` | matched white matter ≥ 6 mm |
| `FakeControl_m10` / `Fake_m10` | matched white matter ≥ 10 mm (sensitivity analysis) |
| block A / block B | the two acquisition batches (HC + PPCS, November 2025; RHI, August 2026) |

---

## 0. Preprocessing (outside this repository)

DWI preprocessing (denoising, Gibbs-ringing removal, TOPUP/eddy) and DTI fitting were
run with **EvgBrainPipe** (separate repository). Its outputs under
`<DATA_ROOT>/sub-XXX/pipeline_output/DATASINK/DKISpace/` are the inputs of everything below.
Two segmentations are also upstream of this repository and were produced with methods
developed by co-authors:

- **Perivascular spaces** were segmented on the T1w (RORPO filtering, white-matter
  hyperintensities removed with the FLAIR, SynthSeg parcellation) following
  Duarte Coello R. et al., medRxiv 2025, doi:10.1101/2025.10.22.25337423
  (<https://www.medrxiv.org/content/10.1101/2025.10.22.25337423v4.full>).
- **Parasagittal dural (PSD) spaces** were segmented along the superior sagittal sinus
  on the FLAIR for regional volumetry with SPESIS (Hett K.,
  <https://github.com/hettk/spesis>).

## 1. Registration and model fitting (native diffusion space)

| step | script | command | writes |
|---|---|---|---|
| 1 | `step01_coregister_T1w_and_bsmall.py` | `python step01_coregister_T1w_and_bsmall.py` | b=0 ↔ T1w transforms (ANTs) |
| 2 | `step02_run_fwdki_on_cohort.py` | `python step02_run_fwdki_on_cohort.py` | free-water-eliminated DKI maps → `DKISpace/fwe_dki/` (uses `fwdki.py`) |
| 3 | `step03_make_lowb_maps.py` | `python step03_make_lowb_maps.py` | low-b (b = 100, 500) tensor maps → `DKISpace/lowb/` |
| 4 | `step04_fit_data_to_DKI_model.py` | `python step04_fit_data_to_DKI_model.py` | full DKI fit, tensor components Dxx…Dyz and V1 → `DKISpace/tensor/` |
| 5 | `step05_fix_V1_3D-to-4D_issue.py` | `python step05_fix_V1_3D-to-4D_issue.py` | V1 stored as a 4-D vector image after warping |
| 6 | `step06_pvs_orientation_analysis.py` | `python step06_pvs_orientation_analysis.py` | first version of the Dpar/Dperp/θ maps (superseded by step 13; kept for the record) |

## 2. Frame diagnostics and the gradient-convention fix

Steps 7–12 are the quality-control chain that found and corrected two acquisition-batch
problems. None of them produces a number reported in the paper, but the reported
orientation results depend on the fix they led to, so they are part of the pipeline.

| step | script | command | purpose |
|---|---|---|---|
| 7 | `step07_check_native_tensor_frame.py` | `python step07_check_native_tensor_frame.py` | coherence of the native tensors with brainstem anatomy, per batch |
| 8 | `step08_diagnose_block_difference.py` | `python step08_diagnose_block_difference.py` | localises the batch difference (b0-vs-T1w header reflection) |
| 9 | `step09_determine_gradient_convention.py` | `python step09_determine_gradient_convention.py --all --out ./gradient_convention` | tractography-based scoring of the gradient sign/permutation conventions per participant: 4 with `--conventions basic` (the default, as run) or 24 with `--conventions extended` |
| 10 | `step10_fix_blockA_bvec_convention.py` | `python step10_fix_blockA_bvec_convention.py` | writes the corrected `DKI_corrected_convfix.bvec` for the block-A participants |
| 11 | `step11_rerun_fullDKI_convfix_gradients.py` | `python step11_rerun_fullDKI_convfix_gradients.py` | re-fits the DKI model with the corrected gradients |
| 12 | `step12_determine_gradient_convention.py` | `python step12_determine_gradient_convention.py --all --grad-stem DKI_corrected` and `python step12_determine_gradient_convention.py --all` | confirms the convention after the fix (corrected and original gradient sets) |

**What was found.** In native space, before any warping, the tensors of the two
acquisition batches differed in how well their principal direction followed brainstem
anatomy: block A scored at chance (coherence 0.50) whereas block B scored 0.72. Scoring
the gradient conventions per participant (the four sign flips of the default set) showed that block A was stored
with the z-component of the gradient directions sign-reversed (the bvecs of the 43 block-A
participants listed in step 10 were corrected; 24/24 block-B participants preferred the stored convention). After
correcting the bvecs and re-fitting (steps 10–11), all participants scored best with the
stored convention (step 12). Independently, step 8 showed that the b=0 and T1w images of
25 participants are stored with frames that differ by a reflection (determinant −1),
almost perfectly aligned with batch. A reflection cannot be represented by a rotation, so
any correction that forces a proper rotation scrambles the PVS-axis-to-tensor angle for
those participants; this is handled in step 13 by mapping the PVS axis through the two
image headers with the determinant preserved (`--rotate-u header`).

## 3. Orientation maps and warping to T1w space

Step 15 was run before steps 13 and 14.

| step | script | command | writes |
|---|---|---|---|
| 15 | `step15_warp_all_maps_to_T1w.py` | `python step15_warp_all_maps_to_T1w.py --apply` | all scalar maps (DKI, low-b, tensor components, V1) warped into T1w space (dry run without `--apply`) |
| 13 | `step13_pvs_orientation_analysis.py` | `python step13_pvs_orientation_analysis.py --rotate-u header --overwrite` | Dpar, Dperp, Dpar−Dperp, θ, cluster linearity and distance-to-PVS maps in `<PVS_ROOT>/sub-XXX/DKI_in_T1w_space/`, with a provenance sidecar. `--rotate-u header` is the mode used for every reported result; the default `off` reproduces the earlier diagnostic run. |
| 14 | `step14_theta_alignment.py` | `python step14_theta_alignment.py --out ./theta_alignment` | cohort-level θ statistics and the convergence check (reflected vs non-reflected participants) |

## 4. Tier extraction and spatial-specificity statistics

| step | script | command | writes |
|---|---|---|---|
| 16 | `step16_pvs_reference_tissue_extraction.py` | `python step16_pvs_reference_tissue_extraction.py > extraction_$(date +%Y%m%d_%H%M).log 2>&1` | per-cluster and per-participant metrics in the PVS, adjacent WM, matched WM (≥ 6 and ≥ 10 mm) and ventricular CSF tiers, cluster geometry and the fibre angle in every tier → `pvs_tier_metrics_v9.xlsx` (the *tier workbook* below) |
| 17 | `step17_checking_registration_confounds.py` | `python step17_checking_registration_confounds.py --root <PVS_ROOT> --q1 --q2 --q3 --metrics-xlsx pvs_tier_metrics_v8.xlsx` | registration-quality covariates per participant → `registration_confounds_subject_covariates.csv` (Q1: is the residual warp rotation PVS-coincident or group-dependent; Q2/Q3: head-size and scaling covariates). The v8 and v9 workbooks hold identical per-participant medians for the scalar metrics; v9 only adds the fibre-angle columns. |

### Group × tier mixed-effects models (step 18)

`value ~ Group × Tier + ICV + |log det J| + scale_x`, random intercept per participant,
likelihood-ratio test of the interaction (ML), cluster-robust Wald F as fallback, and a
participant-level permutation test (labels permuted between participants, tiers kept
together). One Benjamini–Hochberg family across the 14 metrics.

| run | command | role |
|---|---|---|
| 18a — **primary** | `python step18_test_group_by_tier_interaction.py --results pvs_tier_metrics_v9.xlsx --covariates registration_confounds_subject_covariates.csv --tiers PVS NearPVS Fake_m6 --within-refs Fake_m6 NearPVS --include-orientation --n-perm 10000 --out MixedEffects_v12_3TIERS_cov_10k.xlsx` | the three tiers defined in the Methods (PVS, adjacent WM, matched WM ≥ 6 mm); source of the supplementary group × tier table |
| 18b — sensitivity, all tiers | `python step18_test_group_by_tier_interaction.py --results pvs_tier_metrics_v9.xlsx --covariates registration_confounds_subject_covariates.csv --include-orientation --within-refs Fake_m6 NearPVS --n-perm 10000 --out MixedEffects_v11_ALLTIERS_cov_10k.xlsx` | adds the full 3–6 mm shell and the ≥ 10 mm tier |
| 18c — sensitivity, reflected participants | `python step18_test_group_by_tier_interaction.py --results pvs_tier_metrics_v9.xlsx --covariates registration_confounds_subject_covariates.csv --include-orientation --no-covariates --within-refs Fake_m6 NearPVS --subject-list subject_lists/det_neg_reflected_v8.txt --out MixedEffects_v11_detneg_nocov.xlsx` | the 25 participants with a b0/T1w header reflection, to verify that the step-13 correction removed the batch effect; covariates are dropped because they are collinear with group in this subset |

### Figures for the tier analysis (step 19)

| step | command | writes |
|---|---|---|
| 19 | `python step19_plot_pvs_tier_figures.py --mixed-xlsx MixedEffects_v12_3TIERS_cov_10k.xlsx` | main Figure 2 (`Fig_main_PVS_distinct`), the paired PVS-vs-reference supplementary figures (`Fig_suppl_ref-<tier>`), the match-quality and eligible-WM figures, and `Figure_statistics_v16-3.xlsx` |

The reported Figure_statistics_v16-3.xlsx and figures were produced with --mixed-xlsx MixedEffects_v11_ALLTIERS_cov_10k.xlsx; the 3-tier workbook gives identical output.

### Free-water-elimination reliability gate (step 20)

| step | command | result |
|---|---|---|
| 20 | `python step20_fwe_reliability_gate.py --subjects sub-105 sub-125 sub-141 sub-155 sub-160 sub-170 --metrics-xlsx pvs_tier_metrics_v8.xlsx --out ./fwe_gate` | residual-bootstrap precision of each free-water-corrected metric in PVS and matched WM, with a priori pass/fail thresholds (reliability ratio R < 0.20 pass, > 0.50 fail; noise amplification ≤ 4). All seven corrected metrics failed (R = 0.8–1.8): subject-level estimation noise was comparable to the between-participant variance within PVS, reflecting the 1/(1−f) amplification of the deconfounding step in high-f voxels. The FWE metrics are therefore excluded from the paper. Uses `fwdki.py` directly. |

## 5. Regional group comparisons, volumes and supplementary material

| step | script / command | writes |
|---|---|---|
| 21a–b | `step21a_SPM12_eTIV.m` (calls the function in `step21b_run_spm12_TIV_batch.m`; MATLAB + SPM12) | SPM12 tissue-class total intracranial volume per participant → `TIV_spm12_results.csv`, copied by hand into `TIV_spm12_results.xlsx` (one sheet per group) |
| 21c | `step21c_plot_eTIV_three_methods.py` (see `step21c_plot_eTIV_three_methods.bash` for the arguments) | eTIV agreement between SPM12 (`TIV_ml`), the SynthStrip brain mask, and recon-all's atlas-scaling eTIV (`EstimatedTotalIntraCranialVol`, labelled "recon-all" in the plot), read from `FreeSurfer/synthseg/volumetry_scores.csv`. These files were assembled by hand from recon-all output (`asegstats2table`) outside this pipeline; despite their folder and file name they are not SynthSeg output; supplementary figure |
| 22 | `python step22_create_csv_based_on_niftis_and_SPM12ICV.py` | per-participant regional medians of every metric in the three pre-specified regions (centrum semiovale, basal ganglia, midbrain), normalised volumes and PVS volume fraction |
| 23 | `python step23_plot_all_markers_across_cohorts_and_ROIs.py --min-pvs-voxels 10 --max-diffusivity 3.0e-3 --also-plot-absolute` (`step23_plot_all_markers_across_cohorts_and_ROIs.bash`) | the 39-test family of group comparisons (Kruskal–Wallis, pairwise Mann–Whitney with Cliff's Δ, Benjamini–Hochberg) → `stats_all_corrected_v13.xlsx`, the heat map (Figure 3) and Figures 4–6; its statistics are tabulated in Supplementary Tables S2, S3 and S6 by step 25 |
| 24 | `step24_psd_group_stats_and_plots.py` (see `step24_psd_group_stats_and_plots.bash`) | parasagittal dural (PSD) space volumes — total and frontal/parietal/occipital segments along the superior sagittal sinus, from the upstream SPESIS segmentation — compared across groups (Kruskal–Wallis, pairwise Mann–Whitney with Cliff's Δ; supplementary figure and table) |
| 25 | `step25_make_supplementary_tables.bash` (step-23 workbook and the primary step-18 workbook `MixedEffects_v12_3TIERS_cov_10k.xlsx`) | Supplementary Tables S2–S6 (S2 within-PVS metrics and relative volumes, S3 PVS-weighted metrics, S4 within-participant contrasts and group × tier interaction, S5 group sizes and exclusions, S6 SynthSeg structure volumes) from the workbooks above |
| 26 | `python step26_qc_registration_overlay.py` | visual registration check: MD (in T1w space) over the T1w in 12 axial slices per participant, PNGs in `<PVS_ROOT>/QC_MD_overlays_png/`; the co-registration supplementary figure was assembled from them by hand |

---

## Scripts not released

Not released: `step20a_verify_ants_convention.py` (an ANTs-convention correction that was
set aside in favour of the measured-limitation approach: the residual registration rotation
is reported as a limitation, supported by step 17 Q1), `step24_wholeWM_benchmark.py` (an
exploratory whole-white-matter benchmark tier not used in the paper). The S0
extrapolation harness used during the b=0 investigation, which was superseded by step 20, is
included as `validate_s0_extrapolation.py`.
