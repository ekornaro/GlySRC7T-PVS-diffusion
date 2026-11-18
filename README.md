# GlySRC7T: diffusion MRI of perivascular spaces at 7 T

Analysis code for the Brain Communications manuscript on perivascular spaces (PVS) studied
with 7 T diffusion kurtosis imaging in three groups:

| group | n | in the code |
|---|---|---|
| healthy controls | 20 | `HC` |
| persistent post-concussive symptoms | 18 | `PPCS` |
| repetitive head impacts (boxers) | 21 | `RHI`; sheet and folder name `Boxers` |

The pipeline does the following:

1. Fits DKI, low-b and free-water models in native diffusion space.
2. Fixes a gradient-convention error in one acquisition batch.
3. Warps every map to T1w space.
4. Asks whether diffusion inside PVS differs from the white matter around them (the spatial-specificity analysis) and whether that difference depends on group.
5. Compares regional diffusion and volumes across the three groups.

The scripts are numbered in the order they were run. `docs/pipeline_notes.md` is the run log.
[Notes on individual steps](#notes-on-individual-steps) lists corrections and design choices a reader should know
before re-running anything.

No data are included. The participant data are not public.

---

## Setup

### 1. Paths

All site-specific paths are read from `config.yaml`:

```bash
cp config.example.yaml config.yaml      # then edit the paths
python3 pipeline_config.py              # prints every resolved path
```

- `config.yaml` is git-ignored.
- `$GLYSRC7T_CONFIG` can point to a config file elsewhere.
- Command-line flags still override the config defaults.
- Analysis parameters (thresholds, seeds, margins, draw and permutation counts) are not in the config. They stay in the scripts and their flags.
- MATLAB reads the config through `pipeline_config_value.m`; the bash wrappers call `python3 pipeline_config.py <key>`.

| key | what it points to |
|---|---|
| `data_root` | EvgBrainPipe output tree, `sub-XXX/pipeline_output/DATASINK/DKISpace/` |
| `pvs_root` | T1w-space tree: `sub-XXX/{PVS,Segmentation,DKI_in_T1w_space,FreeSurfer}/` |
| `pvs_segmentations_dir` | parent of the `Lund_7T_<case>` segmentation trees (step22) |
| `study_xlsx` | study spreadsheet: sheets `HC`, `PPCS`, `Boxers`; inclusion column `Will we used it?` = `Yes` |
| `synthseg_labels_xlsx` | SynthSeg label table (step22) |
| `freesurfer_root` | tree holding `sub-XXX/FreeSurfer/synthseg/volumetry_scores.csv`, the hand-assembled recon-all eTIV files read by step21c |
| `psd_csv` | SPESIS parasagittal-dural-space volumes (step24) |
| `spm12_dir` | SPM12 installation (step21b) |
| `fwdki_py` | the free-water DKI module, `fwdki.py` in this repository (step20) |
| `tiv_spm12_csv`, `tiv_spm12_xlsx` | SPM12 TIV output of step21b, and the hand-made workbook built from it (see step 21) |
| `workbooks_dir` | tier workbooks: steps 16–20 and 25 read and write here |
| `results_dir` | per-participant spreadsheets (step22) and group statistics (step23, step25) |
| `tier_figures_dir` | output of step19 |

### 2. Software

**Python 3.12.** Install with conda or pip:

```bash
conda env create -f environment.yml && conda activate glysrc7t
# or
pip install -r requirements.txt
```

Versions are pinned to the environment the analysis ran in. The main packages are:

| package | version |
|---|---|
| numpy | 2.2.6 |
| scipy | 1.14.1 |
| pandas | 2.2.2 |
| dipy | 1.12.1 |
| nibabel | 5.2.1 |
| statsmodels | 0.14.6 |
| scikit-image | 0.26.0 |
| matplotlib | 3.9.2 |
| seaborn | 0.13.2 |
| openpyxl | 3.1.5 |
| python-docx | 1.2.0 |
| tabulate | 0.9.0 |

**Software outside Python:**

| tool | version used | needed by |
|---|---|---|
| ANTs (`antsRegistration`, `antsApplyTransforms`) | 2.5.4.post7-g9c98082 | steps 01, 05, 07, 09, 12, 15, 20 (must be on `PATH`) |
| MATLAB | R2025a (25.1.0.2973910, Update 1) | steps 21a–b |
| SPM12 | release 7771 | steps 21a–b |

The following tools ran upstream of this repository. Their outputs are inputs here:

- **EvgBrainPipe:** DWI denoising, Gibbs-ringing removal, TOPUP/eddy and DTI.
- **FreeSurfer `mri_synthseg`:** SynthSeg parcellation.
- **FreeSurfer `recon-all`:** `EstimatedTotalIntraCranialVol`, the third eTIV method of step21c.
- **SynthStrip:** brain mask.
- **The PVS segmentation:** Duarte Coello et al., medRxiv 2025, doi:10.1101/2025.10.22.25337423.
- **SPESIS:** parasagittal dural spaces, <https://github.com/hettk/spesis>.

---

## Pipeline

Run the commands from this folder unless a step says otherwise. Steps 17–20 were run from the
workbooks folder, because their default inputs and outputs are bare file names:

```bash
REPO=$PWD
cd "$(python3 pipeline_config.py workbooks_dir)"
python "$REPO/step18_..."
```

Where a command below differs from the run log, a note under the block says how.

### Native space: registration and model fitting

```bash
python step01_coregister_T1w_and_bsmall.py          # b0 <-> T1w transforms (ANTs); takes no arguments
python step02_run_fwdki_on_cohort.py                # free-water DKI -> DKISpace/fwe_dki/ (uses fwdki.py)
python step03_make_lowb_maps.py                     # low-b (b = 100, 500) tensor maps -> DKISpace/lowb/
python step04_fit_data_to_DKI_model.py              # full DKI, tensor components and V1 -> DKISpace/tensor/
python step05_fix_V1_3D-to-4D_issue.py              # V1 stored as a 4-D vector image
python step06_pvs_orientation_analysis.py           # DIAGNOSTIC: first orientation maps, superseded by step13
```

### Frame diagnostics and the gradient-convention fix (QC chain)

None of these steps produces a number reported in the paper. The reported orientation results
depend on the fix in steps 10–11.

```bash
python step07_check_native_tensor_frame.py
python step08_diagnose_block_difference.py
python step09_determine_gradient_convention.py --all --out ./gradient_convention
python step10_fix_blockA_bvec_convention.py         # writes DKI_corrected_convfix.bvec for the 43 block-A participants
python step11_rerun_fullDKI_convfix_gradients.py
python step12_determine_gradient_convention.py --all --grad-stem DKI_corrected   # run 1: original gradients
python step12_determine_gradient_convention.py --all                             # run 2: corrected (convfix) gradients
```

- step09 scores 4 conventions with `--conventions basic`, the default and the setting that was run, or 24 with `--conventions extended`.
- Steps 09, 12 and 20 apply the T1w→native transforms in the corrected order, `-t [Affine,1] -t InverseWarp`, and were re-run with it in 2026-11. See [Notes on individual steps](#notes-on-individual-steps).

### Warping to T1w space and orientation maps

```bash
python step15_warp_all_maps_to_T1w.py --apply                     # runs first; without --apply it is a dry run
python step13_pvs_orientation_analysis.py --rotate-u header --overwrite
python step14_theta_alignment.py --out ./theta_alignment
python step26_qc_registration_overlay.py                         # visual registration check (PNGs)
```

- step15 was run before step13 and step14. Only the `--apply` run counts.
- step13: the reported results used `--rotate-u header`. The script default, `off`, reproduces the earlier diagnostic behaviour.
- step26 reads `dkifit_MD_toT1w_Warped.nii.gz`. On the current tree, step15 writes `MD_toT1w_Warped.nii.gz` instead, so a re-run today finds only some participants. The PNGs used for the figure date from 2026-08-20.

### Tier extraction and the spatial-specificity analysis

These steps are run from `workbooks_dir`; see the top of this section.

```bash
python $REPO/step16_pvs_reference_tissue_extraction.py > extraction_$(date +%Y%m%d_%H%M).log 2>&1
python $REPO/step17_checking_registration_confounds.py --root <PVS_ROOT> --q1 --q2 --q3 \
       --metrics-xlsx pvs_tier_metrics_v8.xlsx
```

- step16 writes the tier workbook `pvs_tier_metrics_v9.xlsx` to `workbooks_dir`.
- step17's `--root` defaults to `pvs_root`, so it can be left out.

Group × tier mixed-effects models (step18). The permutation seed is fixed at 0, and the
reported runs used `--n-perm 10000`; the script default is 1000.

```bash
# 18a  PRIMARY: the three tiers defined in the Methods
python $REPO/step18_test_group_by_tier_interaction.py --results pvs_tier_metrics_v9.xlsx \
       --covariates registration_confounds_subject_covariates.csv --tiers PVS NearPVS Fake_m6 \
       --within-refs Fake_m6 NearPVS --include-orientation --n-perm 10000 --out MixedEffects_v12_3TIERS_cov_10k.xlsx
# 18b  SENSITIVITY: all tiers (adds the full 3-6 mm shell and the >= 10 mm tier)
python $REPO/step18_test_group_by_tier_interaction.py --results pvs_tier_metrics_v9.xlsx \
       --covariates registration_confounds_subject_covariates.csv --include-orientation \
       --within-refs Fake_m6 NearPVS --n-perm 10000 --out MixedEffects_v11_ALLTIERS_cov_10k.xlsx
# 18c  SENSITIVITY: the 25 participants with a b0/T1w header reflection, no covariates
python $REPO/step18_test_group_by_tier_interaction.py --results pvs_tier_metrics_v9.xlsx \
       --covariates registration_confounds_subject_covariates.csv --include-orientation --no-covariates \
       --within-refs Fake_m6 NearPVS --subject-list subject_lists/det_neg_reflected_v8.txt --out MixedEffects_v11_detneg_nocov.xlsx
```

18c was run without `--n-perm`, so it used the default of 1000 permutations.

```bash
python $REPO/step19_plot_pvs_tier_figures.py --mixed-xlsx MixedEffects_v12_3TIERS_cov_10k.xlsx
```

The reported Figure_statistics_v16-3.xlsx and figures were produced with --mixed-xlsx MixedEffects_v11_ALLTIERS_cov_10k.xlsx; the 3-tier workbook gives identical output.
step19 jitters the data points without a fixed seed, so a re-run reproduces every number and q-value but not the exact positions of the dots.

step19 writes to `tier_figures_dir`:

- main Figure 2;
- the supplementary figures (see the map below);
- `Figure_statistics_v16-3.xlsx`.

### Free-water reliability gate (diagnostic; free-water-corrected metrics excluded from the paper)

```bash
python $REPO/step20_fwe_reliability_gate.py --subjects sub-105 sub-125 sub-141 sub-155 sub-160 sub-170 \
       --metrics-xlsx pvs_tier_metrics_v8.xlsx --out ./fwe_gate
```

Reading the output (the verdict file and R for `f` are under [Notes on individual steps](#notes-on-individual-steps)):

- **Regression check.** The re-run of 2026-11 (old transform order) gives median R = 0.79–1.75 across the six participants, matching the 0.8–1.8 of the original run. The corrected-order run gives 0.84–1.72.
- **`dirsub_roimed_sd` is exactly 0 for MKt and RKt.** This is a property of the corrected fit on this data, not a failed computation:
  - On the unresampled signal, 34–42 % of the ROI voxels come out of the corrected kurtosis fit as exactly 0.
  - A further 27–45 % are clipped at DIPY's lower kurtosis bound of −3/7.
  - So more than half of the ROI sits on these two values, and the ROI median is 0 in every direction-subsample draw.
  - AKt has only about 2 % of voxels at the bound, so its median moves between draws.

### Volumes, regional group comparisons and supplementary material

```bash
matlab -batch "run('step21a_SPM12_eTIV.m')"          # SPM12 TIV (calls step21b_run_spm12_TIV_batch.m) -> TIV_spm12_results.csv
# MANUAL STEP: build TIV_spm12_results.xlsx from the CSV (below)
bash step21c_plot_eTIV_three_methods.bash            # eTIV by three methods (supplementary figure)
python step22_create_csv_based_on_niftis_and_SPM12ICV.py vol
python step22_create_csv_based_on_niftis_and_SPM12ICV.py diff
bash step23_plot_all_markers_across_cohorts_and_ROIs.bash
bash step24_psd_group_stats_and_plots.bash           # PSD volumes; without --apply the script is a dry run
bash step25_make_supplementary_tables.bash           # Supplementary Tables S2-S6 -> ./supplementary_tables/
```

**Manual step between step21b and step21c/step22.** The SPM12 output
`TIV_spm12_results.csv` has one row per participant. Its columns are:

- `subject`, `TIV_ml`, `GM_ml`, `WM_ml`, `CSF_ml`;
- `TIV_masked_ml`, `GM_masked_ml`, `WM_masked_ml`, `CSF_masked_ml`;
- `T1_path`.

The rows were copied by hand into `TIV_spm12_results.xlsx`, one sheet per group: `HC`, `PPCS`,
`Boxers`, and `Miscelleneous` for participants not analysed.

- Step21c and step22 read only the sheets `HC`, `PPCS` and `Boxers`, and from each only the columns `subject` (e.g. `sub-105_ses-01_T1wbrain_restore`) and `TIV_ml`, converted to mm³.
- The columns added by hand (`Subj`, `Stratum`, `Will we used it?`) are not read.
- In the workbook used for the paper, all 60 rows match the CSV exactly.

Notes on the commands:

- **step21c, the three eTIV methods:**
  - SPM12 `TIV_ml`;
  - the SynthStrip brain-mask volume;
  - recon-all's atlas-scaling eTIV (`EstimatedTotalIntraCranialVol`), labelled "recon-all" in the plot. step21c reads it from the first row, second column of `<freesurfer_root>/sub-XXX/FreeSurfer/synthseg/volumetry_scores.csv`.
  - Those files were assembled by hand from recon-all output (`asegstats2table`) outside this pipeline, so that step21c finds all its inputs in one place. Despite their folder and file name, they are not SynthSeg output.
- **step22** was run twice: `vol` (volumetry) and `diff` (diffusion). The default is `diff`; `both` does both in one run. It reads `synthseg_labels_xlsx` and `tiv_spm12_xlsx` when it loads, so even `--help` needs those two files.
- **step23:** the thresholds `--min-pvs-voxels 10 --max-diffusivity 3.0e-3` are required flags and are set in the wrapper. `--also-plot-absolute` adds plots of absolute volumes that belong to no FDR family and are never tested.
- **step25:** the wrapper uses the step23 workbook `stats_all_corrected_v13.xlsx` and the primary step18 workbook `MixedEffects_v12_3TIERS_cov_10k.xlsx`.

---

## Notes on individual steps

None of the points below affects a result reported in the paper.

**Transform order (steps 09, 12, 20).** `antsApplyTransforms` applies the listed transforms
last-to-first, so taking a T1w-space image into native diffusion space needs the inverse of the
b0→T1w registration as `-t [0GenericAffine.mat,1] -t 1InverseWarp.nii.gz` (as in step07).
Steps 09, 12 and 20 use this order in this release and were re-run with it in 2026-11. Runs
before that listed the two transforms the other way round. The re-runs showed:

- step12 on the original gradients (`--grad-stem DKI_corrected`): all 41 scored block-A
  participants score best with `flip_z`, all 24 scored block-B participants with the stored
  convention (`none`). This is the split that led to the fix in step10.
- step12 on the corrected gradients: all 65 scored participants score best with `none`.
- step09: all 65 score best with `none`, because it reads the corrected gradients (next item).
- sub-101, sub-127, sub-130 and sub-173 were not scored in any of the three runs, because an
  input file is missing.
- step20: see "step20 verdicts" below.

**step09 reads the corrected gradients when they exist.** step09 uses
`DKI_corrected_convfix.bvec/.bval` whenever that pair exists, and falls back to the first
`*.bvec`/`*.bval` in the folder only when it does not. After step10 the corrected pair exists,
so a re-run of step09 scores the corrected gradients. The check before the correction is step12
with `--grad-stem DKI_corrected`.

**step20 verdicts.** step20's numbers are not reported in the paper.

- In the corrected-order run, `fwe_gate_verdict.json` gives `FAIL (f sanity)` for every metric.
  The f sanity check (f ≥ 0.7 in the ventricles, ≤ 0.25 in WM, PVS > WM) failed in sub-155 and
  sub-160, and the rule lets that failure override R for all metrics.
- Median R across the six participants is 0.84–1.72 for the seven corrected metrics. Without the
  override, all seven still fail on R > 0.50; MDt, MKt and RKt also exceed the noise-amplification
  limit of 4.
- R for `f` is not interpretable; read it as n/a. It divides the SD of a fraction by the
  between-participant SD of `lowb_ADC`, a diffusivity (`COUNTERPART` in the script). The rule was
  left unchanged.

**Run 18c subject list.** `subject_lists/det_neg_reflected_v8.txt` is not in the repository. It
lists, one ID per line, the 25 participants whose b0 and T1w frames differ by a reflection. To
rebuild it, take the participants with `u_rel_det` = −1 in the `orientation_maps_qc.csv` that
step13 writes with `--rotate-u header` (to `--out`, default `./pvs_orientation_maps_qc`):

```bash
mkdir -p subject_lists
python -c "import pandas as pd; d = pd.read_csv('pvs_orientation_maps_qc/orientation_maps_qc.csv'); \
d.loc[d.u_rel_det < 0, 'Subject'].to_csv('subject_lists/det_neg_reflected_v8.txt', index=False, header=False)"
```

The list should hold 25 IDs. step18 stops if an ID in it is not in the tier workbook.

**step18 permutation design.**

- The seed is fixed at 0. The reported runs used `--n-perm 10000`; the script default is 1000,
  and run 18c used that default.
- The HC-versus-RHI permutation test draws group labels from all 59 participants, so its null
  hypothesis is full exchangeability across the three groups. A variant permuting only within the
  41 HC and RHI participants was not run.

**step06 versus step13.** step06 is the earlier, diagnostic orientation analysis, without the
header correction. Every reported orientation result comes from step13 with `--rotate-u header`,
which maps the PVS axis through the two image headers with the determinant kept, so the 25
participants whose b0 and T1w frames differ by a reflection are handled correctly. step13's
default is still `off`, which reproduces the earlier behaviour.

---

## Figures and tables

| item | source |
|---|---|
| Figure 1, overview | made by hand |
| Figure 2 (PVS distinct from reference tissue; fibre angle θ) | step19 |
| Figure 3, heat map | step23 |
| Figures 4–6, panels A–C (images) | made by hand |
| Figures 4–6, panels D–L (metrics) | step23 |
| Table 1, demographics | compiled by hand |
| Suppl. Table S1, acquisition parameters | compiled by hand |
| Suppl. Tables S2–S6 | step25, from the step23 and step18a workbooks |
| Suppl. figure: eTIV by three methods | step21c |
| Suppl. figure: PSD volumetry | step24 |
| Suppl. figure: matched white-matter comparability (match quality and eligible WM) | step19 |
| Suppl. figures: the two paired PVS-versus-reference figures | step19 |
| Suppl. figures: CSF volumetry; tissue and PVS volumetry; PVS-weighted metrics; FA within PVS | step23 |
| Suppl. figure: co-registration check | assembled by hand from the step26 PNGs |
| Supplement text, handedness check (median θ 44.0° in 25 vs 43.6° in 34 participants, Mann–Whitney P = 0.57) | not step13 or step14 (see below) |

θ has no separate supplementary figure; it is shown in Figure 2.

**The handedness check** comes from neither step13 nor step14:

- step13 and step14 print a similar reflected-versus-non-reflected comparison. On their saved outputs (58 participants) it gives 45.3° vs 44.1° (P = 0.81) and 44.8° vs 46.1° (P = 0.35) respectively.
- The reported numbers are the per-participant median of `NearPVS_Theta_deg`, the fibre angle in the adjacent-WM tier. They come from the step16 tier workbook (`Case_A_Near_PVS`), split by `u_rel_det` in its `Subject_QC` sheet: 25 participants with det = −1 and 34 with det = +1.
- A 13-line script outside this repository, `check_T1w_handedness.py`, computed them.
- The same numbers can be reproduced from `pvs_tier_metrics_v9.xlsx` alone.

---

## Primary, sensitivity and diagnostic runs

| run | status |
|---|---|
| step13 `--rotate-u header` | **primary** orientation maps |
| step06; step13 with the default `off` | diagnostic, superseded |
| steps 07–12 | quality-control chain; only the fix from steps 10–11 enters the results |
| step12 | run twice, on the original (`--grad-stem DKI_corrected`) and the corrected gradients |
| step14 | cohort-level θ statistics and the reflection convergence check |
| step16 tiers `Original`, `NearPVS`, `Fake_m6` | **primary** tiers |
| step16 tiers `NearPVSfull`, `Fake_m10` | sensitivity |
| step18a, 3 tiers, `MixedEffects_v12_3TIERS_cov_10k.xlsx` | **primary**; Suppl. Table S4 |
| step18b, all tiers, `MixedEffects_v11_ALLTIERS_cov_10k.xlsx` | sensitivity |
| step18c, 25 reflected participants, `MixedEffects_v11_detneg_nocov.xlsx` | sensitivity |
| step20 | diagnostic: the free-water-corrected metrics failed and are excluded |
| step23 families A, E, F | **primary** group comparisons |
| step23 `ABS` absolute volumes | display only, not tested |
| `validate_s0_extrapolation.py` | b=0 investigation, superseded by step20; steps 03 and 11 import its file finders |

---

## Glossary

| term | meaning |
|---|---|
| `Original` | the PVS tier |
| `NearPVS` | adjacent white matter, 3–6 mm from the PVS in the same bundle (size-matched draws) |
| `NearPVSfull` | the full 3–6 mm shell (sensitivity only) |
| `FakeControl_m6` / `Fake_m6` | matched white matter ≥ 6 mm from any PVS |
| `FakeControl_m10` / `Fake_m10` | matched white matter ≥ 10 mm (sensitivity analysis) |
| block A / block B | the two acquisition batches: HC + PPCS (November 2025) and RHI (August 2026) |
| `u_rel_det`, reflected participants | det = −1 between the b0 and T1w image frames (25 participants); handled by `--rotate-u header` |
| `convfix` | the corrected block-A gradient table, `DKI_corrected_convfix.bvec` (step10) |
| families A, E, F | the FDR families of step23: A = diffusion within PVS (including low-b), PVS volume / ROI volume and ROI volume / eTIV in the three PVS ROIs; E = PVS-weighted diffusion; F = SynthSeg structure volume / eTIV |
| `_mf` (`MK_mf`, …) | kurtosis maps after median filtering |
| `lowb_*` | tensor metrics from the low-b shells (b = 100, 500) |
| `PVSROI` / `synthseg` | the study ROI (centrum semiovale, basal ganglia, midbrain) / the whole SynthSeg structure |

**Version tags in file names.** Tags such as `v8`, `v9`, `v11`, `v12`, `v13` and `v16-3` record
the revision-era version of the script, or of the run, that wrote a file. They are kept so that
the files match the manuscript's provenance records. The numbers count separately for each
script, so `v12` in two file names does not mean the same run.

| file | written by |
|---|---|
| `pvs_tier_metrics_v8.xlsx`, `pvs_tier_metrics_v9.xlsx` | step16, versions 8 and 9. They hold identical per-participant medians for the scalar metrics; v9 adds the fibre-angle columns. step17 and step20 read v8; steps 18 and 19 read v9. |
| `MixedEffects_v11_…`, `MixedEffects_v12_…` | runs of step18 |
| `…_v13` (`stats_all_corrected_v13.xlsx`, `conservative_PVS_segmentation_v13/`) | step23, version 13 |
| `Figure_statistics_v16-3.xlsx` | step19, version 16-3 |
| `det_neg_reflected_v8.txt` | the list of reflected participants from the v8 extraction |

---

## Citation and licence

Citation metadata are in `CITATION.cff`. The code is released under the MIT licence; see `LICENSE`.
