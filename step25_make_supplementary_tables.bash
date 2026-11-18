#! /bin/bash
# Supplementary Tables from the step23 statistics workbook (v13 output) and the
# primary step18 workbook (3 tiers: PVS, NearPVS, Fake_m6). Run from this folder.
STATS_DIR="$(python3 pipeline_config.py results_dir)/plots_and_stats/tissues-PVS_VOL-DIFF_nifti-based_and_SPM12-eTIV-based/conservative_PVS_segmentation_v13/stats_families_AEF/based_on_basic_info_created_by_Evgenios_updated_on_2026-08-28_PVS_segmentation_xlsx"
python step25_make_supplementary_tables.py --stats "$STATS_DIR/stats_all_corrected_v13.xlsx" --mixed-effects "$(python3 pipeline_config.py workbooks_dir)/MixedEffects_v12_3TIERS_cov_10k.xlsx" --out ./supplementary_tables
