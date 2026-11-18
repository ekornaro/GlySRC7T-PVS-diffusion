#! /bin/bash
# Group statistics and figures (families A, E, F) as run for the paper.
# The two inclusion thresholds are explicit here; the script defaults are not used.
python step23_plot_all_markers_across_cohorts_and_ROIs.py --min-pvs-voxels 10 --max-diffusivity 3.0e-3 --also-plot-absolute
