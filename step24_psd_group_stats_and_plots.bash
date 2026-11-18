#! /bin/bash

python3 step24_psd_group_stats_and_plots.py --psd-csv "$(python3 pipeline_config.py psd_csv)" --study-xlsx "$(python3 pipeline_config.py study_xlsx)" --out-dir . --apply
