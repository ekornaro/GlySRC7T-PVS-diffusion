% ============================================================
% Wrapper script to compute TIV for all subjects using SPM12
% Calls: step21b_run_spm12_TIV_batch.m
% ============================================================

clear; clc;

% ---------- SETTINGS (paths from config.yaml) ----------
% Path to scripts (this repository)
scripts_dir = fileparts(mfilename('fullpath'));
addpath(scripts_dir);

% Folder containing T1w NIfTI files (can be nested)
t1_root_dir = pipeline_config_value('pvs_root');
% ^ example: contains sub-XXX folders or flat .nii files

% Output CSV file
out_csv = pipeline_config_value('tiv_spm12_csv');

% ---------- SANITY CHECKS ----------
if ~isfolder(t1_root_dir)
    error('T1 root directory not found: %s', t1_root_dir);
end

fprintf('T1 root directory : %s\n', t1_root_dir);
fprintf('Output CSV        : %s\n', out_csv);

% ---------- RUN TIV PIPELINE ----------
step21b_run_spm12_TIV_batch(t1_root_dir, out_csv);

fprintf('\nAll done.\nResults saved to:\n%s\n', out_csv);
