function step21b_run_spm12_TIV_batch(t1_dir, out_csv)
% step21b_run_spm12_TIV_batch
% - Runs SPM12 segmentation on each T1w
% - Computes TIV from modulated warped tissue maps:
%     TIV = sum(mwc1 + mwc2 + mwc3) * voxel_volume
% - Optionally restricts computation using an MNI intracranial mask
%
% Usage example:
%   step21b_run_spm12_TIV_batch('/path/to/T1w_folder', '/path/to/TIV_results.csv')

    % --------- SPM12 path (spm12_dir in config.yaml) ----------
    spm_dir = pipeline_config_value('spm12_dir');
    addpath(spm_dir);
    spm('Defaults','fMRI');
    spm_jobman('initcfg');

    if nargin < 2
        error('Need t1_dir and out_csv');
    end
    if ~isfolder(t1_dir)
        error('t1_dir not found: %s', t1_dir);
    end

    % Find T1w NIfTIs
    t1s = dir(fullfile(t1_dir, '**', 'Segmentation','*T1wbrain_restore.nii.gz'));
    masks = dir(fullfile(t1_dir, '**', 'Segmentation','*_ICV_synthstrip.nii.gz'));
    if isempty(t1s)
        error('T1w NIFTI files not found...');
    end

    if isempty(masks)
        error('ICV masks NIFTI files not found...');
    end

    % ---------- Create (or reuse) an intracranial mask in MNI space ----------
    % This is a pragmatic approximation. Malone et al. used a manually-corrected ICV mask.
    % You can replace this mask with your own if you have one.
    %mask_path = fullfile(spm_dir, 'tpm', 'ICV_mask_from_TPM.nii');
    %if ~exist(mask_path)
    %    fprintf('Creating ICV mask: %s\n', mask_path);
    %    make_icv_mask_from_spm_tpm(spm_dir, mask_path);
    %end

    % Prepare output
    fid = fopen(out_csv, 'w');
    if fid < 0
        error('Cannot write: %s', out_csv);
    end
    fprintf(fid, 'subject,TIV_ml,GM_ml,WM_ml,CSF_ml,TIV_masked_ml,GM_masked_ml,WM_masked_ml,CSF_masked_ml,T1_path\n');

    % Process each file
    for i = 1:numel(t1s)
        t1_path = fullfile(t1s(i).folder, t1s(i).name);
        gunzip(t1_path);
        t1_path = erase(t1_path, '.gz');  % now ends with .nii
        mask_path = fullfile(masks(i).folder, masks(i).name);
        gunzip(mask_path);
        mask_path = erase(mask_path, '.gz');  % now ends with .nii
        
        [~, subjname, ~] = fileparts(t1_path);

        fprintf('\n[%d/%d] Segmenting: %s\n', i, numel(t1s), t1_path);

        % --- Run segmentation (New Segment) ---
        matlabbatch = {};
        %matlabbatch{1}.spm.spatial.preproc.channel.vols = {t1_path};
        matlabbatch{1}.spm.spatial.preproc.channel.vols = { [t1_path ',1'] };
        matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001;
        matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1]; % save bias-corrected + bias field

        TPM = { fullfile(spm_dir,'tpm','TPM.nii') };

        % Tissue classes 1..6 as per SPM12 defaults; we write modulated+warped for GM/WM/CSF
        for tc = 1:6
            matlabbatch{1}.spm.spatial.preproc.tissue(tc).tpm = {sprintf('%s,%d', TPM{1}, tc)};
            %matlabbatch{1}.spm.spatial.preproc.tissue(tc).ngaus = [1 1 2 3 4 2](tc);
            ngaus_vals = [1 1 2 3 4 2];
            matlabbatch{1}.spm.spatial.preproc.tissue(tc).ngaus = ngaus_vals(tc);
            if tc <= 3
                matlabbatch{1}.spm.spatial.preproc.tissue(tc).native = [1 0];   % c1/c2/c3
                matlabbatch{1}.spm.spatial.preproc.tissue(tc).warped = [1 1];   % wc* and mwc*
            else
                matlabbatch{1}.spm.spatial.preproc.tissue(tc).native = [0 0];
                matlabbatch{1}.spm.spatial.preproc.tissue(tc).warped = [0 0];
            end
        end

        matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
        matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
        matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
        matlabbatch{1}.spm.spatial.preproc.warp.samp = 3;
        matlabbatch{1}.spm.spatial.preproc.warp.write = [1 1]; % y_ and iy_

        spm_jobman('run', matlabbatch);

        % --- Locate outputs in the same folder as T1 ---
        % mwc1/mwc2/mwc3 are modulated, warped tissues in MNI space
        mwc1 = fullfile(t1s(i).folder, ['mwc1' subjname '.nii']);
        mwc2 = fullfile(t1s(i).folder, ['mwc2' subjname '.nii']);
        mwc3 = fullfile(t1s(i).folder, ['mwc3' subjname '.nii']);

        if ~(isfile(mwc1) && isfile(mwc2) && isfile(mwc3))
            warning('Missing mwc files for %s (expected %s etc). Skipping.', t1_path, mwc1);
            continue;
        end

        % --- Compute volumes (whole image) ---
        ref = mwc1; % reference grid
        mask_resliced = fullfile(t1s(i).folder, 'mask_in_mwc1_space.nii');
        mask_path = reslice_mask_to_ref(mask_path, ref, mask_resliced);
        [GM_ml, WM_ml, CSF_ml] = tissue_volume_ml(mwc1, mwc2, mwc3, '');
        TIV_ml = GM_ml + WM_ml + CSF_ml;

        % --- Compute volumes within ICV mask (recommended) ---
        [GMm_ml, WMm_ml, CSFm_ml] = tissue_volume_ml(mwc1, mwc2, mwc3, mask_path);
        TIV_masked_ml = GMm_ml + WMm_ml + CSFm_ml;

        fprintf('  TIV (ml): %.2f | Masked TIV (ml): %.2f\n', TIV_ml, TIV_masked_ml);

        fprintf(fid, '%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%s\n', ...
            subjname, TIV_ml, GM_ml, WM_ml, CSF_ml, TIV_masked_ml, GMm_ml, WMm_ml, CSFm_ml, t1_path);
    end

    fclose(fid);
    fprintf('\nDone. Wrote: %s\n', out_csv);
end

function make_icv_mask_from_spm_tpm(spm_dir, out_mask)
% Creates a simple intracranial mask in MNI space from SPM TPM priors.
% This is NOT the manually-corrected mask used in Malone et al.; it's a practical default.

    tpm = fullfile(spm_dir, 'tpm', 'TPM.nii');
    V1 = spm_vol([tpm ',1']); % GM prior
    V2 = spm_vol([tpm ',2']); % WM prior
    V3 = spm_vol([tpm ',3']); % CSF prior

    Y1 = spm_read_vols(V1);
    Y2 = spm_read_vols(V2);
    Y3 = spm_read_vols(V3);

    % Conservative-ish mask: where GM+WM+CSF prior > 0.2
    M = (Y1 + Y2 + Y3) > 0.2;

    % Basic cleanup
    M = spm_bwlabel(M, 26) > 0; % keep labeled components
    % Keep largest component
    L = bwlabeln(M, 26);
    if max(L(:)) > 0
        counts = accumarray(L(L>0), 1);
        [~, idx] = max(counts);
        M = (L == idx);
    end

    Vout = V1;
    Vout.fname = out_mask;
    Vout.dt = [spm_type('uint8') 0];
    spm_write_vol(Vout, uint8(M));
end

function [GM_ml, WM_ml, CSF_ml] = tissue_volume_ml(mwc1, mwc2, mwc3, mask_path)
% Volumes in ml from modulated warped tissue maps.
% If mask_path is provided, integrates only within mask.

    V1 = spm_vol(mwc1); Y1 = spm_read_vols(V1);
    V2 = spm_vol(mwc2); Y2 = spm_read_vols(V2);
    V3 = spm_vol(mwc3); Y3 = spm_read_vols(V3);

    if ~isempty(mask_path)
        Vm = spm_vol(mask_path);
        M = spm_read_vols(Vm) > 0;  % binary mask expected
        Y1(~M) = 0; Y2(~M) = 0; Y3(~M) = 0;
    end

    % voxel volume (mm^3)
    vx_mm3 = abs(det(V1.mat));
    % tissue maps are "amount of tissue" per voxel after modulation.
    GM_ml  = sum(Y1(:)) * vx_mm3 / 1000.0;
    WM_ml  = sum(Y2(:)) * vx_mm3 / 1000.0;
    CSF_ml = sum(Y3(:)) * vx_mm3 / 1000.0;
end

function mask_r = reslice_mask_to_ref(mask_path, ref_path, out_path)
    % Reslice mask to the grid of ref_path using nearest-neighbor
    flags = struct('interp',0,'wrap',[0 0 0],'mask',0,'which',1,'mean',0);
    spm_reslice(char(ref_path, mask_path), flags);  % writes r<maskname>

    [p,n,e] = fileparts(mask_path);
    if strcmp(e,'.gz'); error('Gunzip mask first'); end

    rmask = fullfile(p, ['r' n e]); % SPM output
    movefile(rmask, out_path);
    mask_r = out_path;
end