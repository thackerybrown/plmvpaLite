function [subj num_runs num_TP_per_run] = JR_generate_mvpa_workspace_mat_file_no_detrend(subj_id, roi_name, data_imgs_to_use, mvpa_dir);

exp_name = 'PAST';
roi_file = ['/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/' roi_name '.nii'];
num_TP_per_run = 203;

% load user-created filename and onsets lists into workspace

load([mvpa_dir '/' data_imgs_to_use]); %loads predefined cell array called raw_filenames into memory
num_runs = length(raw_filenames)/num_TP_per_run; %calculate the number of runs (allows script to work flexibly for subjects with missing runs)
[subj] = JR_mvpa_load_and_preprocess_raw_data_no_detrend(subj_id, exp_name, roi_name, roi_file, raw_filenames, num_runs, num_TP_per_run);
% voxel_means = mean(subj.patterns{end}.mat,2);
% zero_voxel_inds = find(voxel_means==0);
% subj.patterns{end}.mat(zero_voxel_inds,:)=[]; % remove voxels with no data
% mask_vox = find(subj.masks{end}.mat);
% subj.masks{end}.mat(mask_vox(zero_voxel_inds))=0; % remove these voxels from mask as well
% display([num2str(length(zero_voxel_inds)) ' zero-value voxels removed from pattern and mask structs'])
% clear voxel_means zero_voxel_inds;
save_cmd = ['save ' mvpa_dir '/' subj_id '_' roi_name '_NO_DETREND_' data_imgs_to_use(15:end)];
eval(save_cmd);
