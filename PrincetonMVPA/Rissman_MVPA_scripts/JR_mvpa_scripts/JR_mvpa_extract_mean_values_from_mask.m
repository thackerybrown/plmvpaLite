function []=JR_mvpa_extract_mean_values_from_mask

%subj_array = [1 2 4 6 7 9 12 13 14];
subj_array = [1 2 4 6 7 8 9 11 12 13 14];
%subj_array = [1 2 4:9 11:14];


topdir = pwd;
PAST_dir = '/Users/Jesse/fMRI/data/PAST/fMRI';

mvpa_results_dir = [PAST_dir '/mvpa_results/searchlight_maps_unsmoothed'];
mvpa_results_dir = '/Users/Jesse/fMRI/data/PAST/fMRI/mvpa_results/improved_importance_maps/HC_hits_vs_LC_hits_unbalanced_bins_TR_3_4';

%mvpa_results_dir = [PAST_dir '/mvpa_results/weight_maps/R_hits_vs_HC_hits_unbalanced_bins_TR_3_4'];
%mvpa_results_dir = [PAST_dir '/mvpa_results/weight_maps/HC_hits_vs_LC_hits_unbalanced_bins_TR_3_4'];

maskfile = ['/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/spm5_AAL_ROIs/rMNI_Hippocampus_R.img']
M = spm_vol(maskfile);
m_data = spm_read_vols(M);
m_voxels = find(m_data);


for i = 1:length(subj_array)
    subj_id = ['s1' prepend(num2str(subj_array(i)))];
    %datafile = [mvpa_results_dir '/s' subj_id '_R_hits_vs_HC_hits_3vox_radius_searchlight.img'];
    datafile = [mvpa_results_dir '/' subj_id '_HC_hits.img'];
    V=spm_vol(datafile);
    v_data = spm_read_vols(V);
    data_at_mask_voxels(i) = mean(v_data(m_voxels));
end

data_at_mask_voxels
mean(data_at_mask_voxels)

    
    