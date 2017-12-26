function [] = JR_batch_generate_mvpa_workspace_mats(subj_array)

roi_name = 'merged_AAL_ROIs_FIXED_HOLES';  % name of anatomical mask
data_imgs_to_use = 'raw_filenames_wa.mat'; % .mat file containing names of all functional images
expt_dir = '/Users/Jesse/fMRI/data/PAST/fMRI';


for b=subj_array
    tic
    if b<10
        subj_id=strcat('s10', num2str(b))
    else
        subj_id=strcat('s1', num2str(b))
    end
    
    
    mvpa_dir = [expt_dir '/' subj_id '/mvpa'];
    [subj num_runs num_TP_per_run]= JR_generate_mvpa_workspace_mat_file(subj_id, roi_name, data_imgs_to_use, mvpa_dir);
    clear subj;
end