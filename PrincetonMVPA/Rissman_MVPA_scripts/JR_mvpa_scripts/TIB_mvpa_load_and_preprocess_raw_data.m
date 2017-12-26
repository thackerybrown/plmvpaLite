function [subj] = TIB_mvpa_load_and_preprocess_raw_data(subj_id, exp_name, roi_name, roi_file, raw_filenames, num_runs, num_TP_per_run)

    % initialize subj structure
    subj = init_subj(exp_name,subj_id);

    % load mask file
    subj = load_spm_mask(subj,roi_name,roi_file);
    
    % load functional data
    subj = load_analyze_pattern_PMVPA(subj,'epi',roi_name,raw_filenames,'single',true); %use single precision format to save RAM

    % move pattern to hard drive to save RAM (optional)
    %subj = move_pattern_to_hd(subj, 'epi');

    % make runs vector
    subj = init_object(subj,'selector','runs');

    %%builds run index numbers into something called "runs" such that each
    %%trial of run 2 has a value of 2, run 3 a value of 3, etc etc.
    trial_idx = 0;
    %variable run length version (assumes num_TP_per_run is a vector of run lengths):
    for r = 1:num_runs
        runs(trial_idx+1:(num_TP_per_run(r)+trial_idx)) = r;%for each run, assign the run # to each element of a list num_TP long for that run
        trial_idx = trial_idx + num_TP_per_run(r);
    end
    
    %fixed run length version:
%     for r = 1:num_runs
%         runs(trial_idx:num_TP_per_run*r) = r;
%         trial_idx = trial_idx + num_TP_per_run;
%     end
    
    %% fix runs labels for subject 01 from COUNTERMEASURES study
%     if strcmp(subj_id,'CM001p')
%         load('eleven_run_labels.mat');  %includes s01's Run 5b
%     end

    subj = set_mat(subj,'selector','runs',runs);

    % detrend the timeseries data
    subj = detrend_runs(subj,'epi','runs');  % not in mvpa tutorial, but seems important to do

    % move pattern to hard drive to save RAM (optional)
    %subj = move_pattern_to_hd(subj, 'epi_d');
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','epi');
    
    % high-pass filter the timeseries data
    subj = hpfilter_runs(subj,'epi_d','runs',128,2); % remove frequencies below .01 Hz (adjust as desired)

     % clean up workspace
    subj = remove_mat(subj,'pattern','epi_d');
    
    % move pattern to hard drive to save RAM (optional)
    %subj = move_pattern_to_hd(subj, 'epi_d_hp');
    
    % zscore the data from each run
    subj = zscore_runs(subj,'epi_d_hp','runs'); % gives each voxel a mean of 1 and variance of 0 across all timepoints of each run
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','epi_d_hp');
    
    %save final pattern in single precision form (8 sig figs) to save RAM and HD space    
    subj.patterns{end}.mat = single(subj.patterns{end}.mat);