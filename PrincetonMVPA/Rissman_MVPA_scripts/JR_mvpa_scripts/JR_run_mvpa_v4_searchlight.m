function [subj results] = JR_run_mvpa_v4_searchlight(subj_id, mvpa_workspace, cond1,cond2, balanced_bins)

% runs the mvpa analysis, start to finish
% example usage:  [subj results] = JR_run_mvpa_v4('s106', 's106_merged_AAL_ROIs_8mm_smoothing.mat')

tic % start the clock

%%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flags.save_workspace = 1; % 1 = yes please, 0 = no thanks

% unless a previously created workspace is specified, load in the preprocessed data
if ~exist('mvpa_workspace')

    exp_name = 'PAST';
    roi_file = '/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/merged_AAL_ROIs.nii';
    roi_name = 'merged_AAL_ROIs'
    num_TP_per_run = 203;

    % load user-created filename and onsets lists into workspace
    load('raw_filenames_wa.mat') %loads predefined cell array called raw_filenames into memory
    num_runs = length(raw_filenames)/num_TP_per_run; %calculate the number of runs (allows script to work flexibly for subjects with missing runs)
    load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']);
    load('/Users/Jesse/fMRI/data/PAST/fMRI/vol_info.mat'); %get functional data resolution info for spm .img writing
    [subj] = JR_mvpa_load_and_preprocess_raw_data(subj_id, exp_name, roi_name, roi_file, raw_filenames, num_runs, num_TP_per_run);

    if flags.save_workspace == 1
        save_cmd = ['save ' subj_id '_' roi_name '_unsmoothed.mat'];
        eval(save_cmd);
    end
else
    eval(['load ' mvpa_workspace])
    load('/Users/Jesse/fMRI/data/PAST/fMRI/vol_info.mat'); %get functional data resolution info for spm .img writing
end

% Set flags (% 1 = yes please, 0 = no thanks)
flags.equate_number_of_old_new_trials_per_subjective_bin = balanced_bins; % balance the number of trials for each subjective response
flags.equate_number_of_trials_in_cond_1_and_2 = 1; % balance the number of trials in conditions 1 & 2 (this is a good thing and prevents biasing the classifier to pi
flags.plot_mean_timecourses = 0;
flags.plot_ROC_curve = 0;
flags.display_performance_breakdown = 0;
flags.generate_importance_maps = 0;
flags.remove_outlier_trials = 3;


% specify which conditions to use for classification (must correspond to the names of conditions specified below)
condnames =  {cond1,cond2};

TRs_to_average_over = [3 4]; %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier

num_results_iter = 1; % number of times to run the cross validation process

anova_p_thresh = 1;  %p threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)


searchlight_path = '/Users/Jesse/fMRI/data/PAST/fMRI/mvpa_results/searchlight_maps_unsmoothed_round5';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract info about conditions from onsets file
num_conds = size(onsets,2);
all_regs = zeros(num_conds,num_runs*num_TP_per_run); % initialize regs matrix as conditions x timepoints

for cond = 1: num_conds-1 %(exclude last condition ("no_response" trials)
    for trial = 1: length(onsets{cond})
        time_idx = onsets{cond}(trial)/2+1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
        all_regs(cond,time_idx) = 1;
    end
end

% SPECIAL FIX FOR s104 and s105 because reg #6 (NEW_recollect) is a fake
% placeholder trial
if ~isempty(find(strcmp(subj_id,{'s104','s105','s108','s113'})))
    all_regs(6,:)=0;
end

% Artificially balance the number of trials given each subjective memory response
if flags.equate_number_of_old_new_trials_per_subjective_bin == 1
    for j = 1:5
        OLD_trials = find(all_regs(j,:));
        NEW_trials = find(all_regs(j+5,:));

        num_OLD = length(OLD_trials);
        num_NEW = length(NEW_trials);

        if num_OLD > num_NEW
            rand_array = rand(1,num_OLD);
            [sorted inds]= sort(rand_array);
            trials_to_cut = OLD_trials(inds(1:num_OLD-num_NEW));
            all_regs(j,trials_to_cut)=0;
        elseif num_OLD < num_NEW
            rand_array = rand(1,num_NEW);
            [sorted inds]= sort(rand_array);
            trials_to_cut = NEW_trials(inds(1:num_NEW-num_OLD));
            all_regs(j+5,trials_to_cut)=0;
        end
    end
end


% define conditions of interest
% specify names 'OLD_recollect'    'OLD_hc_old'    'OLD_lc_old'    'OLD_lc_new'    'OLD_hc_new'    'NEW_recollect' 'NEW_hc_old'    'NEW_lc_old'    'NEW_lc_new'    'NEW_hc_new'    'no_resp'

Objective_old = sum(all_regs(1:5,:));
Objective_new = sum(all_regs(6:10,:));

Subjective_old = sum(all_regs([1 2 3 6 7 8],:));
Subjective_new = sum(all_regs([4 5 9 10],:));

Subjective_old_HC_only = sum(all_regs([1 2 6 7],:));
Subjective_new_HC_only = sum(all_regs([5 10],:));

Hits = sum(all_regs([1 2 3],:));
Misses = sum(all_regs([4 5],:));
CRs = sum(all_regs([9 10],:));
FAs = sum(all_regs(6:8,:));

R_hits = all_regs(1,:);
HC_hits = all_regs(2,:);
LC_hits = all_regs(3,:);
LC_misses = all_regs(4,:);
HC_misses = all_regs(5,:);
R_FAs = all_regs(6,:);
HC_FAs = all_regs(7,:);
LC_FAs = all_regs(8,:);
LC_CRs = all_regs(9,:);
HC_CRs = all_regs(10,:);

R_and_HC_hits = sum(all_regs([1 2],:));

%no_resp = all_regs(11,:); %excluded from analysis

%assign conditions to train/test classifier on
regs_of_interest = [];
eval(['regs_of_interest(1,:) = ' condnames{1} ';'])
eval(['regs_of_interest(2,:) = ' condnames{2} ';'])


if flags.equate_number_of_trials_in_cond_1_and_2 == 1

    cond1_trials = find(regs_of_interest(1,:));
    cond2_trials = find(regs_of_interest(2,:));
    num_cond1 = length(cond1_trials);
    num_cond2 = length(cond2_trials);

    if num_cond1 > num_cond2
        rand_array = rand(1,num_cond1);
        [sorted inds]= sort(rand_array);
        trials_to_cut = cond1_trials(inds(1:num_cond1-num_cond2));
        regs_of_interest(1,trials_to_cut) = 0;
        display([num2str(length(trials_to_cut)) ' trials cut from ' condnames{1}]);
    elseif num_cond1 < num_cond2
        rand_array = rand(1,num_cond2);
        [sorted inds]= sort(rand_array);
        trials_to_cut = cond2_trials(inds(1:num_cond2-num_cond1));
        regs_of_interest(2,trials_to_cut) = 0;
        display([num2str(length(trials_to_cut)) ' trials cut from ' condnames{2}]);
    else
        display('Trial numbers are already balanced');
    end
end

display([num2str(count(regs_of_interest(1,:)==1)) ' trials in condition ' condnames{1}])
display([num2str(count(regs_of_interest(2,:)==1)) ' trials in condition ' condnames{2}])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select TRs of interest (to correspond with peak post-stim BOLD response)

all_trials = sum(all_regs,1); % vector of all trial   (patterns{4} contains fully preprocessed data)

data_by_TR(1,:,:) = subj.patterns{end}.mat(:,find(all_trials)+0); % 1st TR (0-2 sec)
data_by_TR(2,:,:) = subj.patterns{end}.mat(:,find(all_trials)+1); % 2nd TR (2-4 sec)
data_by_TR(3,:,:) = subj.patterns{end}.mat(:,find(all_trials)+2); % 3rd TR (4-6 sec)
data_by_TR(4,:,:) = subj.patterns{end}.mat(:,find(all_trials)+3); % 4th TR (6-8 sec)
data_by_TR(5,:,:) = subj.patterns{end}.mat(:,find(all_trials)+4); % 5th TR (8-10 sec)

if length(TRs_to_average_over)==1
    temporally_condensed_data = squeeze(data_by_TR(TRs_to_average_over,:,:)); %remove singleton dimension
else
    temporally_condensed_data = squeeze(mean(data_by_TR(TRs_to_average_over,:,:),1));
end

clear data_by_TR; %clean up matlab workspace to save memory

if flags.plot_mean_timecourses == 1 % plot mean timecourses for whole ROI (useful check of data quality & overall activity pattern)

    c1_ind = find(regs_of_interest(1,:));
    for t = 1: length(c1_ind)
        c1_ts(t,:) = mean(subj.patterns{end}.mat(:,c1_ind(t):c1_ind(t)+6),1);
    end

    c2_ind = find(regs_of_interest(2,:));
    for t = 1: length(c2_ind)
        c2_ts(t,:) = mean(subj.patterns{end}.mat(:,c2_ind(t):c2_ind(t)+6),1);
    end

    figure;
    plot(mean(c1_ts,1),'b');
    hold on;
    plot(mean(c2_ts,1),'r');
end


% condense regs by removing zeros
trial_counter = 1;
for i = 1: size(all_regs,2)
    if ~isempty(find(all_regs(:,i))) % if not a rest timepoint
        condensed_regs_of_interest(:,trial_counter) = regs_of_interest(:,i);
        condensed_regs_all(:,trial_counter) = all_regs(:,i);
        condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
        trial_counter = trial_counter + 1;
    end
end

% initialize regressors object
subj = init_object(subj,'regressors','conds');
subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
subj = set_objfield(subj,'regressors','conds','condnames',condnames);

% add new condensed activation pattern
subj = duplicate_object(subj,'pattern','spiral_d_hp_z','spiral_d_hp_z_condensed');
subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed',temporally_condensed_data);

zhist = sprintf('Pattern ''%s'' created by JR custom code','spiral_d_hp_z_condensed');
subj = add_history(subj,'pattern','spiral_d_hp_z_condensed',zhist,true);

% clean up workspace
subj = remove_mat(subj,'pattern','spiral_d_hp_z');
clear mean_data;

% update run vector to condensed format
subj.selectors{1}.mat = condensed_runs;
subj.selectors{1}.matsize = size(condensed_runs);

% "activate" only those trials of interest (from regs_of_interest) before
% creating cross-validation indices

active_trials = find(sum(condensed_regs_of_interest));


if flags.remove_outlier_trials ~= 0
    % remove outlier trials (timepoints)
    mean_across_voxels = mean(subj.patterns{end}.mat(:,active_trials),1);
    z_mean_across_voxels = zscore(mean_across_voxels);
    upper_outliers = find(z_mean_across_voxels> flags.remove_outlier_trials);
    lower_outliers = find(z_mean_across_voxels< -1 * flags.remove_outlier_trials);
    all_outliers = union(upper_outliers,lower_outliers)
    active_trials(all_outliers) = [];
end

actives_selector = zeros(1,size(condensed_regs_all,2)); % intialize vector of all zeros
actives_selector(active_trials) = 1; % remove all non-"regs_of_interst" trials (set to one)
subj = init_object(subj,'selector','conditions_of_interest'); %initialize selector object
subj = set_mat(subj,'selector','conditions_of_interest',actives_selector);
subj = JR_create_xvalid_indices_searchlight(subj,'runs','actives_selname','conditions_of_interest');

% zscore temporally-condensed data; active trials only (second round of z-scoring)
subj.patterns{end}.mat(:,active_trials) = zscore(subj.patterns{end}.mat(:,active_trials)')';


searchlight_radius = 3;
subj.adj_sphere = create_adj_list(subj,roi_name,'radius',searchlight_radius);

%specify training/testing functions for within-sphere classification
class_args.train_funct_name = 'train_gnb';
class_args.test_funct_name = 'test_gnb';

scratch.class_args = class_args;
scratch.perfmet_funct = 'perfmet_maxclass';
scratch.perfmet_args = struct([]);

statmap_srch_arg.adj_list = subj.adj_sphere;
statmap_srch_arg.obj_funct = 'statmap_classify';
statmap_srch_arg.scratch = scratch;

subj = feature_select( ...
    subj, ...
    'spiral_d_hp_z_condensed', ... % data
    'conds', ... % binary regs (for GNB)
    'runs_xval', ... % selector
    'statmap_funct','statmap_searchlight', ... % function
    'statmap_arg',statmap_srch_arg, ...
    'new_map_patname','spiral_d_hp_z_condensed_srch', ...
    'thresh',[]);

% create masks from the statmaps, by picking the best N values (e.g., 500) in each
% statmap to use for classification on the remaining run of test trials
% NOTE: keeps top N voxels; not top N spheres
% subj = create_sorted_mask( ...
%     subj,'spiral_d_hp_z_condensed_srch', ...
%     'spiral_d_hp_z_condensed_srch_500',500, ...
%     'descending',true);


% % run classifier (No hidden layer NN Toolbox backprop algorithm)
% class_args.train_funct_name = 'train_bp';
% class_args.test_funct_name = 'test_bp';
% class_args.nHidden = 0;
%
% [subj results] = cross_validation(subj,'spiral_d_hp_z_condensed','conds','runs_xval_sl','spiral_d_hp_z_condensed_srch_500',class_args);

%% WRITE OUT MEAN SEARCHLIGHT MAP TO .IMG FILE

% average searchlight performance maps across runs
for r=1:num_runs
    sl_voxel_values(:,r)=subj.patterns{end-num_runs+r}.mat;
end
sl_mean_voxel_values = mean(sl_voxel_values,2);

vol_info.fname = [searchlight_path '/' subj_id '_' condnames{1} '_vs_' condnames{2} '_3vox_radius_searchlight.img'];

sl_map = zeros(vol_info.dim);
included_voxels = find(subj.masks{1}.mat);
sl_map(included_voxels) = sl_mean_voxel_values.*100; %multiply by 100 to improve scaling for visualization

spm_write_vol(vol_info, sl_map);



time2finish = toc/60;
display(['Finished in ' num2str(time2finish) ' minutes']);