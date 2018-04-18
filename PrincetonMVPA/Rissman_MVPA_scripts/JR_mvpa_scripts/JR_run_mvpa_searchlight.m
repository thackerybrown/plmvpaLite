function [subj results] = JR_run_mvpa_searchlight(mvpa_workspace)

tic

if ~exist('mvpa_workspace')

    %%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subj_id = 's102';
    exp_name = 'PAST';
    roi_file = 'omnibus_F_10_no_motor.img';
    %anova_p_thresh = .005;  %p threshold for feature selection ANOVA
    generate_importance_maps = 1; % 1 = yes please, 2 = no thanks
    equate_number_of_trials_per_bin = 0; % 1 = yes please, 2 = no thanks
    equate_subjective_status = 0; % 1 = yes please, 2 = no thanks
    vol_info = spm_vol(['/Users/Jesse/fMRI/data/PAST/fMRI/' subj_id '/anatomical/' subj_id '_mean_func.img']); %get functional data resolution info for spm .img writing

    num_runs = 10;
    num_TP_per_run = 203;
    total_TP = num_runs * num_TP_per_run;

    % load user-created filename and onsets lists into workspace
    load('raw_filename_list_SRA.mat')
    load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']);

    num_conds = size(onsets,2);
    all_regs = zeros(num_conds,total_TP); % initialize regs matrix as conditions x timepoints

    for cond = 1: num_conds-1 %(exclude last condition ("no_response" trials)
        for trial = 1: length(onsets{cond})
            time_idx = onsets{cond}(trial)/2+1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
            all_regs(cond,time_idx) = 1;
        end
    end

    
    % SPECIAL FIX FOR s104 and s105 because reg #6 (NEW_recollect) is a fake
    % placeholder trial
    if ~isempty(find(strcmp(subj_id,{'s104','s105'})))
        all_regs(6,:)=0;
    end

    
    if equate_number_of_trials_per_bin == 1
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
    
    if equate_subjective_status == 1
        sOld_trials = find(sum(all_regs([1:3,6:8],:)));
        sNew_trials = find(sum(all_regs([4:5,9:10],:)));

        num_sOld = length(sOld_trials);
        num_sNew = length(sNew_trials);

        if num_sOld > num_sNew
            rand_array = rand(1,num_sOld);
            [sorted inds]= sort(rand_array);
            trials_to_cut = sOld_trials(inds(1:num_sOld-num_sNew));
            all_regs([1:3,6:8],trials_to_cut)=0;
        elseif num_sOld < num_sNew
            rand_array = rand(1,num_sNew);
            [sorted inds]= sort(rand_array);
            trials_to_cut = sNew_trials(inds(1:num_sNew-num_sOld));
            all_regs([4:6,9:10],trials_to_cut)=0;
        end
    end
                        
    

    % condnames = names; %specify names 'OLD_recollect'    'OLD_hc_old'    'OLD_lc_old'    'OLD_lc_new'    'OLD_hc_new'    'NEW_recollect' 'NEW_hc_old'    'NEW_lc_old'    'NEW_lc_new'    'NEW_hc_new'    'no_resp'

    % define conditions of interest
    objective_old = sum(all_regs(1:5,:)); 
    objective_new = sum(all_regs(6:10,:)); 
    
    subjective_old_all = sum(all_regs([1 2 3 6 7 8],:)); 
    subjective_new_all = sum(all_regs([4 5 9 10],:));
    
    subjective_old_hc_only = sum(all_regs([1 2 6 7],:));
    subjective_new_hc_only = sum(all_regs([5 10],:));
    
    correct_recollect = all_regs(1,:);
    correct_hc_old = all_regs(2,:);
    
    hits = sum(all_regs([1 2 3],:));
    crs = sum(all_regs([9 10],:));
    
    r_and_hc_hits = sum(all_regs([1 2],:));
    lc_hits = all_regs(3,:); 
    lc_crs = all_regs(9,:);
    lc_misses = all_regs(4,:);
    fas = sum(all_regs(6:8,:)); 
    
    %no_resp = all_regs(11,:); %excluded from analysis
    
     
    %choose conditions to train/test classifier on
    regs_of_interest = [];
    regs_of_interest(1,:) = objective_old;
    regs_of_interest(2,:) = objective_new;
   
    condnames =  {'objective_old','objective_new'}; %{'old','new'};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % initialize subj structure
    subj = init_subj(exp_name,subj_id);

    % load mask file
    subj = load_analyze_mask(subj,roi_file,roi_file);

    
    % load functional data
    subj = load_analyze_pattern(subj,'spiral',roi_file,raw_filenames);

    % move pattern to hard drive
    % subj = move_pattern_to_hd(subj, 'spiral');

    % make runs vector
    subj = init_object(subj,'selector','runs');

    trial_idx = 1;
    for r = 1:num_runs

        runs(trial_idx:203*r) = r;
        trial_idx = trial_idx + num_TP_per_run;
    end

    subj = set_mat(subj,'selector','runs',runs);

    % detrend the timeseries data
    subj = detrend_runs(subj,'spiral','runs');  % not in mvpa tutorial, but seems important to do

    % move pattern to hard drive
    %subj = move_pattern_to_hd(subj, 'spiral_d');
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','spiral');
    
    % high-pass filter the timeseries data
    subj = hpfilter_runs(subj,'spiral_d','runs',100,2); % remove frequencies below .01 Hz

     % clean up workspace
    subj = remove_mat(subj,'pattern','spiral_d');
    
    % move pattern to hard drive
    %subj = move_pattern_to_hd(subj, 'spiral_d_hp');
    
    % zscore the data from each run
    subj = zscore_runs(subj,'spiral_d_hp','runs'); % Z-score the data
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','spiral_d_hp');
    
    %save_cmd = ['save ' subj_id '_spiral_d_hp_z_' roi_file(1:end-4) '_' mvpa_datetime '.mat'];
    %eval(save_cmd); 

else
    load(mvpa_workspace); %load in pre-saved workspace files and continue
end

% select TRs of interest (to correspond with peak post-stim BOLD response)

all_trials = sum(all_regs,1); % vector of all trial   (patterns{4} contains fully preprocessed data)

TR0 = subj.patterns{4}.mat(:,find(all_trials)-1); % -1 TR (-2-0 sec)
TR1 = subj.patterns{4}.mat(:,find(all_trials)+0); % 1st TR (0-2 sec)
TR2 = subj.patterns{4}.mat(:,find(all_trials)+1); % 2nd TR (2-4 sec)
TR3 = subj.patterns{4}.mat(:,find(all_trials)+2); % 3rd TR (4-6 sec)
TR4 = subj.patterns{4}.mat(:,find(all_trials)+3); % 4th TR (6-8 sec)
TR5 = subj.patterns{4}.mat(:,find(all_trials)+4); % 5th TR (8-10 sec)

mean_data = (TR3+TR4)/2; %take uniform average of peak signal TRs

% BASELINE CORRECT MEAN DATA (optional)
%mean_baseline = TR1;
%mean_data = mean_data - mean_baseline;



clear TR0 TR1 TR2 TR3 TR4 TR5; %clean up matlab workspace to save memory

% plot mean timecourses for whole ROI (useful check of data quality & overall activity pattern)

new_ind = find(objective_new);
for t = 1: length(new_ind)
     new_ts(t,:) = mean(subj.patterns{4}.mat(:,new_ind(t):new_ind(t)+6),1);
end

old_ind = find(objective_old);
for t = 1: length(old_ind)
     old_ts(t,:) = mean(subj.patterns{4}.mat(:,old_ind(t):old_ind(t)+6),1);
end

figure;
plot(mean(old_ts,1),'b');
hold on;
plot(mean(new_ts,1),'r');


% condense regs by removing zeros

trial_counter = 1;
for i = 1: size(all_regs,2)
    if ~isempty(find(all_regs(:,i))) % if not a rest timepoint
        condensed_regs_of_interest(:,trial_counter) = regs_of_interest(:,i);
        condensed_regs_all(:,trial_counter) = all_regs(:,i); 
        condensed_runs(1,trial_counter) = runs(i);
        trial_counter = trial_counter + 1;
    end
end

% initialize regressors object
subj = init_object(subj,'regressors','conds');
subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
subj = set_objfield(subj,'regressors','conds','condnames',condnames);

% add new condensed activation pattern
subj = duplicate_object(subj,'pattern','spiral_d_hp_z','spiral_d_hp_z_condensed');
subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed',mean_data);

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

%temp_sel = ones(1,size(condensed_regs_all,2)); % intialize vector of all ones
%temp_sel(find(sum(condensed_regs_of_interest)==0)) = 0; % remove all non-"regs_of_interst" trials (set to zero)
%subj = init_object(subj,'selector','conditions_of_interest'); %initialize selector object
%subj = set_mat(subj,'selector','conditions_of_interest',temp_sel);
%subj = create_xvalid_indices(subj,'runs','actives_selname','conditions_of_interest');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                       SPHERICAL SEARCHLIGHT                                %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Create group of cross-validation selectors. There will be one for each
% iteration of the cross-validation process. For each iteration, we're going to withhold one run for final-test data,
% we'll withhold one run for searchlight-generalization data, and we'll use the remainder for training data
runs = get_mat(subj,'selector','runs');
nRuns = max(runs);
nTimepoints = length(runs);

%  Assign runs to be used for training testing data (give value of 1)
% NOTE:  We'll initially give this value to every timepoint, and then we'll
% modify this matrix to include 2's and 3's (see below)

runs_xval_sl = ones(nRuns,nTimepoints);

%  Assign runs to be used for final testing data (give value of 2)
for r=1:nRuns
  cur_final_test_run = find(runs==r);
  runs_xval_sl(r, cur_final_test_run) = 2;
end

%  Assign runs to be used for searchlight generalization (give value of 3)
for r=1:nRuns
  cur_searchlight_gen_run = find(runs== nRuns-r+1);
  runs_xval_sl(r, cur_searchlight_gen_run) = 3;
end

for r=1:nRuns
  cur_name = sprintf('runs_xval_sl_%i',r);
  subj = initset_object(subj, 'selector', cur_name, ...
                        runs_xval_sl(r,:), ...
                        'group_name', 'runs_xval_sl' );
end

subj.adj_sphere = create_adj_list(subj,roi_file,'radius',3);

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
    'runs_xval_sl', ... % selector
    'statmap_funct','statmap_searchlight', ... % function
    'statmap_arg',statmap_srch_arg, ...
    'new_map_patname','spiral_d_hp_z_condensed_srch', ...
    'thresh',[]);

% create masks from the statmaps, by picking the best N values (e.g., 500) in each
% statmap to use for classification on the remaining run of test trials
% NOTE: keeps top N voxels; not top N spheres
subj = create_sorted_mask( ...
    subj,'spiral_d_hp_z_condensed_srch', ...
    'spiral_d_hp_z_condensed_srch_500',500, ...
    'descending',true);    


% run classifier (No hidden layer NN Toolbox backprop algorithm)
class_args.train_funct_name = 'train_bp';
class_args.test_funct_name = 'test_bp';
class_args.nHidden = 0;

[subj results] = cross_validation(subj,'spiral_d_hp_z_condensed','conds','runs_xval_sl','spiral_d_hp_z_condensed_srch_500',class_args);

%% WRITE OUT MEAN SEARCHLIGHT MAP TO .IMG FILE

% average searchlight performance maps across runs
for r=1:nRuns 
    sl_voxel_values(:,r)=subj.patterns{end-nRuns+i}.mat;
end
sl_mean_voxel_values = mean(sl_voxel_values,2);

vol_info.fname = [subj_id '_' condnames{1} '_vs_' condnames{2} '_searchlight_map.img'];

sl_map = zeros(vol_info.dim);
included_voxels = find(subj.masks{1}.mat);
sl_map(included_voxels) = sl_mean_voxel_values.*100; %multiply by 100 to improve scaling for visualization

spm_write_vol(vol_info, sl_map);






time2finish = toc/60;
disp(['time2finish = ' num2str(time2finish) ' minutes']);

display('all done');