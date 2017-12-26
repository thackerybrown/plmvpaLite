function []= MU4_run_mvpa(subj_array, condition1, condition2, nVox, penalty)

% Workflow script for running binary MVPA classification on Jesse Rissman's
% face recognition memory dataset. Requires Princeton MVPA toolbox
% functions to be located in your path.

% example usage:  LTC_run_mvpa({'LTC001'},'OLD','NEW',0,10)

% Inputs
% 
% subj_array: array of subject id #'s (e.g., [1 2 4 5 8]) to run analysis on; these #'s get converted to strings below (e.g., 's111', 's112', 's114', etc.)
% 
% condition1: name of condition1 (e.g., 'Hits')
% condition2: name of condition2 (e.g., 'CRs')


%%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for b=1:length(subj_array)  % run loop for each subject in subj array (subject ID #'s specified in short format)
    tic % record the start time
    
    subj_id = subj_array{b};
    
    expt_name = 'Circmaze';  
    expt_dir = '/Users/thackery/Work/Circmaze/cm03'; % top-level expt. directory where individual subject folders live
    mvpa_dir = [expt_dir '/' subj_id '/model_mux_er_MNI/mvpa_results_mux_er']; % mvpa directory within each subject's folder (will be created below if it doesn't exist)
    %spm_onsets_dir = [expt_dir '/' subj_id '/mvpa/']]; % directory where each subject's SPM-formatted onsets.mat files lives
    
    % change mask back to SEPT09_MVPA_MASK ?????
    roi_name = 'mux1mask_minusmotor';  % name of mask to be used for voxel selection (can be a small ROI, a whole-brain mask, or anywhere in between)
    roi_file = [expt_dir '/Masks/' roi_name '.img']; % specific path to mask file
    data_imgs_to_use = 'raw_filenames.mat'; % .mat file containing names of all functional images (must exist for each subject; can be created by running cellstr(SPM.xY.P) on subject's SPM.mat file)
    
    num_TP_per_run = 108; % number of TRs per scanning run (coded for fixed value; adjust code structure if variable TR counts across runs)
    
    load([mvpa_dir '/onsets_mvpa.mat']); % load in your SPM-formatted onsets file (specifies onsets of each event time in seconds)
    load([mvpa_dir '/' data_imgs_to_use]); %loads predefined cell array called raw_filenames into memory
    num_runs = length(raw_filenames)/num_TP_per_run; %calculate the number of runs (allows script to work flexibly for subjects with missing runs)
    
    % specify previously saved mvpa workspace to bypass time-consuming data extraction and preprocessing (if not yet created, this will specify it's name and location)
    mvpa_workspace = [mvpa_dir '/' subj_id '_' roi_name '_s4mm.mat'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set flags (% unless otherwise indicated: 1 = yes please, 0 = no thanks)
    
    flags.num_full_iter = 1; % number of times to run the entire classification process, including feature selection
    flags.num_results_iter = 1; % number of times to run the post-feature selection classification process (select subset of the data and train/test classifier)
    flags.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data
    flags.equate_number_of_trials_in_cond_1_and_2 = 1; % equate number of trials in conditions 1 and 2 (RECOMMENDED)
    flags.anova_p_thresh = 1;  % p-value threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)
    flags.anova_nVox_thresh = nVox; % alternative to specifying p-value threshold; uses top N voxels (0 = DON'T PERFORM ANY FEATURE SELECTION)
    flags.perform_second_round_of_zscoring = 1;  % z-score data again immediately prior to classification (because trial selection and TR selection undoes initial z-scoring)
    flags.remove_artdetect_outliers = 0; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
    flags.artdetect_motion_thresh = 0; % specify ArtDetect bin for motion outliers (requires that custom ArtDetect scripts have already been run to flag outlier trials)
    flags.artdetect_global_signal_thresh = 0; % specify ArtDetect bin for global signal outliers (requires that custom ArtDetect scripts have already been run to flag outlier trials)
    flags.remove_outlier_trials = 3;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)
    flags.generate_importance_maps = 1; % 1=generate importance maps based on classification weights (scaled by the mean univariate activity of each condition)
    flags.write_data_log_to_text_file=1; % save a .txt file that logs critical classification performance data
    flags.save_data_log_as_mat_file =0; % save a .mat file that logs critical classification performance data
    flags.lambda = penalty; % penalty parameter for logistic regression classifier
    flags.optimize_penalty_param = 0; % 1 = peformed nested cross-validation for empirical penalty optimization
    
    % specify which conditions to use for classification (must correspond to the names of conditions specified below)
    condnames =  {condition1, condition2};
    
%     if strcmp(TrialPortion,'ScOb12'); %use Scene-Obj1 and Scene-Obj2 portions of trial
      TRs_to_average_over = [4 5]; %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier
      TR_weights = [0 0 0 25 75 0]; % should sum to 1 (e.g., [0 0 .5 .5 0 0] indicates to average the 3rd and 4th post-stimulus TRs)
%     %       TR_weights = [.5 .5 0 0 0]; % should sum to 1 (e.g., [0 0 .5 .5 0 0] indicates to average the 3rd and 4th post-stimulus TRs)
%         %TR_weights = [.0072 .2168 .3781 .2742 .1237];  % from SPM canonical values at 1,3,5,7,and 9 sec post-stimulus
%     elseif strcmp(TrialPortion,'ScOnly'); %use SceneOnly portion of trial (but jittered CSI, so only use 1st TR)
%         TRs_to_average_over = [1 2]; %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier
%         TR_weights = [0.5 0.5]; % should sum to 1 (e.g., [0 0 .5 .5 0 0] indicates to average the 3rd and 4th post-stimulus TRs)
%     end
    
    % specify classifier
    class_args.train_funct_name = 'train_pLR';
    class_args.test_funct_name = 'test_pLR';
    class_args.penalty = flags.lambda;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     traintest = {'EXP_only','CM_only','EXP>CM','EXP>CM_first_run','EXP>CM_last_run','CM>EXP'};
    
    
    % specify directory names, create new directories if necessary
    if flags.optimize_penalty_param == 1
        dir_str = [condition1 '_vs_' condition2 '_' class_args.train_funct_name(6:end) 'OPTIMAL_pen_' num2str(flags.anova_nVox_thresh) 'vox'];
    else
        dir_str = [condition1 '_vs_' condition2 '_' class_args.train_funct_name(6:end) '_pen' num2str(class_args.penalty) '_' num2str(flags.anova_nVox_thresh) 'vox'];
    end
    
    if flags.generate_importance_maps == 1
        importance_maps_dir=[expt_dir '/mvpa_results/within_subj_class/importance_maps/' dir_str];
        if ~exist(importance_maps_dir,'dir')
            mkdir(importance_maps_dir);
        end
    end    
       
    if flags.write_data_log_to_text_file == 1
        results_data_logs_txt_dir=[expt_dir '/mvpa_results/within_subj_class/class_perf/' dir_str];
        if ~exist(results_data_logs_txt_dir, 'dir')
            mkdir(results_data_logs_txt_dir)
        end
    end
    
    if flags.save_data_log_as_mat_file == 1
        results_data_logs_mat_dir=[expt_dir '/mvpa_results/within_subj_class/data_logs/' dir_str];
        if ~exist(results_data_logs_mat_dir, 'dir')
            mkdir(results_data_logs_mat_dir)
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~exist(mvpa_workspace,'file') %if a workspace .mat file has not yet been created for this subject
        
        % call script to run time-consuming data preprocessing routines
        % (load pattern, detrend, high-pass filter, and z-score)
        [subj] = JR_mvpa_load_and_preprocess_raw_data(subj_id, expt_name, roi_name, roi_file, raw_filenames, num_runs, num_TP_per_run);
        save(mvpa_workspace, 'subj');
    else
        load(mvpa_workspace)  %load workspace
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Extract info about conditions from onsets file (uses SPM-based coding
    % of the onsets (in seconds) of trials in each condition of interest)
    num_conds = size(onsets,2);
    
    all_onsets = zeros(num_conds,size(subj.patterns{end}.mat,2)); % initialize all_onsets matrix as conditions x timepoints
    
    for cond = 1: num_conds
        for trial = 1: length(onsets{cond})
            time_idx = round(onsets{cond}(trial)/2)+1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
            all_onsets(cond,time_idx) = 1;
        end
    end
    
    % condense all_onsets matrix from having one column per TR to having
    % once column per trial (remove TRs that have 0's for all
    % conditions, i.e. rest timepoints)
    
    condensed_regs_all = [];
    condensed_runs = [];
    trial_counter = 1;
    for i = 1: size(all_onsets,2)
        if ~isempty(find(all_onsets(:,i))) % if not a rest timepoint
            condensed_regs_all(:,trial_counter) = all_onsets(:,i);
            condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
            trial_counter = trial_counter + 1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % select TRs of interest (to correspond with peak post-stim BOLD response)    
    all_trials = sum(all_onsets,1); % vector of all trials
            data_by_TR(1,:,:) = TR_weights(1)*subj.patterns{end}.mat(:,find(all_trials)+0); % 1st TR (0-2 sec)
            data_by_TR(2,:,:) = TR_weights(2)*subj.patterns{end}.mat(:,find(all_trials)+1); % 2nd TR (2-4 sec)
            data_by_TR(3,:,:) = TR_weights(3)*subj.patterns{end}.mat(:,find(all_trials)+2); % 3rd TR (4-6 sec)
            data_by_TR(4,:,:) = TR_weights(4)*subj.patterns{end}.mat(:,find(all_trials)+3); % 4th TR (6-8 sec)
            data_by_TR(5,:,:) = TR_weights(5)*subj.patterns{end}.mat(:,find(all_trials)+4); % 5th TR (8-10 sec)
            data_by_TR(6,:,:) = TR_weights(5)*subj.patterns{end}.mat(:,find(all_trials)+5); % 6th TR (10-12 sec)

    temporally_condensed_data = squeeze(sum(data_by_TR(TRs_to_average_over,:,:),1));
    clear data_by_TR; % clean up matlab workspace to save memory
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Exclude trials determined to be outliers by custom ArtDetect script
    % Guide to outlier file cell arrays...
    % Movement thresholds: .2 .25 .3 .35 .4 .4 .5
    % Global signal thresholds: 2 2.5 3 3.5 4 4.5 5
    
    if flags.remove_artdetect_outliers == 1
        load([expt_dir '/outlier_indices/' subj_id '_outlier_indices']); %load outlier indices
        
        m_outliers = movement_outlier_trials{flags.artdetect_motion_thresh};  % remove trials with more than .35mm/TR of movement
        gs_outliers = global_signal_outlier_trials{flags.artdetect_global_signal_thresh}; % remove trials with global signal change of +/- 3.5 SD from mean
        combined_outliers = union(m_outliers,gs_outliers);
        
        condensed_regs_all(:,combined_outliers) = 0;
        
        display([num2str(length(m_outliers)) ' movement outlier trials flagged']);
        display([num2str(length(gs_outliers)) ' global signal outlier trials flagged']);
        display([num2str(length(combined_outliers)) ' total outlier trials excluded']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % on-the-fly outlier detection/removal;
    if flags.remove_outlier_trials ~= 0
        mean_across_voxels = mean(temporally_condensed_data,1);
        z_mean_across_voxels = zscore(mean_across_voxels);
        upper_outliers = find(z_mean_across_voxels> flags.remove_outlier_trials);
        lower_outliers = find(z_mean_across_voxels< -1 * flags.remove_outlier_trials);
        all_outliers = union(upper_outliers,lower_outliers)
        condensed_regs_all(:,all_outliers) = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    subj_original = subj; % backup the subj structure before entering k-loop
    condensed_regs_all_original = condensed_regs_all; % backup condensed_regs_all before entering k-loop
    x = 0;  %initialize the counter x (gets incremented during each classification run-through)
    for k = 1:flags.num_full_iter
        
        subj = subj_original;
        condensed_regs_all = condensed_regs_all_original;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define conditions of interest (EXPERIMENT-SPECIFIC CODE NEEDED TO
        % INDEX TRIALS FROM EACH CONDITION OF INTEREST)
        %
        
       
        % define mnemonic states based on individual conditions or combinations of conditions
        
        %CONDITION NUMBERS FOR MODEL1:
%             1= face
%             2= scram
%             3= scene
%             4= objects
%             5= blank


        face = condensed_regs_all([1],:); %"OLD" = objective old status
        scram = condensed_regs_all([2],:);
        scenes = condensed_regs_all([3],:);
        objects = condensed_regs_all([4],:); %"OLD" = objective old status

        
%         New_obj_ScOb12 = sum(condensed_regs_all([26,27,29,30],:));
% 
% %         old = sum(condensed_regs_all([1:5 9:13],:)); %"old" = subjective old status
% %         new = sum(condensed_regs_all([7:8 15:16],:));
%         
%         Hits = sum(condensed_regs_all([1:3,7:9],:));
%         Misses = sum(condensed_regs_all([4:6,10:12],:));
%         FAs = sum(condensed_regs_all([16:18,22:24,28:30],:));
%         CRs = sum(condensed_regs_all([13:15,19:21,25:27],:));
% 
%         Hits_ScOb12 = sum(condensed_regs_all([2:3,8:9],:));
%         Misses_ScOb12 = sum(condensed_regs_all([5:6,11:12],:));
%         FAs_ScOb12 = sum(condensed_regs_all([17:18,23:24,29:30],:));
%         CRs_ScOb12 = sum(condensed_regs_all([14:15,20:21,26:27],:));
% 
%         ObjPres_ScOnly = sum(condensed_regs_all([1,4,7,10,13,16,19,22],:));
%         ObjAbs_ScOnly = sum(condensed_regs_all([25,28],:));

%         Recoll=sum(condensed_regs_all([1:2],:));
%         StRecoll=sum(condensed_regs_all([1],:));
%         MedRecoll=sum(condensed_regs_all([2],:));
%         FalseRecoll=sum(condensed_regs_all([9:10],:));
%         All_R_resps=sum(condensed_regs_all([1:2,9:10],:)); %all Recollect responses, regardless of accuracy
% 
%         Famil=sum(condensed_regs_all([3:4],:));
%         StFamil=sum(condensed_regs_all([3],:));
%         MedFamil=sum(condensed_regs_all([4],:));
%         FalseFamil=sum(condensed_regs_all([11:12],:));
%         All_F_resps=sum(condensed_regs_all([3:4,11:12],:));
% 
%         Know=sum(condensed_regs_all([5],:));
%         FalseKnow=sum(condensed_regs_all([13],:));
%         All_K_resps=sum(condensed_regs_all([5,13],:));
        
        
        
        
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%% Specify which training/testing scheme you wish you to use %%%%%%%%%%%%%%%%%
         
%           if which_traintest == 1  %% Train and Test on EXPLICIT runs only
%               condensed_runs(condensed_runs>0)=4;  % give all countermeasures runs a run label of 4; since there will be no active_trials from these runs, cross-validation will be 4-fold
%           elseif which_traintest == 2  %% Train and Test on COUNTERMEASURES runs only
%               condensed_runs(condensed_runs<=4)=5;  % give all explicit memory runs a run label of 5; since there will be no active_trials from these runs, cross-validation will be 4-fold
%               condensed_runs = condensed_runs - 4;  % subtract 4 from each run label, turning run5 into run1, run6 into run2, etc.
%           elseif which_traintest == 3  %% Train on EXPLICIT runs and Test on COUNTERMEASURES runs
%               condensed_runs(condensed_runs==1)=2; % data from runs 1-4 will be the first training set
%               condensed_runs(condensed_runs==2)=1; % data from runs 5-8 will be the first testing set
%           elseif which_traintest == 4  %% Train on EXPLICIT runs and Test on FIRST COUNTERMEASURES run only
%               condensed_runs(condensed_runs<=4)=2; % data from runs 1-4 will be the first training set
%               condensed_runs(condensed_runs==5)=1; % data from run 5 will be the first testing set
%               condensed_runs(condensed_runs>5)=0; % data from other CM runs ignored
%           elseif which_traintest == 5  %% Train on EXPLICIT runs and Test on LAST COUNTERMEASURES run only
%               condensed_runs(condensed_runs<=4)=2; % data from runs 1-4 will be the first training set
%               condensed_runs(condensed_runs==8)=1; % data from run 8 will be the first testing set
%               condensed_runs(intersect(find(condensed_runs>4), find(condensed_runs<8)))=0; % data from other CM runs ignored
%           elseif which_traintest == 6  %% Train on COUNTERMEASURES runs and Test on EXPLICIT runs
%               condensed_runs(condensed_runs<=4)=1; % data from runs 5-8 will be the first training set
%               condensed_runs(condensed_runs>4)=2; % data from runs 1-4 will be the first testing set
%           end
%                            
%           if which_traintest==3 || which_traintest==4 || which_traintest==5
%               % code the CM_OLD trials as examples of condition 1 (e.g., Hits) and the CM_NEW trials as examples of condition 2 (e.g., CRs)
%               eval([condnames{1} ' = ' condnames{1} ' + CM_hits;']);
%               eval([condnames{2} ' = ' condnames{2} ' + CM_FAs;']);
%           elseif which_traintest==6
%               eval([condnames{1} ' = ' condnames{1} ' + EX_hits;']);
%               eval([condnames{2} ' = ' condnames{2} ' + EX_CRs;']);
%           end
         
         
        num_runs = length(unique(nonzeros(condensed_runs)));  %(MU) set num_runs as the length of unique run numbers (i.e., if condensed_runs has values 1 1 1 2 2 2 3 3 3, num_runs=3)
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
        
        %assign conditions to train/test classifier on
        condensed_regs_of_interest = [];
        eval(['condensed_regs_of_interest(1,:) = ' condnames{1} ';'])
        eval(['condensed_regs_of_interest(2,:) = ' condnames{2} ';'])
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Set up regressors and selectors
        
        % initialize regressors object
        subj = init_object(subj,'regressors','conds');
        subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
        subj = set_objfield(subj,'regressors','conds','condnames',condnames);
        
        % add new condensed activation pattern
        subj = duplicate_object(subj,'pattern','epi_d_hp_z','epi_d_hp_z_condensed');
        subj = set_mat(subj,'pattern','epi_d_hp_z_condensed',temporally_condensed_data,'ignore_diff_size',true);
        zhist = sprintf('Pattern ''%s'' created by JR custom code','epi_d_hp_z_condensed');
        subj = add_history(subj,'pattern','epi_d_hp_z_condensed',zhist,true);
        
        % clean up workspace to save RAM
        subj = remove_mat(subj,'pattern','epi_d_hp_z'); % remove original uncondensed pattern (full timeseries)
        
        % update run vector to condensed format
        subj.selectors{1}.mat = condensed_runs;
        subj.selectors{1}.matsize = size(condensed_runs);
        
        % "activate" only those trials of interest (from regs_of_interest) before creating cross-validation indices
        active_trials = find(sum(condensed_regs_of_interest)); % find those trials where the sum (across columns) is not zeros
        
        actives_selector = zeros(1,size(condensed_regs_all,2)); % create new selector object; begin by intializing vector of all zeros
        actives_selector(active_trials) = 1; % remove all non-"regs_of_interst" trials (set to one)
        
        subj = init_object(subj,'selector','conditions_of_interest'); %initialize new selector object called 'conditions_of_interest'
        subj = set_mat(subj,'selector','conditions_of_interest',actives_selector); % set this selector object to have the values specified by actives_selector
        
        
        subj = create_xvalid_indices(subj,'runs','actives_selname','conditions_of_interest', 'ignore_jumbled_runs', true,'ignore_runs_zeros',true); % create cross-validation indices (new selector group), using only the the trials specified by actives_selector
        
        
        
        % run feature selection ANOVA: specify pvalue (if desired)
        statmap_arg.use_mvpa_ver = true;
        if flags.anova_p_thresh ~= 1
            
            subj = JR_feature_select(subj,'epi_d_hp_z_condensed','conds','runs_xval','thresh',flags.anova_p_thresh, 'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
        else
            classifier_mask = subj.masks{1}.name; % use original mask
        end
        
        % run feature selection ANOVA: specify #of voxels (if desired)
        if flags.anova_nVox_thresh ~=0
            
            subj = JR_feature_select_top_N_vox(subj,'epi_d_hp_z_condensed','conds','runs_xval','nVox_thresh',flags.anova_nVox_thresh,'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
        else
            classifier_mask = subj.masks{1}.name; % use original mask
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        subj_prebalancing = subj; % backup the subj structure before entering n-loop
        active_trials_prebalancing = active_trials; % backup condensed_regs_all before entering n-loop
        
        for n = 1: flags.num_results_iter
            
            % restore original versions before continuing with next loop iteration
            active_trials = active_trials_prebalancing;
            subj = subj_prebalancing;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % balance the number of Class A and Class B trials within each run 
            % (prevents classifier from becoming biased to predict one
            % condition much more frequently than the other)
            
            subj = create_balanced_xvalid_selectors(subj,'conds','runs_xval'); % creates a new 'runs_xval_bal' selector group
            
            % after running the create_balanced_xvalid_selectors function,
            % which excludes random subsets of trials from each run to
            % balance the classes, we need to get a new list of the 'active trials'
            new_active_trials =[];
            for rr = 1:num_runs
                new_active_trials = horzcat(new_active_trials, find(subj.selectors{end-num_runs+rr}.mat==2));
            end
            new_actives_selector= zeros(1,size(condensed_regs_all,2));
            new_actives_selector(new_active_trials)=1;
            subj = init_object(subj,'selector','conditions_of_interest_bal_within_runs'); %initialize selector object
            subj = set_mat(subj,'selector','conditions_of_interest_bal_within_runs',new_actives_selector);
            active_trials = new_active_trials;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if flags.perform_second_round_of_zscoring == 1  % z-score the data prior to classification (active trials only)
                subj.patterns{5}.mat(:,active_trials) = zscore(subj.patterns{5}.mat(:,active_trials)')';
                display('Performing second round of z-scoring')
            end
            
            if flags.optimize_penalty_param == 1 && x == 0 % find empirically optimal penalty via nested cross-validation
                %only run this computationally-intensive function during the first pass through this script (i.e., when x==0), and use the resulting penalty value for all subsequent classifications
                [subj best_penalties penalty_iteration_results] = optimal_pLR_penalty(subj,'epi_d_hp_z_condensed','conds','runs_final_xval','runs_final',classifier_mask,'conditions_of_interest_final','use_iteration_perf',false,'perform_final_classification',false);
                class_args.penalty = best_penalties; % because 'use_iteration_perf' is set to 'false', best_penalties will be only a single value (averaged across the nested cross-validation iterations)
            end
            
            % OPTIONAL: To ensure that there is no bias present in the analysis, scramble the condition labels and then measure classification performance (should be at chance)
            %[subj] =  JR_scramble_regressors(subj,'conds','runs','conditions_of_interest_bal_within_runs','conds_scrambled')
             

            
            %%%%%%%%%%%%%%%%%%%%%% RUN THE CLASSIFIER (CROSS-VALIDATION)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for p = 1:flags.num_iter_with_same_data  % can run this loop multiple times if using a stochastic (i.e. non-deterministic) classification algorithm like backpropagation neural nets
                x=x+1; % increment results iteration counter
                
                [subj results{x}] = cross_validation(subj,'epi_d_hp_z_condensed','conds','runs_xval_bal',classifier_mask,class_args);
                %[subj results{x}] = cross_validation(subj,'epi_d_hp_z_condensed','conds_scrambled','runs_xval_bal',classifier_mask,class_args); % if using scrambled regressors
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % do some important RAM clean-up and archive the feature (voxel) weights
                for y = 1:num_runs
                    if flags.generate_importance_maps == 1
                        % save weights to pass to JR_interpret_weights
                        if strcmp(class_args.train_funct_name,'train_bp')
                            results_IW{x}.iterations(y).scratchpad.net.IW{1} = results{x}.iterations(y).scratchpad.net.IW{1};
                        elseif strcmp(class_args.train_funct_name,'train_pLR')
                            results_IW{x}.iterations(y).scratchpad.net.IW{1} = results{x}.iterations(y).scratchpad.weights';
                        elseif strcmp(class_args.train_funct_name,'train_svdlr')
                            results_IW{x}.iterations(y).scratchpad.net.IW{1} = results{x}.iterations(y).scratchpad.W';
                        else
                            results_IW{x}.iterations(y).scratchpad.net.IW{1} = results{x}.iterations(y).scratchpad.w';
                        end
                    end
                    results{x}.iterations(y).scratchpad.net.inputs{1}.exampleInput=[]; % delete huge data object from results scratchpad to free up RAM
                end
                
                
                % analyze the results in more detail
                correct_vector = []; desireds_vector = []; guesses_vector = []; prob_class_A_vector = [];
                
%                 if which_traintest >=3
%                 	num_testing_sets = 1;  %% only summarize and record data from COUNTERMEASURES runs (1st cross-validation testing set)
%                 else
%                     num_testing_sets = num_runs;
%                 end
                
                
                for a = 1:num_runs % concatenate the results across all cross-validation testing sets
                    correct_vector = horzcat(correct_vector,results{x}.iterations(a).perfmet.corrects); % create a vector of the whether or not the classifier's guess on each trial was correct (1=correct; 2=incorrect)
                    desireds_vector = horzcat(desireds_vector,results{x}.iterations(a).perfmet.desireds); % create a vector of the true class of each trial (i.e., the desired outcome) (1=Class A; 2=Class B)
                    guesses_vector = horzcat(guesses_vector,results{x}.iterations(a).perfmet.guesses); % create a vector of the classifier's guesses on each trial (1=Class A; 2=Class B)
                    prob_class_A_vector = horzcat(prob_class_A_vector, results{x}.iterations(a).acts(1,:)); % create a vector of the classifier's scalar probability estimates that each trial is a trial of Class A
                end
                
                overall_accuracy = mean(correct_vector); %accuracy measure based on performance across all trials (rather than averaging the performance of each run)
                overall_hit_rate = mean(correct_vector(desireds_vector==1)); %probability of labeling examples of Cond1 as Cond1
                overall_fa_rate = 1-mean(correct_vector(desireds_vector==2)); %probability of labeling examples of Cond2 as Cond1
                overall_d_prime = norminv(overall_hit_rate)-norminv(overall_fa_rate); %measure of classifier sensitivity              
                    
                % sort by absolute value of classifier "confidence"
                [abs_sorted_diffs abs_ind] = sort(abs(prob_class_A_vector-0.5),2,'descend');
                abs_correct_sorted = correct_vector(abs_ind);
                abs_desireds_sorted = desireds_vector(abs_ind);
                
                
                % compute accuracy for top N % of trials (sorted by classifier 'confidence')
                num_trials(x) = length(abs_correct_sorted);
                display(['Number of testing trials per bin: ' num2str(num_trials(x)/2)])
                
                bin_intervals =1:-.05:.05; % twenty bins, ranging from 1.0 (performance across all trials) to 0.5 (performance for Top 5% of trials, ranked based on classifier confidence)
                for acc_bin = 1:20
                    included_trial_inds = 1:ceil(num_trials(x)*bin_intervals(acc_bin));
                    class_A_inds = find(abs_desireds_sorted==1);
                    class_B_inds = find(abs_desireds_sorted==2);
                    
                    acc_percentiles(acc_bin)= mean(abs_correct_sorted(included_trial_inds));
                end
                
                % print the top 100%, 75%, 50%, and 25% accuracy to command window
                display(acc_percentiles([1 6 11 16]))
                
                % record number of trials in each class at each classification confidence quartile
                class_counts_by_quartile_rank(1,1) = count(abs_desireds_sorted((1:ceil(num_trials(x)*1.0)))==1);
                class_counts_by_quartile_rank(1,2) = count(abs_desireds_sorted((1:ceil(num_trials(x)*1.0)))==2);
                class_counts_by_quartile_rank(2,1) = count(abs_desireds_sorted((1:ceil(num_trials(x)*.75)))==1);
                class_counts_by_quartile_rank(2,2) = count(abs_desireds_sorted((1:ceil(num_trials(x)*.75)))==2);
                class_counts_by_quartile_rank(3,1) = count(abs_desireds_sorted((1:ceil(num_trials(x)*.50)))==1);
                class_counts_by_quartile_rank(3,2) = count(abs_desireds_sorted((1:ceil(num_trials(x)*.50)))==2);
                class_counts_by_quartile_rank(4,1) = count(abs_desireds_sorted((1:ceil(num_trials(x)*.25)))==1);
                class_counts_by_quartile_rank(4,2) = count(abs_desireds_sorted((1:ceil(num_trials(x)*.25)))==2);                          
                              
                % sort by signed classifier "confidence" (for ROI curves)
                [sorted_diffs ind] = sort(prob_class_A_vector,2,'descend');
                correct_sorted = correct_vector(ind);
                desireds_sorted = desireds_vector(ind);
                
                % create continuous ROC function
                for i = 1:length(sorted_diffs);
                    hit_rate(i) = length(intersect(find(desireds_sorted == 1),[1:i])) / length(find(desireds_sorted == 1));
                    fa_rate(i) = length(intersect(find(desireds_sorted == 2),[1:i])) / length(find(desireds_sorted == 2));
                end
                
                % plot ROC curve as Matlab figure
                %                     figure
                %                     plot(fa_rate,hit_rate,'.-')
                %                     hold on
                %                     plot([0 1],[0 1],'r')
                %                     xlabel('P(Old|New)')
                %                     ylabel('P(Old|Old)')
                
                
                %compute and display area under the ROC curve (AUC); with Hit/FA rates computed across full data continuum
                auc_overall = auroc(hit_rate',fa_rate')
                
                % create ROC function with 80 discrete bins (this allows the ROC curves to be more easily averaged across individuals)
                roc_bin_intervals = .975:-.025:-1;
                for bin_num = 1:80
                    hits_80(bin_num)=length(intersect(find(desireds_sorted == 1),find(sorted_diffs>roc_bin_intervals(bin_num)))) / length(find(desireds_sorted == 1));
                    fas_80(bin_num)=length(intersect(find(desireds_sorted == 2),find(sorted_diffs>roc_bin_intervals(bin_num)))) / length(find(desireds_sorted == 2));
                end
                auc_80_bins = auroc(hits_80',fas_80');
                
                
                if flags.write_data_log_to_text_file==1
                    
                    data_log.overall_acc(x)=overall_accuracy;
                    data_log.hits(x)=overall_hit_rate;
                    data_log.FAs(x)=overall_fa_rate;
                    data_log.d_prime(x)=overall_d_prime;
                    data_log.acc_percentiles(x,:) = acc_percentiles;
                    data_log.penalty_param(x) = class_args.penalty;
                    data_log.class_counts_by_quartile_rank(:,:,x) = class_counts_by_quartile_rank;
                    data_log.auc_overall(x) = auc_overall;
                    data_log.auc_80_bins(x) = auc_80_bins;
                    data_log.roc_80_bin_hits(x,:)= hits_80;
                    data_log.roc_80_bin_fas(x,:)= fas_80;
                end
            end
        end
    end
    
    
    if flags.save_data_log_as_mat_file ==1;
        save_cmd = ['save ' results_data_logs_mat_dir '/' subj_id '_' condnames{1} '_vs_' condnames{2} '.mat data_log flags'];
        eval(save_cmd);
    end
    
    if flags.write_data_log_to_text_file==1
        
        filename= [results_data_logs_txt_dir '/' subj_id '_' condnames{1} '_vs_' condnames{2} '.txt'];
        fid=fopen(filename, 'wt');
        fprintf(fid, '%s\r\n', ['subj_id = ' subj_id]);
        fprintf(fid, '%s\r\n', ['ROI_name = ' roi_name]);
        fprintf(fid, '%s\r\n', ['data_imgs_to_use =' data_imgs_to_use]);
        fprintf(fid, '%s\r\n', ['TR_weights = ' num2str(TR_weights)]);
        fprintf(fid, '%s\r\n', ['classification:' condnames{1} ' vs. ' condnames{2}]);
        fprintf(fid, '%s\r\n', ['Number of Trials Per Bin: ' num2str(mean(num_trials)/2)]);
        fprintf(fid, '%s\r\n', ['flags.perform_second_round_of_zscoring = ' num2str(flags.perform_second_round_of_zscoring)]);
        fprintf(fid, '%s\r\n', ['flags.remove_mvpa_outlier_trials (std dev) = ' num2str(flags.remove_outlier_trials)]);
        fprintf(fid, '%s\r\n', ['flags.remove_artdetect_outlier_trials (std dev) = ' num2str(flags.remove_artdetect_outliers)]);
        fprintf(fid, '%s\r\n', ['flags.artdetect_motion_thresh = ' num2str(flags.artdetect_motion_thresh)]);
        fprintf(fid, '%s\r\n', ['flags.artdetect_global_signal_thresh = ' num2str(flags.artdetect_global_signal_thresh)]);

        if isfield(class_args, 'penalty')
            fprintf(fid, '%s\r\n', ['penalty param = ' num2str(class_args.penalty)]);
        end
        fprintf(fid, '\n\n');
        
        for q=1:x
            fprintf(fid, '%4.4f\t', q);           
            fprintf(fid, '%4.4f\t', data_log.auc_overall(q));
            fprintf(fid, '%4.0f\t', (num_trials(q)/2));
            fprintf(fid, '%4.4f\t', data_log.acc_percentiles(q,:));
            fprintf(fid, '\t');
            fprintf(fid, '%4.4f\t', data_log.roc_80_bin_hits(q,:));
            fprintf(fid, '\t');
            fprintf(fid, '%4.4f\t', data_log.roc_80_bin_fas(q,:));
            fprintf(fid, '\n');                     
        end
        
        fprintf(fid, '%s\t', 'mean');        
        fprintf(fid, '%4.4f\t', mean(data_log.auc_overall));
        fprintf(fid, '%4.0f\t', (mean(num_trials(q))/2));
        fprintf(fid, '%4.4f\t', mean(data_log.acc_percentiles,1));
        fprintf(fid, '\t');
        fprintf(fid, '%4.4f\t', mean(data_log.roc_80_bin_hits,1));
        fprintf(fid, '\t');
        fprintf(fid, '%4.4f\t', mean(data_log.roc_80_bin_fas,1));
        fprintf(fid, '\n');
        fclose(fid);
    end
    
    
    
    if flags.generate_importance_maps == 1;
        
        subj = JR_extract_voxel_importance_values(subj, results,results_IW);
               
        load([expt_dir '/Masks/vol_info.mat']); %get functional data resolution info for spm .img writing
        % To create vol_info.mat, run the following command in the matlab workspace:
        % vol_info = spm_vol('name_of_image_file_with_same_dimensions_as_data_being_used_by_the_classifier.img'); save vol_info.mat vol_info
        
        impmap1 = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 1
        impmap2 = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 2
        
        if flags.anova_p_thresh == 1 % NO ANOVA VERSION
            
            voxel_inds = find(subj.masks{end}.mat); %get mask voxel indices
            for j = 1:num_runs
                temp1 = zeros(vol_info.dim); %initialize appropriately sized matrix
                temp2 = zeros(vol_info.dim); %initialize appropriately sized matrix
                temp1(voxel_inds)=subj.patterns{end-num_runs+j}.mat(:,1); %store impmap values at appropriate voxel indices
                temp2(voxel_inds)=subj.patterns{end-num_runs+j}.mat(:,2); %store impmap values at appropriate voxel indices
                impmap1 = impmap1+temp1; %add values cumulatively across iterations
                impmap2 = impmap2+temp2;
            end
            
            impmap1 = impmap1/num_runs*1000; %compute average and multiply by 1000 for scaling
            impmap2 = impmap2/num_runs*1000; %compute average and multiply by 1000 for scaling
            
            vol_info.fname = [importance_maps_dir '/' subj_id '_' condnames{1} '.img'];
            spm_write_vol(vol_info,impmap1);
            vol_info.fname = [importance_maps_dir '/' subj_id '_' condnames{2} '.img'];
            spm_write_vol(vol_info,impmap2);
            
        else % ANOVA VERSION
            
            for j = 1:num_runs
                temp1 = zeros(vol_info.dim); %initialize appropriately sized matrix
                temp2 = zeros(vol_info.dim); %initialize appropriately sized matrix
                voxel_inds{j} = find(subj.masks{end-num_runs+j}.mat); %get mask voxel indices
                temp1(voxel_inds{j})=subj.patterns{end-num_runs+j}.mat(:,1); %store impmap values at appropriate voxel indices
                temp2(voxel_inds{j})=subj.patterns{end-num_runs+j}.mat(:,2); %store impmap values at appropriate voxel indices
                impmap1 = impmap1+temp1; %add values cumulatively across iterations
                impmap2 = impmap2+temp2;
            end
            
            %sum across masks to get composite mask (where value of each voxel = number of runs for which that voxel was included)
            composite_mask = zeros(vol_info.dim);
            for i = 2:size(subj.masks,2)  %exclude first mask (it's the starting ROI)
                composite_mask = composite_mask+subj.masks{i}.mat;
            end
            voxels_to_exclude = find(composite_mask<=5);  % exclude voxels that exist for fewer than 5 of the ANOVA masks
            impmap1(voxels_to_exclude)=0;
            impmap2(voxels_to_exclude)=0;
            
            impmap1_avg = impmap1./composite_mask * 1000;  % divide by number of observations contributing to each sum (to get avg) and multiply by 1000 for scaling
            impmap2_avg = impmap2./composite_mask * 1000;  % divide by number of observations contributing to each sum (to get avg) and multiply by 1000 for scaling
            
            vol_info.fname = [importance_maps_dir '/' subj_id '_' condnames{1} '_p' num2str(flags.anova_p_thresh) '.img'];
            spm_write_vol(vol_info,impmap1_avg);
            vol_info.fname = [importance_maps_dir '/' subj_id '_' condnames{2} '_p' num2str(flags.anova_p_thresh) '.img'];
            spm_write_vol(vol_info,impmap2_avg);
        end
        
    end
    
    time2finish = toc/60;
    display(['Finished ' subj_id ' in ' num2str(time2finish) ' minutes']);
    
    keep b subj_array condition1 condition2 nVox penalty   % use keep.m function to clear all variables from memory except for these four
    
end


