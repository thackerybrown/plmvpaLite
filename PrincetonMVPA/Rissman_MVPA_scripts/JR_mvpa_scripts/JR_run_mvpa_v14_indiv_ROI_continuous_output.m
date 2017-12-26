function []= JR_run_mvpa_v14_indiv_ROI_continuous_output(subj_array, condition1, condition2, nVox)

% example usage:  JR_run_mvpa_general([1 2 4 5 8 12 13 16], 'Hits', 'CRs')


%%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub_num = 0;

for b=subj_array  % allows subject ID #'s to be specified in short format
    
    sub_num = sub_num + 1;
    
    tic
    if b<10
        subj_id=strcat('s10', num2str(b))
    else
        subj_id=strcat('s1', num2str(b))
    end
    
    %s
    expt_dir = '/Users/Jesse/fMRI/data/PAST/fMRI';
    mvpa_dir = [expt_dir '/' subj_id '/mvpa'];
    behav_dir = ['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id]; 
    
    
    roi_name = 'SEPT09_MVPA_MASK';
    data_imgs_to_use = 'raw_filenames_s8mm_wa.mat'; % .mat file containing names of all functional images
    
    
    
    % load some .mat files into memory
    load([expt_dir '/vol_info.mat']); %get functional data resolution info for spm .img writing
    load([behav_dir '/onsets.mat']); % load in your SPM-formatted onsets file
    
    % OPTIONAL:  specify previously saved mvpa workspace to bypass time-consuming data extraction and preprocessing
    %mvpa_workspace = [expt_dir '/' subj_id '/mvpa/' subj_id '_merged_AAL_ROIs_FIXED_HOLES_wa.mat'];
    
    mvpa_workspace = [expt_dir '/' subj_id '/mvpa/' subj_id '_SEPT09_MVPA_MASK_s8mm_wa.mat'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set flags (% unless otherwise indicated: 1 = yes please, 0 = no thanks)
    flags.num_full_iter = 50;
    flags.num_results_iter = 1; % number of times to run the entire classification process (select subset of the data and train/test classifier)
    flags.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data
    
    flags.equate_number_of_old_new_trials_per_subjective_bin = 0; % equate number of trials per subjective bin
    flags.equate_number_of_trials_in_cond_1_and_2 = 0; % equate number of trials in conditions 1 and 2 (RECOMMENDED)
    flags.anova_p_thresh = 1;  % p-value threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)
    flags.GLM_anova_p_thresh = 0;
    flags.anova_nVox_thresh = nVox; % alternative to specifying p-value threshold; uses top N voxels (0 = DON'T PERFORM ANY FEATURE SELECTION)
    flags.perform_second_round_of_zscoring = 0;  % z-score data again immediately prior to classification (NOT RECOMMENDED)
    flags.remove_artdetect_outliers = 0; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
    flags.artdetect_motion_thresh = 7; % specify ArtDetect bin for motion outliers
    flags.artdetect_global_signal_thresh = 5; % specify ArtDetect bin for global signal outliers
    flags.remove_outlier_trials = 0;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)
    flags.plot_ROC_curve = 1;
    flags.display_performance_breakdown = 1;
    flags.generate_importance_maps = 0;
    flags.generate_weight_maps = 0;
    flags.write_data_log_to_text_file=0;
    flags.save_data_log_as_mat_file =0;
    flags.optimize_penalty_param = 0;
    
    
    % specify which conditions to use for classification (must correspond to the names of conditions specified below)
    condnames =  {condition1, condition2};
    
    TRs_to_average_over = [1 2 3 4 5 6]; %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier
    TR_weights = [0 0 .5 .5 0 0]; % should sum to 1
    %TR_weights = [.0072 .2168 .3781 .2742 .1237 0];  % from SPM canonical
    %values at 1,3,5,7,and 9 sec post-stimulus
    
    % classifier parameters
    %     class_args.train_funct_name = 'train_bp';
    %     class_args.test_funct_name = 'test_bp';
    %     class_args.nHidden = 0;
    
    class_args.train_funct_name = 'train_pLR';
    class_args.test_funct_name = 'test_pLR';
    class_args.penalty = 1;
    %class_args.lambda = class_args.penalty;
    %class_args.mode = 'l2';
    
    %      class_args.train_funct_name = 'train_smlr';
    %      class_args.test_funct_name = 'test_smlr';
    
    if flags.optimize_penalty_param == 1
        dir_str = ['AAL_' condition1 '_vs_' condition2 '_' class_args.train_funct_name(6:end) 'OPTIMAL_pen_' num2str(flags.anova_nVox_thresh) 'vox'];
    else
        dir_str = ['AAL_' condition1 '_vs_' condition2 '_' class_args.train_funct_name(6:end) 'pen' num2str(class_args.penalty) '_' num2str(flags.anova_nVox_thresh) 'vox'];
    end
    
    if flags.generate_importance_maps == 1
        importance_maps_dir=[expt_dir '/mvpa_results/FINAL_ROC_data/within_subj_class/FINAL_IMP_MAPS/' dir_str];
        if ~exist(importance_maps_dir,'dir')
            mkdir(importance_maps_dir);
        end
    end
    
    if flags.generate_weight_maps == 1
        weight_maps_dir=[expt_dir '/mvpa_results/FINAL_ROC_data/within_subj_class/FINAL_WEIGHT_MAPS/' dir_str];
        if ~exist(weight_maps_dir,'dir')
            mkdir(weight_maps_dir);
        end
    end
    
    
    if flags.write_data_log_to_text_file == 1
        xls_results_data_logs_txt_dir=[expt_dir '/mvpa_results/FINAL_ROC_data/within_subj_class/FINAL_CLASS_PERF/' dir_str];
        if ~exist(xls_results_data_logs_txt_dir, 'dir')
            mkdir(xls_results_data_logs_txt_dir)
        end
    end
    
    if flags.save_data_log_as_mat_file == 1
        xls_results_data_logs_mat_dir=[expt_dir '/mvpa_results/FINAL_ROC_data/within_subj_class/data_logs/' dir_str];
        if ~exist(xls_results_data_logs_mat_dir, 'dir')
            mkdir(xls_results_data_logs_mat_dir)
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~exist(mvpa_workspace,'file')
        [subj num_runs num_TP_per_run]= JR_generate_mvpa_workspace_mat_file(subj_id, roi_name, data_imgs_to_use, mvpa_dir); % generate and save workspace
    else
        eval(['load ' mvpa_workspace])  %load workspace
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Load ROI file names and store in ROI_file_names%
    %     roi_img_with_labels = '/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/spm5_AAL_ROIs/rAAL_cluster_labeled.nii'; % for writing classification performance results to .img
    %
    %     roi_folder_name = '/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/spm5_AAL_ROIs/';
    %     load([roi_folder_name 'AAL_ROI_names.mat']);
    %     ROI_file_names=filename;
    number_ROI_files=1; %size(ROI_file_names,2);
    %
    %     %%%Add each ROI file specified in ROI_file_names to subj structure%%%
    %     for i=1:number_ROI_files
    %
    %         current_mask_name = ROI_file_names{i};
    %         current_mask_file = strcat(roi_folder_name, current_mask_name);
    subj = load_spm_mask(subj, 'bilateral_hippocampus', '/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/spm5_AAL_ROIs/rMNI_Hippocampus_B.img');
    
    %    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    % Extract info about conditions from onsets file
    num_conds = size(onsets,2);
    all_regs = zeros(num_conds,num_runs*num_TP_per_run); % initialize regs matrix as conditions x timepoints
    
    for cond = 1: num_conds
        for trial = 1: length(onsets{cond})
            time_idx = onsets{cond}(trial)/2+1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
            all_regs(cond,time_idx) = 1;
        end
    end
    
    % Get rid of any fake "No Response" regressors (condition 11)
    if length(find(all_regs(11,:)))==1 && find(all_regs(11,:)) == (num_runs*num_TP_per_run - 2); %if there is only one value, and it's onset is 2 TRs seconds before the end of the experiment
        all_regs(11,:)=0;
    end
    
    % SPECIAL FIX because reg #6 (NEW_recollect) is a fake placeholder trial for some subjects
    if ~isempty(find(strcmp(subj_id,{'s104','s105','s108','s113','s117','s119','s120','s121','s203','s204','s205','s207','s208'})))
        all_regs(6,:)=0;
    end
    
    % SPECIAL FIX because reg #7 (NEW_HC_old) is a fake placeholder trial for some subjects
    if ~isempty(find(strcmp(subj_id,{'s117','s207'})))
        all_regs(7,:)=0;
    end
    
    % SPECIAL FIX because reg #5 (OLD_HC_new) is a fake placeholder trial for some subjects
    if ~isempty(find(strcmp(subj_id,{'s118','s205','s207'})))
        all_regs(5,:)=0;
    end
    
    
    % condense regs by removing zeros
    condensed_regs_all = [];
    condensed_runs = [];
    trial_counter = 1;
    for i = 1: size(all_regs,2)
        if ~isempty(find(all_regs(:,i))) % if not a rest timepoint
            %condensed_regs_of_interest(:,trial_counter) = regs_of_interest(:,i);
            condensed_regs_all(:,trial_counter) = all_regs(:,i);
            condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
            trial_counter = trial_counter + 1;
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select TRs of interest (to correspond with peak post-stim BOLD response)
    
    %% convert data pattern to single precision to save RAM
    %subj.patterns{end}.mat = single(subj.patterns{end}.mat);
    
    all_trials = sum(all_regs,1); % vector of all trials
    data_by_TR(1,:,:) = TR_weights(1)*subj.patterns{end}.mat(:,find(all_trials)+0); % 1st TR (0-2 sec)
    data_by_TR(2,:,:) = TR_weights(2)*subj.patterns{end}.mat(:,find(all_trials)+1); % 2nd TR (2-4 sec)
    data_by_TR(3,:,:) = TR_weights(3)*subj.patterns{end}.mat(:,find(all_trials)+2); % 3rd TR (4-6 sec)
    data_by_TR(4,:,:) = TR_weights(4)*subj.patterns{end}.mat(:,find(all_trials)+3); % 4th TR (6-8 sec)
    data_by_TR(5,:,:) = TR_weights(5)*subj.patterns{end}.mat(:,find(all_trials)+4); % 5th TR (8-10 sec)
    data_by_TR(6,:,:) = TR_weights(6)*subj.patterns{end}.mat(:,find(all_trials)+5); % 5th TR (8-10 sec)
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
    
    if flags.remove_outlier_trials ~= 0
        % remove outlier trials (timepoints)
        mean_across_voxels = mean(temporally_condensed_data,1);
        z_mean_across_voxels = zscore(mean_across_voxels);
        upper_outliers = find(z_mean_across_voxels> flags.remove_outlier_trials);
        lower_outliers = find(z_mean_across_voxels< -1 * flags.remove_outlier_trials);
        all_outliers = union(upper_outliers,lower_outliers)
        %active_trials(all_outliers) = [];
        condensed_regs_all(:,all_outliers) = 0;
    end
    
    
    
    subj_original = subj;
    condensed_regs_all_original = condensed_regs_all;
    x = 0;  %initialize the counter x (gets incremented during each classification run-through)
    num_biased_iterations = 0;
    for k = 1:flags.num_full_iter
        
        subj = subj_original;
        condensed_regs_all = condensed_regs_all_original;
        
        %            condensed_runs = shuffle(condensed_runs);
        %            display('shuffling runs vector')
        
        % Artificially balance the number of trials given each subjective memory response
        if flags.equate_number_of_old_new_trials_per_subjective_bin == 1
            for j = 1:5
                OLD_trials = find(condensed_regs_all(j,:));
                NEW_trials = find(condensed_regs_all(j+5,:));
                
                num_OLD = length(OLD_trials);
                num_NEW = length(NEW_trials);
                
                if num_OLD > num_NEW
                    rand_array = rand(1,num_OLD);
                    [sorted inds]= sort(rand_array);
                    trials_to_cut = OLD_trials(inds(1:num_OLD-num_NEW));
                    condensed_regs_all(j,trials_to_cut)=0;
                elseif num_OLD < num_NEW
                    rand_array = rand(1,num_NEW);
                    [sorted inds]= sort(rand_array);
                    trials_to_cut = NEW_trials(inds(1:num_NEW-num_OLD));
                    condensed_regs_all(j+5,trials_to_cut)=0;
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % define conditions of interest
        % specify names 'OLD_recollect'    'OLD_hc_old'    'OLD_lc_old'    'OLD_lc_new'    'OLD_hc_new'    'NEW_recollect' 'NEW_hc_old'    'NEW_lc_old'    'NEW_lc_new'    'NEW_hc_new'    'no_resp'
        
        Objective_old = sum(condensed_regs_all(1:5,:));
        Objective_new = sum(condensed_regs_all(6:10,:));
        
        Objective_old_LC_only = sum(condensed_regs_all(3:4,:));
        Objective_new_LC_only = sum(condensed_regs_all(8:9,:));
        
        
        Subjective_old = sum(condensed_regs_all([1 2 3 6 7 8],:));
        Subjective_new = sum(condensed_regs_all([4 5 9 10],:));
        
        Subjective_old_HC_only = sum(condensed_regs_all([1 2 6 7],:));
        Subjective_new_HC_only = sum(condensed_regs_all([5 10],:));
        
        Hits = sum(condensed_regs_all([1 2 3],:));
        Misses = sum(condensed_regs_all([4 5],:));
        CRs = sum(condensed_regs_all([9 10],:));
        FAs = sum(condensed_regs_all(6:8,:));
        
        R_hits = condensed_regs_all(1,:);
        HC_hits = condensed_regs_all(2,:);
        LC_hits = condensed_regs_all(3,:);
        LC_misses = condensed_regs_all(4,:);
        HC_misses = condensed_regs_all(5,:);
        R_FAs = condensed_regs_all(6,:);
        HC_FAs = condensed_regs_all(7,:);
        LC_FAs = condensed_regs_all(8,:);
        LC_CRs = condensed_regs_all(9,:);
        HC_CRs = condensed_regs_all(10,:);
        
        R_and_HC_hits = sum(condensed_regs_all([1 2],:));
        
        HC_resp = sum(condensed_regs_all([2 5 7 10],:));
        LC_resp = sum(condensed_regs_all([3 4 6 8 9],:));
        
        no_resp = condensed_regs_all(11,:); %excluded from analysis
        
        %assign conditions to train/test classifier on
        condensed_regs_of_interest = [];
        eval(['condensed_regs_of_interest(1,:) = ' condnames{1} ';'])
        eval(['condensed_regs_of_interest(2,:) = ' condnames{2} ';'])
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set up regressors and selectors
        
        
        % initialize regressors object
        subj = init_object(subj,'regressors','conds');
        subj = set_mat(subj,'regressors','conds',condensed_regs_of_interest);
        subj = set_objfield(subj,'regressors','conds','condnames',condnames);
        
        
        % add new condensed activation pattern
        subj = duplicate_object(subj,'pattern','spiral_d_hp_z','spiral_d_hp_z_condensed');
        subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed',temporally_condensed_data,'ignore_diff_size',true);
        zhist = sprintf('Pattern ''%s'' created by JR custom code','spiral_d_hp_z_condensed');
        subj = add_history(subj,'pattern','spiral_d_hp_z_condensed',zhist,true);
        
        % clean up workspace to save RAM
        subj = remove_mat(subj,'pattern','spiral_d_hp_z');
        
        % update run vector to condensed format
        subj.selectors{1}.mat = condensed_runs;
        subj.selectors{1}.matsize = size(condensed_runs);
        
        % "activate" only those trials of interest (from regs_of_interest) before creating cross-validation indices
        active_trials = find(sum(condensed_regs_of_interest));
        
        %         if flags.remove_outlier_trials ~= 0
        %             % remove outlier trials (timepoints)
        %             mean_across_voxels = mean(subj.patterns{end}.mat(:,active_trials),1);
        %             z_mean_across_voxels = zscore(mean_across_voxels);
        %             upper_outliers = find(z_mean_across_voxels> flags.remove_outlier_trials);
        %             lower_outliers = find(z_mean_across_voxels< -1 * flags.remove_outlier_trials);
        %             all_outliers = union(upper_outliers,lower_outliers)
        %             active_trials(all_outliers) = [];
        %         end
        
        actives_selector = zeros(1,size(condensed_regs_all,2)); % intialize vector of all zeros
        actives_selector(active_trials) = 1; % remove all non-"regs_of_interst" trials (set to one)
        
        
        subj = init_object(subj,'selector','conditions_of_interest'); %initialize selector object
        subj = set_mat(subj,'selector','conditions_of_interest',actives_selector);
        subj = create_xvalid_indices(subj,'runs','actives_selname','conditions_of_interest');
        
        % run feature selection ANOVA: specify pvalue (if desired)
        statmap_arg.use_mvpa_ver = true;
        if flags.anova_p_thresh ~= 1
            
            subj = JR_feature_select(subj,'spiral_d_hp_z_condensed','conds','runs_xval','thresh',flags.anova_p_thresh, 'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            %subj = JR_feature_select(subj,'spiral_d_hp_z_condensed','conds','runs_xval','thresh',flags.anova_p_thresh, 'greater_than', true,'statmap_funct','statmap_RFE');
            classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
        else
            classifier_mask = subj.masks{1}.name; % use original mask
        end
        
        % run feature selection ANOVA: specify #of voxels (if desired)
        if flags.anova_nVox_thresh ~=0
            
            subj = JR_feature_select_top_N_vox(subj,'spiral_d_hp_z_condensed','conds','runs_xval','nVox_thresh',flags.anova_nVox_thresh,'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            %subj = JR_feature_select_top_N_vox_iterative(subj,'spiral_d_hp_z_condensed','conds','runs_xval','nVox_thresh',flags.anova_nVox_thresh,'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
        else
            classifier_mask = subj.masks{1}.name; % use original mask
        end
        
        % copy important variables before looping, since they'll get modified
        % by the subsequent code
        
        
        subj_orig = subj;
        active_trials_orig = active_trials;
        
        
        
        for n = 1: flags.num_results_iter
            
            % restored original versions before continuing with next loop iteration
            active_trials = active_trials_orig;
            %             [rand_actives rand_idx] = sort(rand(1,length(active_trials)));
            %             active_trials(rand_idx(1:80))=[];
            
            subj = subj_orig; % overwrite subj struct w/ original
            
            
            if flags.equate_number_of_trials_in_cond_1_and_2 == 1
                
                cond1_trials = find(condensed_regs_of_interest(1,active_trials));
                cond2_trials = find(condensed_regs_of_interest(2,active_trials));
                num_cond1 = length(cond1_trials);
                num_cond2 = length(cond2_trials);
                
                if num_cond1 > num_cond2
                    rand_array = rand(1,num_cond1);
                    [sorted inds]= sort(rand_array);
                    trials_to_cut = cond1_trials(inds(1:num_cond1-num_cond2));
                    %condensed_regs_of_interest(1,active_trials(trials_to_cut)) = 0;
                    active_trials(trials_to_cut) = [];
                    
                    display([num2str(length(trials_to_cut)) ' trials cut from ' condnames{1}]);
                elseif num_cond1 < num_cond2
                    rand_array = rand(1,num_cond2);
                    [sorted inds]= sort(rand_array);
                    trials_to_cut = cond2_trials(inds(1:num_cond2-num_cond1));
                    %condensed_regs_of_interest(2,active_trials(trials_to_cut)) = 0;
                    active_trials(trials_to_cut) = [];
                    display([num2str(length(trials_to_cut)) ' trials cut from ' condnames{2}]);
                else
                    display('Trial numbers are already balanced');
                end
            end
            
            actives_selector = zeros(1,size(condensed_regs_all,2)); % intialize vector of all zeros
            actives_selector(active_trials) = 1; % remove all non-"regs_of_interst" trials (set to one)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % update the number of trials in each condition and record the
            % new set of indices
            
            
            cond1_trials = intersect(find(condensed_regs_of_interest(1,:)),active_trials);
            cond2_trials = intersect(find(condensed_regs_of_interest(2,:)),active_trials);
            
            num_cond1 = length(cond1_trials);
            num_cond2 = length(cond2_trials);
            
            
            display([num2str(num_cond1) ' trials in condition ' condnames{1}])
            display([num2str(num_cond2) ' trials in condition ' condnames{2}])
            
            num_xvalid_bins = 10;
            num_runs = num_xvalid_bins;
            
            run_labels1 = repmat([1:num_xvalid_bins],1,ceil(num_cond1/num_xvalid_bins));
            run_labels1 = run_labels1(1:num_cond1); %truncate the length of this vector
            
            run_labels2 = repmat([1:num_xvalid_bins],1,ceil(num_cond2/num_xvalid_bins));
            run_labels2 = run_labels2(1:num_cond2); %truncate the length of this vector
            
            
            shuffled_run_inds_for_cond1_trials = shuffle(run_labels1);
            shuffled_run_inds_for_cond2_trials = shuffle(run_labels2);
            
            
            subj = duplicate_object(subj,'selector','runs','runs_final');
            subj.selectors{end}.mat(cond1_trials)=shuffled_run_inds_for_cond1_trials;
            subj.selectors{end}.mat(cond2_trials)=shuffled_run_inds_for_cond2_trials;
            
            subj.selectors{end}.mat(subj.selectors{end}.mat>num_xvalid_bins) = 1;  %replace any extra xvalid labels with 1
            
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %subj = duplicate_object(subj,'selector','runs','runs_final');
            subj = init_object(subj,'selector','conditions_of_interest_final'); %initialize selector object
            subj = set_mat(subj,'selector','conditions_of_interest_final',actives_selector);
            
            subj = create_xvalid_indices(subj,'runs_final','actives_selname','conditions_of_interest_final');
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % create a new set of xvalid sectors that balance the number of
            % Class A and Class B trials within the training set and within
            % the testing set
            subj = create_balanced_xvalid_selectors(subj,'conds','runs_final_xval');
            
            % First, compile a vector that indexes all training trials (1's in the
            % x-valid selectors)
            
                                    
            
            % For this analysis, we don't care if the testing
            % set data are balanced for the two classes, so
            % restore the testing trials to those from the
            % original xvalid selector set (runs_final_xval)
            
            for kk = 1:num_runs
                unbalanced_testing_set_indices = find(subj.selectors{end-num_runs*2+kk}.mat==2);
                subj.selectors{end-num_runs+kk}.mat(unbalanced_testing_set_indices)=2;
            end
            
            compiled_test_indices = [];
            for jj = 1:num_runs
                compiled_test_indices = [compiled_test_indices find(subj.selectors{end-num_runs+jj}.mat==2)];
            end
            
            
            % creat new actives_selector (critical if subsequently scrambling regressors
            % with JR_scramble_regressors function)
            %                         new_active_trials =[];
            %                         for rr = 1:num_runs
            %                             new_active_trials = horzcat(new_active_trials, find(subj.selectors{end-num_runs+rr}.mat==2));
            %                         end
            %                         new_actives_selector= zeros(1,size(condensed_regs_all,2));
            %                         new_actives_selector(new_active_trials)=1;
            %                             subj = init_object(subj,'selector','conditions_of_interest_final_bal_within_runs'); %initialize selector object
            %                         subj = set_mat(subj,'selector','conditions_of_interest_final_bal_within_runs',new_actives_selector);
            %
            %
            %                         active_trials = new_active_trials;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            %
            if flags.perform_second_round_of_zscoring == 1
                subj.patterns{5}.mat(:,active_trials) = zscore(subj.patterns{5}.mat(:,active_trials)')';
                display('Performing second round of z-scoring')
            end
            %
            if flags.optimize_penalty_param == 1 && x == 0 % find optimal averaged penalty first time only (to save time)
                [subj best_penalties penalty_iteration_results] = optimal_pLR_penalty(subj,'spiral_d_hp_z_condensed','conds','runs_final_xval','runs_final',classifier_mask,'conditions_of_interest_final','use_iteration_perf',false,'perform_final_classification',false);
                class_args.penalty = best_penalties
            end
            
            % [subj] =  JR_scramble_regressors(subj,'conds','runs_final','conditions_of_interest_final_bal_within_runs','conds_scrambled')
            
            
            %%%%%%%%%%%%%%%%%%%%%% RUN THE CLASSIFIER (CROSS-VALIDATION)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            
            for r = 1:number_ROI_files
                
                classifier_mask = subj.masks{r+1}.name;
                
                for p = 1:flags.num_iter_with_same_data
                    x=x+1; % increment results iteration counter
                    
                    [subj results{x}] = cross_validation(subj,'spiral_d_hp_z_condensed','conds','runs_final_xval_bal',classifier_mask,class_args);
                    % [subj results{x}] = cross_validation(subj,'spiral_d_hp_z_condensed','conds_scrambled','runs_final_xval_bal',classifier_mask,class_args);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    % do some important RAM clean-up and data archiving
                    for y = 1:num_runs
                        if flags.generate_importance_maps == 1
                            if strcmp(class_args.train_funct_name,'train_bp')
                                results_IW{x}.iterations(y).scratchpad.net.IW{1} = results{x}.iterations(y).scratchpad.net.IW{1}; % save weights to pass to JR_interpret_weights
                            elseif strcmp(class_args.train_funct_name,'train_pLR')
                                results_IW{x}.iterations(y).scratchpad.net.IW{1} = results{x}.iterations(y).scratchpad.weights';
                            elseif strcmp(class_args.train_funct_name,'train_svdlr')
                                results_IW{x}.iterations(y).scratchpad.net.IW{1} = results{x}.iterations(y).scratchpad.W';
                            else
                                results_IW{x}.iterations(y).scratchpad.net.IW{1} = results{x}.iterations(y).scratchpad.w'; % save weights to pass to JR_interpret_weights
                            end
                        end
                        results{x}.iterations(y).scratchpad.net.inputs{1}.exampleInput=[]; % delete huge data object from results scratchpad to free up RAM
                    end
                    
                    if flags.display_performance_breakdown == 1
                        % analyze the results in more detail
                        correct_vector = [];
                        desireds_vector = [];
                        guesses_vector = [];
                        acts_vector = [];
                        for a = 1:num_runs
                            correct_vector = horzcat(correct_vector,results{x}.iterations(a).perfmet.corrects);
                            desireds_vector = horzcat(desireds_vector,results{x}.iterations(a).perfmet.desireds);
                            guesses_vector = horzcat(guesses_vector,results{x}.iterations(a).perfmet.guesses);
                            acts_vector = horzcat(acts_vector, results{x}.iterations(a).acts(1,:));
                        end
                        
                        overall_accuracy = mean(correct_vector);
                        overall_hit_rate = mean(correct_vector(desireds_vector==1));
                        overall_fa_rate = 1-mean(correct_vector(desireds_vector==2));
                        overall_d_prime = norminv(overall_hit_rate)-norminv(overall_fa_rate);
                        
                        bias_index = count(guesses_vector==1)/length(guesses_vector);
                        if bias_index <= .25 || bias_index>=.75
                            %x = x - 1; % reset x counter and jump back to beginning of FOR loop
                            num_biased_iterations = num_biased_iterations+1
                            %break;
                        end
                        
                        [Y,I] = sort(compiled_test_indices);
                                                
                        data_log{r}.accuracy(k)=overall_accuracy;
                        data_log{r}.acts(k,:)=acts_vector(I); % sort by order of compiled_test_indices
                        data_log{r}.active_trials = active_trials;
                        data_log{r}.class_A_acts(k,:) = data_log{r}.acts(k,desireds_vector(I)==1);
                        data_log{r}.class_B_acts(k,:) = data_log{r}.acts(k,desireds_vector(I)==2);
                        
                    end
                end  
            end
        end
    end
    % avg over all k rounds of classification (each with different training
    % set randomization; script not yet adapted for running multiple ROIs
    R_hits_pmod_vector = mean(data_log{r}.class_A_acts,1);
    HC_hits_pmod_vector = mean(data_log{r}.class_B_acts,1);
    
    cd(behav_dir);
    save R_hippo_evidence_pmod.mat R_hits_pmod_vector;
    save HC_hippo_evidence_pmod.mat HC_hits_pmod_vector;
    cd(expt_dir)
   
    
    results_array(sub_num,1) = mean(data_log{r}.accuracy);
    results_array(sub_num,2) = length(R_hits_pmod_vector);
    results_array(sub_num,3) = length(HC_hits_pmod_vector);
    
    
    time2finish = toc/60;
    display(['Finished ' subj_id ' in ' num2str(time2finish) ' minutes']);
    
    keep b subj_array condition1 condition2 nVox sub_num results_array
    
end
csvwrite([condition1 '_vs_' condition2 '_bilat_hippo_class_accuracy.csv'],results_array);

keyboard;



