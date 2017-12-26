function JR_run_mvpa_across_subj_class_AUC(subj_array,condition1, condition2)


concat.pattern = [];
concat.regressors = [];
concat.actives = [];
concat.subj_id = [];

expt_dir = '/Users/Jesse/fMRI/data/PAST/fMRI';
cd(expt_dir);

% load some .mat files into memory
load([expt_dir '/vol_info.mat']); %get functional data resolution info for spm .img writing



importance_maps_dir=[expt_dir '/mvpa_results/FINAL_ROC_data/across_subj_class/FINAL_IMP_MAPS/pen1/'];
if ~exist(importance_maps_dir,'dir')
    mkdir(importance_maps_dir);
end
xls_results_data_logs_txt_dir=[expt_dir '/mvpa_results/FINAL_ROC_data/across_subj_class/FINAL_CLASS_PERF/pen1/'];
if ~exist(xls_results_data_logs_txt_dir, 'dir')
    mkdir(xls_results_data_logs_txt_dir)
end


xls_results_data_logs_mat_dir=[expt_dir '/mvpa_results/FINAL_ROC_data/across_subj_class/data_logs/pen1/'];
if ~exist(xls_results_data_logs_mat_dir, 'dir')
    mkdir(xls_results_data_logs_mat_dir)
end





for sub_idx = 1:length(subj_array)
    
    subj_id = ['s1' prepend(num2str(subj_array(sub_idx)))]
    
    
    mvpa_workspace = [expt_dir '/' subj_id '/mvpa/' subj_id '_SEPT09_MVPA_MASK_s8mm_wa.mat'];
    load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']); % load in your SPM-formatted onsets file
    
    cd([subj_id '/mvpa'])
    
    
    
    %%%%%%% specify user-defined variables
    %%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set flags (% unless otherwise indicated: 1 = yes please, 0 = no thanks)
    flags.num_full_iter = 1;
    flags.num_results_iter = 1; % number of times to run the entire classification process (select subset of the data and train/test classifier)
    flags.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data
    
    flags.equate_number_of_old_new_trials_per_subjective_bin = 0; % equate number of trials per subjective bin
    flags.equate_number_of_trials_in_cond_1_and_2 = 1; % equate number of trials in conditions 1 and 2 (RECOMMENDED)
    flags.anova_p_thresh = 1;  % p-value threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)
    flags.GLM_anova_p_thresh = 0;
    flags.anova_nVox_thresh = 0; % alternative to specifying p-value threshold; uses top N voxels (0 = DON'T PERFORM ANY FEATURE SELECTION)
    flags.perform_second_round_of_zscoring = 1;  % z-score data again immediately prior to classification (NOT RECOMMENDED)
    flags.remove_artdetect_outliers = 1; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
    flags.artdetect_motion_thresh = 7; % specify ArtDetect bin for motion outliers
    flags.artdetect_global_signal_thresh = 5; % specify ArtDetect bin for global signal outliers
    flags.remove_outlier_trials = 0;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)
    flags.plot_ROC_curve = 1;
    flags.display_performance_breakdown = 1;
    flags.generate_importance_maps = 0;
    flags.write_data_log_to_text_file=1;
    flags.save_data_log_as_mat_file =1;
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
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~exist(mvpa_workspace,'file')
        [subj num_runs num_TP_per_run]= JR_generate_mvpa_workspace_mat_file(subj_id, roi_name, data_imgs_to_use, mvpa_dir); % generate and save workspace
    else
        eval(['load ' mvpa_workspace])  %load workspace
    end
    
    
    
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
    for k = 1:flags.num_full_iter
        
        subj = subj_original;
        condensed_regs_all = condensed_regs_all_original;
        
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
        
        % "activate" only those trials of interest (from regs_of_interest) before creating cross-validation indices
        active_trials = find(sum(condensed_regs_of_interest));
        
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
        
        display([num2str(count(condensed_regs_of_interest(1,active_trials)==1)) ' trials in condition ' condnames{1}])
        display([num2str(count(condensed_regs_of_interest(2,active_trials)==1)) ' trials in condition ' condnames{2}])
        
        actives_selector = zeros(1,size(condensed_regs_all,2)); % intialize vector of all zeros
        actives_selector(active_trials) = 1; % remove all non-"regs_of_interst" trials (set to one)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        subj = init_object(subj,'selector','conditions_of_interest_final'); %initialize selector object
        subj = set_mat(subj,'selector','conditions_of_interest_final',actives_selector);
        
        % zscore temporally-condensed data; active trials only (second round of z-scoring)
        subj.patterns{end}.mat(:,active_trials) = zscore(subj.patterns{end}.mat(:,active_trials)')';
        
        voxel_inds{sub_idx} = find(subj.masks{end}.mat);
        
        if sub_idx == 1
            common_voxel_inds = voxel_inds{1};
            voxels_to_include_current = find(ismember(common_voxel_inds,voxel_inds{sub_idx}));
            voxels_to_include_archive = find(ismember(common_voxel_inds,voxel_inds{sub_idx}));
            concat.pattern = subj.patterns{end}.mat(voxels_to_include_current,active_trials);
            
        else
            archive_voxel_inds = common_voxel_inds;  % keep record common voxel inds up to this point
            common_voxel_inds = intersect(common_voxel_inds,voxel_inds{sub_idx}); % find which voxels of the continuously incrementing array exist for the current subject
            %voxels_to_include_current = find(ismember(common_voxel_inds,voxel_inds{sub_idx}));
            %voxels_to_include_archive = find(ismember(common_voxel_inds,archive_voxel_inds));
            
            voxels_to_include_current = find(ismember(voxel_inds{sub_idx},common_voxel_inds));  % find which of the present subject's voxels are in common with the previous collection
            voxels_to_include_archive = find(ismember(archive_voxel_inds,common_voxel_inds)); % find which of the previous collection of voxels are in common with the current conjunction
            
            display([num2str(length(common_voxel_inds)) ' common voxels remaining in the mask'])
            concat.pattern = horzcat(concat.pattern(voxels_to_include_archive,:),subj.patterns{end}.mat(voxels_to_include_current,active_trials));
        end
        
        concat.subj_id = horzcat(concat.subj_id, repmat(sub_idx,1,length(active_trials)));
        concat.regressors = horzcat(concat.regressors,subj.regressors{1}.mat(:,find(subj.selectors{2}.mat)));  %get regressors only for active trials
        %concat.actives = horzcat(concat.actives, subj.selectors{2}.mat);
        
        number_of_trials(sub_idx)= length(subj.selectors{2}.mat);
        
        clear condensed_regs_of_interest condensed_regs_all condensed_runs
        
        
        cd(expt_dir);
        
    end
end



subj.patterns{end}.mat = concat.pattern;
subj.patterns{end}.matsize = size(concat.pattern);
clear concat.pattern

subj.regressors{1}.mat = concat.regressors;
subj.regressors{1}.matsize = size(concat.regressors);
clear concat.regressors


% subj.selectors{2}.mat = concat.actives;
% subj.selectors{2}.matsize = size(concat.actives);

% trial_counter = 1;
%
% for r = 1:sub_idx  % do for number of subjects
%
%     new_selectors{r} = subj.selectors{3};  % copy selector template
%     new_selectors{r}.name = ['runs_xval_' num2str(r)];
%     new_selectors{r}.mat = ones(1,length(concat.actives)).*concat.actives;
%     new_selectors{r}.mat(1,trial_counter:trial_counter+number_of_trials(r)-1) = new_selectors{r}.mat(1,trial_counter:trial_counter+number_of_trials(r)-1)*2;
%     new_selectors{r}.matsize = size(new_selectors{r}.mat);
%
%     trial_counter = trial_counter+number_of_trials(r);
%
% end
%
% subj.selectors = new_selectors;
% clear new_selectors


subj.selectors{1}.mat = concat.subj_id;

subj = create_xvalid_indices(subj,'runs')

% need to replace original whole brain mask with new reduced sized mask so
% that feature selection doesn't complain

subj.masks{1}.mat = zeros(subj.masks{1}.matsize);
subj.masks{1}.mat(common_voxel_inds) = 1;


% save uber_subj.mat subj

%% RUN CLASSIFICATION BY TRAINING ON ALL BUT ONE SUBJECT (CONCATENATED
%% DATA) AND THEN TEST ON EACH HELD OUT SUBJECT

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

subj.patterns{5}.mat = single(subj.patterns{5}.mat);  % make pattern a 'single' matrix to save ram and speed up classification (Note: doesn't work with backprop, though)

[subj results{1}] = cross_validation(subj,'spiral_d_hp_z_condensed','conds','runs_xval',classifier_mask,class_args)



for i = 1:length(subj_array)
    subj_id = ['s1' prepend(num2str(subj_array(i)))]
    
    subj_labels_subj_i = concat.subj_id(concat.subj_id == i);
    
    correct_vector = results{1}.iterations(i).perfmet.corrects(find(subj_labels_subj_i==i));
    desireds_vector = results{1}.iterations(i).perfmet.desireds(find(subj_labels_subj_i==i));
    acts_diff_vector = results{1}.iterations(i).acts(1,find(subj_labels_subj_i==i)) - results{1}.iterations(i).acts(2,find(subj_labels_subj_i==i));
    
    [abs_sorted_diffs abs_ind] = sort(abs(acts_diff_vector),2,'descend');
    abs_correct_sorted = correct_vector(abs_ind);
    num_trials = length(abs_correct_sorted);
    
    group_trained_perf_by_quartile_rank(1,i) = mean(abs_correct_sorted(1:ceil(num_trials*1.0)));
    group_trained_perf_by_quartile_rank(2,i) = mean(abs_correct_sorted(1:ceil(num_trials*.75)));
    group_trained_perf_by_quartile_rank(3,i) = mean(abs_correct_sorted(1:ceil(num_trials*.50)));
    group_trained_perf_by_quartile_rank(4,i) = mean(abs_correct_sorted(1:ceil(num_trials*.25)));
    
    
    bin_intervals =1:-.05:.05;
    for acc_bin = 1:20
        acc_percentiles(acc_bin)= mean(abs_correct_sorted(1:ceil(num_trials*bin_intervals(acc_bin))));
    end
    
    % print the top 100%, 75%, 50%, and 25% accuracy to the screen
    display(acc_percentiles([1 6 11 16]))
    
    
    
    if flags.plot_ROC_curve == 1
        
        % sort by signed classifier "confidence" (for ROI curves)
        [sorted_diffs ind] = sort(acts_diff_vector,2,'descend');
        correct_sorted = correct_vector(ind);
        desireds_sorted = desireds_vector(ind);
        
        % create continuous ROC function
        hit_rate = []; fa_rate = [];
        for j = 1:length(sorted_diffs);
            hit_rate(j) = length(correct_sorted(intersect(find(desireds_sorted == 1),[1:j]))) / length(find(desireds_sorted == 1));
            fa_rate(j) = length(correct_sorted(intersect(find(desireds_sorted == 2),[1:j]))) / length(find(desireds_sorted == 2));
        end
        
        %                     figure
        %                     plot(fa_rate,hit_rate,'.-')
        %                     hold on
        %                     plot([0 1],[0 1],'r')
        %                     xlabel('P(Old|New)')
        %                     ylabel('P(Old|Old)')
        
        auc_overall = auroc(hit_rate',fa_rate')
        
        % create ROC function with 80 bins, based on
        roc_bin_intervals = .975:-.025:-1;
        for bin_num = 1:80
            hits_80(bin_num)=length(correct_sorted(intersect(find(desireds_sorted == 1),find(sorted_diffs>roc_bin_intervals(bin_num))))) / length(find(desireds_sorted == 1));
            fas_80(bin_num)=length(correct_sorted(intersect(find(desireds_sorted == 2),find(sorted_diffs>roc_bin_intervals(bin_num))))) / length(find(desireds_sorted == 2));
        end
        auc_80_bins = auroc(hits_80',fas_80');
        
        
        
    end
    
    if flags.write_data_log_to_text_file==1
        
        %         data_log.overall_acc(i)=overall_accuracy;
        %         data_log.hits(i)=overall_hit_rate;
        %         data_log.FAs(i)=overall_fa_rate;
        %         data_log.d_prime(i)=overall_d_prime;
        %data_log.classification_accuracy_by_resp(i,:)=classification_accuracy_by_resp;
        %data_log.number_trials_per_bin(i,:)=number_of_trials_per_bin;
        %data_log.acc_sorted_by_classifier_confidence(i,:)=acc_sorted_by_classifier_confidence;
        data_log.acc_percentiles(i,:) = acc_percentiles;
        data_log.penalty_param(i) = class_args.penalty;
        %data_log.class_counts_by_quartile_rank(:,:,x) = class_counts_by_quartile_rank;
        %data_log.number_of_trials_per_bin_by_quartile_rank(:,:,x) = number_of_trials_per_bin_by_quartile_rank;
        data_log.auc_overall(i) = auc_overall;
        data_log.auc_80_bins(i) = auc_80_bins;
        %         data_log.roc_continuous_hits(i,:)= hit_rate;
        %         data_log.roc_continuous_fas(i,:)= fa_rate;
        data_log.roc_80_bin_hits(i,:)= hits_80;
        data_log.roc_80_bin_fas(i,:)= fas_80;
        data_log.subj_id{i}=subj_id;
    end
    
    
end
if flags.save_data_log_as_mat_file ==1;
    save_cmd = ['save ' xls_results_data_logs_mat_dir '/'  condnames{1} '_vs_' condnames{2} '.mat data_log flags'];
    eval(save_cmd);
end

if flags.write_data_log_to_text_file==1
    
    filename= [xls_results_data_logs_txt_dir '/' condnames{1} '_vs_' condnames{2} '.txt'];
    fid=fopen(filename, 'wt');
    
    for q=1:i
        fprintf(fid, '%s\t',data_log.subj_id{q});
        %fprintf(fid, '%4.4f\t', data_log.overall_acc(q));
        %fprintf(fid, '%4.4f\t', data_log.hits(q));
        %fprintf(fid, '%4.4f\t', data_log.FAs(q));
        %fprintf(fid, '%4.4f\t', data_log.d_prime(q));
        %fprintf(fid, '%4.4f\t', data_log.classification_accuracy_by_resp(q,:)); %[data_log.classification_accuracy_by_resp(1,q) '\t' data_log.classification_accuracy_by_resp(q,1) '\t' data_log.classification_accuracy_by_resp(q,1) '\t' data_log.classification_accuracy_by_resp(q,4) '\t' data_log.classification_accuracy_by_resp(q,5) '\t' data_log.classification_accuracy_by_resp(q,6) '\t' data_log.classification_accuracy_by_resp(q,7) '\t' data_log.classification_accuracy_by_resp(q,8) '\t' data_log.classification_accuracy_by_resp(q,9) '\t' data_log.classification_accuracy_by_resp(q,10) '\t']);
        %fprintf(fid, '%4.4f\t', data_log.acc_sorted_by_classifier_confidence(q,:));
        fprintf(fid, '%4.4f\t', data_log.auc_overall(q));
        fprintf(fid, '\t');
        fprintf(fid, '%4.4f\t', data_log.acc_percentiles(q,:));
        fprintf(fid, '\t');
        fprintf(fid, '%4.4f\t', data_log.roc_80_bin_hits(q,:));
        fprintf(fid, '\t');
        fprintf(fid, '%4.4f\t', data_log.roc_80_bin_fas(q,:));
        
        %fprintf(fid, '%3.0f\t', reshape(data_log.class_counts_by_quartile_rank(:,:,q)',1,8));
        fprintf(fid, '\n');
        %for p=1:size(data_log.acc_sorted_by_classifier_confidence, 1)
        %    fprintf(fid, '%4.4f\t', data_log.acc_sorted_by_classifier_confidence(q,p));
        %end
        
    end
    fprintf(fid, '%s\t', 'mean');
    %fprintf(fid, '%4.4f\t', mean(data_log.overall_acc));
    %fprintf(fid, '%4.4f\t', mean(data_log.hits));
    %fprintf(fid, '%4.4f\t', mean(data_log.FAs));
    %fprintf(fid, '%4.4f\t', mean(data_log.d_prime));
    %fprintf(fid, '%4.4f\t', mean(data_log.classification_accuracy_by_resp,1));
    %fprintf(fid, '%4.4f\t', mean(data_log.acc_sorted_by_classifier_confidence,1));
    fprintf(fid, '%4.4f\t', mean(data_log.auc_overall));
    fprintf(fid, '\t');
    fprintf(fid, '%4.4f\t', mean(data_log.acc_percentiles,1));
    fprintf(fid, '\t');
    fprintf(fid, '%4.4f\t', mean(data_log.roc_80_bin_hits,1));
    fprintf(fid, '\t');
    fprintf(fid, '%4.4f\t', mean(data_log.roc_80_bin_fas,1));
    fprintf(fid, '\n');
    fclose(fid);
end















%
%
%
%
%     %% DO IT AGAIN WITH THE FLIPPED SELECTORS TO TRAIN ON ONE SUBJECT AND TEST
%     %% ON EACH INDIVIDUAL SUBJECT
%
%     subj = duplicate_object(subj,'selector','runs','runs_flipped');
%     subj = create_xvalid_indices_flip_train_test(subj,'runs_flipped');
%
%     subj = duplicate_object(subj,'pattern','spiral_d_hp_z_condensed','spiral_d_hp_z_condensed_copy');
%
%
%     % run feature selection ANOVA: specify pvalue (if desired)
%     statmap_arg.use_mvpa_ver = true;
%     if flags.anova_p_thresh ~= 1
%
%         subj = JR_feature_select(subj,'spiral_d_hp_z_condensed_copy','conds','runs_flipped_xval','thresh',flags.anova_p_thresh, 'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
%         %subj = JR_feature_select(subj,'spiral_d_hp_z_condensed','conds','runs_xval','thresh',flags.anova_p_thresh, 'greater_than', true,'statmap_funct','statmap_RFE');
%         classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
%     else
%         classifier_mask = subj.masks{1}.name; % use original mask
%     end
%
%     % run feature selection ANOVA: specify #of voxels (if desired)
%     if flags.anova_nVox_thresh ~=0
%
%         subj = JR_feature_select_top_N_vox(subj,'spiral_d_hp_z_condensed_copy','conds','runs_flipped_xval','nVox_thresh',flags.anova_nVox_thresh,'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
%         %subj = JR_feature_select_top_N_vox_iterative(subj,'spiral_d_hp_z_condensed','conds','runs_xval','nVox_thresh',flags.anova_nVox_thresh,'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
%         classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
%     else
%         classifier_mask = subj.masks{1}.name; % use original mask
%     end
%
%     subj.patterns{5}.mat = single(subj.patterns{5}.mat);  % make pattern a 'single' matrix to save ram and speed up classification (Note: doesn't work with backprop, though)
%
%     [subj results{2}] = cross_validation(subj,'spiral_d_hp_z_condensed_copy','conds','runs_flipped_xval',classifier_mask,class_args)
%
%
%
%     for i = 1:length(subj_array)
%
%         subj_labels_without_subj_i = concat.subj_id(concat.subj_id ~= i);
%
%         for j = 1:length(subj_array)
%
%             correct_vector = results{2}.iterations(i).perfmet.corrects(find(subj_labels_without_subj_i==j));
%             acts_diff_vector = results{2}.iterations(i).acts(1,find(subj_labels_without_subj_i==j)) - results{2}.iterations(i).acts(2,find(subj_labels_without_subj_i==j));
%
%             [abs_sorted_diffs abs_ind] = sort(abs(acts_diff_vector),2,'descend');
%             abs_correct_sorted = correct_vector(abs_ind);
%             num_trials = length(abs_correct_sorted);
%
%             top100_pct_perf_matrix(i,j) = mean(correct_vector);
%
%             top10_pct_perf_matrix(i,j)=mean(abs_correct_sorted(1:ceil(num_trials*.10)));
%             top25_pct_perf_matrix(i,j)=mean(abs_correct_sorted(1:ceil(num_trials*.25)));
%             top50_pct_perf_matrix(i,j)=mean(abs_correct_sorted(1:ceil(num_trials*.50)));
%             top75_pct_perf_matrix(i,j)=mean(abs_correct_sorted(1:ceil(num_trials*.75)));
%
%         end
%     end
%
%     keyboard;











