function JR_mvpa_compile_indiv_subj_data_for_across_subj_class(subj_array)



  expt_dir = '/Users/Jesse/fMRI/data/PAST/fMRI';
  cd(expt_dir);
  

 concat.pattern = [];
    concat.regressors = [];
    concat.actives = [];
    concat.subj_id = [];

    
for sub_idx = 1:length(subj_array)

    subj_id = ['s1' prepend(num2str(subj_array(sub_idx)))]
    
       
    mvpa_workspace = [expt_dir '/' subj_id '/mvpa/' subj_id '_NEW_MVPA_MASK_s8mm_wa.mat'];
    
    cd([subj_id '/mvpa'])



    %%%%%%% specify user-defined variables
    %%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    flags.save_workspace = 0; % 1 = yes please, 0 = no thanks

    % unless a previously created workspace is specified, load in the preprocessed data
    if ~exist('mvpa_workspace')

        exp_name = 'PAST';
        roi_file = '/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/merged_AAL_ROIs.nii';
        roi_name = 'merged_AAL_ROIs'
        num_TP_per_run = 203;

        % load user-created filename and onsets lists into workspace
        load('raw_filenames_s8mm_wa.mat') %loads predefined cell array called raw_filenames into memory
        num_runs = length(raw_filenames)/num_TP_per_run; %calculate the number of runs (allows script to work flexibly for subjects with missing runs)
        load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']);
        vol_info = spm_vol(raw_filenames{1}); %get functional data resolution info for spm .img writing
        [subj] = JR_mvpa_load_and_preprocess_raw_data(subj_id, exp_name, roi_name, roi_file, raw_filenames, num_runs, num_TP_per_run);

        if flags.save_workspace == 1
            save_cmd = ['save ' subj_id '_' roi_name '_8mm_smoothing.mat'];
            eval(save_cmd);
        end
    else
        eval(['load ' mvpa_workspace])
    end

     load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']); % load in your SPM-formatted onsets file

    % Set flags (% 1 = yes please, 0 = no thanks)
    flags.equate_number_of_old_new_trials_per_subjective_bin = 0; % balance the number of trials for each subjective response
    flags.equate_number_of_trials_in_cond_1_and_2 = 1; % balance the number of trials in conditions 1 & 2 (this is a good thing and prevents biasing the classifier to pi
    flags.plot_mean_timecourses = 0;
    flags.plot_ROC_curve = 0;
    flags.display_performance_breakdown = 0;
    flags.generate_importance_maps = 0;
    flags.anova_p_thresh = 1;
    flags.anova_nVox_thresh = 0;
     flags.remove_artdetect_outliers = 1; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
    flags.artdetect_motion_thresh = 7;
    flags.artdetect_global_signal_thresh = 5;


    % specify which conditions to use for classification (must correspond to the names of conditions specified below)
    condnames =  {'Hits','CRs'};

    TRs_to_average_over = [3 4]; %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier

    %num_results_iter = 2; % number of times to run the cross validation process

    %anova_p_thresh = .05;  %p threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)


    % classifier parameters
%     class_args.train_funct_name = 'train_bp';
%     class_args.test_funct_name = 'test_bp';
%     class_args.nHidden = 0;
    
    class_args.train_funct_name = 'train_pLR';
    class_args.test_funct_name = 'test_pLR';
    class_args.penalty = 10000;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Extract info about conditions from onsets file
    num_conds = size(onsets,2);
    all_regs = zeros(num_conds,num_runs*num_TP_per_run); % initialize regs matrix as conditions x timepoints

    for cond = 1: num_conds %(exclude last condition ("no_response" trials)
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
    if ~isempty(find(strcmp(subj_id,{'s104','s105','s108','s113','s117','119','120','121','203','204','205','207','208'})))
        all_regs(6,:)=0;
    end
    
    % SPECIAL FIX because reg #7 (NEW_HC_old) is a fake placeholder trial for some subjects
    if ~isempty(find(strcmp(subj_id,{'s117','207'})))
        all_regs(7,:)=0;
    end
    
    % SPECIAL FIX because reg #5 (OLD_HC_new) is a fake placeholder trial for some subjects
    if ~isempty(find(strcmp(subj_id,{'118','205','207'})))
        all_regs(5,:)=0;
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

    % save final pattern in single precision form (8 sig figs) to save RAM and HD space
    %temporally_condensed_data = single(temporally_condensed_data); %

    clear data_by_TR; %clean up matlab workspace to save memory
    
    
        if flags.remove_artdetect_outliers == 1
            load([expt_dir '/outlier_indices/' subj_id '_outlier_indices']); %load outlier indices

            m_outliers = movement_outlier_trials{flags.artdetect_motion_thresh};  % remove trials with more than .35mm/TR of movement
            gs_outliers = global_signal_outlier_trials{flags.artdetect_global_signal_thresh}; % remove trials with global signal change of +/- 3.5 SD from mean
            combined_outliers = union(m_outliers,gs_outliers);

            %temporally_condensed_data(:,combined_outliers)=[]; % remove outlier trials from fmri data
            %condensed_regs_of_interest(:,combined_outliers) = [];
            condensed_regs_all(:,combined_outliers) = 0;
            %condensed_runs(combined_outliers) = [];

            display([num2str(length(m_outliers)) ' movement outlier trials flagged']);
            display([num2str(length(gs_outliers)) ' global signal outlier trials flagged']);
            display([num2str(length(combined_outliers)) ' total outlier trials excluded']);
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

    nonactive_trials = find(sum(condensed_regs_of_interest)==0);
    active_trials = find(sum(condensed_regs_of_interest));
    actives_selector = ones(1,size(condensed_regs_all,2)); % intialize vector of all ones
    actives_selector(nonactive_trials) = 0; % remove all non-"regs_of_interst" trials (set to zero)
    subj = init_object(subj,'selector','conditions_of_interest'); %initialize selector object
    subj = set_mat(subj,'selector','conditions_of_interest',actives_selector);
    %subj = create_xvalid_indices(subj,'runs','actives_selname','conditions_of_interest');
    
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
        common_voxel_inds = intersect(common_voxel_inds,voxel_inds{sub_idx}); % find which voxels of the continuously incrementing array existe for the current subject
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
    
        subj_labels_subj_i = concat.subj_id(concat.subj_id == i);  
 
        correct_vector = results{1}.iterations(i).perfmet.corrects(find(subj_labels_subj_i==i));
        acts_diff_vector = results{1}.iterations(i).acts(1,find(subj_labels_subj_i==i)) - results{1}.iterations(i).acts(2,find(subj_labels_subj_i==i));
        
        [abs_sorted_diffs abs_ind] = sort(abs(acts_diff_vector),2,'descend');
        abs_correct_sorted = correct_vector(abs_ind);
        num_trials = length(abs_correct_sorted);
        
        group_trained_perf_by_quartile_rank(1,i) = mean(abs_correct_sorted(1:ceil(num_trials*1.0)));
        group_trained_perf_by_quartile_rank(2,i) = mean(abs_correct_sorted(1:ceil(num_trials*.75)));
        group_trained_perf_by_quartile_rank(3,i) = mean(abs_correct_sorted(1:ceil(num_trials*.50)));
        group_trained_perf_by_quartile_rank(4,i) = mean(abs_correct_sorted(1:ceil(num_trials*.25)));         
end





%% DO IT AGAIN WITH THE FLIPPED SELECTORS TO TRAIN ON ONE SUBJECT AND TEST
%% ON EACH INDIVIDUAL SUBJECT

subj = duplicate_object(subj,'selector','runs','runs_flipped');
subj = create_xvalid_indices_flip_train_test(subj,'runs_flipped');

subj = duplicate_object(subj,'pattern','spiral_d_hp_z_condensed','spiral_d_hp_z_condensed_copy');


% run feature selection ANOVA: specify pvalue (if desired)
        statmap_arg.use_mvpa_ver = true;
        if flags.anova_p_thresh ~= 1

            subj = JR_feature_select(subj,'spiral_d_hp_z_condensed_copy','conds','runs_flipped_xval','thresh',flags.anova_p_thresh, 'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            %subj = JR_feature_select(subj,'spiral_d_hp_z_condensed','conds','runs_xval','thresh',flags.anova_p_thresh, 'greater_than', true,'statmap_funct','statmap_RFE');
            classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
        else
            classifier_mask = subj.masks{1}.name; % use original mask
        end

        % run feature selection ANOVA: specify #of voxels (if desired)
        if flags.anova_nVox_thresh ~=0

            subj = JR_feature_select_top_N_vox(subj,'spiral_d_hp_z_condensed_copy','conds','runs_flipped_xval','nVox_thresh',flags.anova_nVox_thresh,'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            %subj = JR_feature_select_top_N_vox_iterative(subj,'spiral_d_hp_z_condensed','conds','runs_xval','nVox_thresh',flags.anova_nVox_thresh,'statmap_funct','statmap_anova','statmap_arg',statmap_arg);
            classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
        else
            classifier_mask = subj.masks{1}.name; % use original mask
        end
        
subj.patterns{5}.mat = single(subj.patterns{5}.mat);  % make pattern a 'single' matrix to save ram and speed up classification (Note: doesn't work with backprop, though)

[subj results{2}] = cross_validation(subj,'spiral_d_hp_z_condensed_copy','conds','runs_flipped_xval',classifier_mask,class_args)



for i = 1:length(subj_array)
    
    subj_labels_without_subj_i = concat.subj_id(concat.subj_id ~= i);
    
    for j = 1:length(subj_array)
        
        correct_vector = results{2}.iterations(i).perfmet.corrects(find(subj_labels_without_subj_i==j));
        acts_diff_vector = results{2}.iterations(i).acts(1,find(subj_labels_without_subj_i==j)) - results{2}.iterations(i).acts(2,find(subj_labels_without_subj_i==j));
        
        [abs_sorted_diffs abs_ind] = sort(abs(acts_diff_vector),2,'descend');
        abs_correct_sorted = correct_vector(abs_ind);
        num_trials = length(abs_correct_sorted);
        
        top100_pct_perf_matrix(i,j) = mean(correct_vector);
        
        top10_pct_perf_matrix(i,j)=mean(abs_correct_sorted(1:ceil(num_trials*.10)));
        top25_pct_perf_matrix(i,j)=mean(abs_correct_sorted(1:ceil(num_trials*.25)));
        top50_pct_perf_matrix(i,j)=mean(abs_correct_sorted(1:ceil(num_trials*.50)));
        top75_pct_perf_matrix(i,j)=mean(abs_correct_sorted(1:ceil(num_trials*.75)));       
        
     end
end

keyboard;






    




