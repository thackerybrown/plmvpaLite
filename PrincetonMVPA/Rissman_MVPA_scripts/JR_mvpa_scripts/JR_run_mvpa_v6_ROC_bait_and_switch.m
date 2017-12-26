function []= JR_run_mvpa_v6_ROC_bait_and_switch(subj_array, condition1, condition2, balance_per_subj_bin);

%%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for b=subj_array
    tic
    if b<10
        subj_id=strcat('s10', num2str(b))
    else
        subj_id=strcat('s1', num2str(b))
    end

    if balance_per_subj_bin==1
        balanced_or_unbal='balanced_bins';
    else
        balanced_or_unbal='unbalanced_bins';
    end

    PAST_dir = '/Users/Jesse/fMRI/data/PAST/fMRI';
    mvpa_dir = [PAST_dir '/' subj_id '/mvpa'];
    importance_maps_dir=[PAST_dir '/mvpa_results/improved_importance_maps/' condition1 '_vs_' condition2 '_' balanced_or_unbal '_FIXED_AAL_4mm'];
    xls_results_data_logs_dir=[PAST_dir '/mvpa_results/xls_results_data_logs/' condition1 '_vs_' condition2 '_' balanced_or_unbal '_FIXED_AAL_4mm'];

    if ~exist(importance_maps_dir,'dir')
        mkdir(importance_maps_dir);
    end
    if ~exist(xls_results_data_logs_dir, 'dir')
        mkdir(xls_results_data_logs_dir)
    end

    load([PAST_dir '/vol_info.mat']); %get functional data resolution info for spm .img writing
    load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']);

    roi_name = 'merged_AAL_ROIs_FIXED_HOLES';
    data_imgs_to_use = 'raw_filenames_s4mm_wa.mat';
    %data_imgs_to_use = 'raw_filenames_wa.mat';

    mvpa_workspace = [PAST_dir '/' subj_id '/mvpa/' subj_id '_merged_AAL_ROIs_FIXED_HOLES_wa.m.mat'];

    num_results_iter = 1; % number of times to run the cross validation process
    num_iterations_with_same_classifier = 1;

    %anova_nVox_thresh = 1000;

    condnames_bait =  {'R_hits', 'HC_hits'};
    condnames_switch = {'Objective_old','Objective_new'};




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set flags (% 1 = yes please, 0 = no thanks)
    flags.equate_number_of_old_new_trials_per_subjective_bin = balance_per_subj_bin; % equate_per_subjective_bin;
    flags.equate_number_of_trials_in_cond_1_and_2 = 1;
    flags.anova_p_thresh = 1;  %p threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)
    flags.perform_second_round_of_zscoring = 0;
    flags.remove_artdetect_outliers = 1; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
    flags.artdetect_motion_thresh = 4;
    flags.artdetect_global_signal_thresh = 4;
    flags.remove_outlier_trials = 4;  % how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)
    flags.use_AAL_ROIs = 0; % 1 = use specified AAL ROIs from composite image; 0 = use individually-specified ROIs (specify below)
    flags.plot_ROC_curve = 1;
    flags.display_performance_breakdown = 1;
    flags.generate_importance_maps = 0;
    flags.write_data_log_to_text_file=0;

    % specify which conditions to use for classification (must correspond to the names of conditions specified below)
    condnames =  {condition1, condition2};

    TRs_to_average_over = [1 2 3 4 5]; %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier
    TR_weights = [0 0 .5 .5 0];
    %TR_weights = [.0072 .2168 .3781 .2742 .1237];  % from SPM canonical values at 1,3,5,7,and 9 sec post-stimulus

    % classifier parameters
    class_args.train_funct_name = 'train_logreg';
    class_args.test_funct_name = 'test_logreg';
    class_args.penalty = 40;
    class_args.scale_penalty = false;

    %     class_args.train_funct_name = 'train_bp';
    %     class_args.test_funct_name = 'test_bp';
    %     class_args.nHidden = 0;



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~exist('mvpa_workspace','var')
        [subj num_runs num_TP_per_run]= JR_generate_mvpa_workspace_mat_file(subj_id, roi_name, data_imgs_to_use, mvpa_dir); % generate and save workspace
    else
        eval(['load ' mvpa_workspace])  %load workspace
    end


    %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %Load ROI file names and store in ROI_file_names%

    if flags.use_AAL_ROIs == 1

        roi_img_with_labels = '/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/spm5_AAL_ROIs/rAAL_cluster_labeled.nii'; % for writing classification performance results to .img
        roi_folder_name = '/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/spm5_AAL_ROIs/';
        load([roi_folder_name 'AAL_ROI_names.mat']);
        ROI_file_names=filename;
        number_ROI_files=size(ROI_file_names,2);
        %ROIs_to_use = [59 61 3 39 40]
        ROIs_to_use = [59 60 61 62 3 4 39 40 35 36]
        %ROIs_to_use = [59 61 3 77 39]

        %%%Add each ROI file specified in ROI_file_names to subj structure%%%
        for h=1:length(ROIs_to_use)
            i = ROIs_to_use(h);
            current_mask_name = ROI_file_names{i};
            current_mask_file = strcat(roi_folder_name, current_mask_name);
            subj = load_spm_mask(subj, current_mask_name, current_mask_file);
        end

    else

        ROI_path = '/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/from_ben/';
        ROIs_to_use = {'L_ANG','L_IPS'};

        for h=1:length(ROIs_to_use)
            current_mask_name = ROIs_to_use{h};
            current_mask_file = [ROI_path ROIs_to_use{h} '.img'];
            subj = load_spm_mask(subj, current_mask_name, current_mask_file)
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subj_orig = subj; % save copy of original subj struct

    for r=2:length(ROIs_to_use)+1 % add one b/c whole brain mask is the first mask

        classifier_mask = subj.masks{r}.name;
        ROC{r}.mask = classifier_mask;

        x = 0;
        for n = 1: num_results_iter

            subj = subj_orig; % overwrite subj struct w/ original

            % Extract info about conditions from onsets file
            num_conds = size(onsets,2);
            all_regs = zeros(num_conds,num_runs*num_TP_per_run); % initialize regs matrix as conditions x timepoints

            for cond = 1: num_conds-1 %(exclude last condition ("no_response" trials)
                for trial = 1: length(onsets{cond})
                    time_idx = onsets{cond}(trial)/2+1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
                    all_regs(cond,time_idx) = 1;
                end
            end

            % SPECIAL FIX FOR because reg #6 (NEW_recollect) is a fake
            % placeholder trial for some subjects
            if ~isempty(find(strcmp(subj_id,{'s104','s105','s108','s113','s117'})))
                all_regs(6,:)=0;
            end

            % SPECIAL FIX FOR because reg #7 (NEW_HC_old) is a fake
            % placeholder trial for some subjects
            if ~isempty(find(strcmp(subj_id,{'s117'})))
                all_regs(7,:)=0;
            end




            % condense regs by removing zeros
            condensed_regs_all = [];
            condensed_runs = [];
            condensed_regs_bait = [];
            condensed_regs_switch = [];

            trial_counter = 1;
            for i = 1: size(all_regs,2)
                if ~isempty(find(all_regs(:,i))) % if not a rest trial
                    condensed_regs_bait(:,trial_counter) = regs_bait(:,i);
                    condensed_regs_switch(:,trial_counter) = regs_switch(:,i);
                    condensed_regs_all(:,trial_counter) = all_regs(:,i); %save
                    condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
                    trial_counter = trial_counter + 1;
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % select TRs of interest (to correspond with peak post-stim BOLD response)

            all_trials = sum(all_regs,1); % vector of all trial   (patterns{4} contains fully preprocessed data)
            data_by_TR(1,:,:) = TR_weights(1)*subj.patterns{end}.mat(:,find(all_trials)+0); % 1st TR (0-2 sec)
            data_by_TR(2,:,:) = TR_weights(2)*subj.patterns{end}.mat(:,find(all_trials)+1); % 2nd TR (2-4 sec)
            data_by_TR(3,:,:) = TR_weights(3)*subj.patterns{end}.mat(:,find(all_trials)+2); % 3rd TR (4-6 sec)
            data_by_TR(4,:,:) = TR_weights(4)*subj.patterns{end}.mat(:,find(all_trials)+3); % 4th TR (6-8 sec)
            data_by_TR(5,:,:) = TR_weights(5)*subj.patterns{end}.mat(:,find(all_trials)+4); % 5th TR (8-10 sec)
            temporally_condensed_data = squeeze(sum(data_by_TR(TRs_to_average_over,:,:),1));

            clear data_by_TR; %clean up matlab workspace to save memory
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Exclude trials determined to be outliers by ArtDetect script
            % Guide to outlier file cell arrays...
            % Movement thresholds: .2 .25 .3 .35 .4 .4 .5
            % Global signal thresholds: 2 2.5 3 3.5 4 4.5 5

            if flags.remove_artdetect_outliers == 1
                load([PAST_dir '/outlier_indices/' subj_id '_outlier_indices']); %load outlier indices

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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

            %no_resp = condensed_regs_all(11,:); %excluded from analysis

            %assign conditions to train/test classifier on
            condensed_regs_of_interest = [];
            eval(['condensed_regs_of_interest(1,:) = ' condnames{1} ';'])
            eval(['condensed_regs_of_interest(2,:) = ' condnames{2} ';'])

            if flags.equate_number_of_trials_in_cond_1_and_2 == 1

                cond1_trials = find(condensed_regs_of_interest(1,:));
                cond2_trials = find(condensed_regs_of_interest(2,:));
                num_cond1 = length(cond1_trials);
                num_cond2 = length(cond2_trials);

                if num_cond1 > num_cond2
                    rand_array = rand(1,num_cond1);
                    [sorted inds]= sort(rand_array);
                    trials_to_cut = cond1_trials(inds(1:num_cond1-num_cond2));
                    condensed_regs_of_interest(1,trials_to_cut) = 0;
                    display([num2str(length(trials_to_cut)) ' trials cut from ' condnames{1}]);
                elseif num_cond1 < num_cond2
                    rand_array = rand(1,num_cond2);
                    [sorted inds]= sort(rand_array);
                    trials_to_cut = cond2_trials(inds(1:num_cond2-num_cond1));
                    condensed_regs_of_interest(2,trials_to_cut) = 0;
                    display([num2str(length(trials_to_cut)) ' trials cut from ' condnames{2}]);
                else
                    display('Trial numbers are already balanced');
                end
            end

            display([num2str(count(condensed_regs_of_interest(1,:)==1)) ' trials in condition ' condnames{1}])
            display([num2str(count(condensed_regs_of_interest(2,:)==1)) ' trials in condition ' condnames{2}])



















            % initialize regressors object (bait round)
            subj = init_object(subj,'regressors','conds_bait');
            subj = set_mat(subj,'regressors','conds_bait',condensed_regs_bait);
            subj = set_objfield(subj,'regressors','conds_bait','condnames',condnames_bait);

            % initialize regressors object (switch round)
            subj = init_object(subj,'regressors','conds_switch');
            subj = set_mat(subj,'regressors','conds_switch',condensed_regs_switch);
            subj = set_objfield(subj,'regressors','conds_switch','condnames',condnames_switch);



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

            % "activate" only those trials of interest (bait) before
            % creating cross-validation indices

            nonactive_trials = find(sum(condensed_regs_bait)==0);
            active_trials = find(sum(condensed_regs_bait));
            actives_selector = ones(1,size(condensed_regs_bait,2)); % intialize vector of all ones
            actives_selector(nonactive_trials) = 0; % remove all non-"regs_of_interst" trials (set to zero)
            subj = init_object(subj,'selector','conditions_of_interest_bait'); %initialize selector object
            subj = set_mat(subj,'selector','conditions_of_interest_bait',actives_selector);
            subj = create_xvalid_indices(subj,'runs','actives_selname','conditions_of_interest_bait');

            % "activate" only those trials of interest (switch) before
            % creating cross-validation indices

            nonactive_trials = find(sum(condensed_regs_switch)==0);
            active_trials = find(sum(condensed_regs_switch));
            actives_selector = ones(1,size(condensed_regs_switch,2)); % intialize vector of all ones
            actives_selector(nonactive_trials) = 0; % remove all non-"regs_of_interst" trials (set to zero)
            subj = init_object(subj,'selector','conditions_of_interest_switch'); %initialize selector object
            subj = set_mat(subj,'selector','conditions_of_interest_switch',actives_selector);

            subj = duplicate_object(subj,'selector','runs','runs_switch');
            subj = create_xvalid_indices(subj,'runs_switch','actives_selname','conditions_of_interest_switch');



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
            subj = create_xvalid_indices(subj,'runs','actives_selname','conditions_of_interest');

            if flags.perform_second_round_of_zscoring == 1
                % zscore temporally-condensed data; active trials only (second round of z-scoring)
                %             for r=1:length(active_trials)
                %                 subj.patterns{end}.mat(:,active_trials(r)) = subj.patterns{end}.mat(:,active_trials(r))-mean(subj.patterns{end}.mat(:,active_trials),2);
                %             end
                subj.patterns{end}.mat(:,active_trials) = zscore(subj.patterns{end}.mat(:,active_trials)')';
                display('Performing second round of z-scoring')
            end


            %             % run feature selection ANOVA: specify pvalue (if desired)
            %             if flags.anova_p_thresh ~= 1
            %                 subj = feature_select(subj,'spiral_d_hp_z_condensed','conds','runs_xval','thresh',flags.anova_p_thresh);
            %                 classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
            %             else
            %                 classifier_mask = subj.masks{1}.name; % use original mask
            %             end
            %
            %             % run feature selection ANOVA: specify #of voxels (if desired)
            %             if exist('anova_nVox_thresh','var')
            %                 subj = JR_feature_select_top_N_vox(subj,'spiral_d_hp_z_condensed','conds','runs_xval','nVox_thresh',anova_nVox_thresh);
            %                 classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
            %             else
            %                 classifier_mask = subj.masks{1}.name; % use original mask
            %             end


            
            
            %%%%%%%%%%%%%%%%%%%%%% RUN THE CLASSIFIER (CROSS-VALIDATION)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            correct_vector = [];
            desireds_vector = [];
            guesses_vector = [];
            acts_diff_vector = [];
            compiled_acts_diff = [];
            compiled_guesses_vector = [];
            compiled_correct_vector = [];
            compiled_desireds_vector = [];


            for p = 1:num_iterations_with_same_classifier
                x=x+1;

                [subj results{x}] = cross_validation(subj,'spiral_d_hp_z_condensed','conds','runs_xval',classifier_mask,class_args);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % do some important RAM clean-up and data archiving
                for y = 1:num_runs
                    if flags.generate_importance_maps == 1
                        results_IW{x}.iterations(y).scratchpad.net.IW{1} = results{x}.iterations(y).scratchpad.net.IW{1}; % save weights to pass to JR_interpret_weights
                    end
                    results{x}.iterations(y).scratchpad.net.inputs{1}.exampleInput=[]; % delete huge data object from results scratchpad to free up RAM
                end

                if flags.display_performance_breakdown == 1
                    % analyze the results in more detail

                    for a = 1:num_runs
                        correct_vector = horzcat(correct_vector,results{x}.iterations(a).perfmet.corrects);
                        desireds_vector = horzcat(desireds_vector,results{x}.iterations(a).perfmet.desireds);
                        guesses_vector = horzcat(guesses_vector,results{x}.iterations(a).perfmet.guesses);
                        acts_diff_vector = horzcat(acts_diff_vector, results{x}.iterations(a).acts(1,:)-results{x}.iterations(a).acts(2,:));
                    end

                    compiled_acts_diff(p,:) = acts_diff_vector;
                    compiled_guesses_vector(p,:) = guesses_vector;
                    compiled_correct_vector(p,:) = correct_vector;
                    compiled_desireds_vector(p,:) = desireds_vector;

                    overall_accuracy = mean(correct_vector)
                    overall_hit_rate = mean(correct_vector(desireds_vector==1));
                    overall_fa_rate = 1-mean(correct_vector(desireds_vector==2));
                    overall_d_prime = norminv(overall_hit_rate)-norminv(overall_fa_rate);

                    correct_vector = [];
                    desireds_vector = [];
                    guesses_vector = [];
                    acts_diff_vector = [];

                    %                     % classification_accuracy_by_resp (for full classification only)
                    %                     for b = 1:10
                    %                         classification_accuracy_by_resp(b) = mean(correct_vector(find(condensed_regs_all(b,active_trials))));
                    %                         mean_acts_diffs_by_resp(b) = mean(acts_diff_vector(find(condensed_regs_all(b,active_trials))));
                    %                         number_of_trials_per_bin(b) = length(find(condensed_regs_all(b,active_trials)));
                    %                     end
                end

                
                
        % APPLY TRAINED CLASSIFIER ON NEW DATA (SWITCH PATTERNS)
        [switch_results switch_acts_all switch_acts_test_idx] = apply_trained_classifier(subj,results{x},'conds_switch','runs_switch_xval');

        switch_correct_vector = [];
        switch_desireds_vector = [];
        switch_guesses_vector = [];
        switch_acts_diff_vector = [];
        for a = 1:num_runs
            switch_correct_vector = horzcat(switch_correct_vector,switch_results.iterations(a).perfmet.corrects);
            switch_desireds_vector = horzcat(switch_desireds_vector,switch_results.iterations(a).perfmet.desireds);
            switch_guesses_vector = horzcat(switch_guesses_vector,switch_results.iterations(a).perfmet.guesses);
            switch_acts_diff_vector = horzcat(switch_acts_diff_vector,switch_results.iterations(a).acts(1,:)-switch_results.iterations(a).acts(2,:));
        end
        
            switch_overall_accuracy = mean(correct_vector)
            switch_overall_hit_rate = mean(correct_vector(desireds_vector==1))
            switch_overall_fa_rate = 1-mean(correct_vector(desireds_vector==2))
            switch_overall_d_prime = norminv(overall_hit_rate)-norminv(overall_fa_rate)          

            % sort by absolute value of classifier "confidence"
            [switch_abs_sorted_diffs switch_abs_ind] = sort(abs(switch_acts_diff_vector),2,'descend');
            switch_abs_correct_sorted = correct_vector(switch_abs_ind);
            switch_abs_desireds_sorted = desireds_vector(switch_abs_ind);
                
                
                
                
                
                
                
                
                
                
                
             
            end

            if flags.plot_ROC_curve == 1


                if p>1
                    mean_acts_diff_vector = mean(compiled_acts_diff,1);
                else
                    mean_acts_diff_vector = compiled_acts_diff;
                end

                desireds_vector = compiled_desireds_vector(1,:); % take first b/c all desireds vectors are the same
                guesses_vector = sign(mean_acts_diff_vector);
                guesses_vector(guesses_vector==-1)=2;
                ROC{r}.compiled_accuracy(n) = mean(guesses_vector==desireds_vector)

                % sort by signed classifier "confidence" (for ROI curves)
                [sorted_diffs ind] = sort(mean_acts_diff_vector,2,'descend');
                desireds_sorted = desireds_vector(ind);

                %atanh_sorted_diffs = atanh(sorted_diffs);
                %                 atanh_cutoff_bins = [9 8 7 6 5 6 4 2 1 .75 .5 .25 0 -.25 -.5 .75 1 2 3 4 5 6 7 8 9];
                %                 for bin = 1:length(atanh_cutoff_bins)
                %                 hit_count(bin) = length(find(desireds_sorted(atanh_sorted_diffs>atanh_cutoff_bins(bin))==1));
                %                 fa_count(bin) = length(find(desireds_sorted(atanh_sorted_diffs>atanh_cutoff_bins(bin))==2));
                %                 end


                num_trials = length(sorted_diffs);

                %                 bin_size = 30;
                %                 for j= 1:floor(num_trials/bin_size)
                %
                %                     hit_count(j) = length(find(desireds_sorted(1:bin_size*j)==1));
                %                     fa_count(j) = length(find(desireds_sorted(1:bin_size*j)==2));
                %                 end


                num_bins = 8;
                for j= 1:num_bins

                    hit_count(j) = length(find(desireds_sorted(1:ceil(j*num_trials/num_bins))==1));
                    fa_count(j) = length(find(desireds_sorted(1:ceil(j*num_trials/num_bins))==2));
                end


                hit_rate = hit_count / length(find(desireds_sorted == 1));
                fa_rate = fa_count / length(find(desireds_sorted == 2));

                [Ro,Rn,dp,crit,SSE,pts,flag]=dpsd(hit_rate,fa_rate,'trim');

                ROC{r}.dpsd_recollection_8_bin(n) = Ro;
                ROC{r}.dpsd_familiarity_8_bin(n) = dp;
                ROC{r}.hit_rate_8_bin(n,:)= hit_rate;
                ROC{r}.fa_rate_8_bin(n,:)=fa_rate;

                clear hit_count fa_count;

                num_bins = 5;
                for j= 1:num_bins

                    hit_count(j) = length(find(desireds_sorted(1:ceil(j*num_trials/num_bins))==1));
                    fa_count(j) = length(find(desireds_sorted(1:ceil(j*num_trials/num_bins))==2));
                end


                hit_rate = hit_count / length(find(desireds_sorted == 1));
                fa_rate = fa_count / length(find(desireds_sorted == 2));

                [Ro,Rn,dp,crit,SSE,pts,flag]=dpsd(hit_rate,fa_rate,'trim');

                ROC{r}.dpsd_recollection_5_bin(n) = Ro;
                ROC{r}.dpsd_familiarity_5_bin(n) = dp;
                ROC{r}.hit_rate_5_bin(n,:)= hit_rate;
                ROC{r}.fa_rate_5_bin(n,:)=fa_rate;

                clear hit_count fa_count;

                %
                %                 % create continuous ROC function
                %                 for i = 1:length(sorted_diffs);
                %                     hit_rate(i) = length(correct_sorted(intersect(find(desireds_sorted == 1),[1:i]))) / length(find(desireds_sorted == 1));
                %                     fa_rate(i) = length(correct_sorted(intersect(find(desireds_sorted == 2),[1:i]))) / length(find(desireds_sorted == 2));
                %                 end
                %
                %                 figure
                %                 plot(fa_rate,hit_rate,'.-')
                %                 hold on
                %                 plot([0 1],[0 1],'r')
                %                 xlabel('P(Old|New)')
                %                 ylabel('P(Old|Old)')
            end



            if flags.write_data_log_to_text_file==1

                filename= [xls_results_data_logs_dir '/' subj_id '_' condnames{1} '_vs_' condnames{2} '.txt'];
                fid=fopen(filename, 'wt');
                fprintf(fid, '%s\r\n', ['subj_id = ' subj_id]);
                fprintf(fid, '%s\r\n', ['ROI_name = ' roi_name]);
                fprintf(fid, '%s\r\n', ['data_imgs_to_use =' data_imgs_to_use]);
                fprintf(fid, '%s\r\n', ['TR_weights = ' num2str(TR_weights)]);
                fprintf(fid, '%s\r\n', ['classification:' condnames{1} ' vs. ' condnames{2}]);
                fprintf(fid, '%s\r\n', ['flags.equate_number_of_trials_in_cond_1_and_2 = ' num2str(flags.equate_number_of_trials_in_cond_1_and_2)]);
                fprintf(fid, '%s\r\n', ['flags.equate_number_of_old_new_trials_per_subjective_bin = ' num2str(flags.equate_number_of_old_new_trials_per_subjective_bin)]);
                fprintf(fid, '%s\r\n', ['flags.perform_second_round_of_zscoring = ' num2str(flags.perform_second_round_of_zscoring)]);
                fprintf(fid, '%s\r\n', ['flags.remove_mvpa_outlier_trials (std dev) = ' num2str(flags.remove_outlier_trials)]);
                fprintf(fid, '%s\r\n', ['flags.remove_artdetect_outlier_trials (std dev) = ' num2str(flags.remove_artdetect_outliers)]);
                fprintf(fid, '%s\r\n', ['flags.artdetect_motion_thresh = ' num2str(flags.artdetect_motion_thresh)]);
                fprintf(fid, '%s\r\n', ['flags.artdetect_global_signal_thresh = ' num2str(flags.artdetect_global_signal_thresh)]);
                fprintf(fid, '\n\n');

                fprintf(fid, 'results_iter\toverall_acc\tHits\tFAs\td_prime');
                fprintf(fid, '#trials\t\t\t\t\t');
                fprintf(fid, '%3.0f\t', number_of_trials_per_bin(1:10));
                fprintf(fid, '\n');

                for q=1:x
                    fprintf(fid, '%4.4f\t', q);
                    fprintf(fid, '%4.4f\t', data_log.overall_acc(q));
                    fprintf(fid, '%4.4f\t', data_log.hits(q));
                    fprintf(fid, '%4.4f\t', data_log.FAs(q));
                    fprintf(fid, '%4.4f\t', data_log.d_prime(q));
                    %fprintf(fid, '%4.4f\t', data_log.classification_accuracy_by_resp(q,:)); %[data_log.classification_accuracy_by_resp(1,q) '\t' data_log.classification_accuracy_by_resp(q,1) '\t' data_log.classification_accuracy_by_resp(q,1) '\t' data_log.classification_accuracy_by_resp(q,4) '\t' data_log.classification_accuracy_by_resp(q,5) '\t' data_log.classification_accuracy_by_resp(q,6) '\t' data_log.classification_accuracy_by_resp(q,7) '\t' data_log.classification_accuracy_by_resp(q,8) '\t' data_log.classification_accuracy_by_resp(q,9) '\t' data_log.classification_accuracy_by_resp(q,10) '\t']);
                    %fprintf(fid, '%4.4f\t', data_log.acc_sorted_by_classifier_confidence(q,:));
                    fprintf(fid, '\n');
                    %for p=1:size(data_log.acc_sorted_by_classifier_confidence, 1)
                    %    fprintf(fid, '%4.4f\t', data_log.acc_sorted_by_classifier_confidence(q,p));
                    %end

                end
                fprintf(fid, '%s\t', 'mean');
                fprintf(fid, '%4.4f\t', mean(data_log.overall_acc));
                fprintf(fid, '%4.4f\t', mean(data_log.hits));
                fprintf(fid, '%4.4f\t', mean(data_log.FAs));
                fprintf(fid, '%4.4f\t', mean(data_log.d_prime));
                fprintf(fid, '%4.4f\t', mean(data_log.classification_accuracy_by_resp,1));
                fprintf(fid, '%4.4f\t', mean(data_log.acc_sorted_by_classifier_confidence,1));
                fprintf(fid, '\n');
                fclose(fid);
            end
        end

    end
    save_cmd = ['save ' subj_id '_ROI_ROC_data_logreg.mat ROC'];
    eval(save_cmd);
    clear ROC;


    time2finish = toc/60;
    display(['Finished ' subj_id ' in ' num2str(time2finish) ' minutes']);
end
