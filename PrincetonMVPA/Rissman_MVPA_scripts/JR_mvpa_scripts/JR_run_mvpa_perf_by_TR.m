function [subj results total_perf] = JR_run_mvpa_v5_batch(subj_arr, cond1, cond2, balance_per_subj_bin)


for b=subj_arr
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

    importance_maps_dir=[PAST_dir '/mvpa_results/importance_maps/' cond1 '_vs_' cond2 '_' balanced_or_unbal ];
    xls_results_logs_dir=[PAST_dir '/mvpa_results/xls_results_logs/' cond1 '_vs_' cond2 '_' balanced_or_unbal ];

    if ~exist(importance_maps_dir,'dir')
        mkdir(importance_maps_dir);
    end
    if ~exist(xls_results_logs_dir, 'dir')
        mkdir(xls_results_logs_dir)
    end




    mvpa_workspace = [PAST_dir '/' subj_id '/mvpa/' subj_id '_merged_AAL_ROIs_8mm_smoothing.mat'];
    %%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    flags.save_workspace = 1; % 1 = yes please, 0 = no thanks

    % unless a previously created workspace is specified, load in the preprocessed data
    if ~exist('mvpa_workspace')


        exp_name = 'PAST';
        roi_file = '/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/merged_AAL_ROIs.nii';
        roi_name = 'merged_AAL_ROIs'
        num_TP_per_run = 203;

        % load user-created filename and onsets lists into workspace
        data_imgs_to_use = 'raw_filenames_s8mm_wa.mat';
        load(data_imgs_to_use); %loads predefined cell array called raw_filenames into memory
        num_runs = length(raw_filenames)/num_TP_per_run; %calculate the number of runs (allows script to work flexibly for subjects with missing runs)
        load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']);
        load('/Users/Jesse/fMRI/data/PAST/fMRI/vol_info.mat'); %get functional data resolution info for spm .img writing
        [subj] = JR_mvpa_load_and_preprocess_raw_data(subj_id, exp_name, roi_name, roi_file, raw_filenames, num_runs, num_TP_per_run);

        if flags.save_workspace == 1
            save_cmd = ['save ' subj_id '_' roi_name '_8mm_smoothing_' datetime '.mat'];
            eval(save_cmd);
        end
    else
        eval(['load ' mvpa_workspace])
        if exist('log')
            clear log;  % just in case an old 'log' variable exists in memory
            display('log variable cleared from memory');
        end
        load('/Users/Jesse/fMRI/data/PAST/fMRI/vol_info.mat'); %get functional data resolution info for spm .img writing

    end

    % Set flags (% 1 = yes please, 0 = no thanks)
    flags.equate_number_of_old_new_trials_per_subjective_bin =balance_per_subj_bin; % equate_per_subjective_bin;
    flags.equate_number_of_trials_in_cond_1_and_2 = 1;
    flags.plot_mean_timecourses = 0;
    flags.plot_ROC_curve = 0;
    flags.display_performance_breakdown = 1;
    flags.generate_importance_maps = 1;
    flags.write_log_to_text_file=1;


    % specify which conditions to use for classification (must correspond to the names of conditions specified below)
    condnames =  { cond1 , cond2};

    data_imgs_to_use = 'raw_filenames_s8mm_wa.mat';

    TRs_to_average_over = [4]; %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier

    num_results_iter = 4; % number of times to run the cross validation process

    anova_p_thresh = 1;  %p threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)


    % classifier parameters
    class_args.train_funct_name = 'train_bp';
    class_args.test_funct_name = 'test_bp';
    class_args.nHidden = 0;

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
    if ~isempty(find(strcmp(subj_id,{'s104','s105'})))
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

    % condense regs by removing zeros
    all_trials = sum(all_regs,1); % vector of all trial   (patterns{4} contains fully preprocessed data)

    TR0 = subj.patterns{end}.mat(:,find(all_trials)-1); % -1 TR (-2-0 sec)
    TR1 = subj.patterns{end}.mat(:,find(all_trials)+0); % 1st TR (0-2 sec)
    TR2 = subj.patterns{end}.mat(:,find(all_trials)+1); % 2nd TR (2-4 sec)
    TR3 = subj.patterns{end}.mat(:,find(all_trials)+2); % 3rd TR (4-6 sec)
    TR4 = subj.patterns{end}.mat(:,find(all_trials)+3); % 4th TR (6-8 sec)
    TR5 = subj.patterns{end}.mat(:,find(all_trials)+4); % 5th TR (8-10 sec)
    TR6 = subj.patterns{end}.mat(:,find(all_trials)+5); % 6th TR (10-12 sec)


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
    subj = duplicate_object(subj,'pattern','spiral_d_hp_z','spiral_d_hp_z_condensed_TR0');
    subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed_TR0',TR0);

    % clean up workspace
    subj = remove_mat(subj,'pattern','spiral_d_hp_z');
    clear mean_data; clear TR0;

    subj = duplicate_object(subj,'pattern','spiral_d_hp_z_condensed_TR0','spiral_d_hp_z_condensed_TR1');
    subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed_TR1',TR1);

    clear TR1;

    subj = duplicate_object(subj,'pattern','spiral_d_hp_z_condensed_TR0','spiral_d_hp_z_condensed_TR2');
    subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed_TR2',TR2);

    clear TR2;

    subj = duplicate_object(subj,'pattern','spiral_d_hp_z_condensed_TR0','spiral_d_hp_z_condensed_TR3');
    subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed_TR3',TR3);

    clear TR3;

    subj = duplicate_object(subj,'pattern','spiral_d_hp_z_condensed_TR0','spiral_d_hp_z_condensed_TR4');
    subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed_TR4',TR4);

    clear TR4;

    subj = duplicate_object(subj,'pattern','spiral_d_hp_z_condensed_TR0','spiral_d_hp_z_condensed_TR5');
    subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed_TR5',TR5);

    clear TR5;

    subj = duplicate_object(subj,'pattern','spiral_d_hp_z_condensed_TR0','spiral_d_hp_z_condensed_TR6');
    subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed_TR6',TR6);

    clear TR6;


    % update run vector to condensed format
    subj.selectors{1}.mat = condensed_runs;
    subj.selectors{1}.matsize = size(condensed_runs);

    % "activate" only those trials of interest (from regs_of_interest) before
    % creating cross-validation indices

    temp_sel = ones(1,size(condensed_regs_all,2)); % intialize vector of all ones
    temp_sel(find(sum(condensed_regs_of_interest)==0)) = 0; % remove all non-"regs_of_interst" trials (set to zero)
    subj = init_object(subj,'selector','conditions_of_interest'); %initialize selector object
    subj = set_mat(subj,'selector','conditions_of_interest',temp_sel);
    subj = create_xvalid_indices(subj,'runs','actives_selname','conditions_of_interest');

    % run feature selection ANOVA
    %subj = feature_select(subj,'spiral_d_hp_z_condensed','conds','runs_xval','thresh',anova_p_thresh);

    % run classifier (hidden layer netlab backprop algorithm)
    % class_args.train_funct_name = 'train_bp_netlab';
    % class_args.test_funct_name = 'test_bp_netlab';
    % class_args.nHidden = 10;
    % [subj results] = cross_validation(subj,'spiral_d_hp_z_condensed','conds','runs_xval',subj.masks{2}.group_name,class_args);

    % run classifier (No hidden layer NN Toolbox backprop algorithm)
    class_args.train_funct_name = 'train_bp';
    class_args.test_funct_name = 'test_bp';
    class_args.nHidden = 0;

    for r = 1:4
        for v= 1:7
            pattern_name = ['spiral_d_hp_z_condensed_TR' num2str(v-1)];

            [subj results] = cross_validation(subj,pattern_name,'conds','runs_xval',subj.masks{1}.name,class_args);

            correct_vector = [];
            desireds_vector = [];
            guesses_vector = [];
            acts_diff_vector = [];
            for a = 1:num_runs
                correct_vector = horzcat(correct_vector,results.iterations(a).perfmet.corrects);

            end

            overall_accuracy = mean(correct_vector)



            total_perf(r,v) = overall_accuracy;
           
        end

    end
end

