function [subj results] = JR_run_mvpa_v4_bait_and_switch

tic

subjects_to_include = [6:14];

for a=1:length(subjects_to_include)
    b = subjects_to_include(a);
    if b<10
        subj_id=strcat('s10', num2str(b))
    else
        subj_id=strcat('s1', num2str(b))
    end

    mvpa_workspace=strcat('/Users/Jesse/fMRI/data/PAST/fMRI/' , subj_id, '/mvpa/', subj_id, '_merged_AAL_ROIs_8mm_smoothing.mat');
    %%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    flags.save_workspace = 1; % 1 = yes please, 0 = no thanks

    % unless a previously created workspace is specified, load in the preprocessed data
    if ~exist('mvpa_workspace')

        %subj_id = 's102';
        exp_name = 'PAST';
        roi_file = '/Users/Jesse/fMRI/data/PAST/fMRI/4x4x4/merged_AAL_ROIs.nii';
        roi_name = 'merged_AAL_ROIs'
        num_TP_per_run = 203;

        % load user-created filename and onsets lists into workspace
        data_imgs_to_use = 'raw_filenames_s8mm_wa.mat';
        load(data_imgs_to_use); %loads predefined cell array called raw_filenames into memory
        num_runs = length(raw_filenames)/num_TP_per_run; %calculate the number of runs (allows script to work flexibly for subjects with missing runs)
        load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']);
        vol_info = spm_vol(raw_filenames{1}); %get functional data resolution info for spm .img writing
        [subj] = JR_mvpa_load_and_preprocess_raw_data(subj_id, exp_name, roi_name, roi_file, raw_filenames, num_runs, num_TP_per_run);

        if flags.save_workspace == 1
            save_cmd = ['save ' subj_id '_' roi_name '_4mm_smoothing_' datetime '.mat'];
            eval(save_cmd);
        end
    else
        eval(['load ' mvpa_workspace])
    end

    % Set flags (% 1 = yes please, 0 = no thanks)
    flags.equate_number_of_old_new_trials_per_subjective_bin =1 % equate_per_subjective_bin;
    flags.equate_number_of_trials_in_cond_1_and_2 = 0;
    flags.plot_mean_timecourses = 0;
    flags.plot_ROC_curve = 0;
    flags.display_performance_breakdown = 1;
    flags.generate_importance_maps = 0;
    flags.save_data_in_log=1;
    flags.write_log_to_text_file=1;


    % specify which conditions to use for classification (must correspond to the names of conditions specified below)
    condnames_bait =  {'Hits', 'FAs'};
    condnames_switch = {'Misses','CRs'};

    data_imgs_to_use = 'raw_filenames_s8mm_wa.mat';

    TRs_to_average_over = [4]; %which post-stimulus TRs should be used (and if more than one, averaged across) before feeding data to the classifier

    num_results_iter = 2; % number of times to run the cross validation process

    anova_p_thresh = 1;  %p threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)


    % classifier parameters
    class_args.train_funct_name = 'train_bp_netlab';
    class_args.test_funct_name = 'test_bp_netlab';
    class_args.nHidden = 10;

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
    regs_bait = [];
    eval(['regs_bait(1,:) = ' condnames_bait{1} ';'])
    eval(['regs_bait(2,:) = ' condnames_bait{2} ';'])

    regs_switch = [];
    eval(['regs_switch(1,:) = ' condnames_switch{1} ';'])
    eval(['regs_switch(2,:) = ' condnames_switch{2} ';'])

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
            display([num2str(length(trials_to_cut)) ' trials cut from ' condnames_bait{1}]);
        elseif num_cond1 < num_cond2
            rand_array = rand(1,num_cond2);
            [sorted inds]= sort(rand_array);
            trials_to_cut = cond2_trials(inds(1:num_cond2-num_cond1));
            regs_of_interest(2,trials_to_cut) = 0;
            display([num2str(length(trials_to_cut)) ' trials cut from ' condnames_bait{2}]);
        else
            display('Trial numbers are already balanced');
        end
    end

    display([num2str(count(regs_bait(1,:)==1)) ' baiting trials in condition ' condnames_bait{1}])
    display([num2str(count(regs_bait(2,:)==1)) ' baiting trials in condition ' condnames_bait{2}])
    display([num2str(count(regs_switch(1,:)==1)) ' switch trials in condition ' condnames_switch{1}])
    display([num2str(count(regs_switch(2,:)==1)) ' switch trials in condition ' condnames_switch{2}])
    
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
    if ~isempty(find(all_regs(:,i))) % if not a rest trial
        condensed_regs_bait(:,trial_counter) = regs_bait(:,i);
        condensed_regs_switch(:,trial_counter) = regs_switch(:,i);
        condensed_regs_all(:,trial_counter) = all_regs(:,i); %save
        condensed_runs(1,trial_counter) = subj.selectors{1}.mat(i);
        trial_counter = trial_counter + 1;
    end
end

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
    
    

    % zscore temporally-condensed data; active trials only (second round of z-scoring)
    % subj.patterns{end}.mat(:,active_trials) = zscore(subj.patterns{end}.mat(:,active_trials)')';

    % run feature selection ANOVA (if desired)
    if anova_p_thresh ~= 1
        subj = feature_select(subj,'spiral_d_hp_z_condensed','conds_bait','runs_xval','thresh',anova_p_thresh);
        classifier_mask = subj.masks{end}.group_name; % use group of masks created by ANOVA
    else
        classifier_mask = subj.masks{1}.name; % use original mask
    end

    %%%%%%%%%%%%%%%%%%%%%% RUN THE CLASSIFIER (CROSS-VALIDATION)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for x = 1: num_results_iter

        [subj results{x}] = cross_validation(subj,'spiral_d_hp_z_condensed','conds_bait','runs_xval',classifier_mask,class_args);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % do some important RAM clean-up and data archiving
        for y = 1:num_runs
            if flags.generate_importance_maps == 1
                results_IW{x}.iterations(y).scratchpad.net.IW{1} = single(results{x}.iterations(y).scratchpad.net.IW{1}); % save weights to pass to JR_interpret_weights
            end
            results{x}.iterations(y).scratchpad.net.inputs{1}.exampleInput=[]; % delete huge data object from results scratchpad to free up RAM
        end

        if flags.display_performance_breakdown == 1
            % analyze the results in more detail
            correct_vector = [];
            desireds_vector = [];
            guesses_vector = [];
            acts_diff_vector = [];
            for a = 1:num_runs
                correct_vector = horzcat(correct_vector,results{x}.iterations(a).perfmet.corrects);
                desireds_vector = horzcat(desireds_vector,results{x}.iterations(a).perfmet.desireds);
                guesses_vector = horzcat(guesses_vector,results{x}.iterations(a).perfmet.guesses);
                acts_diff_vector = horzcat(acts_diff_vector, results{x}.iterations(a).acts(1,:)-results{x}.iterations(a).acts(2,:));
            end

            overall_accuracy = mean(correct_vector)
            overall_hit_rate = mean(correct_vector(desireds_vector==1))
            overall_fa_rate = 1-mean(correct_vector(desireds_vector==2))
            overall_d_prime = norminv(overall_hit_rate)-norminv(overall_fa_rate)


            % classification_accuracy_by_resp (for full classification only)
            for b = 1:10
                classification_accuracy_by_resp(b) = mean(correct_vector(find(condensed_regs_all(b,active_trials))));
                mean_acts_diffs_by_resp(b) = mean(acts_diff_vector(find(condensed_regs_all(b,active_trials))));
                number_of_trials_per_bin(b) = length(find(condensed_regs_all(b,active_trials)));
            end

          
            % sort by absolute value of classifier "confidence"
            [abs_sorted_diffs abs_ind] = sort(abs(acts_diff_vector),2,'descend');
            abs_correct_sorted = correct_vector(abs_ind);
            abs_desireds_sorted = desireds_vector(abs_ind);

            % print results of top 50,100, 150, etc. trials, sorted by
            % classifier "confidence"
            num_trials = length(abs_correct_sorted);
            if num_trials>50

                for j= 1:floor(num_trials/50)
                    acc_sorted_by_classifier_confidence(j)=mean(abs_correct_sorted(1:50*j));
                end
                acc_sorted_by_classifier_confidence(j+1)=mean(abs_correct_sorted(1:end));

              
            end
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

            % print results of top 50,100, 150, etc. trials, sorted by
            % classifier "confidence"
            
            acc_sorted_by_classifier_confidence = zeros(1,8);
            
            num_trials = length(switch_abs_correct_sorted);
            if num_trials>50

                for j= 1:floor(num_trials/50)
                    switch_acc_sorted_by_classifier_confidence(j)=mean(switch_abs_correct_sorted(1:50*j));
                end
                switch_acc_sorted_by_classifier_confidence(j+1)=mean(switch_abs_correct_sorted(1:end));

               
            end
        
        
        

        if flags.plot_ROC_curve == 1

            % sort by signed classifier "confidence" (for ROI curves)
            [sorted_diffs ind] = sort(acts_diff_vector,2,'descend');
            correct_sorted = correct_vector(ind);
            desireds_sorted = desireds_vector(ind);

            % create continuous ROC function
            for i = 1:length(sorted_diffs);
                hit_rate(i) = length(correct_sorted(intersect(find(desireds_sorted == 1),[1:i]))) / length(find(desireds_sorted == 1));
                fa_rate(i) = length(correct_sorted(intersect(find(desireds_sorted == 2),[1:i]))) / length(find(desireds_sorted == 2));
            end

            figure
            plot(fa_rate,hit_rate,'.-')
            hold on
            plot([0 1],[0 1],'r')
            xlabel('P(Old|New)')
            ylabel('P(Old|Old)')
        end

        if flags.save_data_in_log==1
            log.overall_acc(x)=switch_overall_accuracy;
            log.hits(x)=switch_overall_hit_rate;
            log.FAs(x)=switch_overall_fa_rate;
            log.d_prime(x)=switch_overall_d_prime;
            log.classification_accuracy_by_resp(x,:)=classification_accuracy_by_resp;  %print values from original classification
            log.number_trials_per_bin(x,:)=number_of_trials_per_bin;
            log.acc_sorted_by_classifier_confidence(x,:)=switch_acc_sorted_by_classifier_confidence;
        end
    end


    if flags.write_log_to_text_file==1
        filename=strcat(subj_id, condnames_bait{1}, '_vs_', condnames_bait{2}, '_as_bait.xls');
        fid=fopen(filename, 'wt');
        fprintf(filename, 'wt');
        fprintf(fid, '%s\r\n', ['subj_id = ' subj_id]);
        fprintf(fid, '%s\r\n', ['ROI_name = ' roi_name]);
        fprintf(fid, '%s\r\n\', ['data_imgs_to_use =' data_imgs_to_use]);
        fprintf(fid, '%s\r\n', ['TRs_to_average_over = ' num2str(TRs_to_average_over)]);
        fprintf(fid, '%s\r\n\', ['bait classification:' condnames_bait{1} ' vs. ' condnames_bait{2}]);
        fprintf(fid, '%s\r\n\', ['switch classification:' condnames_switch{1} ' vs. ' condnames_switch{2}]);        
        fprintf(fid, '%s\r\n\', ['flags.equate_number_of_trials_in_cond_1_and_2 = ' num2str(flags.equate_number_of_trials_in_cond_1_and_2)]);
        fprintf(fid, '%s\r\n\', ['flags.equate_number_of_old_new_trials_per_subjective_bin = ' num2str(flags.equate_number_of_old_new_trials_per_subjective_bin)]);


        fprintf(fid, '\n\n\n\n');
        fprintf(fid, 'results_iter\toverall_acc\tHits\tFAs\td_prime\tOLD_recollect\tOLD_hc_old\tOLD_lc_old\tOLD_lc_new\tOLD_hc_new\tNEW_recollect\tNEW_hc_old\tNEW_lc_old\tNEW_lc_new\tNEW_hc_new\tno_resp\tTop50\tTop100\tTop150\tTop200\tTop250\n');
        fprintf(fid, '0\t0\t0\t0\t0\t');
        fprintf(fid, '%3.0f\t', [number_of_trials_per_bin(1),  number_of_trials_per_bin(2),  number_of_trials_per_bin(3),  number_of_trials_per_bin(4), number_of_trials_per_bin(5),  number_of_trials_per_bin(6), number_of_trials_per_bin(7),  number_of_trials_per_bin(8), number_of_trials_per_bin(9), number_of_trials_per_bin(10) '/n']);
        fprintf(fid, '\n');

        for q=1:x
            fprintf(fid, '%3.3f\t', q);
            fprintf(fid, '%3.3f\t', log.overall_acc(q));
            fprintf(fid, '%3.3f\t', log.hits(q));
            fprintf(fid, '%3.3f\t', log.FAs(q));
            fprintf(fid, '%3.3f\t', log.d_prime(q));
            fprintf(fid, '%3.3f\t', log.classification_accuracy_by_resp(q,:)); %[log.classification_accuracy_by_resp(1,q) '\t' log.classification_accuracy_by_resp(q,1) '\t' log.classification_accuracy_by_resp(q,1) '\t' log.classification_accuracy_by_resp(q,4) '\t' log.classification_accuracy_by_resp(q,5) '\t' log.classification_accuracy_by_resp(q,6) '\t' log.classification_accuracy_by_resp(q,7) '\t' log.classification_accuracy_by_resp(q,8) '\t' log.classification_accuracy_by_resp(q,9) '\t' log.classification_accuracy_by_resp(q,10) '\t']);
            fprintf(fid, '%3.3f\t', log.acc_sorted_by_classifier_confidence(q,:));
            fprintf(fid, '\n');
            %for p=1:size(log.acc_sorted_by_classifier_confidence, 1)
            %    fprintf(fid, '%3.3f\t', log.acc_sorted_by_classifier_confidence(q,p));
            %end

        end
        fprintf(fid, '%s\t', 'mean');
        fprintf(fid, '%3.3f\t', mean(log.overall_acc));
        fprintf(fid, '%3.3f\t', mean(log.hits));
        fprintf(fid, '%3.3f\t', mean(log.FAs));
        fprintf(fid, '%3.3f\t', mean(log.d_prime));
        fprintf(fid, '%3.3f\t', mean(log.classification_accuracy_by_resp,1));
        fprintf(fid, '%3.3f\t', mean(log.acc_sorted_by_classifier_confidence,1));
        fprintf(fid, '\n');
    end



    if flags.generate_importance_maps == 1;

        subj = JR_interpret_weights(subj, results,results_IW);

        impmap1 = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 1
        impmap2 = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 2

        if anova_p_thresh == 1 % NO ANOVA VERSION

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

            vol_info.fname = [subj_id '_' condnames_bait{1} '_no_anova_' datetime '.img'];
            spm_write_vol(vol_info,impmap1);
            vol_info.fname = [subj_id '_' condnames_bait{2} '_no_anova_' datetime '.img'];
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

            %sum across masks to get composite mask (where value of each voxel =
            %number of runs for which that voxel was included)
            composite_mask = zeros(vol_info.dim);
            for i = 2:size(subj.masks,2)  %exclude first mask (it's the starting ROI)
                composite_mask = composite_mask+subj.masks{i}.mat;
            end
            voxels_to_exclude = find(composite_mask<6);  % exclude voxels that exist for less than 6 of the ANOVA masks
            impmap1(voxels_to_exclude)=0;
            impmap2(voxels_to_exclude)=0;

            impmap1_avg = impmap1./composite_mask * 1000;  % divide by number of observations contributing to each sum (to get avg) and multiply by 1000 for scaling
            impmap2_avg = impmap2./composite_mask * 1000;  % divide by number of observations contributing to each sum (to get avg) and multiply by 1000 for scaling

            vol_info.fname = [subj_id '_' condnames_bait{1} '_p' num2str(anova_p_thresh) '_' datetime '.img'];
            spm_write_vol(vol_info,impmap1_avg);
            vol_info.fname = [subj_id '_' condnames_bait{2} '_p' num2str(anova_p_thresh) '_' datetime '.img'];
            spm_write_vol(vol_info,impmap2_avg);
        end

    end

    clear log;

end
time2finish = toc/60;
display(['Finished in ' num2str(time2finish) ' minutes']);