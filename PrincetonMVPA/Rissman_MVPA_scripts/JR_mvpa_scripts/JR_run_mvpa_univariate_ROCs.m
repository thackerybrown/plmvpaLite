function []= JR_run_mvpa_univariate_ROCs(subj_array, condition1, condition2, nVox,penalty)

% example usage:  JR_run_mvpa_general([1 2 4 5 8 12 13 16], 'Hits', 'CRs')


%%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for b=subj_array  % allows subject ID #'s to be specified in short format
    tic
    if b<10
        subj_id=strcat('s10', num2str(b))
    else
        subj_id=strcat('s1', num2str(b))
    end
    
    %s
    expt_dir = '/Users/Jesse/fMRI/data/PAST/fMRI';
    mvpa_dir = [expt_dir '/' subj_id '/mvpa'];
    
    
    roi_name = 'SEPT09_MVPA_MASK';
    data_imgs_to_use = 'raw_filenames_s8mm_wa.mat'; % .mat file containing names of all functional images
    
    
    
    % load some .mat files into memory
    load([expt_dir '/vol_info.mat']); %get functional data resolution info for spm .img writing
    load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']); % load in your SPM-formatted onsets file
    
    % OPTIONAL:  specify previously saved mvpa workspace to bypass time-consuming data extraction and preprocessing
    %mvpa_workspace = [expt_dir '/' subj_id '/mvpa/' subj_id '_merged_AAL_ROIs_FIXED_HOLES_wa.mat'];
    
    mvpa_workspace = [expt_dir '/' subj_id '/mvpa/' subj_id '_SEPT09_MVPA_MASK_s8mm_wa.mat'];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set flags (% unless otherwise indicated: 1 = yes please, 0 = no thanks)
    flags.num_full_iter = 1;
    flags.num_results_iter = 1; % number of times to run the entire classification process (select subset of the data and train/test classifier)
    flags.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data
    
    flags.equate_number_of_old_new_trials_per_subjective_bin = 0; % equate number of trials per subjective bin
    flags.equate_number_of_trials_in_cond_1_and_2 = 0; % equate number of trials in conditions 1 and 2 (RECOMMENDED)
    flags.anova_p_thresh = 1;  % p-value threshold for feature selection ANOVA (1 = DON'T PERFORM ANY FEATURE SELECTION)
    flags.GLM_anova_p_thresh = 0;
    flags.anova_nVox_thresh = nVox; % alternative to specifying p-value threshold; uses top N voxels (0 = DON'T PERFORM ANY FEATURE SELECTION)
    flags.perform_second_round_of_zscoring = 1;  % z-score data again immediately prior to classification (NOT RECOMMENDED)
    flags.remove_artdetect_outliers = 1; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
    flags.artdetect_motion_thresh = 7; % specify ArtDetect bin for motion outliers
    flags.artdetect_global_signal_thresh = 5; % specify ArtDetect bin for global signal outliers
    flags.remove_outlier_trials = 0;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)
    flags.plot_ROC_curve = 1;
    flags.display_performance_breakdown = 1;
    flags.generate_importance_maps = 1;
    flags.generate_weight_maps = 1;
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
    class_args.penalty = penalty;
    
    
    dir_str = [condition1 '_vs_' condition2 '_5bin_ROCs'];
    
    
    if flags.generate_importance_maps == 1
        rec_fam_dpsd_maps_dir=[expt_dir '/mvpa_results/REC_FAM_dpsd_maps/' dir_str];
        if ~exist( rec_fam_dpsd_maps_dir,'dir')
            mkdir( rec_fam_dpsd_maps_dir);
        end
    end
    
    
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
        
        
        
        %% compute univariate voxel-by-voxel ROC curves %%
        
        
        active_data_pat = subj.patterns{end}.mat(:,active_trials);
        num_voxels = size(active_data_pat,1);
        num_trials = size(active_data_pat,2);
        auc_est = zeros(1,num_voxels);
        rec_est = zeros(1,num_voxels);
        fam_est = zeros(1,num_voxels);
        num_bins = 5;
        bin_size = ceil(num_trials/num_bins);
        bin_array = [bin_size bin_size*2 bin_size*3 bin_size*4 num_trials];
        
        OLD_items = find(subj.regressors{1}.mat(1,active_trials));
        NEW_items = find(subj.regressors{1}.mat(2,active_trials));
        
        old_new_labels = zeros (1,length(active_trials));
        old_new_labels(OLD_items) = 1;  % old items coded with 1, new items coded with 2
        old_new_labels(NEW_items) = 2;  % old items coded with 1, new items coded with 2
        
        
        for vv = 1:size(active_data_pat,1)
            current_vox = active_data_pat(vv,:);
            [ranked_values ind] = sort(current_vox,2,'descend');
            old_new_labels_sorted = old_new_labels(ind);
            
            %remove outlier data points%
            mean_value = mean(ranked_values);
            std_value = std(ranked_values);
            upper_outliers = find(ranked_values>mean_value+3*std_value);
            lower_outliers = find(ranked_values<mean_value-3*std_value);
            all_outliers = union(upper_outliers,lower_outliers);
            num_outliers = length(all_outliers);
            ranked_values(all_outliers)=[];
            old_new_labels_sorted(all_outliers)=[];
            num_OLD = count(old_new_labels_sorted == 1);
            num_NEW = count(old_new_labels_sorted == 2);
            
            
            
            % create continuous ROC curve
            %            for i = 1:length(ranked_values);
            %                hit_rate(i) = length(intersect(find(old_new_labels_sorted == 1),[1:i])) / length(find(old_new_labels_sorted == 1));
            %                fa_rate(i) = length(intersect(find(old_new_labels_sorted == 2),[1:i])) / length(find(old_new_labels_sorted == 2));
            %            end
            
            % create 5-bin ROC curve
            for i = 1:num_bins-1;
                hit_rate(i) = length(find(old_new_labels_sorted(1:i*bin_size)==1)) / num_OLD;
                fa_rate(i) =  length(find(old_new_labels_sorted(1:i*bin_size)==2)) / num_NEW;
            end
            hit_rate(i+1)=1;
            fa_rate(i+1)=1;
            
            auc_est(vv) = auroc(hit_rate',fa_rate');
            
            if auc_est(vv)>= 0.55
                [Ro,Rn,dp,crit,SSE,pts,flag]=dpsd(hit_rate,fa_rate,'trim');
                rec_est(vv) = Ro;
                fam_est(vv) = dp;
            end
        end
        
        
        
        auc_map = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 1
        fam_map = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 2
        rec_map = zeros(vol_info.dim);
        
        voxel_inds = find(subj.masks{end}.mat); %get mask voxel indices
        
        temp = zeros(vol_info.dim); %initialize appropriately sized matrix        
        temp(voxel_inds)=auc_est;
        
        vol_info.fname = [rec_fam_dpsd_maps_dir '/' subj_id '_AUC.img'];
        spm_write_vol(vol_info,temp);
        
        temp = zeros(vol_info.dim); %initialize appropriately sized matrix
        temp(voxel_inds)=rec_est;
        
        vol_info.fname = [rec_fam_dpsd_maps_dir '/' subj_id '_REC.img'];
        spm_write_vol(vol_info,temp);
        
        temp = zeros(vol_info.dim); %initialize appropriately sized matrix
        temp(voxel_inds)=fam_est;
        
        vol_info.fname = [rec_fam_dpsd_maps_dir '/' subj_id '_FAM.img'];
        spm_write_vol(vol_info,temp);       
      
        
    end
    
    
    time2finish = toc/60;
    display(['Finished ' subj_id ' in ' num2str(time2finish) ' minutes']);
    
    keep b subj_array condition1 condition2 nVox penalty
    
end


