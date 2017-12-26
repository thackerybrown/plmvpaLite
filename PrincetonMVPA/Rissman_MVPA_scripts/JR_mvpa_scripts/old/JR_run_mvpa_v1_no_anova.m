function [subj results] = JR_run_mvpa_v1_no_anova(mvpa_workspace)

tic

if ~exist('mvpa_workspace')

    %%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subj_id = 's101';
    exp_name = 'PAST';
    roi_file = 'whole_brain_no_motor.img';
    %anova_p_thresh = .005;  %p threshold for feature selection ANOVA
    generate_importance_maps = 0; % 1 = yes please, 2 = no thanks
    vol_info = spm_vol(['/Users/Jesse/fMRI/data/PAST/fMRI/' subj_id '/anatomical/' subj_id '_mean_func.img']); %get functional data resolution info for spm .img writing

    num_runs = 8;
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

    
    % SPECIAL FIX FOR s104 because reg #6 (NEW_recollect) is a fake
    % placeholder trial
    if strcmp(subj_id,'s104')
        all_regs(6,:)=0;
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
    
    lc_hits = all_regs(3,:); 
    fas = sum(all_regs(6:8,:)); 
    
    %no_resp = all_regs(11,:); %excluded from analysis
    
     
    %choose conditions to train/test classifier on
    regs_of_interest(1,:) = correct_recollect;
    regs_of_interest(2,:) = correct_hc_old;
   
    condnames =  {'correct_recollect','correct_hc_old'}; %{'old','new'};

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
    
    save_cmd = ['save ' subj_id '_spiral_d_hp_z_' roi_file(1:end-4) '_' datetime '.mat'];
    eval(save_cmd); 

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

mean_data = (TR3+TR4+TR5)/3; %take uniform average of peak signal TRs

clear TR0 TR1 TR2 TR3 TR4 TR5; %clean up matlab workspace to save memory


% BASELINE CORRECT MEAN DATA (optional)
% mean_baseline = (t0+t1)/2;
% mean_data = mean_data - mean_baseline;


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
    if ~isempty(find(all_regs(:,i))) % if not a rest trial
        condensed_regs_of_interest(:,trial_counter) = regs_of_interest(:,i);
        condensed_regs_all(:,trial_counter) = all_regs(:,i); %save
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
% [subj results] = cross_validation(subj,'spiral_d_hp_z','conds','runs_xval',subj.masks{2}.group_name,class_args);

% run classifier (No hidden layer NN Toolbox backprop algorithm)
class_args.train_funct_name = 'train_bp';
class_args.test_funct_name = 'test_bp';
class_args.nHidden = 0;

[subj results] = cross_validation(subj,'spiral_d_hp_z_condensed','conds','runs_xval',subj.masks{1}.name,class_args);


correct_vector = [];
desireds_vector = [];
guesses_vector = [];
acts_diff_vector = [];
for a = 1:num_runs
    correct_vector = horzcat(correct_vector,results.iterations(a).perfmet.corrects);
    desireds_vector = horzcat(desireds_vector,results.iterations(a).perfmet.desireds);
    guesses_vector = horzcat(guesses_vector,results.iterations(a).perfmet.guesses);
    acts_diff_vector = horzcat(acts_diff_vector, results.iterations(a).acts(1,:)-results.iterations(a).acts(2,:));
end

% classification_accuracy_by_resp = [];
for b = 1:10
    classification_accuracy_by_resp(b) = mean(correct_vector(find(condensed_regs_all(b,:))));
    mean_acts_diffs_by_resp(b) = mean(acts_diff_vector(find(condensed_regs_all(b,:))));
end

classification_accuracy_by_resp
mean_acts_diffs_by_resp

[sorted_diffs ind] = sort(acts_diff_vector,2,'descend');
correct_sorted = correct_vector(ind);
desireds_sorted = desireds_vector(ind);

% create continuous ROC function 
 for i = 1:length(sorted_diffs);
 hit_rate(i) = length(correct_sorted(intersect(find(desireds_sorted == 1),[1:i]))) / length(find(desireds_sorted == 1));
 fa_rate(i) = length(correct_sorted(intersect(find(desireds_sorted == 2),[1:i]))) / length(find(desireds_sorted == 2));
 end

% created binned ROC
counter = 1;
for i = 1:8
hit(i) = length(correct_sorted(intersect(find(desireds_sorted == 1),[1:counter+45]))) / length(find(desireds_sorted == 1));
fa(i) = length(correct_sorted(intersect(find(desireds_sorted == 2),[1:counter+45]))) / length(find(desireds_sorted == 2));
counter = counter+46;
end
figure
plot(fa_rate,hit_rate,'.-')
hold on
plot([0 1],[0 1],'r')
xlabel('P(Old|New)')
ylabel('P(Old|Old)')

if generate_importance_maps == 1;  % NO ANOVA VERSION
           
    subj = interpret_weights(subj, results);
        
    immap1 = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 1
    immap2 = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 2
    
          
        voxel_inds = find(subj.masks{end}.mat); %get mask voxel indices
        immap1(voxel_inds)=subj.patterns{end}.mat(:,1); %store immap values at appropriate voxel indices
        immap2(voxel_inds)=subj.patterns{end}.mat(:,2); %store immap values at appropriate voxel indices
        
        immap1 = immap1*1000;
        immap2 = immap2*1000;       
    
    vol_info.fname = [condnames{1} '_no_anova_' datetime '.img'];
    spm_write_vol(vol_info,immap1);
    vol_info.fname = [condnames{2} '_no_anova_' datetime '.img'];
    spm_write_vol(vol_info,immap2);
end
    
time2finish = toc/60;
disp(['time2finish = ' num2str(time2finish) ' minutes']);

display('all done');