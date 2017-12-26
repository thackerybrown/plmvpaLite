function [subj results] = JR_run_mvpa_v1(mvpa_workspace)

if ~exist('mvpa_workspace')

    %%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subj_id = 's104';
    exp_name = 'PAST';
    roi_file = 'whole_brain.img';
    anova_p_thresh = .0001;  %p threshold for feature selection ANOVA

    num_runs = 9;
    num_TP_per_run = 203;
    total_TP = num_runs * num_TP_per_run;

    % load user-created filename and onsets lists into workspace
    load('raw_filename_list_SRA.mat')
    load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']);

    num_conds = size(onsets,2);
    all_regs = zeros(num_conds,total_TP); % initialize regs matrix as conditions x timepoints

    for cond = 1: num_conds
        for trial = 1: length(onsets{cond})
            time_idx = onsets{cond}(trial)/2+1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
            all_regs(cond,time_idx) = 1;
        end
    end
    % condnames = names; %specify names 'OLD_recollect'    'OLD_hc_old'    'OLD_lc_old'    'OLD_lc_new'    'OLD_hc_new'    'NEW_recollect' 'NEW_hc_old'    'NEW_lc_old'    'NEW_lc_new'    'NEW_hc_new'    'no_resp'

    % define conditions of interest
    objective_old = sum(all_regs(1:5,:)); 
    objective_new = sum(all_regs(7:10,:)); % changed for s104 from objective_new = sum(all_regs(6:10,:));
    
    subjective_old_all = sum(all_regs([1 2 3 6 7 8],:));
    subjective_new_all = sum(all_regs([4 5 9 10],:));
    
    subjective_old_hc_only = sum(all_regs([1 2 6 7],:));
    subjective_old_hc_only = sum(all_regs([5 10],:));
    
    correct_recollect = all_regs(1,:);
    correct_hc_old = all_regs(2,:);    
    
    no_resp = all_regs(11,:); %excluded from analysis

    %choose conditions to train/test classifier on
    regs(1,:) = objective_old;
    regs(2,:) = objective_new;
   
    condnames = {'old','new'};

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
    
    save_cmd = ['save ' subj_id '_spiral_d_hp_z_' date '.mat'];
    eval(save_cmd); 

else
    load(mvpa_workspace); %load in pre-saved workspace files and continue
end

subj = apply_to_runs(subj, 'spiral_d_hp_z', 'runs', 'apply_filt');

% select TRs of interest (to correspond with peak post-stim BOLD response)

all_trials = sum(regs,1); % vector of all trial   (patterns{4} contains fully preprocessed data)

t0 = subj.patterns{5}.mat(:,find(all_trials)-1); % -1 TR (-2-0 sec)
t1 = subj.patterns{5}.mat(:,find(all_trials)+0); % 1st TR (0-2 sec)
t2 = subj.patterns{5}.mat(:,find(all_trials)+1); % 2nd TR (2-4 sec)
t3 = subj.patterns{5}.mat(:,find(all_trials)+2); % 3rd TR (4-6 sec)
t4 = subj.patterns{5}.mat(:,find(all_trials)+3); % 4th TR (6-8 sec)
t5 = subj.patterns{5}.mat(:,find(all_trials)+4); % 5th TR (8-10 sec)

mean_data = (t3+t4+t5)/3; %take uniform average of 3 TRs

% BASELINE CORRECT MEAN DATA
% mean_baseline = (t0+t1)/2;
% mean_data = mean_data - mean_baseline;


% plot mean timecourses for whole ROI (useful check of data quality & overall activity pattern)

new_ind = find(objective_new);
for t = 1: length(new_ind)
     new_ts(t,:) = mean(subj.patterns{5}.mat(:,new_ind(t):new_ind(t)+6),1);
end

old_ind = find(objective_old);
for t = 1: length(old_ind)
     old_ts(t,:) = mean(subj.patterns{5}.mat(:,old_ind(t):old_ind(t)+6),1);
end

figure;
plot(mean(old_ts,1),'b');
hold on;
plot(mean(new_ts,1),'r');


% condense regs by removing zeros

trial_counter = 1;
for i = 1: size(regs,2)
    if ~isempty(find(regs(:,i)))
        condensed_regs(:,trial_counter) = regs(:,i);
        condensed_all(:,trial_counter) = all_regs(:,i); %save
        condensed_runs(1,trial_counter) = runs(i);
        trial_counter = trial_counter + 1;
    end
end

subj = init_object(subj,'regressors','conds');
subj = set_mat(subj,'regressors','conds',condensed_regs);
subj = set_objfield(subj,'regressors','conds','condnames',condnames);

% add new condensed activation pattern
subj = duplicate_object(subj,'pattern','spiral_d_hp_z','spiral_d_hp_z_condensed');
subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed',mean_data);

zhist = sprintf('Pattern ''%s'' created by custom code','spiral_d_hp_z_condensed');
subj = add_history(subj,'pattern','spiral_d_hp_z_condensed',zhist,true);

% update run vector to condensed format
subj.selectors{1}.mat = condensed_runs;
subj.selectors{1}.matsize = size(condensed_runs);

% create cross-validation indices
subj = create_xvalid_indices(subj,'runs');

% % remove rest timepoints
%
% temp_sel = ones(1,size(regs,2));
% temp_sel(find(sum(regs)==0)) = 0; % remove all rest timepoints
% subj = init_object(subj,'selector','no_rest');
% subj = set_mat(subj,'selector','no_rest',temp_sel);
% subj = create_xvalid_indices(subj,'runs','actives_selname','no_rest');

% run feature selection ANOVA
subj = feature_select(subj,'spiral_d_hp_z_condensed','conds','runs_xval','thresh',anova_p_thresh);

% run classifier (hidden layer netlab backprop algorithm)
% class_args.train_funct_name = 'train_bp_netlab';
% class_args.test_funct_name = 'test_bp_netlab';
% class_args.nHidden = 10;
% [subj results] = cross_validation(subj,'spiral_d_hp_z','conds','runs_xval',subj.masks{2}.group_name,class_args);

% run classifier (No hidden layer NN Toolbox backprop algorithm)
class_args.train_funct_name = 'train_bp';
class_args.test_funct_name = 'test_bp';
class_args.nHidden = 0;
[subj results] = cross_validation(subj,'spiral_d_hp_z_condensed','conds','runs_xval',subj.masks{2}.group_name,class_args);

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

classification_accuracy_by_resp = [];
for b = 1:10
    classification_accuracy_by_resp(b) = mean(correct_vector(find(condensed_all(b,:))));
end

classification_accuracy_by_resp

[sorted_diffs ind] = sort(acts_diff_vector,2,'descend');
correct_sorted = correct_vector(ind);
desireds_sorted = desireds_vector(ind);

% for i = 1:length(sorted_diffs);
% hit(i) = length(correct_sorted(intersect(find(desireds_sorted == 1),[1:i]))) / length(find(desireds_sorted == 1));
% fa(i) = length(correct_sorted(intersect(find(desireds_sorted == 2),[1:i]))) / length(find(desireds_sorted == 2));
% end

% created binned ROC
counter = 1;
for i = 1:8
hit(i) = length(correct_sorted(intersect(find(desireds_sorted == 1),[1:counter+45]))) / length(find(desireds_sorted == 1));
fa(i) = length(correct_sorted(intersect(find(desireds_sorted == 2),[1:counter+45]))) / length(find(desireds_sorted == 2));
counter = counter+46;
end
figure
plot(fa,hit,'.-')
hold on
plot([0 1],[0 1],'r')
xlabel('P(Old|New)')
ylabel('P(Old|Old)')

display('all done');