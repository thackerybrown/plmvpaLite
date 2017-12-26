function [subj results] = JR_run_mvpa_TR_averaging

%%%%%%% user-defined variables %%%%%%%%%%%%%%%%%%

subj_id = 's102';
exp_name = 'PAST';
roi_filename = 'bilat_IT_H-FA_p0001.nii';
%roi_filename = 'all_hits-cr_p01.nii'; %'new>old_t>3.img';
anova_p_thresh = .1;


num_runs = 10;
num_TP_per_run = 203;
total_TP = num_runs * num_TP_per_run;

% load user-created filename and onsets lists into workspace
load('raw_filename_list.mat')
load('/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/s102/onsets.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initialize subj structure
subj = init_subj(exp_name,subj_id);

% load mask file
subj = load_analyze_mask(subj,roi_filename,roi_filename);

% load functional data

subj = load_analyze_pattern(subj,'epi',roi_filename,raw_filenames);


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
old = sum(all_regs(1:5,:));
new = sum(all_regs(6:10,:));
no_resp = all_regs(11,:);

regs(1,:) = old;
regs(2,:) = new;
%  regs(3,:) = no_resp;  % treat no_resp trials as rest, which will later get excluded from analysis

condnames = {'old','new'};

% make runs vector
subj = init_object(subj,'selector','runs'); 

trial_idx = 1;
for r = 1:num_runs
 
    runs(trial_idx:203*r) = r;
    trial_idx = trial_idx + num_TP_per_run;
end

subj = set_mat(subj,'selector','runs',runs);

% filter the timeseries data
subj = detrend_runs(subj,'epi','runs');  % not in mvpa tutorial, but seems important to do
subj = hpfilter_runs(subj,'epi_d','runs',100,2); % remove frequencies below .01 Hz
subj = zscore_runs(subj,'epi_d_hp','runs'); % Z-score the data


% TR averaging
all_trials = sum(regs,1); % vector of all trial onsets

t1 = subj.patterns{4}.mat(:,find(all_trials)+1); % 2nd TR (2-4 sec)
t2 = subj.patterns{4}.mat(:,find(all_trials)+2); % 3rd TR (4-6 sec)
t3 = subj.patterns{4}.mat(:,find(all_trials)+3); % 4th TR (6-8 sec)
t4 = subj.patterns{4}.mat(:,find(all_trials)+4); % 4th TR (8-10 sec)

%mean_data = (t1+t2+t3)/3;
mean_data = (t2+t3+t4)/3;

% plot mean timecourses

new_ind = find(new);
for t = 1: length(new_ind)
     new_ts(t,:) = mean(subj.patterns{4}.mat(:,new_ind(t):new_ind(t)+6),1);
end
    
old_ind = find(old);
for t = 1: length(old_ind)
     old_ts(t,:) = mean(subj.patterns{4}.mat(:,old_ind(t):old_ind(t)+6),1);
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
        condensed_runs(1,trial_counter) = runs(i);
        trial_counter = trial_counter + 1;
    end
end

subj = init_object(subj,'regressors','conds');
subj = set_mat(subj,'regressors','conds',condensed_regs);
subj = set_objfield(subj,'regressors','conds','condnames',condnames);

% add new condensed activation pattern
subj = duplicate_object(subj,'pattern','epi_d_hp_z','epi_d_hp_z_condensed');
subj = set_mat(subj,'pattern','epi_d_hp_z_condensed',mean_data);

zhist = sprintf('Pattern ''%s'' created by zscore_runs','epi_d_hp_z_condensed');
subj = add_history(subj,'pattern','epi_d_hp_z_condensed',zhist,true);

% created.function = 'condensed_data';
% created.use_mvpa_ver = args.use_mvpa_ver;
% created.patname = patname;
% created.selname = selname;
% created.actives_selname = args.actives_selname;
% subj = add_created(subj,'pattern',args.new_patname,created);


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
subj = feature_select(subj,'epi_d_hp_z_condensed','conds','runs_xval','thresh',anova_p_thresh);
%subj = feature_select(subj,'epi_z','conds','runs_xval','thresh',.001);

% run classifier (hidden layer netlab backprop algorithm)
% class_args.train_funct_name = 'train_bp_netlab';
% class_args.test_funct_name = 'test_bp_netlab';
% class_args.nHidden = 10;
% [subj results] = cross_validation(subj,'epi_d_hp_z','conds','runs_xval',subj.masks{2}.group_name,class_args);

% run classifier (No hidden layer NN Toolbox backprop algorithm)
 class_args.train_funct_name = 'train_bp';
 class_args.test_funct_name = 'test_bp';
 class_args.nHidden = 0;
 [subj results] = cross_validation(subj,'epi_d_hp_z_condensed','conds','runs_xval',subj.masks{2}.group_name,class_args);