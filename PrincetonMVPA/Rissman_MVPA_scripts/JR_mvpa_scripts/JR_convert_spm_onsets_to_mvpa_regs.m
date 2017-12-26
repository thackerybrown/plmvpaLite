function [subj, regs, runs] = convert_spm_onsets_to_mvpa_regs(subj)

num_runs = 8;
num_TP_per_run = 203;
total_TP = num_runs * num_TP_per_run;

load('/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/s101/onsets.mat')

num_conds = size(onsets,2);
all_regs = zeros(num_conds,total_TP); % initialize regs matrix as conditions x timepoints

for cond = 1: num_conds
    for trial = 1: length(onsets{cond})
        time_idx = onsets{cond}(trial)/2+1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
        all_regs(cond,time_idx+2) = 1;  % use data from 1 TRs (4 sec post-stimulus)
        %all_regs(cond,time_idx+2:time_idx+3) = 1;  % use data from 2 TRs (4 and 6 sec post-stimulus)
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


subj = init_object(subj,'regressors','conds');
subj = set_mat(subj,'regressors','conds',regs);
subj = set_objfield(subj,'regressors','conds','condnames',condnames);


% Create 'runs' vector. This is stored simply as a row-vector with the same number of timepoints as the regressors, containing 1s, 2s, 3s, etc. 
% This goes in a selector object (think of these as a set of 'timelabels'), which we'll call 'runs':

subj = init_object(subj,'selector','runs'); 

trial_idx = 1;
for r = 1:num_runs
 
    runs(trial_idx:203*r) = r;
    trial_idx = trial_idx + num_TP_per_run;
end

subj = set_mat(subj,'selector','runs',runs);




