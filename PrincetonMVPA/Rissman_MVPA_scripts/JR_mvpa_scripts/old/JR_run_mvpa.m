function [subj results] = JR_run_mvpa()

%%%%%%% user-defined variables %%%%%%%%%%%%%%%%%%

subj_id = 's101';
exp_name = 'PAST';
%roi_filename = 'whole_brain.img';
roi_filename = 'all_hits-cr_p01.nii';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% initialize subj structure
subj = init_subj(exp_name,subj_id);

% load mask file
subj = load_analyze_mask(subj,roi_filename,roi_filename);

% load user-created filename list into workspace

load('raw_filename_list.mat')

% load functional data

subj = load_analyze_pattern(subj,'epi',roi_filename,raw_filenames);


% get trial/run onset information from spm.mat
[subj, regs, runs] = JR_convert_spm_onsets_to_mvpa_regs(subj);


subj = detrend_runs(subj,'epi','runs');  % not in mvpa tutorial, but seems important to do
subj = hpfilter_runs(subj,'epi_d','runs',100,2); % remove frequencies below .01 Hz
subj = zscore_runs(subj,'epi_d_hp','runs'); % Z-score the data


% create cross-validation indices
%subj = create_xvalid_indices(subj,'runs');

% remove rest timepoints

temp_sel = ones(1,size(regs,2));
temp_sel(find(sum(regs)==0)) = 0; % remove all rest timepoints
subj = init_object(subj,'selector','no_rest');
subj = set_mat(subj,'selector','no_rest',temp_sel);
subj = create_xvalid_indices(subj,'runs','actives_selname','no_rest'); 

% run feature selection ANOVA
subj = feature_select(subj,'epi_d_hp_z','conds','runs_xval');
%subj = feature_select(subj,'epi_z','conds','runs_xval','thresh',.001);

% run classifier (hidden layer netlab backprop algorithm)
class_args.train_funct_name = 'train_bp_netlab';
class_args.test_funct_name = 'test_bp_netlab';
class_args.nHidden = 1;
[subj results] = cross_validation(subj,'epi_d_hp_z','conds','runs_xval',subj.masks{2}.group_name,class_args);

% % run classifier (No hidden layer NN Toolbox backprop algorithm)
% class_args.train_funct_name = 'train_bp';
% class_args.test_funct_name = 'test_bp';
% class_args.nHidden = 0;
% [subj results] = cross_validation(subj,'epi_z_d','conds','runs_xval',subj.masks{2}.group_name,class_args);