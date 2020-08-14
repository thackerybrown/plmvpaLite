function [S idxTr idxTe par]= TIB_mvpa_params_general_CIT(subj_id, task, TRsperRun, imgtype)
% Created for PSYCH 8080, Spr 2018

% establish parameters for mvpa analysis
% <subj_id> - identifier for the given subject. Can be numerical or a
% string
% <task> 'goals' or 'plan'
%
% Study paths, scan parameters, and subject parameters are set up here for
% the analysis in several dat structures (e.g., S.xX; par.xX)

% Important - 09/2018 - existpatmat scenario has a variety of quirks and
% hardcodes. If using this function, check the comments and relevant
% scripts to ensure proper use

% Most variables are assigned to either par or S. At this time, it's not well
% organized what info goes into what.

%% EDIT - You must establish these general parameters

% ~~~ WHAT IS YOUR *Study name* (code looks for this specific folder) ~~~
S.exp_name = 'CIT'; %change this to flexibly redirect the script to different studies in subdirectories

% ~~~ WHAT IS YOUR *subject ID/number PREFIX* (code appends this to the
% front of the sub number provided with the run_mvpa_general function call)
subprefix = 'anlz_';

% ~~~ WHAT IS YOUR *scan TR in the units of your model file (usually seconds)*
% NOTE: if working with Betas - if your "onsets" file is really a numerical list of beta numbers instead of onsets, you must set TR = 1.
% If your onsets file is in seconds, set this as the TR ***even for Betas***.
par.TR = 2;

% ~~~ WHAT ARE YOUR *image dimensions* (4D or 3D files?)
% As of 12/31/17, code only supports 3D images.
% It is planned to support use 4D nifti files for raw BOLD data ('4').
% If you have split them out into TR-by-TR, enter '3'
par.ImgDims = 3;

% ~~~ DO YOU WANT *to read in the BOLD image file names on-the-fly or do you
% have an existing list of image names? On-the-fly adds a little time; but
% if you have it off, and make a change to the images, do must update your
% image list file... otherwise you'll run the old classification problem
% again
par.readimglist = 0; %1=yes, read an existing list; 0 = no, generate on the fly please (slower but recommended).

% ~~~ DO YOU WANT *to read in an existing extracted workspace (e.g., BOLD
% patterns). Not the same as having a simple pre-defined matrix ready for
% classification with one column per pattern (see existpatmat flag on next line)
S.use_premade_workspace = 0;

% ~~~ DO YOU WANT *to read in a .mat file with existing patterns in it
% already (instead of fMRI image files)?
S.existpatmat = 0; %1=yes - skip trying to load image files using SPM. We've already got all patterns in a matrix. Currently (09/2018) the existing pattern matrix is hardcoded in load_matrix_pattern_2D to be a matrix (double) named 'testmat'. \\
% Because so much of the Princeton MVPA toolbox assumes a mask volume is used, a dummy mask file is now included in the PLMVPA_Lite toolkit and called in TIB_run_MVPA_general
% What's that existing pattern .mat file called?
S.datafile = 'testsubmat_fixed.mat'; % added for ADNI neuropsych study. Replaces MRI image (e.g., nii) files with an existing matrix of pattern data

% ~~~ HOW MANY *single trial betas* are there? If data input type is 'raw', number is not used (can be ignored)
S.stbetacount = 100; % NOTE: the code assumes all single trial betas of potential interest are contiguous in the model. If multi-event betas are inteleaved in the .mat model structure, this will need more editing.

%% EDIT - You must establish parameters for this SPECIFIC classification scenario involving the data described in the preceding section

% ~~~ what study conditions or phases do you want to *TRAIN* on?
S.trainTask = 'TvsI';%Circmaze - 'goals' or 'plan'

% ~~~ WHAT is the ROI we're analyzing
S.roi_name = 'all3rois_nowm';%'grp_hippovis'; %S.roi_name = 'HVisCtx_1.nii'; %S.roi_name = 'NativeGM_BOLDres.nii';

% ~~~ what study conditions or phases do you want to *TEST* on?
% NOTE: if the string here is not the same as S.trainTask, the classifier
% will switch to a 'tr1teo' procedure (train one phase, test on the other).
% When would tr1teo be useful? e.g., if you want to train on a functional localizer, and test on a memory retrieval task
S.testTask = 'CMpcmvsCMcmi';%Circamze - 'goals' or 'plan'

% ~~~ what cross-validation procedure do you want? *ignored if S.testTask
% and S.trainTask are not the same
% 'loo' = leave one scan run / data collection bin out for testing. HIGHLY
% RECOMMENDED
% 'nf' = random nfold - generates "pseudo runs" if you want to do cross
% validation but don't want to use existing run structure (e.g., you only
% have one run, there's no such thing as a "run" in your data
S.xvaltype = 'nf';
S.nFolds = 100; % number of cross validation iterations - only used for nFold (as opposed to run-by-run loo)

% ~~~ WHAT IS THE NAME *of your model.mat file* for this analysis?
% NOTE: be default JUST put the end of the name; code will append the
% subject ID in front of this. Edit as appropriate
S.trainonsfname = '_localizer_onsets_test';
S.testonsfname = '_localizer_onsets_test';

S.trainonsfnamebetas = '_CIT_ids';
S.testonsfnamebetas = '_CIT_ids';

% ~~~ WHAT is the preprocessing level of your BOLDs (if input type 'raw')
par.preproc_lvl = ''; % 'a' for slice-time-only, 'u' for realigned-only, 'ua' for realign+unwarped, 'swua' for smoothed, normalized, and... you get the picture. Modify as needed if you changed SPM's prefix append defaults
par.boldnames = [par.preproc_lvl 'run']; %name of image files with preprocessing level prefix
par.runnames = 'run_*'; %what are your run/scan session folders called? If none or only 1 run, set to '' or whatever may be appropriate.
par.imageextension = '.nii'; %are your images .nii or .img?

% ~~~ WHAT is a reference BOLD image we can look to for image dimenstions
% (EVEN if doing Betas)
par.refrun = 'modfold'; %just the run number - we'll fill in the name prefix below automatically
par.ref_funcimage = ['con_0001_01' par.imageextension];

%% EDIT - Path customization (can be changed more below, but not recommended - try to maintain the directory structure and just change these paths for consistency

% ~~~ WHAT IS YOUR *computer base path* (where your study and its subfolders
% live
S.sbasepath = '/mnt/hgfs/Work/mvpa_sample_data/';

% ~~~ WHAT IS THE NAME *of your model folder's directory*?
S.modfold = 'modfold';
S.modfold_singlebetas = 'modfold'; % if your single trial (or person, ...
%...for between-sub classification) betas are in their own folder

% ~~~ WHAT IS THE NAME *of your Masks folder* where your ROIs live
S.maskdir = 'Masks';

% ~~~ WHAT IS THE NAME *of your BOLDs folder*
S.boldsdir = '';

% ~~~ WHAT IS THE NAME *of your BOLD run folder prefix*
% NOTE: be default JUST put the start of the name; code will append the
% specific number to the end of this. Edit as appropriate
par.boldrundirpfx = '';

%% EDIT - Classifier tuning and various settings

% ~~~ Iteration Parameters
S.num_results_iter = 3; % number of times to run the entire classification process (select subset of the data and train/test classifier)
S.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data - useful for non-deterministic cases.

% ~~~ Balancing Parameters
S.equate_number_of_trials_in_groups = 0; % equate number of trials in conditions
S.numBalancedParams = 1; % number of parameters to balance across (e.g., both goal location AND cue in Circmaze data). The code currently (12/29/17) only handles two options - 1 (standard; main class type), or 2 (main class type plus a second parameter, specified in a second file).
S.numBalancedIts = 10; % number of iterations to run, with different randomization for the balancing

% ~~~ Secondary Z-Scoring
% note: by default all raw BOLDs and existpatmats are subjected to a FIRST
% round of z-scoring across all timepoints. Betas are by default not. 
S.perform_second_round_of_zscoring = 0;  % z-score data again immediately prior to classification

% ~~~ Noising/Denoising
S.addnoise = 0; %Overly sparse data (lots of features that often have zeroes)? Add gaussian white noise to pattern matrix to help with overfitting. >0 = yes. Value specified = target SNR.
%S.denoise = 0; %undergo denoising?
%S.denoiseOpt.denoisespec = '10001'; %which parts of the glm output do we want to save?

% ~~~ Signal intensity analysis
S.thisSigIntenseSelector = 'randomNFold_xval'; %which selector to use for signal intensity analysis
S.zscoreIntensityVals = 1; % zscore the intensity values?

% ~~~ Mean Signal Extraction Params
% parameters for selecting the mean signal from a class-specific ROI for each pattern.
S.extractMeanSignal = 0; %1 - do signal extraction. 0 = don't do this.
S.defineROIsFromANOVAFS = 0; % define ROIs using ANOVA-based feature selection, instead of pre-defining them.
%S.logreg_2Features = 0; %perform a logistic regression, using the two extracted intensity vectors
%**********************************************************************************come back to this section once up and running......%%%%
% S.ROI1PatName = [S.preprocPatCondensedName '_ROI1'];
% S.ROI1_name = [ 'occipitoTemporal_faceVsScene_500vox.img'];
% S.ROI1_file  = [par.subdir '/analysis_loc_mnem/' S.ROI1_name];
%
% S.ROI2PatName = [S.preprocPatCondensedName '_ROI2'];
% S.ROI2_name  = ['occipitoTemporal_sceneVsFace_500vox.img'];
% S.ROI2_file   = [par.subdir '/analysis_loc_mnem/' S.ROI2_name];

% ~~~ Importance Maps
S.generate_importance_maps = 1; %visualize classifier weights
S.generateBetaMaps = 1; %use betas, instead of importance values
S.impType = {'pos' 'neg' 'both' 'raw'}; %importance map types
S.regNames = {'CondA' 'CondB'}; % should match number of classes

% ~~~ TR Weighting
%which post-stimulus TRs should be used (and if more than one, averaged
%across) before feeding data to the classifier?  S.TR_weights_set{1} -
%training weights, S.TR_weights_set{2} = testing weights

%NOTE: for beta analyses, we don't average over multiple images because
%different images = different events
S.inputformat = imgtype; %assign input from function call. Either 'raw' for raw bold images or 'betas' for beta images. Selection here automatically changes some params below.
if strcmp(S.inputformat, 'raw')
    %S.TR_weights_set = {{[.0072 .2168 .3781 .2742 .1237] [.0072 .2168 .3781 .2742 .1237]}}; %approximates the canonical haemodynamic response
    S.TR_weights_set = {{[0 0.25 0.5 0.25] [0 0.25 0.5 0.25]}};%use double-bracket structure in case want to set code up to run a sliding window across multiple TR bins
elseif strcmp(S.inputformat, 'betas')
    S.TR_weights_set = {{[1] [1]}};%give full weighting to the 1 and only image corresponding to each event
end

% ~~~ Special types of analysis
S.searchlightAnalysis = 0; % run a searchlight analysis
%S.linReg = 0; % run an analysis with a continuous outcome variable
S.scrambleregs = 1; % run an anlysis with the class labels scrambled on a run-by-run basis.

% ~~~ classifier parameters
S.class_args.train_funct_name = 'train_liblinear_multiclass';%'train_pLR';   %training function
S.class_args.test_funct_name = 'test_liblinear_multiclass';%'test_pLR';      %testing function
S.class_args.classType = 'libLin';
S.perfmet_functs = 'perfmet_maxclass'; % performance metric
S.statmap_funct = 'statmap_anova';%'AG_statmap_anova'; % performance metric

S.class_args.nVox = 0; % number of voxels to select with feature selection e.g. [1000 5000 10000]
S.class_args.fseltype = 'topn'; % feature selection format: top N vox (topn) or random N vox (rand)?
S.class_args.libLin = '-q -s 0 -B 1'; %arguments for liblinear; -s 0 = L2; -s 6 = L1; -s 5 = L1 with L2 loss; -s 3 L2 with L1 loss
S.class_args.constant = true; % include a constant term?
S.class_args.prefitWeights = true;

S.class_args.chooseOptimalPenalty = 0; % 1 = yes. cycle through cost parameters in the training set, and choose the optimal one. Note, this only makes sense in context of loo with >2 runs or for nf with >2 folds, because it subdivides training set into additional 'runs' and performs nested xvalidation.
S.class_args.penaltyRange = [.001 .005 .01 .05 .1 .5 1 5 10 50 100 500 1000 50000]; % a vector "[]" of cost parameters to cycle through
S.class_args.nFoldsPenaltySelection = 10; % number of cross validation folds for penalty parameter selection.

S.class_args.penalty = 1; %uncomment if not using optimal penalty. Typical value is 1. If using sample data provided with plmvpaLite, start with 0.000001 to see how minimal regularization harms performance.
%establishment

%% Default and auto-generated parameters. **Only change if you must for your specific use case to work**
%Functional image scan selectors ***if ALL FILENAMES corresponding to ALL RUNS OF INTEREST are stored in ONE cell...
% ...of raw_filenames.mat (i.e., not broken up by run), set index to 1 or 1:1. Otherwise, create indexing for elements...
% ...of cell array raw_filenames.mat corresponding to task of interest (i.e. if cells runs 1:4 correspond to task 1, we want...
% ....to reference {1}, {2}... in raw_filenames.mat)
par.scansSelect.(S.exp_name).loc = 1:1;%par.scansSelect.goals.loc = 1:1;%
par.scansSelect.(S.exp_name).loc = 1:1;%par.scansSelect.plan.loc = 1:1;%

% how many Timepoints are in each run
par.TRsperRun = TRsperRun;

% ~~~ WHAT IS YOUR *Subject ID/number* (fills in prefix, code fills in the
% rest from run_mvpa_general function call)
par.substr = [subprefix subj_id{1}];
S.subj_id = par.substr;

% Task type
par.task = task; %assign input from function call. Task phase label. For circmaze, this is 'goals' or 'plan'. For localizer (8080), this is 'CM_Localizer'

%input image info
if strcmp(S.inputformat, 'raw')
    data_imgs_to_use = 'raw_filenames.mat';
elseif strcmp(S.inputformat, 'betas')
    data_imgs_to_use = 'beta_filenames.mat'; % analogous to raw_filenames, but with a boolean index for which betas correspond to which conditions. Created by calling TIB_generate_beta_filenames.mat below
end

S.preprocType = 'spm'; % 'spm' for spm preprocessing, 'knk' for kendrick preprocessing

%%model information - define which timepoints or images correspond to which classes of data
if strcmp(S.inputformat, 'raw')
    S.onsets_filename = [S.subj_id S.trainonsfname];%
    S.onsets_filename_tr = [S.subj_id S.trainonsfname];% added for train on 1 phase, test on another - this assumes the data are actually in the same set of files.
    S.onsets_filename_tst = [S.subj_id S.testonsfname];% added for train on 1 phase, test on another - this assumes the data are actually in the same set of files.
elseif strcmp(S.inputformat, 'betas')
    S.onsets_filename = [S.subj_id S.trainonsfnamebetas];
    S.onsets_filename_tr = [S.subj_id S.trainonsfnamebetas];
    S.onsets_filename_tst = [S.subj_id S.testonsfnamebetas];
    
    S.betaidx_filename = [S.subj_id '_betas_idx'];
    S.betaidx_filename_tr = [S.subj_id '_betas_idx_tr'];
    S.betaidx_filename_te = [S.subj_id '_betas_idx_te'];
end

%% autospecified and generated directories. **Only change if you must for your specific use case to work**
S.expt_dir = [S.sbasepath S.exp_name '/'];%study location

par.subdir =[S.expt_dir S.subj_id];%subject location

par.funcdir =[par.subdir '/' S.boldsdir '/'];%subfolder for 'raw' BOLD data. Assumes BOLDs are stored in subfolders labeled 'run_01', etc)

S.workspace_dir = [par.subdir '/mvpa_workspace'];%temporary files workspace; stores loaded pats that can be used to speed up subsequent classifications if you don't want to re-load from image files

%model file directory (onsets.mat and betas in here) - this is still used
%when working with raw data. We must have some way to tell the classifier
%which images correspond to which classes
if strcmp(S.inputformat, 'raw')
    S.mvpa_dir = [S.expt_dir S.subj_id '/' S.modfold '/'];
elseif strcmp(S.inputformat, 'betas')
    S.mvpa_dir = [S.expt_dir S.subj_id '/' S.modfold_singlebetas '/'];
end

%ROI masks (could be whole-brain mask, but the code wants a mask file
S.anat_dir = [S.expt_dir S.subj_id '/' S.maskdir];

S.importance_maps_dir=[S.expt_dir 'ImpMaps_' date];
S.group_mvpa_dir = [S.expt_dir 'mvpa_output_files'];%results .mat files are spit out in here

%% identify betas
if strcmp(S.inputformat, 'betas')%if we are running a beta analysis, take a moment to create "model" files for the beta images
    TIB_generate_beta_filenames_CIT(S);
end

%% Initialize empty vars to be filled below
idxTr = [];
idxTe = [];

%% create raw_filenames on the fly (as opposed to making it by loading structures from an SPM.mat file [for this method, see comment above])

%if fnames are wrong and you want to standardize, explore 'rename'
%function, ex.: "rename -v s/CM001_localizer02_0/run_02_/ *.nii"

if strcmp(S.inputformat, 'raw')
    if par.readimglist == 1
        load([par.funcdir '/' data_imgs_to_use]); %loads predefined cell array
    else %generate them on the fly
        raw_filenames = gen_raw_boldfnames(par, S);
    end
    
else %if using betas...
    raw_filenames = [];
end

if strcmp(S.inputformat, 'betas')
    load([S.mvpa_dir '/' data_imgs_to_use]); %loads predefined cell array
    %called raw_filenames or beta_filenames into memory
end

%% cross-validation scheme
%this code directs the training/testing toward either a cross validation
%across all trials of a given set (normal cross validation) or training on
%one set of data and testing on the other (for example if you had a
%localizer or encoding task and then wanted to test retrieval patterns
if strcmp(S.trainTask, S.testTask)
    S.xval = 1;
    disp(['Cross-validation type is ' S.xvaltype])
    if strcmp(S.xvaltype, 'nf')
        S.thisSelector =  'randomNFold_xval'; % nfold cross validation selector
    elseif strcmp(S.xvaltype, 'loo')
        S.thisSelector =  'leave_one_out_xval'; % leave-one-out cross validation selector
    end
else
    S.xval = 0;
    S.thisSelector = 'TrainTestOneIterGroup'; % train on one group, test on another
    disp(['Instead of cross-validation, running ' S.thisSelector])
end

%% EDIT - Information about Conditions/Classes specific to the training and testing of specific subsets of data

% S.onsetsTrainDir = location of training onsets

% S.condsTrain = conditions on which to train

% S.dnCondsTrain = conditions which which to denoise, if denoising is used

% S.TrainRuns = runs of data on which to train (this is NOT x-validation -
% this is selecting subsets of runs pertinent to analysis (assuming not all
% runs are. For example if the first half of the runs pertained to a
% different task from the other half of the runs)

% S.durTrain = duration, in SECONDS, of runs used for training - this is
% the total number of usable trs from the runs of interest, multiplied by
% TR time to convert into seconds (e.g. for [165 156 153] + TR =2,
% S.durTrain = 948s)

% S.filenames_train = names of data images to use for training

% idxTr = behavioral indices for training task, used by TIB_run_MVPA_general

if strcmp(S.trainTask,'TvsI')
    S.onsetsTrainDir = [S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTrain = {{'t'}  {'i'}} ;%corresponds to the names in the onsets.mat or betas_idx.mat files. This is used to select what is being compared with what.
    S.TrainRuns = par.scansSelect.(par.task).loc;%pull up indexing, defined above, for RUNS corresponding to task of interest (i.e. if runs 2,4,6 correspond to task 1)
    if strcmp(S.inputformat, 'raw')
        S.filenames_train = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_train = beta_filenames;%
    end
    S.durTrain = numel(S.filenames_train) * par.TR;
    %[~, idxTr] = fMRIBehAnalysis_Loc(par);
    
elseif strcmp(S.trainTask,'EAvsAA')
    S.onsetsTrainDir = [S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTrain = {{'EA'}  {'AA'}} ;%corresponds to the names in the onsets.mat or betas_idx.mat files. This is used to select what is being compared with what.
    S.TrainRuns = par.scansSelect.(par.task).loc;%pull up indexing, defined above, for RUNS corresponding to task of interest (i.e. if runs 2,4,6 correspond to task 1)
    if strcmp(S.inputformat, 'raw')
        S.filenames_train = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_train = beta_filenames;%
    end
    S.durTrain = numel(S.filenames_train) * par.TR;
    
elseif strcmp(S.trainTask,'AAvsScene')
    S.onsetsTrainDir = [S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTrain = {{'AA'}  {'Scene'}} ;%corresponds to the names in the onsets.mat or betas_idx.mat files. This is used to select what is being compared with what.
    S.TrainRuns = par.scansSelect.(par.task).loc;%pull up indexing, defined above, for RUNS corresponding to task of interest (i.e. if runs 2,4,6 correspond to task 1)
    if strcmp(S.inputformat, 'raw')
        S.filenames_train = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_train = beta_filenames;%
    end
    S.durTrain = numel(S.filenames_train) * par.TR;
    
elseif strcmp(S.trainTask,'AAvsObj')
    S.onsetsTrainDir = [S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTrain = {{'AA'}  {'Obj'}} ;%corresponds to the names in the onsets.mat or betas_idx.mat files. This is used to select what is being compared with what.
    S.TrainRuns = par.scansSelect.(par.task).loc;%pull up indexing, defined above, for RUNS corresponding to task of interest (i.e. if runs 2,4,6 correspond to task 1)
    if strcmp(S.inputformat, 'raw')
        S.filenames_train = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_train = beta_filenames;%
    end
    S.durTrain = numel(S.filenames_train) * par.TR;
    
elseif strcmp(S.trainTask,'AAvsScrambled')
    S.onsetsTrainDir = [S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTrain = {{'AA'}  {'AA_scrambled'}} ;%corresponds to the names in the onsets.mat or betas_idx.mat files. This is used to select what is being compared with what.
    S.TrainRuns = par.scansSelect.(par.task).loc;%pull up indexing, defined above, for RUNS corresponding to task of interest (i.e. if runs 2,4,6 correspond to task 1)
    if strcmp(S.inputformat, 'raw')
        S.filenames_train = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_train = beta_filenames;%
    end
    S.durTrain = numel(S.filenames_train) * par.TR;
    
elseif strcmp(S.trainTask,'EAvsAAvsScene')
    S.onsetsTrainDir = [S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTrain = {{'EA'} {'AA'} {'Scene'}} ;%corresponds to the names in the onsets.mat or betas_idx.mat files. This is used to select what is being compared with what.
    S.TrainRuns = par.scansSelect.(par.task).loc;%pull up indexing, defined above, for RUNS corresponding to task of interest (i.e. if runs 2,4,6 correspond to task 1)
    if strcmp(S.inputformat, 'raw')
        S.filenames_train = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_train = beta_filenames;%
    end
    S.durTrain = numel(S.filenames_train) * par.TR;
    
elseif strcmp(S.trainTask,'FacevsScenevsObj')
    S.onsetsTrainDir = [S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTrain = {{'Face'}  {'Scene'} {'Obj'}} ;%corresponds to the names in the onsets.mat or betas_idx.mat files. This is used to select what is being compared with what.
    S.TrainRuns = par.scansSelect.(par.task).loc;%pull up indexing, defined above, for RUNS corresponding to task of interest (i.e. if runs 2,4,6 correspond to task 1)
    if strcmp(S.inputformat, 'raw')
        S.filenames_train = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_train = beta_filenames;%
    end
    S.durTrain = numel(S.filenames_train) * par.TR;
end

% testing - this defines the testing set. The code is set up this way to enable us to step outside xval if desired to test on different set of data (e.g., at retrieval)
if strcmp(S.testTask,'TvsI')
    S.onsetsTestDir =[S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTest = {{'t'} {'i'}};
    S.nwayclass = num2str(numel(S.condsTest));%stores the number classification dimensions just for reference (i.e. is this a 5-way or a 2-way/binary classification?)
    S.TestRuns = par.scansSelect.(par.task).loc;
    if strcmp(S.inputformat, 'raw')
        S.filenames_test = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_test = beta_filenames;%
    end
    S.durTest = numel(S.filenames_test) * par.TR;
    %[~, idxTe] = fMRIBehAnalysis_Loc(par);
    
elseif strcmp(S.testTask,'CMpcmvsCMcmi')
    S.onsetsTestDir =[S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTest = {{'cmpcm'} {'cmcmi'}};
    S.nwayclass = num2str(numel(S.condsTest));%stores the number classification dimensions just for reference (i.e. is this a 5-way or a 2-way/binary classification?)
    S.TestRuns = par.scansSelect.(par.task).loc;
    if strcmp(S.inputformat, 'raw')
        S.filenames_test = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_test = beta_filenames;%
    end
    S.durTest = numel(S.filenames_test) * par.TR;
    
elseif strcmp(S.testTask,'EAvsAA')
    S.onsetsTestDir =[S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTest = {{'EA'} {'AA'}};
    S.nwayclass = num2str(numel(S.condsTest));%stores the number classification dimensions just for reference (i.e. is this a 5-way or a 2-way/binary classification?)
    S.TestRuns = par.scansSelect.(par.task).loc;
    if strcmp(S.inputformat, 'raw')
        S.filenames_test = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_test = beta_filenames;%
    end
    S.durTest = numel(S.filenames_test) * par.TR;
    
elseif strcmp(S.testTask,'AAvsScene')
    S.onsetsTestDir =[S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTest = {{'AA'} {'Scene'}};
    S.nwayclass = num2str(numel(S.condsTest));%stores the number classification dimensions just for reference (i.e. is this a 5-way or a 2-way/binary classification?)
    S.TestRuns = par.scansSelect.(par.task).loc;
    if strcmp(S.inputformat, 'raw')
        S.filenames_test = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_test = beta_filenames;%
    end
    S.durTest = numel(S.filenames_test) * par.TR;
    
elseif strcmp(S.testTask,'AAvsObj')
    S.onsetsTestDir =[S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTest = {{'AA'} {'Obj'}};
    S.nwayclass = num2str(numel(S.condsTest));%stores the number classification dimensions just for reference (i.e. is this a 5-way or a 2-way/binary classification?)
    S.TestRuns = par.scansSelect.(par.task).loc;
    if strcmp(S.inputformat, 'raw')
        S.filenames_test = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_test = beta_filenames;%
    end
    S.durTest = numel(S.filenames_test) * par.TR;
    
elseif strcmp(S.testTask,'AAvsScrambled')
    S.onsetsTestDir =[S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTest = {{'AA'} {'AA_scrambled'}};
    S.nwayclass = num2str(numel(S.condsTest));%stores the number classification dimensions just for reference (i.e. is this a 5-way or a 2-way/binary classification?)
    S.TestRuns = par.scansSelect.(par.task).loc;
    if strcmp(S.inputformat, 'raw')
        S.filenames_test = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_test = beta_filenames;%
    end
    S.durTest = numel(S.filenames_test) * par.TR;
    
elseif strcmp(S.testTask,'EAvsAAvsScene')
    S.onsetsTestDir =[S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTest = {{'EA'} {'AA'} {'Scene'}};
    S.nwayclass = num2str(numel(S.condsTest));%stores the number classification dimensions just for reference (i.e. is this a 5-way or a 2-way/binary classification?)
    S.TestRuns = par.scansSelect.(par.task).loc;
    if strcmp(S.inputformat, 'raw')
        S.filenames_test = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_test = beta_filenames;%
    end
    S.durTest = numel(S.filenames_test) * par.TR;
    
elseif strcmp(S.testTask,'FacevsScenevsObj')
    S.onsetsTestDir =[S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTest = {{'Face'} {'Scene'} {'Obj'}};
    S.nwayclass = num2str(numel(S.condsTest));%stores the number classification dimensions just for reference (i.e. is this a 5-way or a 2-way/binary classification?)
    S.TestRuns = par.scansSelect.(par.task).loc;
    if strcmp(S.inputformat, 'raw')
        S.filenames_test = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_test = beta_filenames;%
    end
    S.durTest = numel(S.filenames_test) * par.TR;
    
end

S.condnames = S.condsTrain;
S.regName = 'conds';


%% Smoothing Parameters - ******************integrate with switch mode above
S.funcType = 3;
if strcmp(S.inputformat, 'betas')%override hard code here if we're playing with betas
    S.funcType = 4;
end

S.smoothTxt = { 'unsmoothed' 'smoothed' 'native' 'betas'};
switch S.funcType
    case 1
        par.filesForPatterns = par.wrascanfiles.all;
    case 2
        par.filesForPatterns = par.swrascanfiles.all;
    case 3
        %par.filesForPatterns = par.rascanfiles.all;
        par.filesForPatterns = raw_filenames;
    case 4
        par.filesForPatterns = beta_filenames;
end

%% specify which files to load for classification
if S.xval%if we are doing nfold xval (automatically set above via a 1)
    S.filenames = S.filenames_train;
else
    S.filenames_h{1} = S.filenames_train;
    S.filenames_h{2} = S.filenames_test;
    S.filenames = S.filenames_h{1};%char(S.filenames_h);
end
S.img_files =  mat2cell(S.filenames, [ones(1,size(S.filenames,1))], [size(S.filenames,2)]);

%% Runs Parameters

S.runs_vector = par.TRsperRun; % number of TRs per scanning run (coded for vector of values)

S.meta_runs = S.runs_vector;
S.num_runs = length(S.runs_vector);
S.num_vols = sum(S.runs_vector);
S.TR = par.TR;

%create an index of from which run each beta comes
if strcmp(S.inputformat, 'betas')
    onsetsmat = [S.mvpa_dir S.onsets_filename];
    load(onsetsmat);
    %mask sessionlength to onsets that are single trial regs in the model.
    %ASSUMES that single trial regressors are continuous - e.g., the first
    %n onsets values are those we will classify over. Modify code further
    %if multi-event regressors are interleaved with single-trial ones.
    if S.existpatmat==1
        idxTr.sess(1:length(onsets))=1;%NOTE: modify if needed, but for ADNI we only have one "run" of data for xval
        %temporary for heaton
        %         idxTr.sess(1:44)=1;
        %         idxTr.sess(45:length(onsets))=2;
        
    else
        idxTr.sess = onsets(1:S.stbetacount);
        
        varswecareabout = length(idxTr.sess);
        onsetscount = 0;
        for r = 1:length(S.runs_vector)
            onsetscount = onsetscount + (S.runs_vector(r)*S.TR);%sums the time (s) from the current run with previously iterated through runs to set an increasing threshold
            %s.(sprintf('x%d', r)) =
            %x1= find(S.idxTr.sess)% < 500)
            runsc.(sprintf('x%d', r)) = find([idxTr.sess{:,1:varswecareabout}] <= onsetscount);
            
        end
        for r = 1:length(S.runs_vector)%creates bidx(r) which has beta numbers corresponding to run 1, beta numbers corresponding to run 2, etc
            if r == 1
                runs.bidx1 = runsc.x1;
                idxTr.sess(runs.bidx1) = {1};
            else
                runs.(sprintf('bidx%d', r)) = setdiff(runsc.(sprintf('x%d', r)), runsc.(sprintf('x%d', r-1)));
                idxTr.sess(runs.(sprintf('bidx%d', r))) = {r};
            end
        end
        idxTr.sess = cell2mat(idxTr.sess);%convert to matrix format
    end
end


%% Volume Parameters 
if S.existpatmat==1
    S.vol_info = 'NA';
    S.roi_name = 'NA';
    S.roi_file = 'NA';
    S.secondaryMask = [];
else
    S.vol_info = spm_vol(fullfile(par.funcdir, [par.boldrundirpfx par.refrun], par.ref_funcimage)); %get functional data resolution info for spm .img writing
    
    S.roi_file = [S.expt_dir S.subj_id '/' S.maskdir '/' S.roi_name par.imageextension]; %this is the large-scale ROI (could be wholebrain) that workspace info is calculated for. Saves time to have this volume include any sub-volumes you are interested in (e.g. MTL if you plan on looking in hippo and phc separately)
    
    %Apply another mask to the primary data loaded in the workspace. [] = no secondary mask.
    %This is useful for a number of potential scenarios. It's particularly
    %useful when working with betas, enabling us to filter a hippocampal ROI
    %further to remove zeros or NaNs from outside the brain space mask (due to implicit masking, dropout, etc).
    if strcmp(S.inputformat, 'raw')
        S.secondaryMask = [S.expt_dir S.subj_id '/' S.maskdir '/' S.roi_name par.imageextension]; % secondary mask (the specific classification mask - e.g. hippocampus within MTL)
    elseif strcmp(S.inputformat, 'betas')
        S.secondaryMask = [];
    end
end
%% Workspace Parameters - these files can be huge. In future versions, consider finding ways to pare down.
S.workspace = fullfile(S.workspace_dir, [S.subj_id '_' S.roi_name '_' S.smoothTxt{S.funcType} '_train_' S.trainTask '_test_' S.testTask S.preprocType '.mat']);

%% Pattern names
S.patternType = S.inputformat; %'raw' or 'betas'
if strcmp(S.inputformat, 'raw')
    S.preprocPatName = 'spiral_hp_z';%stands for 'spiral imaging'_'high-pass filtered'_'z-scored'. Doesn't matter really. Just leave this alone.
elseif strcmp(S.inputformat, 'betas')
    S.preprocPatName = 'betas';%'betas_z';%use betas_z if z-scoring betas
end

S.preprocPatCondensedName = [S.preprocPatName '_condensed'];

if isempty(S.secondaryMask)
    S.preprocPatNameFinalMask = S.preprocPatName;
else
    S.preprocPatNameFinalMask = [S.preprocPatName '_masked'];
end

%% UNSUPPORTED -- Artifacts
%S.artFile = (fullfile(par.artrepdir, ['art_global_modified_' par.substr])); %directory where artifact information is contained
S.inactivateArtifacts = 0; %remove artifact trials? 0 = no.

%% UNSUPPORTED -- outlier detection
S.remove_artdetect_outliers = 0; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
S.artdetect_motion_thresh = 0; % specify ArtDetect bin for motion outliers
S.artdetect_global_signal_thresh = 0; % specify ArtDetect bin for global signal outliers
S.remove_outlier_trials = 3;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)

%% UNSUPPORTED -- Subsample
%S.subsampleToMatch = 0; %subsample trials to match quantities across them.
%S.balanceHiAndLowConf = 0;% match N of hi and low confidence trials?

%% UNSUPPORTED -- voxel interactions
S.includeVoxelInteractions = 0; %include interactions among voxels? 0 = no
S.interactionType = 2;
S.intEst = 1;
S.intConcat = 0;
S.intPThresh = .001;
S.intReportIncrement = 100;
S.intFilePrefix = 'intVox';
S.intMaskName = 'interactionMask';
S.intPatName = 'interactions';
S.intGroupName = 'interactionsGroup';
S.intUseIntsWithUniqueInfo = 1;

%% UNSUPPORTED -- classifier parameters
S.nPlsCompsSet = 0; % number of pls components to include. 0 = do not use pls components.
S.class_args.libsvm = '-q -s 0 -t 2 -d 3'; % arguments for libsvm
S.class_args.radialBasisSelection = [];%[.00001 .0001 .001 .01 .1 1 10];

end


function raw_filenames = gen_raw_boldfnames(par, S)

runfolds = dir(fullfile(par.funcdir, par.runnames));%dir(fullfile(par.funcdir, 'localizer*'));%
if strcmp(par.runnames,'') % add contingency for when all the raw filenames are just dumped in your main funcdir (i.e., there are no runfolds)
    runfolds = runfolds(1); % runfolds(1) = '.', aka, funcdir itself
end
for idxr = 1:length(runfolds)
    if ~strcmp(par.runnames,'') % add contingency for when all the raw filenames are just dumped in your main funcdir (i.e., there are no runfolds)
        
        allrawfilenames{idxr,1} = dir(fullfile(par.funcdir, runfolds(idxr).name, ['/' par.boldnames '*' par.imageextension]));%'/swa*.nii'));%
    else
        allrawfilenames{idxr,1} = dir(fullfile(par.funcdir, ['/' par.boldnames '*' par.imageextension]));%'/swa*.nii'));%
    end
    %if 3D images (not recommended) check if the count matches that
    %specified for other stages of the process
    if par.ImgDims == 3
        if length(allrawfilenames{idxr})~=par.TRsperRun(idxr);
            error('your specified run length does not match 3D file count')
        end
    end
    
    for idxf = 1:length(allrawfilenames{idxr})
        if ~strcmp(par.runnames,'') % add contingency for when all the raw filenames are just dumped in your main funcdir (i.e., there are no runfolds)
            allrawfilepaths{idxr,1}{idxf,1} = runfolds(idxr).name;
        else
            allrawfilepaths{idxr,1}{idxf,1} = '';
        end
    end
end
allrawfilenames = vertcat(allrawfilenames{:});
allrawfilepaths = vertcat(allrawfilepaths{:});
for idx = 1:length(allrawfilenames);
    raw_filenames{idx,1} = [par.funcdir char(allrawfilepaths(idx)) '/' allrawfilenames(idx).name];
end

%files may have been read in out of order. This would be very very bad.
%Here, we try to confirm/fix this with a resort - but you *MUST* double
%check that the final file order is correct before proceeding with
%analysis
for idx = 1:length(raw_filenames)
    %first, identify the image number from its name in full
    %('001' from run_001.nii)
    nifti_indices = strfind(raw_filenames{idx,1}, par.imageextension);
    underscore_indices = strfind(raw_filenames{idx,1}, '_'); %assuming the number is preceded by '_', where are the underscores?
    imnum = str2double(raw_filenames{idx,1}(underscore_indices(end)+1:nifti_indices(end)-1));
    raw_filenames{idx,2} = imnum;
    %if length(raw_filenames{idx,1}) == 100%80
    %    raw_filenames{idx,2} = str2double(raw_filenames{idx,1}(length(raw_filenames{idx,1})-9:length(raw_filenames{idx,1})-9));
    %else
    %    raw_filenames{idx,2} = str2double(raw_filenames{idx,1}(length(raw_filenames{idx,1})-10:length(raw_filenames{idx,1})-9));
    %end
    
end

a = sortrows(raw_filenames, 2);
raw_filenames = a(:,1);

%if the BOLD images are 3D instead of 4D (TR-by-TR; NOT recommended, but currently only option supported [12/31/17]),
%we need to modify indices further to avoid introducing a new sorting error
if par.ImgDims == 3
    if ~strcmp(par.runnames,'') % add contingency for when all the raw filenames are just dumped in your main funcdir (i.e., there are no runfolds)
        for idx = 1:length(raw_filenames)
            %first, identify the RUN number from its name in full
            runref_indices = strfind(raw_filenames{idx,1}, ['/' par.boldrundirpfx]);
            runidxnum = str2double(raw_filenames{idx,1}(runref_indices(1)+5:runref_indices(2)-1));
            raw_filenames{idx,3} = runidxnum;
        end
        
        b = sortrows(raw_filenames, 3);
        raw_filenames = b(:,1);
    end
end

%save raw_filenames for reference
savename_rawfnms=[par.funcdir 'raw_filenames.mat'];
save(savename_rawfnms, 'raw_filenames');

end