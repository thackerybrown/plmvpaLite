function [S idxTr idxTe par]= TIB_mvpa_params_8080_betas(subj_id, task, TRsperRun, imgtype)
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

%% establish general parameters~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
idxTr = [];
idxTe = [];

%Study name
S.exp_name = 'CM_localizer'; %change this to flexibly redirect the script to different studies in subdirectories

%Subject ID/number
par.substr = ['CM' subj_id{1}];
S.subj_id = par.substr;

%Task type
par.task = task; %assign input from function call. Task phase label. For circmaze, this is 'goals' or 'plan'. For localizer (8080), this is 'CM_Localizer'

par.TR = 2; %TR (s). NOTE: if working with Betas - if your "onsets" file is really a numerical list of beta numbers instead of onsets, you must set TR = 1. If your onsets file is in seconds, set this as the TR even for Betas.

ImgDims = 3; %as of 12/31/17, code only supports 3D images. %it is highly recommended that you modify to use 4D nifti files for raw BOLD data ('4'). If you have split them out into TR-by-TR, enter '3'

par.readimglist = 0; %1=yes; 0 = no. Flag specifies whether to generate raw_filenames on the fly or to read in a previously-made file

%Functional image scan selectors
%par.scansSelect.goals.loc = 1:1;%***if ALL FILENAMES corresponding to ALL RUNS OF INTEREST are stored in ONE cell of raw_filenames.mat (i.e., not broken up by run), set index to 1 or 1:1. Otherwise, create indexing for elements of cell array raw_filenames.mat corresponding to task of interest (i.e. if cells runs 1:4 correspond to task 1, we want to reference {1}, {2}... in raw_filenames.mat)
%par.scansSelect.plan.loc = 1:1;%***if ALL FILENAMES corresponding to ALL RUNS OF INTEREST are stored in ONE cell of raw_filenames.mat (i.e., not broken up by run), set index to 1 or 1:1. Otherwise, create indexing for elements of cell array raw_filenames.mat corresponding to task of interest (i.e. if cells runs 5:8 correspond to task 2, we want to reference {5}, {6}... in raw_filenames.mat)
par.scansSelect.CM_localizer.loc = 1:1;%***if ALL FILENAMES corresponding to ALL RUNS OF INTEREST are stored in ONE cell of raw_filenames.mat (i.e., not broken up by run), set index to 1 or 1:1. Otherwise, create indexing for elements of cell array raw_filenames.mat corresponding to task of interest (i.e. if cells runs 1:4 correspond to phase 1, we want to reference {1}, {2}... in raw_filenames.mat)
par.scansSelect.CM_localizer.loc = 1:1;%***if ALL FILENAMES corresponding to ALL RUNS OF INTEREST are stored in ONE cell of raw_filenames.mat (i.e., not broken up by run), set index to 1 or 1:1. Otherwise, create indexing for elements of cell array raw_filenames.mat corresponding to task of interest (i.e. if cells runs 5:8 correspond to phase 2, we want to reference {5}, {6}... in raw_filenames.mat)

%input image info
S.inputformat = imgtype; %assign input from function call. Either 'raw' for raw bold images or 'betas' for beta images. Selection here automatically changes some params below.
if strcmp(S.inputformat, 'raw')
    data_imgs_to_use = 'raw_filenames.mat';
elseif strcmp(S.inputformat, 'betas')
    data_imgs_to_use = 'beta_filenames.mat'; % analogous to raw_filenames, but with a boolean index for which betas correspond to which conditions. Created by calling TIB_generate_beta_filenames.mat below
end

S.preprocType = 'spm'; % 'spm' for spm preprocessing, 'knk' for kendrick preprocessing

%% existing pattern handling inputs
S.existpatmat = 0; %1=yes - skip trying to load image files using SPM. We've already got all patterns in a matrix. Currently (09/2018) the existing pattern matrix is hardcoded in load_matrix_pattern_2D to be a matrix (double) named 'testmat'. \\
%Also, because so much of the Princeton MVPA toolbox assumes a mask volume is used, a dummy mask file is now included in the PLMVPA_Lite toolkit and called in TIB_run_MVPA_general
S.datafile = 'testsubmat_fixed.mat'; % added for ADNI neuropsych study. Replaces MRI image (e.g., nii) files with an existing matrix of pattern data

%% tasks or study phases
%set trainTask and testTask to be the same if you want to train and test on the same set of trials via
%cross-validation (see section below). If you want to train on one set of
%data (e.g., a localizer) and test on another (e.g., a retrieval task),
%then specify different tasks or study phases
S.trainTask = 'EAvsScene';%Circmaze - 'goals' or 'plan'
S.testTask = 'EAvsScene';%Circamze - 'goals' or 'plan'

%x-validation info
S.xvaltype = 'loo'; %set to 'loo' for leave-one-out x-validation or 'nf' for nfold using the S.nFolds defined below.

%how many single trial betas are the? Number is not used (can be ignored)
%for studies using raw bolds
S.stbetacount = 12; % NOTE: the code assumes all single trial betas of potential interest are contiguous in the model. If multi-event betas are inteleaved in the .mat model structure, this will need more editing.

%%model information - define which timepoints or images correspond to which classes of data
if strcmp(S.inputformat, 'raw')
    S.onsets_filename = [S.subj_id '_localizer_onsets_test'];%
    S.onsets_filename_tr = [S.subj_id '_localizer_onsets_test'];% added for train on 1 phase, test on another - this assumes the data are actually in the same set of files.
    S.onsets_filename_tst = [S.subj_id '_localizer_onsets_test'];% added for train on 1 phase, test on another - this assumes the data are actually in the same set of files.
elseif strcmp(S.inputformat, 'betas')
    S.onsets_filename = [S.subj_id '_localizer_onsets_oneperblcknew'];
    S.onsets_filename_tr = [S.subj_id '_localizer_onsets_oneperblcknew'];
    S.onsets_filename_tst = [S.subj_id '_localizer_onsets_oneperblcknew'];
    
    S.betaidx_filename = [S.subj_id '_betas_idx'];
    S.betaidx_filename_tr = [S.subj_id '_betas_idx_tr'];
    S.betaidx_filename_te = [S.subj_id '_betas_idx_te'];
end

%% directories~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
S.expt_dir = ['/home/brain/host/mvpa_sample_data/' S.exp_name '/'];%study location

par.subdir =[S.expt_dir S.subj_id];%subject location

par.funcdir =[par.subdir '/bolds/'];%subfolder for 'raw' BOLD data. Assumes BOLDs are stored in subfolders labeled 'run_01', etc)

S.workspace_dir = [par.subdir '/mvpa_workspace'];%temporary files workspace


%model file directory (onsets.mat and betas in here) - this is still used
%when working with raw data. We must have some way to tell the classifier
%which images correspond to which classes
if strcmp(S.inputformat, 'raw')
    S.mvpa_dir = [S.expt_dir S.subj_id '/results01/'];
elseif strcmp(S.inputformat, 'betas')
    S.mvpa_dir = [S.expt_dir S.subj_id '/results01/simpleLSA/'];
end

%ROI masks (could be whole-brain mask, but the code wants a mask file
S.anat_dir = [S.expt_dir S.subj_id '/Masks'];

S.importance_maps_dir=[S.expt_dir 'ImpMaps_' date];
S.group_mvpa_dir = [S.expt_dir 'mvpa_output_files'];%results .mat files are spit out in here

%% identify betas
if strcmp(S.inputformat, 'betas')%if we are running a beta analysis, take a moment to create "model" files for the beta images
    TIB_generate_beta_filenames(S);
end

%% create raw_filenames on the fly (as opposed to making it by loading structures from an SPM.mat file [for this method, see comment above])

%if fnames are wrong and you want to standardize, explore 'rename'
%function, ex.: "rename -v s/CM001_localizer02_0/run_02_/ *.nii"

%specify preprocessing level of BOLDs
preproc_lvl = ''; % 'a' for slice-time-only, 'u' for realigned-only, 'ua' for realign+unwarped, 'swua' for smoothed, normalized, and... you get the picture. Modify as needed if you changed SPM's prefix append defaults
boldnames = [preproc_lvl 'run']; %name of image files with preprocessing level prefix
runnames = 'run_*'; %what are your run/scan session folders called? If none or only 1 run, set to '' or whatever may be appropriate.

if strcmp(S.inputformat, 'raw')
    if par.readimglist == 1
        load([par.funcdir '/' data_imgs_to_use]); %loads predefined cell array
    else %generate them on the fly
        runfolds = dir(fullfile(par.funcdir, runnames));%dir(fullfile(par.funcdir, 'localizer*'));%
        if strcmp(runnames,'') % add contingency for when all the raw filenames are just dumped in your main funcdir (i.e., there are no runfolds)
            runfolds = runfolds(1); % runfolds(1) = '.', aka, funcdir itself
        end
        for idxr = 1:length(runfolds)
            if ~strcmp(runnames,'') % add contingency for when all the raw filenames are just dumped in your main funcdir (i.e., there are no runfolds)
                
                allrawfilenames{idxr,1} = dir(fullfile(par.funcdir, runfolds(idxr).name, ['/' boldnames '*.nii']));%'/swa*.nii'));%
            else
                allrawfilenames{idxr,1} = dir(fullfile(par.funcdir, ['/' boldnames '*.nii']));%'/swa*.nii'));%
            end
            %if 3D images (not recommended) check if the count matches that
            %specified for other stages of the process
            if ImgDims == 3
                if length(allrawfilenames{idxr})~=TRsperRun(idxr);
                    error('your specified run length does not match 3D file count')
                end
            end
            
            for idxf = 1:length(allrawfilenames{idxr})
                if ~strcmp(runnames,'') % add contingency for when all the raw filenames are just dumped in your main funcdir (i.e., there are no runfolds)
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
            nifti_indices = strfind(raw_filenames{idx,1}, '.nii'); %assuming .nii, where does that fall in the string?
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
        if ImgDims == 3
            if ~strcmp(runnames,'') % add contingency for when all the raw filenames are just dumped in your main funcdir (i.e., there are no runfolds)
                for idx = 1:length(raw_filenames)
                    %first, identify the RUN number from its name in full
                    runref_indices = strfind(raw_filenames{idx,1}, '/run_');
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

%% information specific to the training and testing of specific subsets of data

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

if strcmp(S.trainTask,'EAvsScene')
    S.onsetsTrainDir = [S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTrain = {{'EA'}  {'Scene'}} ;%corresponds to the names in the onsets.mat or betas_idx.mat files. This is used to select what is being compared with what.
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

%% testing - this defines the testing set. The code is set up this way to enable us to step outside xval if desired to test on different set of data (e.g., at retrieval)
if strcmp(S.testTask,'EAvsScene')
    S.onsetsTestDir =[S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTest = {{'EA'} {'Scene'}};
    S.nwayclass = num2str(numel(S.condsTest));%stores the number classification dimensions just for reference (i.e. is this a 5-way or a 2-way/binary classification?)
    S.TestRuns = par.scansSelect.(par.task).loc;
    if strcmp(S.inputformat, 'raw')
        S.filenames_test = raw_filenames;%
    elseif strcmp(S.inputformat, 'betas')
        S.filenames_test = beta_filenames;%
    end
    S.durTest = numel(S.filenames_test) * par.TR;
    %[~, idxTe] = fMRIBehAnalysis_Loc(par);
    
elseif strcmp(S.testTask,'EAvsObj')
    S.onsetsTestDir =[S.mvpa_dir];%directory containing onsets.mat or betas_idx.mat file to be loaded in
    S.condsTest = {{'EA'} {'Obj'}};
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

S.runs_vector = TRsperRun; % number of TRs per scanning run (coded for vector of values)

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
    idxTr.sess = onsets(1:S.stbetacount);%length(onsets)-1);%NOTE: Specific to circmaze, we are skipping the last onsets entry because this is an array of "arrows period" onsets.
    
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


%% Volume Parameters
S.vol_info = spm_vol(fullfile(par.funcdir, 'run_01', 'run_01_006.nii')); %get functional data resolution info for spm .img writing

S.roi_name = 'hvis0p1intensthresh.nii';
%S.roi_name = 'HVisCtx_1.nii';
%S.roi_name = 'NativeGM_BOLDres.nii';
S.roi_file = [S.expt_dir S.subj_id '/Masks/' S.roi_name]; %this is the large-scale ROI (could be wholebrain) that workspace info is calculated for. Saves time to have this volume include any sub-volumes you are interested in (e.g. MTL if you plan on looking in hippo and phc separately)

%Apply another mask to the primary data loaded in the workspace. [] = no secondary mask.
%This is useful for a number of potential scenarios. It's particularly
%useful when working with betas, enabling us to filter a hippocampal ROI
%further to remove zeros or NaNs from outside the brain space mask (due to implicit masking, dropout, etc).
if strcmp(S.inputformat, 'raw')
    S.secondaryMask = [S.expt_dir S.subj_id '/Masks/' S.roi_name]; % secondary mask (the specific classification mask - e.g. hippocampus within MTL)
elseif strcmp(S.inputformat, 'betas')
    S.secondaryMask = [S.mvpa_dir 'mask.nii'];
end

%% Workspace Parameters - these files can be huge. In future versions, consider finding ways to pare down.
S.use_premade_workspace = 0;
S.workspace = fullfile(S.workspace_dir, [S.subj_id '_' S.roi_name '_' S.smoothTxt{S.funcType} '_train_' S.trainTask '_test_' S.testTask S.preprocType '.mat']);

%% Pattern names
S.patternType = S.inputformat; %'raw' or 'betas'
if strcmp(S.inputformat, 'raw')
    S.preprocPatName = 'spiral_hp_z';%stands for 'spiral imaging'_'high-pass filtered'_'z-scored'
elseif strcmp(S.inputformat, 'betas')
    S.preprocPatName = 'betas';%'betas_z';%use betas_z if z-scoring betas
end

S.preprocPatCondensedName = [S.preprocPatName '_condensed'];

if isempty(S.secondaryMask)
    S.preprocPatNameFinalMask = S.preprocPatName;
else
    S.preprocPatNameFinalMask = [S.preprocPatName '_masked'];
end

%% Artifacts
%S.artFile = (fullfile(par.artrepdir, ['art_global_modified_' par.substr])); %directory where artifact information is contained
S.inactivateArtifacts = 0; %remove artifact trials? 0 = no.

%% Iteration Parameters
S.num_results_iter = 1; % number of times to run the entire classification process (select subset of the data and train/test classifier)
S.num_iter_with_same_data = 1; % number of times to run the classfication step for a given subset of data - useful for non-deterministic cases

%% Balancing Parameters
S.equate_number_of_trials_in_groups = 0; % equate number of trials in conditions
S.numBalancedParams = 1; % number of parameters to balance across (e.g., both goal location AND cue in Circmaze data). The code currently (12/29/17) only handles two options - 1 (standard; main class type), or 2 (main class type plus a second parameter, specified in a second file).
S.numBalancedIts = 10; % number of iterations to run, with different randomization for the balancing

%% Z-Scoring and outlier detection
S.perform_second_round_of_zscoring = 0;  % z-score data again immediately prior to classification
S.remove_artdetect_outliers = 0; % 1 = remove trials that exhibited movement or global signal artifacts as determined by ArtDetect
S.artdetect_motion_thresh = 0; % specify ArtDetect bin for motion outliers
S.artdetect_global_signal_thresh = 0; % specify ArtDetect bin for global signal outliers
S.remove_outlier_trials = 3;  % on-the-fly outlier detection/removal; specify how many std dev from whole brain mean to exclude as outliers (0 = don't exclude any trials)

%% Importance Maps
S.generate_importance_maps = 0; %visualize classifier weights
S.generateBetaMaps = 1; %use betas, instead of importance values
S.impType = {'pos' 'neg' 'both' 'raw'}; %importance map types
S.regNames = {'CondA' 'CondB'}; % should match number of classes

%% Special types of analysis
S.searchlightAnalysis = 0; % run a searchlight analysis
%S.linReg = 0; % run an analysis with a continuous outcome variable
S.scrambleregs = 0; % run an anlysis with the class labels scrambled on a run-by-run basis.

%% Subsample %%alan stuff only - don't use
%S.subsampleToMatch = 0; %subsample trials to match quantities across them.
%S.balanceHiAndLowConf = 0;% match N of hi and low confidence trials?

%% voxel interactions
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

%% Signal intensity analysis
S.thisSigIntenseSelector = 'randomNFold_xval'; %which selector to use for signal intensity analysis
S.zscoreIntensityVals = 1; % zscore the intensity values?

%% Noising/Denoising
S.addnoise = 0; %Overly sparse data (lots of features that often have zeroes)? Add gaussian white noise to pattern matrix to help with overfitting. >0 = yes. Value specified = target SNR.
S.denoise = 0; %undergo denoising?
S.denoiseOpt.denoisespec = '10001'; %which parts of the glm output do we want to save?

%% Mean Signal Extraction Params
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

%% TR Weighting
%which post-stimulus TRs should be used (and if more than one, averaged
%across) before feeding data to the classifier?  S.TR_weights_set{1} -
%training weights, S.TR_weights_set{2} = testing weights

%NOTE: for beta analyses, we don't average over multiple images because
%different images = different events
if strcmp(S.inputformat, 'raw')
    %S.TR_weights_set = {{[.0072 .2168 .3781 .2742 .1237] [.0072 .2168 .3781 .2742 .1237]}}; %approximates the canonical haemodynamic response
    S.TR_weights_set = {{[0 0.25 0.5 0.25] [0 0.25 0.5 0.25]}};%use double-bracket structure in case want to set code up to run a sliding window across multiple TR bins
elseif strcmp(S.inputformat, 'betas')
    S.TR_weights_set = {{[1] [1]}};%give full weighting to the 1 and only image corresponding to each event
end


%% classifier parameters
S.class_args.train_funct_name = 'train_liblinear_multiclass';%'train_pLR';   %training function
S.class_args.test_funct_name = 'test_liblinear_multiclass';%'test_pLR';      %testing function
S.class_args.classType = 'libLin';
S.perfmet_functs = 'perfmet_maxclass'; % performance metric
S.statmap_funct = 'statmap_anova';%'AG_statmap_anova'; % performance metric
S.nPlsCompsSet = 0; % number of pls components to include. 0 = do not use pls components.
S.nFolds = 8; % number of cross validation iterations - only used for nFold (as opposed to run-by-run leave-one-out)

S.class_args.nVox = 0; % number of voxels to select with feature selection e.g. [1000 5000 10000]
S.class_args.fseltype = 'topn'; % feature selection format: top N vox (topn) or random N vox (rand)?
S.class_args.libLin = '-q -s 0 -B 1'; %arguments for liblinear; -s 0 = L2; -s 6 = L1; -s 5 = L1 with L2 loss; -s 3 L2 with L1 loss
S.class_args.libsvm = '-q -s 0 -t 2 -d 3'; % arguments for libsvm
S.class_args.constant = true; % include a constant term?
S.class_args.prefitWeights = true;
S.class_args.chooseOptimalPenalty = 0; % 1 = yes. cycle through cost parameters in the training set, and choose the optimal one. Note, this only makes sense in context of loo with >2 runs or for nf with >2 folds, because it subdivides training set into additional 'runs' and performs nested xvalidation.
S.class_args.penaltyRange = [.001 .005 .01 .05 .1 .5 1 5 10 50 100 500 1000 50000]; % a vector "[]" of cost parameters to cycle through
S.class_args.radialBasisSelection = [];%[.00001 .0001 .001 .01 .1 1 10];
S.class_args.nFoldsPenaltySelection = 10; % number of cross validation folds for penalty parameter selection.
S.class_args.penalty = 1000; %uncomment if not using optimal penalty. Typical value is 1. If using sample data provided with plmvpaLite, start with 0.000001 to see how minimal regularization harms performance.
%establishment
end