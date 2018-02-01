function par = sub_params_mvpa(subject)
% function par = sub_params(subject)
% Basic parameter specification directory.  Creates par struct which can be
% used in par_ functions as the subpar argument
%
% Modified by TIB 2/18/14
%
% subject = '01'
%


%%
%dir_base
fprintf('\nEstablishing Parameters for Subject %s...', subject);

%spm_get_defaults; %set-up defaults (also makes sure spm is happily in the path)
%
% par.substr = subject;
% % convert sub string to number to simplify things later-- assumes s## format
% par.subnum = str2double(subject(end-1:end));


%% SPECIFY DIRECTORY NAMES
% ---------Directory names----------
% par.pardir = fullfile([dbase 'MCI'],'mciExpt');
% par.scriptsdir = fullfile(par.pardir,'scripts','spm' );
% addpath(par.scriptsdir);
% par.groupbehavdir = fullfile(par.pardir, 'grp_behav_data'); % pkFLAG- group behav data? for saving?
% par.group = 2;
% switch par.group
%     case 1
% Younger adults
%         par.fmridir = fullfile([dbase 'fMRI_MCI'],'young');
%     case 2
%         % Older adults (healthy controls)
%         par.fmridir = fullfile([dbase 'fMRI_MCI'],'old');
%     case 3
%         % MCI patients
%         par.fmridir = fullfile([dbase 'fMRI_MCI'],'mci');
% end
%---sub specific---
par.subdir = ['/media/sf_host/mvpa_sample_data/CM_localizer/CM' subject];
par.funcdir = fullfile(par.subdir,'/bolds/');
par.behavdir = fullfile(par.subdir, '/results01');
par.logdir = fullfile(par.subdir, 'logfiles');

% %% MODEL SELECTION

%
% %% SPECIFY BASIC SCAN PARAMETERS
% %----------Scan Params----------------
% par.numscans = 8; %number of scans (4 for encoding, 4 for retrieval)
% cd(par.rawdir);
% fprintf([num2str(par.numscans) ' scans for this subject\n']);
par.TR = 2; %TR in seconds
par.numslice = 36;
% par.TA = par.TR-par.TR/par.numslice; %TA
% par.sliceorder = 1:1:par.numslice; %slice order (assumes ascending)
par.refslice = floor(par.numslice-1); %reference slice (assumes middle for interleaved)
%
% % anat info
% par.inplaneimg = matchfiles(fullfile(par.anatdir, '*Inplane*.nii'));
% par.inplaneimg = par.inplaneimg{1};
% par.inplregimg = [par.inplaneimg(1:end-4) '_reg.nii']; % axial co-reg'd to meanfunc
% % par.inplregbrainimg = [par.inplregimg(1:end-4) '_brain.nii']; % axial co-reg'd to meanfunc, BET'd
% par.hiresimg = matchfiles(fullfile(par.anatdir, '*SPGR*.nii'));
% if isempty(par.hiresimg) == 1
%     par.hiresimg = 'no_SPGR_collected';
% else
%     par.hiresimg = par.hiresimg{1};
% end
% par.coronalimg = matchfiles(fullfile(par.anatdir, '*Coronal*.nii'));
% par.coronalimg = par.coronalimg{1};
% par.coronalimgor = [par.coronalimg(1:end-4) '_orient.nii'];
%
%
%% CREATE NAMES OF SCAN FILES
% % --------populate image list names------%switch cases depending on 8 run
% %or 16 run variant
% switch par.modelnum
%     case 1 % enc
par.scans_to_include = [1:16];
%     case 2 % encm
%         par.scans_to_include = [1:8];

% end
%
% par.minvol = 0; % fslsplit creates 3D files starting at 0
% k = 0;

%~~~~~3D realigned nifti files~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% for H = 1:length(par.scans_to_include)
%     fnum = 0;
%     I = par.scans_to_include(H);
%
%      runNum = num2str(I);
%
%       runpath = [par.funcdir '/run_' runNum];
%       cd(runpath);
%       tmp=dir('urun*.nii')
%       %tmp=dir('surun*.nii')
%       boldindx = (length(tmp)+5)
%       for J = 6:boldindx %these are the 3D nifti file numbers for this subject. We start at 6 because we dropped the first 6 TRs from each run
%         K = (J-5);
%         strJ = sprintf('%4.4d', J);
%         %par.scanfiles{H}(K,:) = [par.funcdir '/Study_run' runNum '/swaCM' subject '_study0' runNum '_' strJ '.nii']; %smoothed scan files
%         par.scanfiles{H}{K,1} = [par.funcdir '/run_' runNum '/urun' runNum '_' strJ '.nii']; %unsmoothed scan files
%         %par.scanfiles{H}{K,1} = [par.funcdir '/run_' runNum '/surun' runNum '_' strJ '.nii']; %smoothed scan files
%
%     end
% end
% %
%  par.allscanfiles = vertcat(par.scanfiles{:});

% runfolds = dir(fullfile(par.funcdir, 'run_*'));
% for idxr = 1:length(runfolds)
%     allrawfilenames{idxr,1} = dir(fullfile(par.funcdir, runfolds(idxr).name, '/urun*.nii'));
%     for idxf = 1:length(allrawfilenames{idxr})
%         allrawfilepaths{idxr,1}{idxf,1} = runfolds(idxr).name;
%     end
% end
% allrawfilenames = vertcat(allrawfilenames{:});
% allrawfilepaths = vertcat(allrawfilepaths{:});
% for idx = 1:length(allrawfilenames)%-1%note, we are filling in the beta file names based on how many betas OF INTEREST we have (length(betaidx)). We don't care about the error reg betas for this analysis
%     raw_filenames{idx,1} = [par.funcdir char(allrawfilepaths(idx)) '/' allrawfilenames(idx).name]; %create the analog to "raw_filenames.mat" - i.e. a list of all filenames including the path
% end
% for idx = 1:length(raw_filenames)
%     if length(raw_filenames{idx,1}) == 109%80
%         raw_filenames{idx,2} = str2double(raw_filenames{idx,1}(length(raw_filenames{idx,1})-9:length(raw_filenames{idx,1})-9));
%     else
%         raw_filenames{idx,2} = str2double(raw_filenames{idx,1}(length(raw_filenames{idx,1})-10:length(raw_filenames{idx,1})-9));
%     end
% end
% a = sortrows(raw_filenames, 2);
% raw_filenames = a(:,1);
%
% par.allscanfiles = raw_filenames;



% from rsa_CM_Localizer.m

%specify preprocessing level of BOLDs
preproc_lvl = ''; % 'a' for slice-time-only, 'u' for realigned-only, 'ua' for realign+unwarped, 'swua' for smoothed, normalized, and... you get the picture. Modify as needed if you changed SPM's prefix append defaults
boldnames = [preproc_lvl 'run']; %name of image files with preprocessing level prefix

ImgDims = 3; %if working with timeseries, it is highly recommended that you use 4D nifti files ('4'). If you have split them out into TR-by-TR, or are working with betas, enter '3'

TRsperRun = [114,114];

runfolds = dir(fullfile(par.funcdir, 'run_*'));%dir(fullfile(par.funcdir, 'localizer*'));%
for idxr = 1:length(runfolds)
    allrawfilenames{idxr,1} = dir(fullfile(par.funcdir, runfolds(idxr).name, ['/' boldnames '*.nii']));%'/swa*.nii'));%
    
    %if 3D images (not recommended) check if the count matches that
    %specified for other stages of the process
    if ImgDims == 3
        if length(allrawfilenames{idxr})~=TRsperRun(idxr);
            error('your specified run length does not match 3D file count')
        end
    end
    
    for idxf = 1:length(allrawfilenames{idxr})
        allrawfilepaths{idxr,1}{idxf,1} = runfolds(idxr).name;
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
    
end

a = sortrows(raw_filenames, 2);
raw_filenames = a(:,1);

%if the BOLD images are 3D instead of 4D, we need to modify indices further to avoid introducing a new sorting error

for idx = 1:length(raw_filenames)
    %first, identify the RUN number from its name in full
    runref_indices = strfind(raw_filenames{idx,1}, '/run_');
    runidxnum = str2double(raw_filenames{idx,1}(runref_indices(1)+5:runref_indices(1)+6));%runref_indices(2)-1)); %%Warning - this coding assumes the run numbers do not exceed double digits
    raw_filenames{idx,3} = runidxnum;
end

b = sortrows(raw_filenames, 3);
raw_filenames = b(:,1);

par.allscanfiles = raw_filenames;



%~~~~~~~read me some 4D slice-timed nifti files~~~~~~~~~~~~~~~~~~~~~~~~~~~
% for H = 1:length(par.scans_to_include)
%     fnum = 0;
%     I = par.scans_to_include(H);
%      runNum = num2str(I);
%
%       runpath = [par.funcdir '/Study_run' runNum];
%
%       cd(runpath);
%
%       flnm = ['^wastudy0' runNum];
%       x = spm_select('ExtFPList', pwd, flnm, 1:300);
%       for j = 1:length(x)
%           par.warunfiles{H}(j,:) = [x(j,:)];
%
%      end
% end
% par.allscanfiles = vertcat(par.warunfiles{:});


% par.numvols = par.maxvol(I) + 1; % numvols the same for all runs, so just take the number from the last run of either enc/ret
%
% % every 10th raw scan file... used for art_movie_nogui quick flag
% par.artfiles = par.allscanfiles(1:10:(length(par.allscanfiles)),:);
%
%
%% SLICE TIMING, REALIGN, and RESLICE PARAMS
% % % variables for slice timing not specified above
% %
% % par.slicetiming(1) = par.TA / (par.numslice -1);%timing var for slice timing
% % par.slicetiming(2) = par.TR - par.TA; %timing var for slice timing
% %
% %
% %
% realign flags...
par.realflag.quality = 1;
par.realflag.fwhm = 5;
par.realflag.sep = 4;
% par.realflag.rtm = 1; % if field exists, realigns to mean image; default realigns to 1st image
% par.realflag.PW;  %if field exists, then weighting is done...
par.realflag.interp = 2;

% reslice flags
par.reslflag.mask = 1;
par.reslflag.mean = 1;
par.reslflag.interp = 4;
par.reslflag.which = 2;
%
%
% %% COREGISTRATION, SEGMENTING, and NORMALIZATION PARAMS
% % coregistration info
% par.cor_meanfunc = [par.funcdir '/mean.nii'];
% % note that this mean scan file is assuming you realign after slice
% % timing...
%
% % segment info
% par.img2bSeg = par.inplregimg;
% % par.img2bSeg = par.inplregbrainimg;
% % need filepartinfo for below
% [segpath, segname, segext] = fileparts(par.img2bSeg);
% par.segopts = struct('biascor',1,'GM',[0 0 1],'WM',[0 0 1],'CSF',[0 0 0],'cleanup',0);
% par.segbiasfwhm = 60; % 60 is default in gui, 75 is default in command line for reasons unexplained
% % see spm_config_preproc and spm_preproc(_write) for details
%
%
% % normalization:
% % gray matter:
% par.graytemp = '/Applications/spm5/apriori/grey.nii';
% par.grsegs(1,:) = fullfile(segpath, ['c1' segname segext]);
% par.grsegs(2,:) = fullfile(segpath, ['c2' segname segext]);
% par.graysrcimg = fullfile(segpath, ['c1' segname segext]);
% par.graywrimg = fullfile(segpath, ['c1' segname segext]);
% par.grflags.smosrc = 8;
% par.grflags.smoref = 0;
% par.grflags.regtype = 'mni';
% par.grflags.cutoff = 25;
% par.grflags.nits = 16;
% par.grflags.reg = 1;
% par.grwrflags.preserve = 0;
% par.grwrflags.bb = [[-78 -112 -50];[78 76 85]]; % where did this come from?
% par.grwrflags.vox        = [3 3 3];
% par.grwrflags.interp     = 1;
% par.grwrflags.wrap       = [0 0 0];
%
%
% %% SMOOTHING and SPECMASK PARAMS
% % smoothing funcs
% par.smoothkernel = [8 8 8];
%
% % specmaskvars
% par.specwrflags.preserve = 0;
% par.specwrflags.bb = [[-78 -112 -50];[78 76 85]];
% par.specwrflags.vox        = [1 1 1];
% par.specwrflags.interp     = 1;
% par.specwrflags.wrap       = [0 0 0];
%
% par.specsmooth = [20 20 20];
% par.maskimg = [par.anatdir '/mask.nii'];

%explicit mask
par.expmaskimg = [par.subdir '/rbilat_hipp_trace.nii'];
% par.tsmaskimg = [par.anatdir '/tsmask.nii'];
% % par.wmaskimg = [par.anatdir '/wmask.img'];
% % par.swmaskimg = [par.anatdir '/swmask.img'];
% % par.tswmaskimg = [par.anatdir '/tswmask.img'];
% par.addsegs = 'i1 + i2';
% par.maskthresh = 'i1 > .2';

%% MODEL PARAMETERS

% model vars...
par.timing.fmri_t0 = par.refslice;  %micro-time resolution stuff (changed from 8)
par.timing.fmri_t = par.numslice; %used to be 16-- changed based on conversation with melina
par.timing.units = 'secs';
% par.bases.hrf.derivs = [0 0]; %defined below according to model
% Melina says no cost to doing time and  dispersion derivs ([0 0] = no derivs)
% NOTE that this will impact how you make contrasts!!
par.volt = 1;%1 - do not model interactions
% par.sess.scans is specified after populating list names...
%switch par.modelnum
%    case 1 % enc
%         par.sess.multi = {fullfile(par.behavdir, 'enc', ['mci_' par.substr '_enc.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'enc', ['mci_' par.substr '_enc.regs.mat'])};
%    case 2 % enc + motion regs
%         par.sess.multi = {fullfile(par.behavdir, 'enc', ['mci_' par.substr '_encm.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'enc', ['mci_' par.substr '_encm.regs.mat'])};
%     case 3 % ret
%         par.sess.multi = {fullfile(par.behavdir, 'ret', ['mci_' par.substr '_ret.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'ret', ['mci_' par.substr '_ret.regs.mat'])};
%     case 4 % ret + motion regs
%         par.sess.multi = {fullfile(par.behavdir, 'ret', ['mci_' par.substr '_retm.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'ret', ['mci_' par.substr '_retm.regs.mat'])};
%     case 5 % enc_nov
%         par.sess.multi = {fullfile(par.behavdir, 'enc_nov', ['mci_' par.substr '_encNov.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'enc_nov', ['mci_' par.substr '_encNov.regs.mat'])};
%     case 6 % enc_nov + motion regs
%         par.sess.multi = {fullfile(par.behavdir, 'enc_nov', ['mci_' par.substr '_encNovm.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'enc_nov', ['mci_' par.substr '_encNovm.regs.mat'])};
%     case 7 % ret_assoc
%         par.sess.multi = {fullfile(par.behavdir, 'ret_assoc', ['mci_' par.substr '_retAssoc.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'ret_assoc', ['mci_' par.substr '_retAssoc.regs.mat'])};
%     case 8 % ret_assoc + motion regs
%         par.sess.multi = {fullfile(par.behavdir, 'ret_assoc', ['mci_' par.substr '_retAssocm.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'ret_assoc', ['mci_' par.substr '_retAssocm.regs.mat'])};
%     case 9 % enc_epoch
%         par.sess.multi = {fullfile(par.behavdir, 'enc', ['mci_' par.substr '_enc_epoch.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'enc', ['mci_' par.substr '_enc_epoch.regs.mat'])};
%     case 10 % ret_epoch
%         par.sess.multi = {fullfile(par.behavdir, 'ret', ['mci_' par.substr '_ret_epoch.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'ret', ['mci_' par.substr '_ret_epoch.regs.mat'])};
%     case 11 % enc_nov_epoch
%         par.sess.multi = {fullfile(par.behavdir, 'enc_nov', ['mci_' par.substr '_encNov_epoch.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'enc_nov', ['mci_' par.substr '_encNov_epoch.regs.mat'])};
%     case 12 % ret_assoc_epoch
%         par.sess.multi = {fullfile(par.behavdir, 'ret_assoc', ['mci_' par.substr '_retAssoc_epoch.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'ret_assoc', ['mci_' par.substr '_retAssoc_epoch.regs.mat'])};
%     case 13 % enc_item
%         par.sess.multi = {fullfile(par.behavdir, 'enc_item', ['mci_' par.substr '_enc_item.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'enc_item', ['mci_' par.substr '_enc_item.regs.mat'])};
%     case 14 % ret_item
%         par.sess.multi = {fullfile(par.behavdir, 'ret_item', ['mci_' par.substr '_ret_item.onsets.mat'])};
%         par.sess.multi_reg = {fullfile(par.behavdir, 'ret_item', ['mci_' par.substr '_ret_item.regs.mat'])};
% end
par.sess.hpf = 128;%128;  % has anyone played around with this AND linear regs?
% par.sess.cond.help = 'Placeholder';
par.cvi = 'AR(1)'; %note that this actually gets changed to AR(0.2) in spm_fmri_spm_ui.
% It looks like you might be able to give it a custom value
% by simply putting a number in here, but I haven't played around with it
par.global = 'None';
% explicit mask - input {par.tsmaskimg} - which points to a mask file
% (defined above). For no explicit mask, enter []
par.mask = [];
%par.mask = [{par.expmaskimg}];

% contrast vars
% par.constat = 'T';

% % more model stuff
par.sess.scans = par.allscanfiles;

% Expected HRF length - aka FIR lag; FIR overrides HRF specification
%par.bases.fir = 1;
%par.hrf_len = 20;
% switch par.modelnum
%     case 1
% Model enc: 6 regs, conds = {'inin','inrec','innov','rec','rep','junk'},
%  no motion regs, no temporal deriv
par.bases.hrf.derivs = [0 0];
%     case 2
%         % Model encm: 6 regs, conds = {'inin','inrec','innov','rec','rep','junk'},
%         %  yes motion regs, no temporal deriv
%         par.bases.hrf.derivs = [0 0];
%     case 3
%         % Model ret: 7 regs, conds = {'inin','inrec','innov','recin','recrec','recnov','novin','novrec','novnov','junk'},
%         %  no motion regs, no temporal deriv
%         par.bases.hrf.derivs = [0 0];
%     case 4
%         % Model retm: 7 regs, conds = {'inin','inrec','innov','recin','recrec','recnov','novin','novrec','novnov','junk'},
%         %  yes motion regs, no temporal deriv
%         par.bases.hrf.derivs = [0 0];
%     case 5
%         % Model enc_nov: 2 regs, conds = {'nov','rep'},
%         %  no motion regs, no temporal deriv
%         par.bases.hrf.derivs = [0 0];
%     case 6
%         % Model enc_novm: 2 regs, conds = {'nov','rep'},
%         %  yes motion regs, no temporal deriv
%         par.bases.hrf.derivs = [0 0];
%     case 7
%         % Model ret_assoc: 5 regs, conds = {'ahit','amiss','miss','FA','CR'},
%         %  no motion regs, no temporal deriv
%         par.bases.hrf.derivs = [0 0];
%     case 8
%         % Model ret_assocm: 5 regs, conds = {'ahit','amiss','miss','FA','CR'},
%         %  yes motion regs, no temporal deriv
%         par.bases.hrf.derivs = [0 0];
%     case 9
%         % Model enc_epoch: 6 regs, conds = {'inin','inrec','innov','rec','rep','junk'},
%         %  no motion regs, no temporal deriv
%         par.bases.hrf.derivs = [0 0];
%     case 10
%         % Model ret_epoch: 7 regs, conds = {'inin','inrec','innov','recin','recrec','recnov','novin','novrec','novnov','junk'},
%         %  no motion regs, no temporal deriv
%         par.bases.hrf.derivs = [0 0];
%     case 11
%         % Model enc_nov_epoch: 2 regs, conds = {'nov','rep'},
%         %  no motion regs, no temporal deriv
%         par.bases.hrf.derivs = [0 0];
%     case 12
%         % Model ret_assoc_epoch: 5 regs, conds = {'ahit','amiss','miss','FA','CR'},
%         %  no motion regs, no temporal deriv
%         par.bases.hrf.derivs = [0 0];
%     case 13
%         % Model enc_item: 5 regs, conds = {'hit','miss','rec','rep','junk'},
%         %  no motion regs, no temporal deriv
%         par.bases.hrf.derivs = [0 0];
%     case 14
%         % Model ret_item: 5 regs, conds = {'hit','miss','FA','CR','junk'},
%         %  no motion regs, no temporal deriv
%         par.bases.hrf.derivs = [0 0];
% end

%% ---EMAIL OPTIONS---
% *PLEASE* if you are not Thackery, then change these settings  (or disable the
% mailing option) so that you don't inform him of every subject you run, thanks!
%
% par.smtpserv = 'smtp-roam.stanford.edu'; % smtp server
% par.fromaddy = 'thackery@stanford.edu'; % from address
% par.toaddy = 'thackery@stanford.edu'; % to address



fprintf('Done\n');
