function [] = rsa_CM_Localizer(Sub, Mask, TRsperRun)
% code for RSA analysis.

%For alternate code demo purposes, this is built around both 4D and 3D
%image file types. The use case for 4D includes scenarios like "raw" BOLD data or residual time-series

%warning: as of 1/3/2018; unconcatenated 4D analysis has not been debugged.
%at this time only use concatenated 3D

%% flags and parameters
S.TR = 2;
theseTRWeights = [0 0.25 0.5 0.25 0];
weights_str = mat2str(theseTRWeights);% assign values to string for custom output file naming

funcftype = '.nii';
runs_concat = 1; %1 = typical SPM analysis; will have continuous onsets concatenated across runs. 0 = But you might not have bothered creating such a file, or in SST case we are using files from FSL. In this case, the onsets are assumed to "reset" for each run ('raw' unconcatenated onsets)
use_exist_workspace = 0; %1=yes. Load existing pattern workspaces and onsets files. Saves time, but turn off if want to manually re-do pattern generation.

gen_onsetsTR = 1; %1=yes. Typically, you'll use an onsets.mat file with tr-by-tr onsets and names (as used for a beta-series). But if you only have a traditional GLM model with one name for multiple onsets, setting this flag to 1 will auto-populate unique but related names (e.g., Face_1; Face_2...)
runzscore = 1;%1=yes. Standard but controversial preprocessing.

%Subject ID/number
par.substr = ['CM' Sub{1}];
S.subj_id = par.substr;

mask = Mask;
S.exp_name = 'CM_Localizer';
study_prefix = 'CM';

S.inputformat = 'raw'; % are we working with BOLDs/timeseries ('raw') or with beta maps ('betas')?

S.onsets_filename = [S.subj_id '_localizer_onsets_test'];

%specify preprocessing level of BOLDs
preproc_lvl = 'a'; % 'a' for slice-time-only, 'u' for realigned-only, 'ua' for realign+unwarped, 'swua' for smoothed, normalized, and... you get the picture. Modify as needed if you changed SPM's prefix append defaults
boldnames = [preproc_lvl 'run']; %name of image files with preprocessing level prefix

ImgDims = 3; %if working with timeseries, it is highly recommended that you use 4D nifti files ('4'). If you have split them out into TR-by-TR, or are working with betas, enter '3'

%% Directories
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
    S.mvpa_dir = [S.expt_dir S.subj_id '/results01/betaseries_rearranged/'];
end

%ROI masks (could be whole-brain mask, but the code wants a mask file
S.anat_dir = [S.expt_dir S.subj_id '/Masks'];
maskname=[S.anat_dir '/' mask '.nii'];

S.group_mvpa_dir = [S.expt_dir 'RSA_output_files'];%results .mat files are spit out in here

if ~exist(S.group_mvpa_dir)
    mkdir(S.group_mvpa_dir)
end

if ~exist([S.mvpa_dir '/RSA_data/'])
    mkdir([S.mvpa_dir '/RSA_data/'])
end


%% extract patterns
if runs_concat == 1
    %load onsets
    load([S.mvpa_dir S.onsets_filename]);
    
    %runs = [1 2]; %indicate number of runs and which ones to target with analysis
    
    rmat_condensed = [];
    onsets_TRs = [];
    names_TRs = [];%filled in if gen_onsetsTR == 1
    
    %before loading and preprocessing all the pattern data, we can see if its
    %already there for analysis
    if use_exist_workspace && exist([S.mvpa_dir '/RSA_data/' S.subj_id '_' mask '_condensedpats.mat']);
        load([S.mvpa_dir '/RSA_data/' S.subj_id '_' mask '_condensedpats.mat']);
        load([S.mvpa_dir '/RSA_data/' S.subj_id '_onsets_expanded.mat']);
    else
        
        %% load in pattern data
        if ImgDims == 4 %in development
            
            a = [];
            
        elseif ImgDims == 3   
            
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
                %if length(raw_filenames{idx,1}) == 100%80
                %    raw_filenames{idx,2} = str2double(raw_filenames{idx,1}(length(raw_filenames{idx,1})-9:length(raw_filenames{idx,1})-9));
                %else
                %    raw_filenames{idx,2} = str2double(raw_filenames{idx,1}(length(raw_filenames{idx,1})-10:length(raw_filenames{idx,1})-9));
                %end
                
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
            run_sel = b(:,3);%store run numbers for reference
            imgslength = length(raw_filenames);
            
            %% iterate through 3D frames to extract all patterns
            for i=1:imgslength
                betamaps{i} = raw_filenames{i};%[tmp.name ',' num2str(i)];
                
                [b,r] = MAP_getROI(maskname, betamaps{i}, 'vox', 0, '');
                bmat_t(:,i) = b{1}; % returns all voxels, whether or not they have NaNs
                rmat_t(:,i) = r; % returns voxels excluding NaNs
                %meanbetas_t = nanmean(bmat_t(:,:));%get mean beta values from the ROI for each regressor
                
            end
            
            %% optional preprocessing
            if runzscore == 1
                run_sel = cell2mat(run_sel);
                
                pat_t = [];
                
                for r = 1:length(TRsperRun) %for each run
                    activepats = rmat_t(:,logical(run_sel==r));%filter patterns to current run
                    if size(activepats,2)~=TRsperRun(r)
                        error('Your pattern count doesnt match current run length');
                    end
                    
                    pat_t = [pat_t zscore_mvpa(activepats,2)];%2 = z-score within rows (within-voxels)
                end
                
                rmat_t = pat_t;
            end
            %% now compute a mean pattern for __ TRs surrounding the onset
            
            for n = 1:length(names)
                
                
                time_idx = floor(onsets{n}/S.TR) + 1;%convert onsets to TRs
                
                onsets_TR{n} = time_idx; %sort(horzcat(onsets_t{n}, onsets_t{n}{idxThisCond}(enoughTRs_h)));%put the onsets for cond{i} into an array,
                
                theseTRWeights2 = theseTRWeights;
                
                %create weighted mean pattern for that onset
                for nvoxr = 1:length(rmat_t)
                    %tempvals = [];
                    for tidx = 1:length(time_idx)
                        tempvals = theseTRWeights2.*rmat_t(nvoxr,time_idx(tidx):time_idx(tidx)+(length(theseTRWeights2)-1));
                        condmat_condensed_t(nvoxr,tidx) = squeeze(sum(tempvals));
                    end
                end
                rmat_condensed = [rmat_condensed condmat_condensed_t];
                
                if gen_onsetsTR == 1
                    for tidx = 1:length(time_idx)
                        tempnames{tidx} = [names{n} '_' num2str(tidx)];
                    end
                end
                
                names_TRs = [names_TRs tempnames];
            end
        end
        
        
        %% write out pattern info
        savename1 = [S.mvpa_dir '/RSA_data/' S.subj_id '_' mask '_condensedpats.mat'];
        save(savename1,'rmat_condensed');
        
        savename2 = [S.mvpa_dir '/RSA_data/' S.subj_id '_onsets_expanded.mat'];
        save(savename2,'names','onsets','names_TRs','durations');
    end
    
else %in debugging stage as of 1/3/2018
    %load onsets
    load([S.mvpa_dir S.onsets_filename]);
    
    %runs = [1 2]; %indicate number of runs and which ones to target with analysis
    
    rmat_condensed = [];
    onsets_TRs = [];
    
    %before loading and preprocessing all the pattern data, we can see if its
    %already there for analysis
    if exist([S.mvpa_dir '/RSA_data/' S.subj_id '_' mask '_condensedpats.mat']);
        load([S.mvpa_dir '/RSA_data/' S.subj_id '_' mask '_condensedpats.mat']);
    else
        
        %%load in pattern data
        if ImgDims == 4
            
            for rnum = 1:length(TRsperRun)
                
                run = num2str(rnum);
                
                path = [par.funcdir '/run_' run '/'];
                
                cd(path);
                
                
                tmp=dir(['*' funcftype]);
                
                runlength = length(spm_vol(strtok(tmp.name)));
                
                %iterate through 3D frames of 4D file to extract all patterns
                for i=1:runlength
                    betamaps{i} = [tmp.name ',' num2str(i)];
                    
                    
                    [b,r] = MAP_getROI(maskname, betamaps{i}, 'vox', 0, '');
                    bmat_t(:,i) = b{1}; % returns all voxels, whether or not they have NaNs
                    rmat_t(:,i) = r; % returns voxels excluding NaNs
                    %meanbetas_t = nanmean(bmat_t(:,:));%get mean beta values from the ROI for each regressor
                    
                    
                    
                    
                end
                %         names_t = names(cell2mat(runs_onsidx)==num2str(run));%filter to current run
                %         onsets_t = onsets(cell2mat(runs_onsidx)==num2str(run));%filter to current run
                
                names_t = names(cell2mat(runs_onsidx)==rnum);%filter to current run
                onsets_t = onsets(cell2mat(runs_onsidx)==rnum);%filter to current run
                %['ASSIGNED_habit_env4_rep' num2str(rnum)]        
                
                %%      now compute a mean pattern for __ TRs surrounding the onset
                
                % convert onsets to TRs
                for n = 1:length(names_t)
                                 
                    time_idx = floor(onsets_t{n}/S.TR) + 1;
                    enoughTRs_h = (time_idx + (length(theseTRWeights)-1)) <= runlength;
                    lengthdiff = runlength - time_idx; %calculate how many more TRs there are beyond the onset
                    %onsets{i} = sort(horzcat(onsets{i}, theseOnsets(enoughTRs_h)));
                    onsets_TR{n} = time_idx; %sort(horzcat(onsets_t{n}, onsets_t{n}{idxThisCond}(enoughTRs_h)));%put the onsets for cond{i} into an array,
                    
                    %if the n TRs is too short, calculate a new weights vector
                    %to accomodate
                    if enoughTRs_h == 1
                        theseTRWeights2 = theseTRWeights;
                    else
                        theseTRWeights2 = theseTRWeights(1:(lengthdiff+1))/(sum(theseTRWeights(1:(lengthdiff+1))));
                    end
                    %restricting to onsets for which we have enough TRs acquired.
                    %NOTE: This will literally REcreate the onsets component of my model files,
                    %BUT is still useful IF we haven't prescreened them based
                    %on how many TRs we acquired (i.e. we might get 31/32
                    %onsets of the 32nd exceeds the number of TRs we acquired,
                    %based on how many post-onset TRs we want to average across
                    %in our MVPA analysis
                    
                    %create weighted mean pattern for that onset
                    for nvoxr = 1:length(rmat_t)         
                        tempvals = theseTRWeights2.*rmat_t(nvoxr,time_idx:time_idx+(length(theseTRWeights2)-1));
                        rmat_condensed_t(nvoxr,n) = squeeze(sum(tempvals));
                    end
                end
                %% now integrate patterns and indices into master matrices for study-level analysis
                rmat_condensed = [rmat_condensed rmat_condensed_t];
                onsets_TRs = [onsets_TRs onsets_TR];
            end
            
            
        elseif ImgDims == 3
            
            a=[];
            
        end
        
        %% write out pattern info
        savename1 = [basepath '/wagner/thackery/SST_temp/sst' sub '/mvpa_orient/sst' sub '_' mask '_condensedpats.mat'];
        save(savename1,'rmat_condensed');
        
        savename2 = [basepath '/wagner/thackery/SST_temp/sst' sub '/mvpa_orient/onsets_nonconcat.mat'];
        save(savename2,'names','onsets','onsets_TRs','durations','runs_onsidx');
    end
end

%% create indices for patterns of interest

if gen_onsetsTR == 1
    names = names_TRs;
end

% first index different categories
for n = 1:length(names)
    if strfind(names{n},'EA')%
        EA_idx(n)=1;
        AA_idx(n)=0;
        Face_idx(n) = 0;
        Scene_idx(n)=0;
        Obj_idx(n)=0;
        othercond_idx(n) = 0;
    elseif strfind(names{n},'AA')%
        EA_idx(n)=0;
        AA_idx(n)=1;
        Face_idx(n) = 0;
        Scene_idx(n)=0;
        Obj_idx(n)=0;
        othercond_idx(n) = 0;
    elseif strfind(names{n},'Face')%
        EA_idx(n)=0;
        AA_idx(n)=0;
        Face_idx(n) = 1;
        Scene_idx(n)=0;
        Obj_idx(n)=0;
        othercond_idx(n) = 0;
    elseif strfind(names{n},'Scene')%
        EA_idx(n)=0;
        AA_idx(n)=0;
        Face_idx(n) = 0;
        Scene_idx(n)=1;
        Obj_idx(n)=0;
        othercond_idx(n) = 0;
    elseif strfind(names{n},'Obj')%
        EA_idx(n)=0;
        AA_idx(n)=0;
        Face_idx(n) = 0;
        Scene_idx(n)=0;
        Obj_idx(n)=1;
        othercond_idx(n) = 0;
    else
        othercond_idx(n) = 1;
        %fprintf('Error! None of categories names found for this idx!\n')
        %return
    end
end

% but some conditions are scrambled faces, not intact. Let's index those
for n = 1:length(names)
    if strfind(names{n},'_scrambled')%
        scrambled_idx(n)=1;
    else
        scrambled_idx(n)=0;
    end
end

%~~~~~~~Find intersections of the instances to examine more specific results
%scrambled faces
EA_scrambled = EA_idx.*scrambled_idx;% ".*" syntax means multiply the corresponding elements of each matrix or vector
AA_scrambled = AA_idx.*scrambled_idx;

%now we can isolate intact faces
EA_intact = EA_idx-EA_scrambled;%
AA_intact = AA_idx-AA_scrambled;


%% sanity checks
% if sum(ea_ex) ~= sum(ea_ex2)
%     disp('Trials from 1st and 2nd do not match');
%     return
% end
%
% if sum(aa_ex) ~= sum(aa_ex2)
%     disp('Trials from 1st and 2nd do not match');
%     return
% end
%
% if sum(ea_ex_corr) ~= sum(ea_ex_corr2)
%     disp('Trials from 1st and 2nd do not match');
%     return
% end
%
% if sum(aa_ex_corr) ~= sum(aa_ex_corr2)
%     disp('Trials from 1st and 2nd do not match');
%     return
% end
%
% if sum(ea_ex_incorr) ~= sum(ea_ex_incorr2)
%     disp('Trials from 1st and 2nd do not match');
%     return
% end
%
% if sum(aa_ex_incorr) ~= sum(aa_ex_incorr2)
%     disp('Trials from 1st and 2nd do not match');
%     return
% end

% create vox x stim matrix of beta values
% for i=1:length(tmp)
%
%     [b,r] = MAP_getROI(maskname, betamaps{i}, 'vox', 0, '');
%     bmat(:,i) = b{1}; % returns all voxels, whether or not they have NaNs
%     rmat(:,i) = r; % returns voxels excluding NaNs
%     meanbetas = nanmean(bmat(:,:));%get mean beta values from the ROI for each regressor
% end



%% create correlation matrix

%optional - toss values below a certain number as well (e.g., maybe a zero
%is equivalent to a NaN for your study. Say you are using a mask on "raw"
%BOLD data that does nothing to account for signal drop-out and extra-brain
%voxels are zeros instead of NaNs - here we can fix that
threshpats = 0; % 1 = YES
if threshpats == 1
    thresh = 0.01*mean(mean(rmat_condensed'));%you come up with your scheme - this example will threshold out anything <99% of the average (quite liberal thresholding)
    x1 =[]; %vector of filtered intensity values
    for p = 1:length(rmat_condensed(1,:))
        x1 = [x1 (rmat_condensed(:,p)>thresh)];
    end
    t = mean(x1');% average across columns - anything less than 1 indicates there are some zeros
    t1 = t'>=1; %threshold once more to voxels that had signal passing our threshold for EVERY pattern
    
    rmat_condensed = rmat_condensed(t1,:);
end

cm = corr(rmat_condensed);
cm(find(~triu(cm,1))) = NaN; % nan out correlations that are redundant
%if calculating similarity off overlapping indices, mask out lower half of
%matrix. Otherwise you will sample the same correlation twice. E.g.,
%"similarity of probes with probes"

%on the other hand, if you are using non-overlapping indices (e.g.,
%remembereds with forgottens) then you'll miss correlations by looking in
%only one half of the matrix. in this case use the full matrix (NaN
%diagonal still)

cm2 = corr(rmat_condensed); %another correlation matrix
cm2(cm2==1)=nan; %nan out the diagonal

%assign to res struct for easy saving
res.cm = cm;
res.cm2 = cm2;

%% global similarity measures
res.EA_w_EA = cm(logical(EA_intact),logical(EA_intact));
res.EA_w_EA_mean = nanmean(res.EA_w_EA(:));

res.AA_w_AA = cm(logical(AA_intact),logical(AA_intact));
res.AA_w_AA_mean = nanmean(res.AA_w_AA(:));

res.EA_w_AA = cm2(logical(EA_intact),logical(AA_intact));
res.EA_w_AA_mean = nanmean(res.EA_w_AA(:));

res.EA_w_Scene = cm2(logical(EA_intact),logical(Scene_idx));
res.EA_w_Scene_mean = nanmean(res.EA_w_Scene(:));


% %% stability item/town specific effects
%
% %within type (probe r1 with probe r2)
% probe_assigned_r1_w_probe_assigned_r2_repst = [cm2(logical(probe_assigned_env1_r1),logical(probe_assigned_env1_r2)) cm2(logical(probe_assigned_env2_r1),logical(probe_assigned_env2_r2)) cm2(logical(probe_assigned_env3_r1),logical(probe_assigned_env3_r2)) cm2(logical(probe_assigned_env4_r1),logical(probe_assigned_env4_r2)) cm2(logical(probe_assigned_env5_r1),logical(probe_assigned_env5_r2)) cm2(logical(probe_assigned_env6_r1),logical(probe_assigned_env6_r2)) cm2(logical(probe_assigned_env7_r1),logical(probe_assigned_env7_r2)) cm2(logical(probe_assigned_env8_r1),logical(probe_assigned_env8_r2)) cm2(logical(probe_assigned_env9_r1),logical(probe_assigned_env9_r2)) cm2(logical(probe_assigned_env10_r1),logical(probe_assigned_env10_r2)) cm2(logical(probe_assigned_env11_r1),logical(probe_assigned_env11_r2)) cm2(logical(probe_assigned_env12_r1),logical(probe_assigned_env12_r2))];
% probe_assigned_r1_w_probe_assigned_r2_repst_mean = nanmean(probe_assigned_r1_w_probe_assigned_r2_repst(:));
%
% probe_assigned_r1_w_probe_assigned_r2_repstcon = [cm2(logical(probe_assigned_env1_r1),logical(probe_assigned_r2-probe_assigned_env1_r2)) cm2(logical(probe_assigned_env2_r1),logical(probe_assigned_r2-probe_assigned_env2_r2)) cm2(logical(probe_assigned_env3_r1),logical(probe_assigned_r2-probe_assigned_env3_r2)) cm2(logical(probe_assigned_env4_r1),logical(probe_assigned_r2-probe_assigned_env4_r2)) cm2(logical(probe_assigned_env5_r1),logical(probe_assigned_r2-probe_assigned_env5_r2)) cm2(logical(probe_assigned_env6_r1),logical(probe_assigned_r2-probe_assigned_env6_r2)) cm2(logical(probe_assigned_env7_r1),logical(probe_assigned_r2-probe_assigned_env7_r2)) cm2(logical(probe_assigned_env8_r1),logical(probe_assigned_r2-probe_assigned_env8_r2)) cm2(logical(probe_assigned_env9_r1),logical(probe_assigned_r2-probe_assigned_env9_r2)) cm2(logical(probe_assigned_env10_r1),logical(probe_assigned_r2-probe_assigned_env10_r2)) cm2(logical(probe_assigned_env11_r1),logical(probe_assigned_r2-probe_assigned_env11_r2)) cm2(logical(probe_assigned_env12_r1),logical(probe_assigned_r2-probe_assigned_env12_r2))];
% probe_assigned_r1_w_probe_assigned_r2_repstcon_mean = nanmean(probe_assigned_r1_w_probe_assigned_r2_repstcon(:));
%
% probe_assigned_r2_w_probe_assigned_r1_repstcon = [cm2(logical(probe_assigned_env1_r2),logical(probe_assigned_r1-probe_assigned_env1_r1)) cm2(logical(probe_assigned_env2_r2),logical(probe_assigned_r1-probe_assigned_env2_r1)) cm2(logical(probe_assigned_env3_r2),logical(probe_assigned_r1-probe_assigned_env3_r1)) cm2(logical(probe_assigned_env4_r2),logical(probe_assigned_r1-probe_assigned_env4_r1)) cm2(logical(probe_assigned_env5_r2),logical(probe_assigned_r1-probe_assigned_env5_r1)) cm2(logical(probe_assigned_env6_r2),logical(probe_assigned_r1-probe_assigned_env6_r1)) cm2(logical(probe_assigned_env7_r2),logical(probe_assigned_r1-probe_assigned_env7_r1)) cm2(logical(probe_assigned_env8_r2),logical(probe_assigned_r1-probe_assigned_env8_r1)) cm2(logical(probe_assigned_env9_r2),logical(probe_assigned_r1-probe_assigned_env9_r1)) cm2(logical(probe_assigned_env10_r2),logical(probe_assigned_r1-probe_assigned_env10_r1)) cm2(logical(probe_assigned_env11_r2),logical(probe_assigned_r1-probe_assigned_env11_r1)) cm2(logical(probe_assigned_env12_r2),logical(probe_assigned_r1-probe_assigned_env12_r1))];
% probe_assigned_r2_w_probe_assigned_r1_repstcon_mean = nanmean(probe_assigned_r2_w_probe_assigned_r1_repstcon(:));
%
%
% %across type (probe with habit), r1 with habit
% probe_assigned_r1_w_habit_assigned_r1_repst = [cm2(logical(probe_assigned_env1_r1),logical(habit_assigned_env1_r1)) cm2(logical(probe_assigned_env2_r1),logical(habit_assigned_env2_r1)) cm2(logical(probe_assigned_env3_r1),logical(habit_assigned_env3_r1)) cm2(logical(probe_assigned_env4_r1),logical(habit_assigned_env4_r1)) cm2(logical(probe_assigned_env5_r1),logical(habit_assigned_env5_r1)) cm2(logical(probe_assigned_env6_r1),logical(habit_assigned_env6_r1)) cm2(logical(probe_assigned_env7_r1),logical(habit_assigned_env7_r1)) cm2(logical(probe_assigned_env8_r1),logical(habit_assigned_env8_r1)) cm2(logical(probe_assigned_env9_r1),logical(habit_assigned_env9_r1)) cm2(logical(probe_assigned_env10_r1),logical(habit_assigned_env10_r1)) cm2(logical(probe_assigned_env11_r1),logical(habit_assigned_env11_r1)) cm2(logical(probe_assigned_env12_r1),logical(habit_assigned_env12_r1))];
% probe_assigned_r1_w_habit_assigned_r1_repst_mean = nanmean(probe_assigned_r1_w_habit_assigned_r1_repst(:));
%
% probe_assigned_r1_w_habit_assigned_r1_repstcon = [cm2(logical(probe_assigned_env1_r1),logical(habit_assigned_r1-habit_assigned_env1_r1)) cm2(logical(probe_assigned_env2_r1),logical(habit_assigned_r1-habit_assigned_env2_r1)) cm2(logical(probe_assigned_env3_r1),logical(habit_assigned_r1-habit_assigned_env3_r1)) cm2(logical(probe_assigned_env4_r1),logical(habit_assigned_r1-habit_assigned_env4_r1)) cm2(logical(probe_assigned_env5_r1),logical(habit_assigned_r1-habit_assigned_env5_r1)) cm2(logical(probe_assigned_env6_r1),logical(habit_assigned_r1-habit_assigned_env6_r1)) cm2(logical(probe_assigned_env7_r1),logical(habit_assigned_r1-habit_assigned_env7_r1)) cm2(logical(probe_assigned_env8_r1),logical(habit_assigned_r1-habit_assigned_env8_r1)) cm2(logical(probe_assigned_env9_r1),logical(habit_assigned_r1-habit_assigned_env9_r1)) cm2(logical(probe_assigned_env10_r1),logical(habit_assigned_r1-habit_assigned_env10_r1)) cm2(logical(probe_assigned_env11_r1),logical(habit_assigned_r1-habit_assigned_env11_r1)) cm2(logical(probe_assigned_env12_r1),logical(habit_assigned_r1-habit_assigned_env12_r1))];
% probe_assigned_r1_w_habit_assigned_r1_repstcon_mean = nanmean(probe_assigned_r1_w_habit_assigned_r1_repstcon(:));
%
% %% reinstatement analysis
% %probe r1 with arrive
% probe_assigned_r1_w_probe_arriv_r1_repst = [cm2(logical(probe_assigned_env1_r1),logical(probe_arriv_env1_r1)) cm2(logical(probe_assigned_env2_r1),logical(probe_arriv_env2_r1)) cm2(logical(probe_assigned_env3_r1),logical(probe_arriv_env3_r1)) cm2(logical(probe_assigned_env4_r1),logical(probe_arriv_env4_r1)) cm2(logical(probe_assigned_env5_r1),logical(probe_arriv_env5_r1)) cm2(logical(probe_assigned_env6_r1),logical(probe_arriv_env6_r1)) cm2(logical(probe_assigned_env7_r1),logical(probe_arriv_env7_r1)) cm2(logical(probe_assigned_env8_r1),logical(probe_arriv_env8_r1)) cm2(logical(probe_assigned_env9_r1),logical(probe_arriv_env9_r1)) cm2(logical(probe_assigned_env10_r1),logical(probe_arriv_env10_r1)) cm2(logical(probe_assigned_env11_r1),logical(probe_arriv_env11_r1)) cm2(logical(probe_assigned_env12_r1),logical(probe_arriv_env12_r1))];
% probe_assigned_r1_w_probe_arriv_r1_repst_mean = nanmean(probe_assigned_r1_w_probe_arriv_r1_repst(:));
%
% probe_assigned_r1_w_habit_arriv_r1_repst = [cm2(logical(probe_assigned_env1_r1),logical(habit_arriv_env1_r1)) cm2(logical(probe_assigned_env2_r1),logical(habit_arriv_env2_r1)) cm2(logical(probe_assigned_env3_r1),logical(habit_arriv_env3_r1)) cm2(logical(probe_assigned_env4_r1),logical(habit_arriv_env4_r1)) cm2(logical(probe_assigned_env5_r1),logical(habit_arriv_env5_r1)) cm2(logical(probe_assigned_env6_r1),logical(habit_arriv_env6_r1)) cm2(logical(probe_assigned_env7_r1),logical(habit_arriv_env7_r1)) cm2(logical(probe_assigned_env8_r1),logical(habit_arriv_env8_r1)) cm2(logical(probe_assigned_env9_r1),logical(habit_arriv_env9_r1)) cm2(logical(probe_assigned_env10_r1),logical(habit_arriv_env10_r1)) cm2(logical(probe_assigned_env11_r1),logical(habit_arriv_env11_r1)) cm2(logical(probe_assigned_env12_r1),logical(habit_arriv_env12_r1))];
% probe_assigned_r1_w_habit_arriv_r1_repst_mean = nanmean(probe_assigned_r1_w_habit_arriv_r1_repst(:));
%
% probe_assigned_r1_w_probe_arriv_r1_repstcon = [cm2(logical(probe_assigned_env1_r1),logical(probe_arriv_r1-probe_arriv_env1_r1)) cm2(logical(probe_assigned_env2_r1),logical(probe_arriv_r1-probe_arriv_env2_r1)) cm2(logical(probe_assigned_env3_r1),logical(probe_arriv_r1-probe_arriv_env3_r1)) cm2(logical(probe_assigned_env4_r1),logical(probe_arriv_r1-probe_arriv_env4_r1)) cm2(logical(probe_assigned_env5_r1),logical(probe_arriv_r1-probe_arriv_env5_r1)) cm2(logical(probe_assigned_env6_r1),logical(probe_arriv_r1-probe_arriv_env6_r1)) cm2(logical(probe_assigned_env7_r1),logical(probe_arriv_r1-probe_arriv_env7_r1)) cm2(logical(probe_assigned_env8_r1),logical(probe_arriv_r1-probe_arriv_env8_r1)) cm2(logical(probe_assigned_env9_r1),logical(probe_arriv_r1-probe_arriv_env9_r1)) cm2(logical(probe_assigned_env10_r1),logical(probe_arriv_r1-probe_arriv_env10_r1)) cm2(logical(probe_assigned_env11_r1),logical(probe_arriv_r1-probe_arriv_env11_r1)) cm2(logical(probe_assigned_env12_r1),logical(probe_arriv_r1-probe_arriv_env12_r1))];
% probe_assigned_r1_w_probe_arriv_r1_repstcon_mean = nanmean(probe_assigned_r1_w_probe_arriv_r1_repstcon(:));

%% univariate control

meanbetas = nanmean(rmat_condensed(:,:));%get mean beta values from the ROI for each regressor


% probe_assigned_r1_mbeta = meanbetas(logical(probe_assigned_r1));
% probe_assigned_r1_mbeta_mean = nanmean(probe_assigned_r1_mbeta(:));
%
% probe_nav_r1_mbeta = meanbetas(logical(probe_nav_r1));
% probe_nav_r1_mbeta_mean = nanmean(probe_nav_r1_mbeta(:));
%
% probe_arriv_r1_mbeta = meanbetas(logical(probe_arriv_r1));
% probe_arriv_r1_mbeta_mean = nanmean(probe_arriv_r1_mbeta(:));
%
% habit_assigned_r1_mbeta = meanbetas(logical(habit_assigned_r1));
% habit_assigned_r1_mbeta_mean = nanmean(habit_assigned_r1_mbeta(:));
%
% habit_nav_r1_mbeta = meanbetas(logical(habit_nav_r1));
% habit_nav_r1_mbeta_mean = nanmean(habit_nav_r1_mbeta(:));
%
% habit_arriv_r1_mbeta = meanbetas(logical(habit_arriv_r1));
% habit_arriv_r1_mbeta_mean = nanmean(habit_arriv_r1_mbeta(:));
%
% %second repetition indices
% probe_assigned_r2_mbeta = meanbetas(logical(probe_assigned_r2));
% probe_assigned_r2_mbeta_mean = nanmean(probe_assigned_r2_mbeta(:));
%
% probe_nav_r2_mbeta = meanbetas(logical(probe_nav_r2));
% probe_nav_r2_mbeta_mean = nanmean(probe_nav_r2_mbeta(:));
%
% probe_arriv_r2_mbeta = meanbetas(logical(probe_arriv_r2));
% probe_arriv_r2_mbeta_mean = nanmean(probe_arriv_r2_mbeta(:));


%probe_assigned_r1_w_probe_assigned_r1 = cm(logical(probe_assigned_r1),logical(probe_assigned_r1));

% %% Activity-correlation control

% %are trial-by-trial rsa scores correlated with univariate signal?
% r_ea_mbetawr = corr(ea_corr_mbeta', nanmean(ea_r_with_ea)');
% r_aa_mbetawr = corr(aa_corr_mbeta', nanmean(aa_r_with_aa)');
% r_ea_mbetawr_inc = corr(ea_incorr_mbeta', nanmean(ea_f_with_ea)');
% r_aa_mbetawr_inc = corr(aa_incorr_mbeta', nanmean(aa_f_with_aa)');
% r_ea_mbetawr2 = corr(ea_corr_mbeta2', nanmean(ea_r_with_ea2)');
% r_aa_mbetawr2 = corr(aa_corr_mbeta2', nanmean(aa_r_with_aa2)');
% r_ea_mbetawr2_inc = corr(ea_incorr_mbeta2', nanmean(ea_f_with_ea2)');
% r_aa_mbetawr2_inc = corr(aa_incorr_mbeta2', nanmean(aa_f_with_aa2)');
%
% %are trial 1-2 stability measurements correlated with trial 1-2 activity
% %differences?
% r_ea_betawstab = corr(ea1_ea2_r_actdiff', ea1_ea2_r);
% r_aa_betawstab = corr(aa1_aa2_r_actdiff', aa1_aa2_r);
% r_ea_betawstab_inc = corr(ea1_ea2_f_actdiff', ea1_ea2_f);
% r_aa_betawstab_inc = corr(aa1_aa2_f_actdiff', aa1_aa2_f);

%fisher's z
%fish = 0.5*log((1+CM2)./(1-CM2))

%% Plots
%plot matrices of interest

subplot(2,2,1), imagesc(cm);
title('Within-cond corrmat')
colormap('hot'); % set the colorscheme
colorbar; % enable colorbar

subplot(2,2,2), imagesc(cm2);
title('Overall corrmat')
colormap('hot'); % set the colorscheme
colorbar; % enable colorbar

subplot(2,2,3), imagesc(cm(logical(EA_intact),logical(EA_intact)));
title('EA with EA across blocks and runs')
colormap('hot'); % set the colorscheme
colorbar; % enable colorbar

subplot(2,2,4), imagesc(cm2(logical(EA_intact),logical(Scene_idx)));
title('EA with Scene across blocks and runs')
colormap('hot'); % set the colorscheme
colorbar; % enable colorbar

% save plot
plot_savename = [S.group_mvpa_dir '/Rcorrs_' S.subj_id '_' mask '_' weights_str '_' S.exp_name '_corrmats.png'];
saveas(gcf,plot_savename);

%% Save data
savename = [S.group_mvpa_dir '/Rcorrs_' S.subj_id '_' mask '_' weights_str '_' S.exp_name '.mat'];
save(savename, 'res');





end
