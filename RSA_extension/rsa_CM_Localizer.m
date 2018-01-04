function [] = rsa_CM_Localizer(Sub, Mask, TRsperRun)
% code for RSA analysis.

%For alternate code demo purposes, this is built around both 4D and 3D
%image file types. The use case for 4D includes scenarios like "raw" BOLD data or residual time-series

%warning: as of 1/3/2018; unconcatenated 4D analysis has not been debugged.
%at this time only use concatenated 3D

%% flags and parameters
S.TR = 2;
theseTRWeights = [0 0.25 0.5 0.25 0];
funcftype = '.nii';
runs_concat = 1; %1 = typical SPM analysis; will have continuous onsets concatenated across runs. 0 = But you might not have bothered creating such a file, or in SST case we are using files from FSL. In this case, the onsets are assumed to "reset" for each run ('raw' unconcatenated onsets)

gen_onsetsTR = 1; %1=yes. Typically, you'll use an onsets.mat file with tr-by-tr onsets and names (as used for a beta-series). But if you only have a traditional GLM model with one name for multile onsets, setting this flag to 1 will auto-populate unique but related names (e.g., Face_1; Face_2...)
%sub = Sub{1};

%Subject ID/number
par.substr = ['CM' Sub{1}];
S.subj_id = par.substr;

mask = Mask;
S.exp_name = 'CM_Localizer';
study_prefix = 'CM';

S.inputformat = 'raw'; % are we working with BOLDs/timeseries ('raw') or with beta maps ('betas')?

S.onsets_filename = [S.subj_id '_localizer_onsets_test'];

%specify preprocessing level of BOLDs
preproc_lvl = ''; % 'a' for slice-time-only, 'u' for realigned-only, 'ua' for realign+unwarped, 'swua' for smoothed, normalized, and... you get the picture. Modify as needed if you changed SPM's prefix append defaults
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
    if exist([S.mvpa_dir '/RSA_data/' S.subj_id '_' mask '_condensedpats.mat']);
        load([S.mvpa_dir '/RSA_data/' S.subj_id '_' mask '_condensedpats.mat']);
        load([S.mvpa_dir '/RSA_data/' S.subj_id '_onsets_expanded.mat']);
    else
        
        %% load in pattern data
        if ImgDims == 4
            
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
                runidxnum = str2double(raw_filenames{idx,1}(runref_indices(1)+5:runref_indices(2)-1));
                raw_filenames{idx,3} = runidxnum;
            end
            
            b = sortrows(raw_filenames, 3);
            raw_filenames = b(:,1);
            
            imgslength = length(raw_filenames);
            
            %% iterate through 3D frames to extract all patterns
            for i=1:imgslength
                betamaps{i} = raw_filenames{i};%[tmp.name ',' num2str(i)];
                
                [b,r] = MAP_getROI(maskname, betamaps{i}, 'vox', 0, '');
                bmat_t(:,i) = b{1}; % returns all voxels, whether or not they have NaNs
                rmat_t(:,i) = r; % returns voxels excluding NaNs
                %meanbetas_t = nanmean(bmat_t(:,:));%get mean beta values from the ROI for each regressor
                
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
%end

%% create indices for patterns of interest

% first index different nav trial stages
for n = 1:length(names)
    if strfind(names{n},'ASSIGNED')%
        assgn_idx(n)=1;
        nav_idx(n)=0;
        arriv_idx(n) = 0;
        othercond_idx(n) = 0;
    elseif strfind(names{n},'NAVIGATE')%
        assgn_idx(n)=0;
        nav_idx(n)=1;
        arriv_idx(n) = 0;
        othercond_idx(n) = 0;
    elseif strfind(names{n},'ARRIVED')%
        assgn_idx(n)=0;
        nav_idx(n)=0;
        arriv_idx(n) = 1;
        othercond_idx(n) = 0;
    else
        othercond_idx(n) = 1;
        %fprintf('Error! Neither Assigned, Navigate, nor Arrived!\n')
        %return
    end
end


% now index run cond
for n = 1:length(names)
    if strfind(names{n},'_habit_')%
        habit_idx(n)=1;
        probe_idx(n)=0;
        
    elseif strfind(names{n},'_shortcut_')%
        habit_idx(n)=0;
        probe_idx(n)=1;
        
        
    else
        %othercond_idx(n) = 1;
        fprintf('Error! Neither habit nor shortcut trials!\n')
        return
    end
end



% now index town
for n = 1:length(names)
    if strfind(names{n},'_env1_')%
        env1_idx(n)=1;
        env2_idx(n)=0;
        env3_idx(n)=0;
        env4_idx(n)=0;
        env5_idx(n)=0;
        env6_idx(n)=0;
        env7_idx(n)=0;
        env8_idx(n)=0;
        env9_idx(n)=0;
        env10_idx(n)=0;
        env11_idx(n)=0;
        env12_idx(n)=0;
    elseif strfind(names{n},'_env2_')%
        env1_idx(n)=0;
        env2_idx(n)=1;
        env3_idx(n)=0;
        env4_idx(n)=0;
        env5_idx(n)=0;
        env6_idx(n)=0;
        env7_idx(n)=0;
        env8_idx(n)=0;
        env9_idx(n)=0;
        env10_idx(n)=0;
        env11_idx(n)=0;
        env12_idx(n)=0;
    elseif strfind(names{n},'_env3_')%
        env1_idx(n)=0;
        env2_idx(n)=0;
        env3_idx(n)=1;
        env4_idx(n)=0;
        env5_idx(n)=0;
        env6_idx(n)=0;
        env7_idx(n)=0;
        env8_idx(n)=0;
        env9_idx(n)=0;
        env10_idx(n)=0;
        env11_idx(n)=0;
        env12_idx(n)=0;
    elseif strfind(names{n},'_env4_')%
        env1_idx(n)=0;
        env2_idx(n)=0;
        env3_idx(n)=0;
        env4_idx(n)=1;
        env5_idx(n)=0;
        env6_idx(n)=0;
        env7_idx(n)=0;
        env8_idx(n)=0;
        env9_idx(n)=0;
        env10_idx(n)=0;
        env11_idx(n)=0;
        env12_idx(n)=0;
    elseif strfind(names{n},'_env5_')%
        env1_idx(n)=0;
        env2_idx(n)=0;
        env3_idx(n)=0;
        env4_idx(n)=0;
        env5_idx(n)=1;
        env6_idx(n)=0;
        env7_idx(n)=0;
        env8_idx(n)=0;
        env9_idx(n)=0;
        env10_idx(n)=0;
        env11_idx(n)=0;
        env12_idx(n)=0;
    elseif strfind(names{n},'_env6_')%
        env1_idx(n)=0;
        env2_idx(n)=0;
        env3_idx(n)=0;
        env4_idx(n)=0;
        env5_idx(n)=0;
        env6_idx(n)=1;
        env7_idx(n)=0;
        env8_idx(n)=0;
        env9_idx(n)=0;
        env10_idx(n)=0;
        env11_idx(n)=0;
        env12_idx(n)=0;
    elseif strfind(names{n},'_env7_')%
        env1_idx(n)=0;
        env2_idx(n)=0;
        env3_idx(n)=0;
        env4_idx(n)=0;
        env5_idx(n)=0;
        env6_idx(n)=0;
        env7_idx(n)=1;
        env8_idx(n)=0;
        env9_idx(n)=0;
        env10_idx(n)=0;
        env11_idx(n)=0;
        env12_idx(n)=0;
    elseif strfind(names{n},'_env8_')%
        env1_idx(n)=0;
        env2_idx(n)=0;
        env3_idx(n)=0;
        env4_idx(n)=0;
        env5_idx(n)=0;
        env6_idx(n)=0;
        env7_idx(n)=0;
        env8_idx(n)=1;
        env9_idx(n)=0;
        env10_idx(n)=0;
        env11_idx(n)=0;
        env12_idx(n)=0;
    elseif strfind(names{n},'_env9_')%
        env1_idx(n)=0;
        env2_idx(n)=0;
        env3_idx(n)=0;
        env4_idx(n)=0;
        env5_idx(n)=0;
        env6_idx(n)=0;
        env7_idx(n)=0;
        env8_idx(n)=0;
        env9_idx(n)=1;
        env10_idx(n)=0;
        env11_idx(n)=0;
        env12_idx(n)=0;
    elseif strfind(names{n},'_env10_')%
        env1_idx(n)=0;
        env2_idx(n)=0;
        env3_idx(n)=0;
        env4_idx(n)=0;
        env5_idx(n)=0;
        env6_idx(n)=0;
        env7_idx(n)=0;
        env8_idx(n)=0;
        env9_idx(n)=0;
        env10_idx(n)=1;
        env11_idx(n)=0;
        env12_idx(n)=0;
    elseif strfind(names{n},'_env11_')%
        env1_idx(n)=0;
        env2_idx(n)=0;
        env3_idx(n)=0;
        env4_idx(n)=0;
        env5_idx(n)=0;
        env6_idx(n)=0;
        env7_idx(n)=0;
        env8_idx(n)=0;
        env9_idx(n)=0;
        env10_idx(n)=0;
        env11_idx(n)=1;
        env12_idx(n)=0;
    elseif strfind(names{n},'_env12_')%
        env1_idx(n)=0;
        env2_idx(n)=0;
        env3_idx(n)=0;
        env4_idx(n)=0;
        env5_idx(n)=0;
        env6_idx(n)=0;
        env7_idx(n)=0;
        env8_idx(n)=0;
        env9_idx(n)=0;
        env10_idx(n)=0;
        env11_idx(n)=0;
        env12_idx(n)=1;
        
    else
        %othercond_idx(n) = 1;
        fprintf('Error! None of the towns are represented!\n')
        return
    end
end


% now index repetition
for n = 1:length(names)
    if strfind(names{n},'_rep1')%
        rep1_idx(n)=1;
        rep2_idx(n)=0;
        
    elseif strfind(names{n},'_rep2')%
        rep1_idx(n)=0;
        rep2_idx(n)=1;
        
        
    else
        %othercond_idx(n) = 1;
        fprintf('Error! Neither repetition is represented!\n')
        return
    end
end

%~~~~~~~Find intersections of the instances to examine more specific results
%first repetition indices
probe_r1 = probe_idx.*rep1_idx;%
probe_assigned_r1 = assgn_idx.*probe_r1;
probe_nav_r1 = nav_idx.*probe_r1;
probe_arriv_r1 = arriv_idx.*probe_r1;

habit_assigned_r1 = assgn_idx.*habit_idx;
habit_nav_r1 = nav_idx.*habit_idx;
habit_arriv_r1 = arriv_idx.*habit_idx;

%second repetition indices
probe_r2 = probe_idx.*rep2_idx;%
probe_assigned_r2 = assgn_idx.*probe_r2;
probe_nav_r2 = nav_idx.*probe_r2;
probe_arriv_r2 = arriv_idx.*probe_r2;


% env1 specific indices
probe_assigned_env1_r1 = probe_assigned_r1.*env1_idx;
probe_nav_env1_r1 = probe_nav_r1.*env1_idx;
probe_arriv_env1_r1 = probe_arriv_r1.*env1_idx;

probe_assigned_env1_r2 = probe_assigned_r2.*env1_idx;
probe_nav_env1_r2 = probe_nav_r2.*env1_idx;
probe_arriv_env1_r2 = probe_arriv_r2.*env1_idx;

habit_assigned_env1_r1 = habit_assigned_r1.*env1_idx;
habit_nav_env1_r1 = habit_nav_r1.*env1_idx;
habit_arriv_env1_r1 = habit_arriv_r1.*env1_idx;

% env2 specific indices
probe_assigned_env2_r1 = probe_assigned_r1.*env2_idx;
probe_nav_env2_r1 = probe_nav_r1.*env2_idx;
probe_arriv_env2_r1 = probe_arriv_r1.*env2_idx;

probe_assigned_env2_r2 = probe_assigned_r2.*env2_idx;
probe_nav_env2_r2 = probe_nav_r2.*env2_idx;
probe_arriv_env2_r2 = probe_arriv_r2.*env2_idx;

habit_assigned_env2_r1 = habit_assigned_r1.*env2_idx;
habit_nav_env2_r1 = habit_nav_r1.*env2_idx;
habit_arriv_env2_r1 = habit_arriv_r1.*env2_idx;

% env3 specific indices
probe_assigned_env3_r1 = probe_assigned_r1.*env3_idx;
probe_nav_env3_r1 = probe_nav_r1.*env3_idx;
probe_arriv_env3_r1 = probe_arriv_r1.*env3_idx;

probe_assigned_env3_r2 = probe_assigned_r2.*env3_idx;
probe_nav_env3_r2 = probe_nav_r2.*env3_idx;
probe_arriv_env3_r2 = probe_arriv_r2.*env3_idx;

habit_assigned_env3_r1 = habit_assigned_r1.*env3_idx;
habit_nav_env3_r1 = habit_nav_r1.*env3_idx;
habit_arriv_env3_r1 = habit_arriv_r1.*env3_idx;

% env4 specific indices
probe_assigned_env4_r1 = probe_assigned_r1.*env4_idx;
probe_nav_env4_r1 = probe_nav_r1.*env4_idx;
probe_arriv_env4_r1 = probe_arriv_r1.*env4_idx;

probe_assigned_env4_r2 = probe_assigned_r2.*env4_idx;
probe_nav_env4_r2 = probe_nav_r2.*env4_idx;
probe_arriv_env4_r2 = probe_arriv_r2.*env4_idx;

habit_assigned_env4_r1 = habit_assigned_r1.*env4_idx;
habit_nav_env4_r1 = habit_nav_r1.*env4_idx;
habit_arriv_env4_r1 = habit_arriv_r1.*env4_idx;

% env5 specific indices
probe_assigned_env5_r1 = probe_assigned_r1.*env5_idx;
probe_nav_env5_r1 = probe_nav_r1.*env5_idx;
probe_arriv_env5_r1 = probe_arriv_r1.*env5_idx;

probe_assigned_env5_r2 = probe_assigned_r2.*env5_idx;
probe_nav_env5_r2 = probe_nav_r2.*env5_idx;
probe_arriv_env5_r2 = probe_arriv_r2.*env5_idx;

habit_assigned_env5_r1 = habit_assigned_r1.*env5_idx;
habit_nav_env5_r1 = habit_nav_r1.*env5_idx;
habit_arriv_env5_r1 = habit_arriv_r1.*env5_idx;

% env6 specific indices
probe_assigned_env6_r1 = probe_assigned_r1.*env6_idx;
probe_nav_env6_r1 = probe_nav_r1.*env6_idx;
probe_arriv_env6_r1 = probe_arriv_r1.*env6_idx;

probe_assigned_env6_r2 = probe_assigned_r2.*env6_idx;
probe_nav_env6_r2 = probe_nav_r2.*env6_idx;
probe_arriv_env6_r2 = probe_arriv_r2.*env6_idx;

habit_assigned_env6_r1 = habit_assigned_r1.*env6_idx;
habit_nav_env6_r1 = habit_nav_r1.*env6_idx;
habit_arriv_env6_r1 = habit_arriv_r1.*env6_idx;

% env7 specific indices
probe_assigned_env7_r1 = probe_assigned_r1.*env7_idx;
probe_nav_env7_r1 = probe_nav_r1.*env7_idx;
probe_arriv_env7_r1 = probe_arriv_r1.*env7_idx;

probe_assigned_env7_r2 = probe_assigned_r2.*env7_idx;
probe_nav_env7_r2 = probe_nav_r2.*env7_idx;
probe_arriv_env7_r2 = probe_arriv_r2.*env7_idx;

habit_assigned_env7_r1 = habit_assigned_r1.*env7_idx;
habit_nav_env7_r1 = habit_nav_r1.*env7_idx;
habit_arriv_env7_r1 = habit_arriv_r1.*env7_idx;

% env8 specific indices
probe_assigned_env8_r1 = probe_assigned_r1.*env8_idx;
probe_nav_env8_r1 = probe_nav_r1.*env8_idx;
probe_arriv_env8_r1 = probe_arriv_r1.*env8_idx;

probe_assigned_env8_r2 = probe_assigned_r2.*env8_idx;
probe_nav_env8_r2 = probe_nav_r2.*env8_idx;
probe_arriv_env8_r2 = probe_arriv_r2.*env8_idx;

habit_assigned_env8_r1 = habit_assigned_r1.*env8_idx;
habit_nav_env8_r1 = habit_nav_r1.*env8_idx;
habit_arriv_env8_r1 = habit_arriv_r1.*env8_idx;

% env9 specific indices
probe_assigned_env9_r1 = probe_assigned_r1.*env9_idx;
probe_nav_env9_r1 = probe_nav_r1.*env9_idx;
probe_arriv_env9_r1 = probe_arriv_r1.*env9_idx;

probe_assigned_env9_r2 = probe_assigned_r2.*env9_idx;
probe_nav_env9_r2 = probe_nav_r2.*env9_idx;
probe_arriv_env9_r2 = probe_arriv_r2.*env9_idx;

habit_assigned_env9_r1 = habit_assigned_r1.*env9_idx;
habit_nav_env9_r1 = habit_nav_r1.*env9_idx;
habit_arriv_env9_r1 = habit_arriv_r1.*env9_idx;

% env10 specific indices
probe_assigned_env10_r1 = probe_assigned_r1.*env10_idx;
probe_nav_env10_r1 = probe_nav_r1.*env10_idx;
probe_arriv_env10_r1 = probe_arriv_r1.*env10_idx;

probe_assigned_env10_r2 = probe_assigned_r2.*env10_idx;
probe_nav_env10_r2 = probe_nav_r2.*env10_idx;
probe_arriv_env10_r2 = probe_arriv_r2.*env10_idx;

habit_assigned_env10_r1 = habit_assigned_r1.*env10_idx;
habit_nav_env10_r1 = habit_nav_r1.*env10_idx;
habit_arriv_env10_r1 = habit_arriv_r1.*env10_idx;

% env11 specific indices
probe_assigned_env11_r1 = probe_assigned_r1.*env11_idx;
probe_nav_env11_r1 = probe_nav_r1.*env11_idx;
probe_arriv_env11_r1 = probe_arriv_r1.*env11_idx;

probe_assigned_env11_r2 = probe_assigned_r2.*env11_idx;
probe_nav_env11_r2 = probe_nav_r2.*env11_idx;
probe_arriv_env11_r2 = probe_arriv_r2.*env11_idx;

habit_assigned_env11_r1 = habit_assigned_r1.*env11_idx;
habit_nav_env11_r1 = habit_nav_r1.*env11_idx;
habit_arriv_env11_r1 = habit_arriv_r1.*env11_idx;

% env12 specific indices
probe_assigned_env12_r1 = probe_assigned_r1.*env12_idx;
probe_nav_env12_r1 = probe_nav_r1.*env12_idx;
probe_arriv_env12_r1 = probe_arriv_r1.*env12_idx;

probe_assigned_env12_r2 = probe_assigned_r2.*env12_idx;
probe_nav_env12_r2 = probe_nav_r2.*env12_idx;
probe_arriv_env12_r2 = probe_arriv_r2.*env12_idx;

habit_assigned_env12_r1 = habit_assigned_r1.*env12_idx;
habit_nav_env12_r1 = habit_nav_r1.*env12_idx;
habit_arriv_env12_r1 = habit_arriv_r1.*env12_idx;




%% convert indices to logicals
% fut=logical(fut);
% curr=logical(curr);
%
%
% probe_r1 = probe_idx.*rep1_idx;%
% probe_assigned_r1 = assgn_idx.*probe_r1;
% probe_nav_r1 = nav_idx.*probe_r1;
% probe_arriv_r1 = arriv_idx.*probe_r1;
%
% habit_assigned_r1 = assgn_idx.*habit_idx;
% habit_nav_r1 = nav_idx.*habit_idx;
% habit_arriv_r1 = arriv_idx.*habit_idx;
%
% %second repetition indices
% probe_r2 = probe_idx.*rep2_idx;%
% probe_assigned_r2 = assgn_idx.*probe_r2;
% probe_nav_r2 = nav_idx.*probe_r2;
% probe_arriv_r2 = arriv_idx.*probe_r2;
%
%
% % env1 specific indices
% probe_assigned_env1_r1 = probe_assigned_r1.*env1_idx;
% probe_nav_env1_r1 = probe_nav_r1.*env1_idx;
% probe_arriv_env1_r1 = probe_arriv_r1.*env1_idx;
%
% probe_assigned_env1_r2 = probe_assigned_r2.*env1_idx;
% probe_nav_env1_r2 = probe_nav_r2.*env1_idx;
% probe_arriv_env1_r2 = probe_arriv_r2.*env1_idx;
%
% habit_assigned_env1_r1 = habit_assigned_r1.*env1_idx;
% habit_nav_env1_r1 = habit_nav_r1.*env1_idx;
% habit_arriv_env1_r1 = habit_arriv_r1.*env1_idx;
%
% % env2 specific indices
% probe_assigned_env2_r1 = probe_assigned_r1.*env2_idx;
% probe_nav_env2_r1 = probe_nav_r1.*env2_idx;
% probe_arriv_env2_r1 = probe_arriv_r1.*env2_idx;
%
% probe_assigned_env2_r2 = probe_assigned_r2.*env2_idx;
% probe_nav_env2_r2 = probe_nav_r2.*env2_idx;
% probe_arriv_env2_r2 = probe_arriv_r2.*env2_idx;
%
% habit_assigned_env2_r1 = habit_assigned_r1.*env2_idx;
% habit_nav_env2_r1 = habit_nav_r1.*env2_idx;
% habit_arriv_env2_r1 = habit_arriv_r1.*env2_idx;
%
% % env3 specific indices
% probe_assigned_env3_r1 = probe_assigned_r1.*env3_idx;
% probe_nav_env3_r1 = probe_nav_r1.*env3_idx;
% probe_arriv_env3_r1 = probe_arriv_r1.*env3_idx;
%
% probe_assigned_env3_r2 = probe_assigned_r2.*env3_idx;
% probe_nav_env3_r2 = probe_nav_r2.*env3_idx;
% probe_arriv_env3_r2 = probe_arriv_r2.*env3_idx;
%
% habit_assigned_env3_r1 = habit_assigned_r1.*env3_idx;
% habit_nav_env3_r1 = habit_nav_r1.*env3_idx;
% habit_arriv_env3_r1 = habit_arriv_r1.*env3_idx;
%
% % env4 specific indices
% probe_assigned_env4_r1 = probe_assigned_r1.*env4_idx;
% probe_nav_env4_r1 = probe_nav_r1.*env4_idx;
% probe_arriv_env4_r1 = probe_arriv_r1.*env4_idx;
%
% probe_assigned_env4_r2 = probe_assigned_r2.*env4_idx;
% probe_nav_env4_r2 = probe_nav_r2.*env4_idx;
% probe_arriv_env4_r2 = probe_arriv_r2.*env4_idx;
%
% habit_assigned_env4_r1 = habit_assigned_r1.*env4_idx;
% habit_nav_env4_r1 = habit_nav_r1.*env4_idx;
% habit_arriv_env4_r1 = habit_arriv_r1.*env4_idx;
%
% % env5 specific indices
% probe_assigned_env5_r1 = probe_assigned_r1.*env5_idx;
% probe_nav_env5_r1 = probe_nav_r1.*env5_idx;
% probe_arriv_env5_r1 = probe_arriv_r1.*env5_idx;
%
% probe_assigned_env5_r2 = probe_assigned_r2.*env5_idx;
% probe_nav_env5_r2 = probe_nav_r2.*env5_idx;
% probe_arriv_env5_r2 = probe_arriv_r2.*env5_idx;
%
% habit_assigned_env5_r1 = habit_assigned_r1.*env5_idx;
% habit_nav_env5_r1 = habit_nav_r1.*env5_idx;
% habit_arriv_env5_r1 = habit_arriv_r1.*env5_idx;
%
% % env6 specific indices
% probe_assigned_env6_r1 = probe_assigned_r1.*env6_idx;
% probe_nav_env6_r1 = probe_nav_r1.*env6_idx;
% probe_arriv_env6_r1 = probe_arriv_r1.*env6_idx;
%
% probe_assigned_env6_r2 = probe_assigned_r2.*env6_idx;
% probe_nav_env6_r2 = probe_nav_r2.*env6_idx;
% probe_arriv_env6_r2 = probe_arriv_r2.*env6_idx;
%
% habit_assigned_env6_r1 = habit_assigned_r1.*env6_idx;
% habit_nav_env6_r1 = habit_nav_r1.*env6_idx;
% habit_arriv_env6_r1 = habit_arriv_r1.*env6_idx;
%
% % env7 specific indices
% probe_assigned_env7_r1 = probe_assigned_r1.*env7_idx;
% probe_nav_env7_r1 = probe_nav_r1.*env7_idx;
% probe_arriv_env7_r1 = probe_arriv_r1.*env7_idx;
%
% probe_assigned_env7_r2 = probe_assigned_r2.*env7_idx;
% probe_nav_env7_r2 = probe_nav_r2.*env7_idx;
% probe_arriv_env7_r2 = probe_arriv_r2.*env7_idx;
%
% habit_assigned_env7_r1 = habit_assigned_r1.*env7_idx;
% habit_nav_env7_r1 = habit_nav_r1.*env7_idx;
% habit_arriv_env7_r1 = habit_arriv_r1.*env7_idx;
%
% % env8 specific indices
% probe_assigned_env8_r1 = probe_assigned_r1.*env8_idx;
% probe_nav_env8_r1 = probe_nav_r1.*env8_idx;
% probe_arriv_env8_r1 = probe_arriv_r1.*env8_idx;
%
% probe_assigned_env8_r2 = probe_assigned_r2.*env8_idx;
% probe_nav_env8_r2 = logical()
% probe_arriv_env8_r2 = logical()
%
% habit_assigned_env8_r1 = logical()
% habit_nav_env8_r1 = logical()
% habit_arriv_env8_r1 = logical()
%
% % env9 specific indices
% probe_assigned_env9_r1 = logical()
% probe_nav_env9_r1 = logical()
% probe_arriv_env9_r1 = logical()
%
% probe_assigned_env9_r2 = logical()
% probe_nav_env9_r2 = logical()
% probe_arriv_env9_r2 = logical()
%
% habit_assigned_env9_r1 = logical()
% habit_nav_env9_r1 = logical()
% habit_arriv_env9_r1 = logical()
%
% % env10 specific indices
% probe_assigned_env10_r1 = logical()
% probe_nav_env10_r1 = logical()
% probe_arriv_env10_r1 = logical()
%
% probe_assigned_env10_r2 = logical()
% probe_nav_env10_r2 = logical()
% probe_arriv_env10_r2 = logical()
%
% habit_assigned_env10_r1 = logical()
% habit_nav_env10_r1 = logical()
% habit_arriv_env10_r1 = logical()
%
% % env11 specific indices
% probe_assigned_env11_r1 = logical()
% probe_nav_env11_r1 = logical()
% probe_arriv_env11_r1 = logical()
%
% probe_assigned_env11_r2 = logical()
% probe_nav_env11_r2 = logical()
% probe_arriv_env11_r2 = logical()
%
% habit_assigned_env11_r1 = logical()
% habit_nav_env11_r1 = logical()
% habit_arriv_env11_r1 = logical()
%
% % env12 specific indices
% probe_assigned_env12_r1 = logical()
% probe_nav_env12_r1 = logical()
% probe_arriv_env12_r1 = logical()
%
% probe_assigned_env12_r2 = logical()
% probe_nav_env12_r2 = logical()
% probe_arriv_env12_r2 = logical()
%
% habit_assigned_env12_r1 = logical()
% habit_nav_env12_r1 = logical()
% habit_arriv_env12_r1 = logical()







%sanity checks
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
%is equivalent to a NaN for your study.
% !! Currently hard-coded for only 2 images !!
% threshpats = 1; % 1 = YES
% if threshpats == 1
%     thresh = 5;
%     x1 = find(rmat(:,1)>thresh);
%     x2 = find(rmat(:,2)>thresh);
%     z = union(x1,x2); %keep only indices shared between patterns
%     rmatthresh = rmat(z,:)
%     %rmatthresh = rmat(x2,:)
%     cm = corr(rmatthresh);
% else





cm = corr(rmat_condensed);
cm(find(~triu(cm,1))) = NaN; % nan out correlations that are redundant
%if calculating similarity off overlapping indices, mask out lower half of
%matrix. Otherwise you will sample the same correlation twice. E.g.,
%"similarity of probes with probes"

%on the other hand, if you are using non-overlapping indices (e.g.,
%remembereds with forgottens) then you'll miss correlations by looking in
%only one half of the matrix. in this case use the full matrix (NaN
%diagonal still)


%cm = (cm-(nanmean(cm(:))))/nanstd(cm(:)); %convert to z-scores against mean

cm2 = corr(rmat_condensed); %another correlation matrix
cm2(cm2==1)=nan; %nan out the diagonal


%% global similarity measures
probe_assigned_r1_w_probe_assigned_r1 = cm(logical(probe_assigned_r1),logical(probe_assigned_r1));
probe_assigned_r1_w_probe_assigned_r1_mean = nanmean(probe_assigned_r1_w_probe_assigned_r1(:));

probe_nav_r1_w_probe_nav_r1 = cm(logical(probe_nav_r1),logical(probe_nav_r1));
probe_nav_r1_w_probe_nav_r1_mean = nanmean(probe_nav_r1_w_probe_nav_r1(:));

probe_arriv_r1_w_probe_arriv_r1 = cm(logical(probe_arriv_r1),logical(probe_arriv_r1));
probe_arriv_r1_w_probe_arriv_r1_mean = nanmean(probe_arriv_r1_w_probe_arriv_r1(:));

habit_assigned_r1_w_habit_assigned_r1 = cm(logical(habit_assigned_r1),logical(habit_assigned_r1));
habit_assigned_r1_w_habit_assigned_r1_mean = nanmean(habit_assigned_r1_w_habit_assigned_r1(:));

habit_nav_r1_w_habit_nav_r1 = cm(logical(habit_nav_r1),logical(habit_nav_r1));
habit_nav_r1_w_habit_nav_r1_mean = nanmean(habit_nav_r1_w_habit_nav_r1(:));

habit_arriv_r1_w_habit_arriv_r1 = cm(logical(habit_arriv_r1),logical(habit_arriv_r1));
habit_arriv_r1_w_habit_arriv_r1_mean = nanmean(habit_arriv_r1_w_habit_arriv_r1(:));

%second repetition indices
probe_assigned_r2_w_probe_assigned_r2 = cm(logical(probe_assigned_r2),logical(probe_assigned_r2));
probe_assigned_r2_w_probe_assigned_r2_mean = nanmean(probe_assigned_r2_w_probe_assigned_r2(:));

probe_nav_r2_w_probe_nav_r2 = cm(logical(probe_nav_r2),logical(probe_nav_r2));
probe_nav_r2_w_probe_nav_r2_mean = nanmean(probe_nav_r2_w_probe_nav_r2(:));

probe_arriv_r2_w_probe_arriv_r2 = cm(logical(probe_arriv_r2),logical(probe_arriv_r2));
probe_arriv_r2_w_probe_arriv_r2_mean = nanmean(probe_arriv_r2_w_probe_arriv_r2(:));

%% stability item/town specific effects

if strcmp(sub,'20') %hard coded in repair for sub 20 who had only one rep of 10-12, and 3 reps of 4-6
    probe_assigned_env10_r1 = probe_assigned_env10_r2;
    probe_assigned_env11_r1 = probe_assigned_env11_r2;
    probe_assigned_env12_r1 = probe_assigned_env12_r2;
    
    probe_nav_env10_r1 = probe_nav_env10_r2;
    probe_nav_env11_r1 = probe_nav_env11_r2;
    probe_nav_env12_r1 = probe_nav_env12_r2;
    
    probe_arriv_env10_r1 = probe_arriv_env10_r2;
    probe_arriv_env11_r1 = probe_arriv_env11_r2;
    probe_arriv_env12_r1 = probe_arriv_env12_r2;
    
    
    testo = zeros(size(probe_assigned_env4_r1));
    [~,i] = unique(probe_assigned_env4_r1, 'first');
    testo(i(2)) = 1;
    probe_assigned_env4_r1 = testo;
    
    testo = zeros(size(probe_assigned_env5_r1));
    [~,i] = unique(probe_assigned_env5_r1, 'first');
    testo(i(2)) = 1;
    probe_assigned_env5_r1 = testo;
    
    testo = zeros(size(probe_assigned_env6_r1));
    [~,i] = unique(probe_assigned_env6_r1, 'first');
    testo(i(2)) = 1;
    probe_assigned_env6_r1 = testo;
    
    
    testo = zeros(size(probe_nav_env4_r1));
    [~,i] = unique(probe_nav_env4_r1, 'first');
    testo(i(2)) = 1;
    probe_nav_env4_r1 = testo;
    
    testo = zeros(size(probe_nav_env5_r1));
    [~,i] = unique(probe_nav_env5_r1, 'first');
    testo(i(2)) = 1;
    probe_nav_env5_r1 = testo;
    
    testo = zeros(size(probe_nav_env6_r1));
    [~,i] = unique(probe_nav_env6_r1, 'first');
    testo(i(2)) = 1;
    probe_nav_env6_r1 = testo;
    
    
    testo = zeros(size(probe_arriv_env4_r1));
    [~,i] = unique(probe_arriv_env4_r1, 'first');
    testo(i(2)) = 1;
    probe_arriv_env4_r1 = testo;
    
    testo = zeros(size(probe_arriv_env5_r1));
    [~,i] = unique(probe_arriv_env5_r1, 'first');
    testo(i(2)) = 1;
    probe_arriv_env5_r1 = testo;
    
    testo = zeros(size(probe_arriv_env6_r1));
    [~,i] = unique(probe_arriv_env6_r1, 'first');
    testo(i(2)) = 1;
    probe_arriv_env6_r1 = testo;
    
    
    
end

%within type (probe r1 with probe r2)
probe_assigned_r1_w_probe_assigned_r2_repst = [cm2(logical(probe_assigned_env1_r1),logical(probe_assigned_env1_r2)) cm2(logical(probe_assigned_env2_r1),logical(probe_assigned_env2_r2)) cm2(logical(probe_assigned_env3_r1),logical(probe_assigned_env3_r2)) cm2(logical(probe_assigned_env4_r1),logical(probe_assigned_env4_r2)) cm2(logical(probe_assigned_env5_r1),logical(probe_assigned_env5_r2)) cm2(logical(probe_assigned_env6_r1),logical(probe_assigned_env6_r2)) cm2(logical(probe_assigned_env7_r1),logical(probe_assigned_env7_r2)) cm2(logical(probe_assigned_env8_r1),logical(probe_assigned_env8_r2)) cm2(logical(probe_assigned_env9_r1),logical(probe_assigned_env9_r2)) cm2(logical(probe_assigned_env10_r1),logical(probe_assigned_env10_r2)) cm2(logical(probe_assigned_env11_r1),logical(probe_assigned_env11_r2)) cm2(logical(probe_assigned_env12_r1),logical(probe_assigned_env12_r2))];
probe_assigned_r1_w_probe_assigned_r2_repst_mean = nanmean(probe_assigned_r1_w_probe_assigned_r2_repst(:));

probe_assigned_r1_w_probe_assigned_r2_repstcon = [cm2(logical(probe_assigned_env1_r1),logical(probe_assigned_r2-probe_assigned_env1_r2)) cm2(logical(probe_assigned_env2_r1),logical(probe_assigned_r2-probe_assigned_env2_r2)) cm2(logical(probe_assigned_env3_r1),logical(probe_assigned_r2-probe_assigned_env3_r2)) cm2(logical(probe_assigned_env4_r1),logical(probe_assigned_r2-probe_assigned_env4_r2)) cm2(logical(probe_assigned_env5_r1),logical(probe_assigned_r2-probe_assigned_env5_r2)) cm2(logical(probe_assigned_env6_r1),logical(probe_assigned_r2-probe_assigned_env6_r2)) cm2(logical(probe_assigned_env7_r1),logical(probe_assigned_r2-probe_assigned_env7_r2)) cm2(logical(probe_assigned_env8_r1),logical(probe_assigned_r2-probe_assigned_env8_r2)) cm2(logical(probe_assigned_env9_r1),logical(probe_assigned_r2-probe_assigned_env9_r2)) cm2(logical(probe_assigned_env10_r1),logical(probe_assigned_r2-probe_assigned_env10_r2)) cm2(logical(probe_assigned_env11_r1),logical(probe_assigned_r2-probe_assigned_env11_r2)) cm2(logical(probe_assigned_env12_r1),logical(probe_assigned_r2-probe_assigned_env12_r2))];
probe_assigned_r1_w_probe_assigned_r2_repstcon_mean = nanmean(probe_assigned_r1_w_probe_assigned_r2_repstcon(:));

probe_assigned_r2_w_probe_assigned_r1_repstcon = [cm2(logical(probe_assigned_env1_r2),logical(probe_assigned_r1-probe_assigned_env1_r1)) cm2(logical(probe_assigned_env2_r2),logical(probe_assigned_r1-probe_assigned_env2_r1)) cm2(logical(probe_assigned_env3_r2),logical(probe_assigned_r1-probe_assigned_env3_r1)) cm2(logical(probe_assigned_env4_r2),logical(probe_assigned_r1-probe_assigned_env4_r1)) cm2(logical(probe_assigned_env5_r2),logical(probe_assigned_r1-probe_assigned_env5_r1)) cm2(logical(probe_assigned_env6_r2),logical(probe_assigned_r1-probe_assigned_env6_r1)) cm2(logical(probe_assigned_env7_r2),logical(probe_assigned_r1-probe_assigned_env7_r1)) cm2(logical(probe_assigned_env8_r2),logical(probe_assigned_r1-probe_assigned_env8_r1)) cm2(logical(probe_assigned_env9_r2),logical(probe_assigned_r1-probe_assigned_env9_r1)) cm2(logical(probe_assigned_env10_r2),logical(probe_assigned_r1-probe_assigned_env10_r1)) cm2(logical(probe_assigned_env11_r2),logical(probe_assigned_r1-probe_assigned_env11_r1)) cm2(logical(probe_assigned_env12_r2),logical(probe_assigned_r1-probe_assigned_env12_r1))];
probe_assigned_r2_w_probe_assigned_r1_repstcon_mean = nanmean(probe_assigned_r2_w_probe_assigned_r1_repstcon(:));



%nav
probe_nav_r1_w_probe_nav_r2_repst = [cm2(logical(probe_nav_env1_r1),logical(probe_nav_env1_r2)) cm2(logical(probe_nav_env2_r1),logical(probe_nav_env2_r2)) cm2(logical(probe_nav_env3_r1),logical(probe_nav_env3_r2)) cm2(logical(probe_nav_env4_r1),logical(probe_nav_env4_r2)) cm2(logical(probe_nav_env5_r1),logical(probe_nav_env5_r2)) cm2(logical(probe_nav_env6_r1),logical(probe_nav_env6_r2)) cm2(logical(probe_nav_env7_r1),logical(probe_nav_env7_r2)) cm2(logical(probe_nav_env8_r1),logical(probe_nav_env8_r2)) cm2(logical(probe_nav_env9_r1),logical(probe_nav_env9_r2)) cm2(logical(probe_nav_env10_r1),logical(probe_nav_env10_r2)) cm2(logical(probe_nav_env11_r1),logical(probe_nav_env11_r2)) cm2(logical(probe_nav_env12_r1),logical(probe_nav_env12_r2))];
probe_nav_r1_w_probe_nav_r2_repst_mean = nanmean(probe_nav_r1_w_probe_nav_r2_repst(:));

probe_nav_r1_w_probe_nav_r2_repstcon = [cm2(logical(probe_nav_env1_r1),logical(probe_nav_r2-probe_nav_env1_r2)) cm2(logical(probe_nav_env2_r1),logical(probe_nav_r2-probe_nav_env2_r2)) cm2(logical(probe_nav_env3_r1),logical(probe_nav_r2-probe_nav_env3_r2)) cm2(logical(probe_nav_env4_r1),logical(probe_nav_r2-probe_nav_env4_r2)) cm2(logical(probe_nav_env5_r1),logical(probe_nav_r2-probe_nav_env5_r2)) cm2(logical(probe_nav_env6_r1),logical(probe_nav_r2-probe_nav_env6_r2)) cm2(logical(probe_nav_env7_r1),logical(probe_nav_r2-probe_nav_env7_r2)) cm2(logical(probe_nav_env8_r1),logical(probe_nav_r2-probe_nav_env8_r2)) cm2(logical(probe_nav_env9_r1),logical(probe_nav_r2-probe_nav_env9_r2)) cm2(logical(probe_nav_env10_r1),logical(probe_nav_r2-probe_nav_env10_r2)) cm2(logical(probe_nav_env11_r1),logical(probe_nav_r2-probe_nav_env11_r2)) cm2(logical(probe_nav_env12_r1),logical(probe_nav_r2-probe_nav_env12_r2))];
probe_nav_r1_w_probe_nav_r2_repstcon_mean = nanmean(probe_nav_r1_w_probe_nav_r2_repstcon(:));


probe_nav_r2_w_probe_nav_r1_repstcon = [cm2(logical(probe_nav_env1_r2),logical(probe_nav_r1-probe_nav_env1_r1)) cm2(logical(probe_nav_env2_r2),logical(probe_nav_r1-probe_nav_env2_r1)) cm2(logical(probe_nav_env3_r2),logical(probe_nav_r1-probe_nav_env3_r1)) cm2(logical(probe_nav_env4_r2),logical(probe_nav_r1-probe_nav_env4_r1)) cm2(logical(probe_nav_env5_r2),logical(probe_nav_r1-probe_nav_env5_r1)) cm2(logical(probe_nav_env6_r2),logical(probe_nav_r1-probe_nav_env6_r1)) cm2(logical(probe_nav_env7_r2),logical(probe_nav_r1-probe_nav_env7_r1)) cm2(logical(probe_nav_env8_r2),logical(probe_nav_r1-probe_nav_env8_r1)) cm2(logical(probe_nav_env9_r2),logical(probe_nav_r1-probe_nav_env9_r1)) cm2(logical(probe_nav_env10_r2),logical(probe_nav_r1-probe_nav_env10_r1)) cm2(logical(probe_nav_env11_r2),logical(probe_nav_r1-probe_nav_env11_r1)) cm2(logical(probe_nav_env12_r2),logical(probe_nav_r1-probe_nav_env12_r1))];
probe_nav_r2_w_probe_nav_r1_repstcon_mean = nanmean(probe_nav_r2_w_probe_nav_r1_repstcon(:));


%arriv
probe_arriv_r1_w_probe_arriv_r2_repst = [cm2(logical(probe_arriv_env1_r1),logical(probe_arriv_env1_r2)) cm2(logical(probe_arriv_env2_r1),logical(probe_arriv_env2_r2)) cm2(logical(probe_arriv_env3_r1),logical(probe_arriv_env3_r2)) cm2(logical(probe_arriv_env4_r1),logical(probe_arriv_env4_r2)) cm2(logical(probe_arriv_env5_r1),logical(probe_arriv_env5_r2)) cm2(logical(probe_arriv_env6_r1),logical(probe_arriv_env6_r2)) cm2(logical(probe_arriv_env7_r1),logical(probe_arriv_env7_r2)) cm2(logical(probe_arriv_env8_r1),logical(probe_arriv_env8_r2)) cm2(logical(probe_arriv_env9_r1),logical(probe_arriv_env9_r2)) cm2(logical(probe_arriv_env10_r1),logical(probe_arriv_env10_r2)) cm2(logical(probe_arriv_env11_r1),logical(probe_arriv_env11_r2)) cm2(logical(probe_arriv_env12_r1),logical(probe_arriv_env12_r2))];
probe_arriv_r1_w_probe_arriv_r2_repst_mean = nanmean(probe_arriv_r1_w_probe_arriv_r2_repst(:));

probe_arriv_r1_w_probe_arriv_r2_repstcon = [cm2(logical(probe_arriv_env1_r1),logical(probe_arriv_r2-probe_arriv_env1_r2)) cm2(logical(probe_arriv_env2_r1),logical(probe_arriv_r2-probe_arriv_env2_r2)) cm2(logical(probe_arriv_env3_r1),logical(probe_arriv_r2-probe_arriv_env3_r2)) cm2(logical(probe_arriv_env4_r1),logical(probe_arriv_r2-probe_arriv_env4_r2)) cm2(logical(probe_arriv_env5_r1),logical(probe_arriv_r2-probe_arriv_env5_r2)) cm2(logical(probe_arriv_env6_r1),logical(probe_arriv_r2-probe_arriv_env6_r2)) cm2(logical(probe_arriv_env7_r1),logical(probe_arriv_r2-probe_arriv_env7_r2)) cm2(logical(probe_arriv_env8_r1),logical(probe_arriv_r2-probe_arriv_env8_r2)) cm2(logical(probe_arriv_env9_r1),logical(probe_arriv_r2-probe_arriv_env9_r2)) cm2(logical(probe_arriv_env10_r1),logical(probe_arriv_r2-probe_arriv_env10_r2)) cm2(logical(probe_arriv_env11_r1),logical(probe_arriv_r2-probe_arriv_env11_r2)) cm2(logical(probe_arriv_env12_r1),logical(probe_arriv_r2-probe_arriv_env12_r2))];
probe_arriv_r1_w_probe_arriv_r2_repstcon_mean = nanmean(probe_arriv_r1_w_probe_arriv_r2_repstcon(:));

probe_arriv_r2_w_probe_arriv_r1_repstcon = [cm2(logical(probe_arriv_env1_r2),logical(probe_arriv_r1-probe_arriv_env1_r1)) cm2(logical(probe_arriv_env2_r2),logical(probe_arriv_r1-probe_arriv_env2_r1)) cm2(logical(probe_arriv_env3_r2),logical(probe_arriv_r1-probe_arriv_env3_r1)) cm2(logical(probe_arriv_env4_r2),logical(probe_arriv_r1-probe_arriv_env4_r1)) cm2(logical(probe_arriv_env5_r2),logical(probe_arriv_r1-probe_arriv_env5_r1)) cm2(logical(probe_arriv_env6_r2),logical(probe_arriv_r1-probe_arriv_env6_r1)) cm2(logical(probe_arriv_env7_r2),logical(probe_arriv_r1-probe_arriv_env7_r1)) cm2(logical(probe_arriv_env8_r2),logical(probe_arriv_r1-probe_arriv_env8_r1)) cm2(logical(probe_arriv_env9_r2),logical(probe_arriv_r1-probe_arriv_env9_r1)) cm2(logical(probe_arriv_env10_r2),logical(probe_arriv_r1-probe_arriv_env10_r1)) cm2(logical(probe_arriv_env11_r2),logical(probe_arriv_r1-probe_arriv_env11_r1)) cm2(logical(probe_arriv_env12_r2),logical(probe_arriv_r1-probe_arriv_env12_r1))];
probe_arriv_r2_w_probe_arriv_r1_repstcon_mean = nanmean(probe_arriv_r2_w_probe_arriv_r1_repstcon(:));


%across type (probe with habit), r1 with habit
probe_assigned_r1_w_habit_assigned_r1_repst = [cm2(logical(probe_assigned_env1_r1),logical(habit_assigned_env1_r1)) cm2(logical(probe_assigned_env2_r1),logical(habit_assigned_env2_r1)) cm2(logical(probe_assigned_env3_r1),logical(habit_assigned_env3_r1)) cm2(logical(probe_assigned_env4_r1),logical(habit_assigned_env4_r1)) cm2(logical(probe_assigned_env5_r1),logical(habit_assigned_env5_r1)) cm2(logical(probe_assigned_env6_r1),logical(habit_assigned_env6_r1)) cm2(logical(probe_assigned_env7_r1),logical(habit_assigned_env7_r1)) cm2(logical(probe_assigned_env8_r1),logical(habit_assigned_env8_r1)) cm2(logical(probe_assigned_env9_r1),logical(habit_assigned_env9_r1)) cm2(logical(probe_assigned_env10_r1),logical(habit_assigned_env10_r1)) cm2(logical(probe_assigned_env11_r1),logical(habit_assigned_env11_r1)) cm2(logical(probe_assigned_env12_r1),logical(habit_assigned_env12_r1))];
probe_assigned_r1_w_habit_assigned_r1_repst_mean = nanmean(probe_assigned_r1_w_habit_assigned_r1_repst(:));

probe_assigned_r1_w_habit_assigned_r1_repstcon = [cm2(logical(probe_assigned_env1_r1),logical(habit_assigned_r1-habit_assigned_env1_r1)) cm2(logical(probe_assigned_env2_r1),logical(habit_assigned_r1-habit_assigned_env2_r1)) cm2(logical(probe_assigned_env3_r1),logical(habit_assigned_r1-habit_assigned_env3_r1)) cm2(logical(probe_assigned_env4_r1),logical(habit_assigned_r1-habit_assigned_env4_r1)) cm2(logical(probe_assigned_env5_r1),logical(habit_assigned_r1-habit_assigned_env5_r1)) cm2(logical(probe_assigned_env6_r1),logical(habit_assigned_r1-habit_assigned_env6_r1)) cm2(logical(probe_assigned_env7_r1),logical(habit_assigned_r1-habit_assigned_env7_r1)) cm2(logical(probe_assigned_env8_r1),logical(habit_assigned_r1-habit_assigned_env8_r1)) cm2(logical(probe_assigned_env9_r1),logical(habit_assigned_r1-habit_assigned_env9_r1)) cm2(logical(probe_assigned_env10_r1),logical(habit_assigned_r1-habit_assigned_env10_r1)) cm2(logical(probe_assigned_env11_r1),logical(habit_assigned_r1-habit_assigned_env11_r1)) cm2(logical(probe_assigned_env12_r1),logical(habit_assigned_r1-habit_assigned_env12_r1))];
probe_assigned_r1_w_habit_assigned_r1_repstcon_mean = nanmean(probe_assigned_r1_w_habit_assigned_r1_repstcon(:));


%nav
probe_nav_r1_w_habit_nav_r1_repst = [cm2(logical(probe_nav_env1_r1),logical(habit_nav_env1_r1)) cm2(logical(probe_nav_env2_r1),logical(habit_nav_env2_r1)) cm2(logical(probe_nav_env3_r1),logical(habit_nav_env3_r1)) cm2(logical(probe_nav_env4_r1),logical(habit_nav_env4_r1)) cm2(logical(probe_nav_env5_r1),logical(habit_nav_env5_r1)) cm2(logical(probe_nav_env6_r1),logical(habit_nav_env6_r1)) cm2(logical(probe_nav_env7_r1),logical(habit_nav_env7_r1)) cm2(logical(probe_nav_env8_r1),logical(habit_nav_env8_r1)) cm2(logical(probe_nav_env9_r1),logical(habit_nav_env9_r1)) cm2(logical(probe_nav_env10_r1),logical(habit_nav_env10_r1)) cm2(logical(probe_nav_env11_r1),logical(habit_nav_env11_r1)) cm2(logical(probe_nav_env12_r1),logical(habit_nav_env12_r1))];
probe_nav_r1_w_habit_nav_r1_repst_mean = nanmean(probe_nav_r1_w_habit_nav_r1_repst(:));

probe_nav_r1_w_habit_nav_r1_repstcon = [cm2(logical(probe_nav_env1_r1),logical(habit_nav_r1-habit_nav_env1_r1)) cm2(logical(probe_nav_env2_r1),logical(habit_nav_r1-habit_nav_env2_r1)) cm2(logical(probe_nav_env3_r1),logical(habit_nav_r1-habit_nav_env3_r1)) cm2(logical(probe_nav_env4_r1),logical(habit_nav_r1-habit_nav_env4_r1)) cm2(logical(probe_nav_env5_r1),logical(habit_nav_r1-habit_nav_env5_r1)) cm2(logical(probe_nav_env6_r1),logical(habit_nav_r1-habit_nav_env6_r1)) cm2(logical(probe_nav_env7_r1),logical(habit_nav_r1-habit_nav_env7_r1)) cm2(logical(probe_nav_env8_r1),logical(habit_nav_r1-habit_nav_env8_r1)) cm2(logical(probe_nav_env9_r1),logical(habit_nav_r1-habit_nav_env9_r1)) cm2(logical(probe_nav_env10_r1),logical(habit_nav_r1-habit_nav_env10_r1)) cm2(logical(probe_nav_env11_r1),logical(habit_nav_r1-habit_nav_env11_r1)) cm2(logical(probe_nav_env12_r1),logical(habit_nav_r1-habit_nav_env12_r1))];
probe_nav_r1_w_habit_nav_r1_repstcon_mean = nanmean(probe_nav_r1_w_habit_nav_r1_repstcon(:));


%arrive
probe_arriv_r1_w_habit_arriv_r1_repst = [cm2(logical(probe_arriv_env1_r1),logical(habit_arriv_env1_r1)) cm2(logical(probe_arriv_env2_r1),logical(habit_arriv_env2_r1)) cm2(logical(probe_arriv_env3_r1),logical(habit_arriv_env3_r1)) cm2(logical(probe_arriv_env4_r1),logical(habit_arriv_env4_r1)) cm2(logical(probe_arriv_env5_r1),logical(habit_arriv_env5_r1)) cm2(logical(probe_arriv_env6_r1),logical(habit_arriv_env6_r1)) cm2(logical(probe_arriv_env7_r1),logical(habit_arriv_env7_r1)) cm2(logical(probe_arriv_env8_r1),logical(habit_arriv_env8_r1)) cm2(logical(probe_arriv_env9_r1),logical(habit_arriv_env9_r1)) cm2(logical(probe_arriv_env10_r1),logical(habit_arriv_env10_r1)) cm2(logical(probe_arriv_env11_r1),logical(habit_arriv_env11_r1)) cm2(logical(probe_arriv_env12_r1),logical(habit_arriv_env12_r1))];
probe_arriv_r1_w_habit_arriv_r1_repst_mean = nanmean(probe_arriv_r1_w_habit_arriv_r1_repst(:));

probe_arriv_r1_w_habit_arriv_r1_repstcon = [cm2(logical(probe_arriv_env1_r1),logical(habit_arriv_r1-habit_arriv_env1_r1)) cm2(logical(probe_arriv_env2_r1),logical(habit_arriv_r1-habit_arriv_env2_r1)) cm2(logical(probe_arriv_env3_r1),logical(habit_arriv_r1-habit_arriv_env3_r1)) cm2(logical(probe_arriv_env4_r1),logical(habit_arriv_r1-habit_arriv_env4_r1)) cm2(logical(probe_arriv_env5_r1),logical(habit_arriv_r1-habit_arriv_env5_r1)) cm2(logical(probe_arriv_env6_r1),logical(habit_arriv_r1-habit_arriv_env6_r1)) cm2(logical(probe_arriv_env7_r1),logical(habit_arriv_r1-habit_arriv_env7_r1)) cm2(logical(probe_arriv_env8_r1),logical(habit_arriv_r1-habit_arriv_env8_r1)) cm2(logical(probe_arriv_env9_r1),logical(habit_arriv_r1-habit_arriv_env9_r1)) cm2(logical(probe_arriv_env10_r1),logical(habit_arriv_r1-habit_arriv_env10_r1)) cm2(logical(probe_arriv_env11_r1),logical(habit_arriv_r1-habit_arriv_env11_r1)) cm2(logical(probe_arriv_env12_r1),logical(habit_arriv_r1-habit_arriv_env12_r1))];
probe_arriv_r1_w_habit_arriv_r1_repstcon_mean = nanmean(probe_arriv_r1_w_habit_arriv_r1_repstcon(:));



%across type (probe with habit), r2 with habit
probe_assigned_r2_w_habit_assigned_r1_repst = [cm2(logical(probe_assigned_env1_r2),logical(habit_assigned_env1_r1)) cm2(logical(probe_assigned_env2_r2),logical(habit_assigned_env2_r1)) cm2(logical(probe_assigned_env3_r2),logical(habit_assigned_env3_r1)) cm2(logical(probe_assigned_env4_r2),logical(habit_assigned_env4_r1)) cm2(logical(probe_assigned_env5_r2),logical(habit_assigned_env5_r1)) cm2(logical(probe_assigned_env6_r2),logical(habit_assigned_env6_r1)) cm2(logical(probe_assigned_env7_r2),logical(habit_assigned_env7_r1)) cm2(logical(probe_assigned_env8_r2),logical(habit_assigned_env8_r1)) cm2(logical(probe_assigned_env9_r2),logical(habit_assigned_env9_r1)) cm2(logical(probe_assigned_env10_r2),logical(habit_assigned_env10_r1)) cm2(logical(probe_assigned_env11_r2),logical(habit_assigned_env11_r1)) cm2(logical(probe_assigned_env12_r2),logical(habit_assigned_env12_r1))];
probe_assigned_r2_w_habit_assigned_r1_repst_mean = nanmean(probe_assigned_r2_w_habit_assigned_r1_repst(:));

probe_assigned_r2_w_habit_assigned_r1_repstcon = [cm2(logical(probe_assigned_env1_r2),logical(habit_assigned_r1-habit_assigned_env1_r1)) cm2(logical(probe_assigned_env2_r2),logical(habit_assigned_r1-habit_assigned_env2_r1)) cm2(logical(probe_assigned_env3_r2),logical(habit_assigned_r1-habit_assigned_env3_r1)) cm2(logical(probe_assigned_env4_r2),logical(habit_assigned_r1-habit_assigned_env4_r1)) cm2(logical(probe_assigned_env5_r2),logical(habit_assigned_r1-habit_assigned_env5_r1)) cm2(logical(probe_assigned_env6_r2),logical(habit_assigned_r1-habit_assigned_env6_r1)) cm2(logical(probe_assigned_env7_r2),logical(habit_assigned_r1-habit_assigned_env7_r1)) cm2(logical(probe_assigned_env8_r2),logical(habit_assigned_r1-habit_assigned_env8_r1)) cm2(logical(probe_assigned_env9_r2),logical(habit_assigned_r1-habit_assigned_env9_r1)) cm2(logical(probe_assigned_env10_r2),logical(habit_assigned_r1-habit_assigned_env10_r1)) cm2(logical(probe_assigned_env11_r2),logical(habit_assigned_r1-habit_assigned_env11_r1)) cm2(logical(probe_assigned_env12_r2),logical(habit_assigned_r1-habit_assigned_env12_r1))];
probe_assigned_r2_w_habit_assigned_r1_repstcon_mean = nanmean(probe_assigned_r2_w_habit_assigned_r1_repstcon(:));


%nav
probe_nav_r2_w_habit_nav_r1_repst = [cm2(logical(probe_nav_env1_r2),logical(habit_nav_env1_r1)) cm2(logical(probe_nav_env2_r2),logical(habit_nav_env2_r1)) cm2(logical(probe_nav_env3_r2),logical(habit_nav_env3_r1)) cm2(logical(probe_nav_env4_r2),logical(habit_nav_env4_r1)) cm2(logical(probe_nav_env5_r2),logical(habit_nav_env5_r1)) cm2(logical(probe_nav_env6_r2),logical(habit_nav_env6_r1)) cm2(logical(probe_nav_env7_r2),logical(habit_nav_env7_r1)) cm2(logical(probe_nav_env8_r2),logical(habit_nav_env8_r1)) cm2(logical(probe_nav_env9_r2),logical(habit_nav_env9_r1)) cm2(logical(probe_nav_env10_r2),logical(habit_nav_env10_r1)) cm2(logical(probe_nav_env11_r2),logical(habit_nav_env11_r1)) cm2(logical(probe_nav_env12_r2),logical(habit_nav_env12_r1))];
probe_nav_r2_w_habit_nav_r1_repst_mean = nanmean(probe_nav_r2_w_habit_nav_r1_repst(:));

probe_nav_r2_w_habit_nav_r1_repstcon = [cm2(logical(probe_nav_env1_r2),logical(habit_nav_r1-habit_nav_env1_r1)) cm2(logical(probe_nav_env2_r2),logical(habit_nav_r1-habit_nav_env2_r1)) cm2(logical(probe_nav_env3_r2),logical(habit_nav_r1-habit_nav_env3_r1)) cm2(logical(probe_nav_env4_r2),logical(habit_nav_r1-habit_nav_env4_r1)) cm2(logical(probe_nav_env5_r2),logical(habit_nav_r1-habit_nav_env5_r1)) cm2(logical(probe_nav_env6_r2),logical(habit_nav_r1-habit_nav_env6_r1)) cm2(logical(probe_nav_env7_r2),logical(habit_nav_r1-habit_nav_env7_r1)) cm2(logical(probe_nav_env8_r2),logical(habit_nav_r1-habit_nav_env8_r1)) cm2(logical(probe_nav_env9_r2),logical(habit_nav_r1-habit_nav_env9_r1)) cm2(logical(probe_nav_env10_r2),logical(habit_nav_r1-habit_nav_env10_r1)) cm2(logical(probe_nav_env11_r2),logical(habit_nav_r1-habit_nav_env11_r1)) cm2(logical(probe_nav_env12_r2),logical(habit_nav_r1-habit_nav_env12_r1))];
probe_nav_r2_w_habit_nav_r1_repstcon_mean = nanmean(probe_nav_r2_w_habit_nav_r1_repstcon(:));


%arrive
probe_arriv_r2_w_habit_arriv_r1_repst = [cm2(logical(probe_arriv_env1_r2),logical(habit_arriv_env1_r1)) cm2(logical(probe_arriv_env2_r2),logical(habit_arriv_env2_r1)) cm2(logical(probe_arriv_env3_r2),logical(habit_arriv_env3_r1)) cm2(logical(probe_arriv_env4_r2),logical(habit_arriv_env4_r1)) cm2(logical(probe_arriv_env5_r2),logical(habit_arriv_env5_r1)) cm2(logical(probe_arriv_env6_r2),logical(habit_arriv_env6_r1)) cm2(logical(probe_arriv_env7_r2),logical(habit_arriv_env7_r1)) cm2(logical(probe_arriv_env8_r2),logical(habit_arriv_env8_r1)) cm2(logical(probe_arriv_env9_r2),logical(habit_arriv_env9_r1)) cm2(logical(probe_arriv_env10_r2),logical(habit_arriv_env10_r1)) cm2(logical(probe_arriv_env11_r2),logical(habit_arriv_env11_r1)) cm2(logical(probe_arriv_env12_r2),logical(habit_arriv_env12_r1))];
probe_arriv_r2_w_habit_arriv_r1_repst_mean = nanmean(probe_arriv_r2_w_habit_arriv_r1_repst(:));

probe_arriv_r2_w_habit_arriv_r1_repstcon = [cm2(logical(probe_arriv_env1_r2),logical(habit_arriv_r1-habit_arriv_env1_r1)) cm2(logical(probe_arriv_env2_r2),logical(habit_arriv_r1-habit_arriv_env2_r1)) cm2(logical(probe_arriv_env3_r2),logical(habit_arriv_r1-habit_arriv_env3_r1)) cm2(logical(probe_arriv_env4_r2),logical(habit_arriv_r1-habit_arriv_env4_r1)) cm2(logical(probe_arriv_env5_r2),logical(habit_arriv_r1-habit_arriv_env5_r1)) cm2(logical(probe_arriv_env6_r2),logical(habit_arriv_r1-habit_arriv_env6_r1)) cm2(logical(probe_arriv_env7_r2),logical(habit_arriv_r1-habit_arriv_env7_r1)) cm2(logical(probe_arriv_env8_r2),logical(habit_arriv_r1-habit_arriv_env8_r1)) cm2(logical(probe_arriv_env9_r2),logical(habit_arriv_r1-habit_arriv_env9_r1)) cm2(logical(probe_arriv_env10_r2),logical(habit_arriv_r1-habit_arriv_env10_r1)) cm2(logical(probe_arriv_env11_r2),logical(habit_arriv_r1-habit_arriv_env11_r1)) cm2(logical(probe_arriv_env12_r2),logical(habit_arriv_r1-habit_arriv_env12_r1))];
probe_arriv_r2_w_habit_arriv_r1_repstcon_mean = nanmean(probe_arriv_r2_w_habit_arriv_r1_repstcon(:));




%% reinstatement analysis
%probe r1 with arrive
probe_assigned_r1_w_probe_arriv_r1_repst = [cm2(logical(probe_assigned_env1_r1),logical(probe_arriv_env1_r1)) cm2(logical(probe_assigned_env2_r1),logical(probe_arriv_env2_r1)) cm2(logical(probe_assigned_env3_r1),logical(probe_arriv_env3_r1)) cm2(logical(probe_assigned_env4_r1),logical(probe_arriv_env4_r1)) cm2(logical(probe_assigned_env5_r1),logical(probe_arriv_env5_r1)) cm2(logical(probe_assigned_env6_r1),logical(probe_arriv_env6_r1)) cm2(logical(probe_assigned_env7_r1),logical(probe_arriv_env7_r1)) cm2(logical(probe_assigned_env8_r1),logical(probe_arriv_env8_r1)) cm2(logical(probe_assigned_env9_r1),logical(probe_arriv_env9_r1)) cm2(logical(probe_assigned_env10_r1),logical(probe_arriv_env10_r1)) cm2(logical(probe_assigned_env11_r1),logical(probe_arriv_env11_r1)) cm2(logical(probe_assigned_env12_r1),logical(probe_arriv_env12_r1))];
probe_assigned_r1_w_probe_arriv_r1_repst_mean = nanmean(probe_assigned_r1_w_probe_arriv_r1_repst(:));


probe_assigned_r1_w_habit_arriv_r1_repst = [cm2(logical(probe_assigned_env1_r1),logical(habit_arriv_env1_r1)) cm2(logical(probe_assigned_env2_r1),logical(habit_arriv_env2_r1)) cm2(logical(probe_assigned_env3_r1),logical(habit_arriv_env3_r1)) cm2(logical(probe_assigned_env4_r1),logical(habit_arriv_env4_r1)) cm2(logical(probe_assigned_env5_r1),logical(habit_arriv_env5_r1)) cm2(logical(probe_assigned_env6_r1),logical(habit_arriv_env6_r1)) cm2(logical(probe_assigned_env7_r1),logical(habit_arriv_env7_r1)) cm2(logical(probe_assigned_env8_r1),logical(habit_arriv_env8_r1)) cm2(logical(probe_assigned_env9_r1),logical(habit_arriv_env9_r1)) cm2(logical(probe_assigned_env10_r1),logical(habit_arriv_env10_r1)) cm2(logical(probe_assigned_env11_r1),logical(habit_arriv_env11_r1)) cm2(logical(probe_assigned_env12_r1),logical(habit_arriv_env12_r1))];
probe_assigned_r1_w_habit_arriv_r1_repst_mean = nanmean(probe_assigned_r1_w_habit_arriv_r1_repst(:));


probe_assigned_r1_w_probe_arriv_r1_repstcon = [cm2(logical(probe_assigned_env1_r1),logical(probe_arriv_r1-probe_arriv_env1_r1)) cm2(logical(probe_assigned_env2_r1),logical(probe_arriv_r1-probe_arriv_env2_r1)) cm2(logical(probe_assigned_env3_r1),logical(probe_arriv_r1-probe_arriv_env3_r1)) cm2(logical(probe_assigned_env4_r1),logical(probe_arriv_r1-probe_arriv_env4_r1)) cm2(logical(probe_assigned_env5_r1),logical(probe_arriv_r1-probe_arriv_env5_r1)) cm2(logical(probe_assigned_env6_r1),logical(probe_arriv_r1-probe_arriv_env6_r1)) cm2(logical(probe_assigned_env7_r1),logical(probe_arriv_r1-probe_arriv_env7_r1)) cm2(logical(probe_assigned_env8_r1),logical(probe_arriv_r1-probe_arriv_env8_r1)) cm2(logical(probe_assigned_env9_r1),logical(probe_arriv_r1-probe_arriv_env9_r1)) cm2(logical(probe_assigned_env10_r1),logical(probe_arriv_r1-probe_arriv_env10_r1)) cm2(logical(probe_assigned_env11_r1),logical(probe_arriv_r1-probe_arriv_env11_r1)) cm2(logical(probe_assigned_env12_r1),logical(probe_arriv_r1-probe_arriv_env12_r1))];
probe_assigned_r1_w_probe_arriv_r1_repstcon_mean = nanmean(probe_assigned_r1_w_probe_arriv_r1_repstcon(:));


%nav
probe_nav_r1_w_probe_arriv_r1_repst = [cm2(logical(probe_nav_env1_r1),logical(probe_arriv_env1_r1)) cm2(logical(probe_nav_env2_r1),logical(probe_arriv_env2_r1)) cm2(logical(probe_nav_env3_r1),logical(probe_arriv_env3_r1)) cm2(logical(probe_nav_env4_r1),logical(probe_arriv_env4_r1)) cm2(logical(probe_nav_env5_r1),logical(probe_arriv_env5_r1)) cm2(logical(probe_nav_env6_r1),logical(probe_arriv_env6_r1)) cm2(logical(probe_nav_env7_r1),logical(probe_arriv_env7_r1)) cm2(logical(probe_nav_env8_r1),logical(probe_arriv_env8_r1)) cm2(logical(probe_nav_env9_r1),logical(probe_arriv_env9_r1)) cm2(logical(probe_nav_env10_r1),logical(probe_arriv_env10_r1)) cm2(logical(probe_nav_env11_r1),logical(probe_arriv_env11_r1)) cm2(logical(probe_nav_env12_r1),logical(probe_arriv_env12_r1))];
probe_nav_r1_w_probe_arriv_r1_repst_mean = nanmean(probe_nav_r1_w_probe_arriv_r1_repst(:));


probe_nav_r1_w_habit_arriv_r1_repst = [cm2(logical(probe_nav_env1_r1),logical(habit_arriv_env1_r1)) cm2(logical(probe_nav_env2_r1),logical(habit_arriv_env2_r1)) cm2(logical(probe_nav_env3_r1),logical(habit_arriv_env3_r1)) cm2(logical(probe_nav_env4_r1),logical(habit_arriv_env4_r1)) cm2(logical(probe_nav_env5_r1),logical(habit_arriv_env5_r1)) cm2(logical(probe_nav_env6_r1),logical(habit_arriv_env6_r1)) cm2(logical(probe_nav_env7_r1),logical(habit_arriv_env7_r1)) cm2(logical(probe_nav_env8_r1),logical(habit_arriv_env8_r1)) cm2(logical(probe_nav_env9_r1),logical(habit_arriv_env9_r1)) cm2(logical(probe_nav_env10_r1),logical(habit_arriv_env10_r1)) cm2(logical(probe_nav_env11_r1),logical(habit_arriv_env11_r1)) cm2(logical(probe_nav_env12_r1),logical(habit_arriv_env12_r1))];
probe_nav_r1_w_habit_arriv_r1_repst_mean = nanmean(probe_nav_r1_w_habit_arriv_r1_repst(:));


probe_nav_r1_w_probe_arriv_r1_repstcon = [cm2(logical(probe_nav_env1_r1),logical(probe_arriv_r1-probe_arriv_env1_r1)) cm2(logical(probe_nav_env2_r1),logical(probe_arriv_r1-probe_arriv_env2_r1)) cm2(logical(probe_nav_env3_r1),logical(probe_arriv_r1-probe_arriv_env3_r1)) cm2(logical(probe_nav_env4_r1),logical(probe_arriv_r1-probe_arriv_env4_r1)) cm2(logical(probe_nav_env5_r1),logical(probe_arriv_r1-probe_arriv_env5_r1)) cm2(logical(probe_nav_env6_r1),logical(probe_arriv_r1-probe_arriv_env6_r1)) cm2(logical(probe_nav_env7_r1),logical(probe_arriv_r1-probe_arriv_env7_r1)) cm2(logical(probe_nav_env8_r1),logical(probe_arriv_r1-probe_arriv_env8_r1)) cm2(logical(probe_nav_env9_r1),logical(probe_arriv_r1-probe_arriv_env9_r1)) cm2(logical(probe_nav_env10_r1),logical(probe_arriv_r1-probe_arriv_env10_r1)) cm2(logical(probe_nav_env11_r1),logical(probe_arriv_r1-probe_arriv_env11_r1)) cm2(logical(probe_nav_env12_r1),logical(probe_arriv_r1-probe_arriv_env12_r1))];
probe_nav_r1_w_probe_arriv_r1_repstcon_mean = nanmean(probe_nav_r1_w_probe_arriv_r1_repstcon(:));

%probe r2 with arrive
probe_assigned_r2_w_probe_arriv_r2_repst = [cm2(logical(probe_assigned_env1_r2),logical(probe_arriv_env1_r2)) cm2(logical(probe_assigned_env2_r2),logical(probe_arriv_env2_r2)) cm2(logical(probe_assigned_env3_r2),logical(probe_arriv_env3_r2)) cm2(logical(probe_assigned_env4_r2),logical(probe_arriv_env4_r2)) cm2(logical(probe_assigned_env5_r2),logical(probe_arriv_env5_r2)) cm2(logical(probe_assigned_env6_r2),logical(probe_arriv_env6_r2)) cm2(logical(probe_assigned_env7_r2),logical(probe_arriv_env7_r2)) cm2(logical(probe_assigned_env8_r2),logical(probe_arriv_env8_r2)) cm2(logical(probe_assigned_env9_r2),logical(probe_arriv_env9_r2)) cm2(logical(probe_assigned_env10_r2),logical(probe_arriv_env10_r2)) cm2(logical(probe_assigned_env11_r2),logical(probe_arriv_env11_r2)) cm2(logical(probe_assigned_env12_r2),logical(probe_arriv_env12_r2))];
probe_assigned_r2_w_probe_arriv_r2_repst_mean = nanmean(probe_assigned_r2_w_probe_arriv_r2_repst(:));

probe_assigned_r2_w_probe_arriv_r2_repstcon = [cm2(logical(probe_assigned_env1_r2),logical(probe_arriv_r2-probe_arriv_env1_r2)) cm2(logical(probe_assigned_env2_r2),logical(probe_arriv_r2-probe_arriv_env2_r2)) cm2(logical(probe_assigned_env3_r2),logical(probe_arriv_r2-probe_arriv_env3_r2)) cm2(logical(probe_assigned_env4_r2),logical(probe_arriv_r2-probe_arriv_env4_r2)) cm2(logical(probe_assigned_env5_r2),logical(probe_arriv_r2-probe_arriv_env5_r2)) cm2(logical(probe_assigned_env6_r2),logical(probe_arriv_r2-probe_arriv_env6_r2)) cm2(logical(probe_assigned_env7_r2),logical(probe_arriv_r2-probe_arriv_env7_r2)) cm2(logical(probe_assigned_env8_r2),logical(probe_arriv_r2-probe_arriv_env8_r2)) cm2(logical(probe_assigned_env9_r2),logical(probe_arriv_r2-probe_arriv_env9_r2)) cm2(logical(probe_assigned_env10_r2),logical(probe_arriv_r2-probe_arriv_env10_r2)) cm2(logical(probe_assigned_env11_r2),logical(probe_arriv_r2-probe_arriv_env11_r2)) cm2(logical(probe_assigned_env12_r2),logical(probe_arriv_r2-probe_arriv_env12_r2))];
probe_assigned_r2_w_probe_arriv_r2_repstcon_mean = nanmean(probe_assigned_r2_w_probe_arriv_r2_repstcon(:));


%nav
probe_nav_r2_w_probe_arriv_r2_repst = [cm2(logical(probe_nav_env1_r2),logical(probe_arriv_env1_r2)) cm2(logical(probe_nav_env2_r2),logical(probe_arriv_env2_r2)) cm2(logical(probe_nav_env3_r2),logical(probe_arriv_env3_r2)) cm2(logical(probe_nav_env4_r2),logical(probe_arriv_env4_r2)) cm2(logical(probe_nav_env5_r2),logical(probe_arriv_env5_r2)) cm2(logical(probe_nav_env6_r2),logical(probe_arriv_env6_r2)) cm2(logical(probe_nav_env7_r2),logical(probe_arriv_env7_r2)) cm2(logical(probe_nav_env8_r2),logical(probe_arriv_env8_r2)) cm2(logical(probe_nav_env9_r2),logical(probe_arriv_env9_r2)) cm2(logical(probe_nav_env10_r2),logical(probe_arriv_env10_r2)) cm2(logical(probe_nav_env11_r2),logical(probe_arriv_env11_r2)) cm2(logical(probe_nav_env12_r2),logical(probe_arriv_env12_r2))];
probe_nav_r2_w_probe_arriv_r2_repst_mean = nanmean(probe_nav_r2_w_probe_arriv_r2_repst(:));

probe_nav_r2_w_probe_arriv_r2_repstcon = [cm2(logical(probe_nav_env1_r2),logical(probe_arriv_r2-probe_arriv_env1_r2)) cm2(logical(probe_nav_env2_r2),logical(probe_arriv_r2-probe_arriv_env2_r2)) cm2(logical(probe_nav_env3_r2),logical(probe_arriv_r2-probe_arriv_env3_r2)) cm2(logical(probe_nav_env4_r2),logical(probe_arriv_r2-probe_arriv_env4_r2)) cm2(logical(probe_nav_env5_r2),logical(probe_arriv_r2-probe_arriv_env5_r2)) cm2(logical(probe_nav_env6_r2),logical(probe_arriv_r2-probe_arriv_env6_r2)) cm2(logical(probe_nav_env7_r2),logical(probe_arriv_r2-probe_arriv_env7_r2)) cm2(logical(probe_nav_env8_r2),logical(probe_arriv_r2-probe_arriv_env8_r2)) cm2(logical(probe_nav_env9_r2),logical(probe_arriv_r2-probe_arriv_env9_r2)) cm2(logical(probe_nav_env10_r2),logical(probe_arriv_r2-probe_arriv_env10_r2)) cm2(logical(probe_nav_env11_r2),logical(probe_arriv_r2-probe_arriv_env11_r2)) cm2(logical(probe_nav_env12_r2),logical(probe_arriv_r2-probe_arriv_env12_r2))];
probe_nav_r2_w_probe_arriv_r2_repstcon_mean = nanmean(probe_nav_r2_w_probe_arriv_r2_repstcon(:));



%% univariate control

meanbetas = nanmean(rmat_condensed(:,:));%get mean beta values from the ROI for each regressor


probe_assigned_r1_mbeta = meanbetas(logical(probe_assigned_r1));
probe_assigned_r1_mbeta_mean = nanmean(probe_assigned_r1_mbeta(:));

probe_nav_r1_mbeta = meanbetas(logical(probe_nav_r1));
probe_nav_r1_mbeta_mean = nanmean(probe_nav_r1_mbeta(:));

probe_arriv_r1_mbeta = meanbetas(logical(probe_arriv_r1));
probe_arriv_r1_mbeta_mean = nanmean(probe_arriv_r1_mbeta(:));

habit_assigned_r1_mbeta = meanbetas(logical(habit_assigned_r1));
habit_assigned_r1_mbeta_mean = nanmean(habit_assigned_r1_mbeta(:));

habit_nav_r1_mbeta = meanbetas(logical(habit_nav_r1));
habit_nav_r1_mbeta_mean = nanmean(habit_nav_r1_mbeta(:));

habit_arriv_r1_mbeta = meanbetas(logical(habit_arriv_r1));
habit_arriv_r1_mbeta_mean = nanmean(habit_arriv_r1_mbeta(:));

%second repetition indices
probe_assigned_r2_mbeta = meanbetas(logical(probe_assigned_r2));
probe_assigned_r2_mbeta_mean = nanmean(probe_assigned_r2_mbeta(:));

probe_nav_r2_mbeta = meanbetas(logical(probe_nav_r2));
probe_nav_r2_mbeta_mean = nanmean(probe_nav_r2_mbeta(:));

probe_arriv_r2_mbeta = meanbetas(logical(probe_arriv_r2));
probe_arriv_r2_mbeta_mean = nanmean(probe_arriv_r2_mbeta(:));


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

%% Plots
%fisher's z
%fish = 0.5*log((1+CM2)./(1-CM2))

%plot matrices of interest
% subplot(3,2,1), imagesc(cm2(logical(probe_assigned_r1),logical(probe_assigned_r2)));
% title('p_a1 with p_a2')
% colormap('jet'); % set the colorscheme
% colorbar; % enable colorbar
% subplot(3,2,2), imagesc(cm2(logical(probe_assigned_r1),logical(habit_assigned_r1)));
% title('p_a1 with h_a1')
% colorbar; % enable colorbar
% %
% subplot(3,2,3), imagesc(cm2(logical(probe_assigned_r2),logical(habit_assigned_r1)));
% title('p_a2 with h_a1')
% colorbar; % enable colorbar
% subplot(3,2,4), imagesc(cm(aa_ex_corr2,aa_ex_corr2));
% title('aa_r with aa_r 2nd rep')
% colorbar; % enable colorbar
%
% subplot(3,2,5), imagesc(cm2(ea_ex2,ea_ex));
% title('ea_1st with ea_2nd')
% colorbar; % enable colorbar
% subplot(3,2,6), imagesc(cm2(aa_ex2,aa_ex));
% title('aa_1st with aa_2nd')
% colorbar; % enable colorbar
%
% figure
% %plot hintons of matrices of interest
% subplot(1,2,1), hintonDiagram(ea_r_with_ea_r, 1, 0);
% subplot(1,2,2), hintonDiagram(aa_r_with_aa_r, 1, 0);


%% Save data
%savename = ['Rcorrs_cm' sub '_' mask '_FutwCurr.mat'];
savename = ['Rcorrs_cm' sub '_' mask '_towncorrs.mat'];

%savename = ['Zcorrs_CM' sub '_' mask '.mat'];
%savename = ['corrs_CM' sub '_COUNTERMEASURES_' mask '.mat']
save(savename, 'probe_assigned_r1_w_probe_assigned_r1', 'probe_assigned_r1_w_probe_assigned_r1_mean', 'probe_nav_r1_w_probe_nav_r1', 'probe_nav_r1_w_probe_nav_r1_mean', 'probe_arriv_r1_w_probe_arriv_r1', 'probe_arriv_r1_w_probe_arriv_r1_mean', 'habit_assigned_r1_w_habit_assigned_r1', 'habit_assigned_r1_w_habit_assigned_r1_mean', 'habit_nav_r1_w_habit_nav_r1', 'habit_nav_r1_w_habit_nav_r1_mean', 'habit_arriv_r1_w_habit_arriv_r1', 'habit_arriv_r1_w_habit_arriv_r1_mean', 'probe_assigned_r2_w_probe_assigned_r2', 'probe_assigned_r2_w_probe_assigned_r2_mean', 'probe_nav_r2_w_probe_nav_r2', 'probe_nav_r2_w_probe_nav_r2_mean', 'probe_arriv_r2_w_probe_arriv_r2', 'probe_arriv_r2_w_probe_arriv_r2_mean', 'probe_assigned_r1_w_probe_assigned_r2_repst', 'probe_assigned_r1_w_probe_assigned_r2_repst_mean', 'probe_assigned_r1_w_probe_assigned_r2_repstcon', 'probe_assigned_r1_w_probe_assigned_r2_repstcon_mean', 'probe_nav_r1_w_probe_nav_r2_repst', 'probe_nav_r1_w_probe_nav_r2_repst_mean', 'probe_nav_r1_w_probe_nav_r2_repstcon', 'probe_nav_r1_w_probe_nav_r2_repstcon_mean', 'probe_arriv_r1_w_probe_arriv_r2_repst', 'probe_arriv_r1_w_probe_arriv_r2_repst_mean', 'probe_arriv_r1_w_probe_arriv_r2_repstcon', 'probe_arriv_r1_w_probe_arriv_r2_repstcon_mean', 'probe_assigned_r1_w_habit_assigned_r1_repst', 'probe_assigned_r1_w_habit_assigned_r1_repst_mean', 'probe_assigned_r1_w_habit_assigned_r1_repstcon', 'probe_assigned_r1_w_habit_assigned_r1_repstcon_mean', 'probe_nav_r1_w_habit_nav_r1_repst', 'probe_nav_r1_w_habit_nav_r1_repst_mean', 'probe_nav_r1_w_habit_nav_r1_repstcon', 'probe_nav_r1_w_habit_nav_r1_repstcon_mean', 'probe_arriv_r1_w_habit_arriv_r1_repst', 'probe_arriv_r1_w_habit_arriv_r1_repst_mean', 'probe_arriv_r1_w_habit_arriv_r1_repstcon', 'probe_arriv_r1_w_habit_arriv_r1_repstcon_mean', 'probe_assigned_r2_w_habit_assigned_r1_repst', 'probe_assigned_r2_w_habit_assigned_r1_repst_mean', 'probe_assigned_r2_w_habit_assigned_r1_repstcon', 'probe_assigned_r2_w_habit_assigned_r1_repstcon_mean', 'probe_nav_r2_w_habit_nav_r1_repst', 'probe_nav_r2_w_habit_nav_r1_repst_mean', 'probe_nav_r2_w_habit_nav_r1_repstcon', 'probe_nav_r2_w_habit_nav_r1_repstcon_mean', 'probe_arriv_r2_w_habit_arriv_r1_repst', 'probe_arriv_r2_w_habit_arriv_r1_repst_mean', 'probe_arriv_r2_w_habit_arriv_r1_repstcon', 'probe_arriv_r2_w_habit_arriv_r1_repstcon_mean', 'probe_assigned_r1_w_probe_arriv_r1_repst', 'probe_assigned_r1_w_probe_arriv_r1_repst_mean', 'probe_assigned_r1_w_probe_arriv_r1_repstcon', 'probe_assigned_r1_w_probe_arriv_r1_repstcon_mean', 'probe_nav_r1_w_probe_arriv_r1_repst', 'probe_nav_r1_w_probe_arriv_r1_repst_mean', 'probe_nav_r1_w_probe_arriv_r1_repstcon', 'probe_nav_r1_w_probe_arriv_r1_repstcon_mean', 'probe_assigned_r2_w_probe_arriv_r2_repst', 'probe_assigned_r2_w_probe_arriv_r2_repst_mean', 'probe_assigned_r2_w_probe_arriv_r2_repstcon', 'probe_assigned_r2_w_probe_arriv_r2_repstcon_mean', 'probe_nav_r2_w_probe_arriv_r2_repst', 'probe_nav_r2_w_probe_arriv_r2_repst_mean', 'probe_nav_r2_w_probe_arriv_r2_repstcon', 'probe_nav_r2_w_probe_arriv_r2_repstcon_mean', 'probe_assigned_r1_mbeta', 'probe_assigned_r1_mbeta_mean', 'probe_nav_r1_mbeta', 'probe_nav_r1_mbeta_mean', 'probe_arriv_r1_mbeta', 'probe_arriv_r1_mbeta_mean', 'habit_assigned_r1_mbeta', 'habit_assigned_r1_mbeta_mean', 'habit_nav_r1_mbeta', 'habit_nav_r1_mbeta_mean', 'habit_arriv_r1_mbeta', 'habit_arriv_r1_mbeta_mean', 'probe_assigned_r2_mbeta', 'probe_assigned_r2_mbeta_mean', 'probe_nav_r2_mbeta', 'probe_nav_r2_mbeta_mean', 'probe_arriv_r2_mbeta', 'probe_arriv_r2_mbeta_mean', 'probe_assigned_r2_w_probe_assigned_r1_repstcon', 'probe_assigned_r2_w_probe_assigned_r1_repstcon_mean', 'probe_nav_r2_w_probe_nav_r1_repstcon', 'probe_nav_r2_w_probe_nav_r1_repstcon_mean', 'probe_arriv_r2_w_probe_arriv_r1_repstcon', 'probe_arriv_r2_w_probe_arriv_r1_repstcon_mean', 'probe_assigned_r1_w_habit_arriv_r1_repst', 'probe_assigned_r1_w_habit_arriv_r1_repst_mean', 'probe_nav_r1_w_habit_arriv_r1_repst', 'probe_nav_r1_w_habit_arriv_r1_repst_mean');





end
