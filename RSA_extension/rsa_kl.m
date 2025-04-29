function [] = rsa_kl(Sub, Mask, TRsperRun)
% code for RSA analysis. Kyoungeun's study

% example call with 'mvpa_sample_data' - rsa_kl({'201'}, 'BA2_Amyg', [1 1 1
% 1 1 1]) %for betas the TRsPerRun isn't used right now, but an entry is
% still needed for the function to run

%For alternate code demo purposes, this is built around both 4D and 3D
%image file types. The use case for 4D includes scenarios like "raw" BOLD data or residual time-series

%warning: as of 1/3/2018; unconcatenated 4D analysis has not been debugged.
%at this time only use concatenated 3D

%% functional data parameters
S.TR = 1; % IF using betas and your "onset times" are not times in seconds but simply rank order numbers, set = 1
theseTRWeights = [0 0 0.2 0.2 0.2 0.2 0.2]; %[0 0.25 0.5 0.25 0]; %ignore for betas
weights_str = mat2str(theseTRWeights);% assign values to string for custom output file naming %ignore for betas

% what extension do your BOLD images have (nii, img)
funcftype = '.nii';

% are we working with BOLDs/timeseries ('raw') or with beta maps ('betas')?
S.inputformat = 'betas';

%specify preprocessing level of BOLDs (if S.inputformat = 'raw')
preproc_lvl = ''; % in SPM's conventions, these would be: 'a' for slice-time-only, 'u' for realigned-only, 'ua' for realign+unwarped, 'swua' for smoothed, normalized, and... you get the picture. Modify as needed if you changed SPM's prefix append defaults
boldnames = [preproc_lvl 'run']; %name of image files with preprocessing level prefix
% if your time-series is instead residuals (e.g., you regressed out
% motion-related signal using a GLM - so the data aren't *truly* "raw"...
%boldnames = ['ResI']; %for residual time-series

%specify beta filename unique identifiers (often simply 'beta'; if S.inputformat = 'betas')
betanames = 'beta'; %name shared across image files to help ensure only those are read. For LSS, often code renames betas according to conditions and events, so this could be set to read in only a specific condition type, or to load in all events ('event' - all patterns as you normally would)
LStype = 'Multievent'; %Multievent, LSS, or LSA will divert code to different accordingly named beta folders

ImgDims = 3; %if working with timeseries, it is highly recommended that you use 4D nifti files ('4'). If you have split them out into TR-by-TR, or are working with betas, enter '3'


%% analysis flags
cutoffonsets = 0; % 1 = yes, cutoff the conditions in the names file you want to analyze (cutoff point specified below)
runs_concat = 1; %1 = typical SPM analysis; will have continuous onsets concatenated across runs. 0 = But you might not have bothered creating such a file, or in SST case we are using files from FSL. In this case, the onsets are assumed to "reset" for each run ('raw' unconcatenated onsets)
use_exist_workspace = 0; %1=yes. Load existing pattern workspaces and onsets files. Saves time, but turn off if want to manually re-do pattern extraction and generation.

gen_onsetsTR = 1; %1=yes. Typically, you'll use an onsets.mat file with tr-by-tr onsets and names (as used for a beta-series). But if you only have a traditional GLM model with one name for multiple onsets, setting this flag to 1 will auto-populate unique but related names (e.g., Face_1; Face_2...)
sortmatbycat = 0; %1=yes. Rearrange names and patterns according to alphabetical order (helps with visualization of corrmats)
runhpfilt = 1;%1=yes. Standard %ignore for betas
runzscore = 1;%1=yes. Standard but controversial preprocessing. %ignore for betas

%optional - toss similarity values from correlation matrices below a certain number as well (e.g., maybe a zero
%is equivalent to a NaN for your study. Say you are using a mask on "raw"
%BOLD data that does nothing to account for signal drop-out and extra-brain
%voxels are zeros instead of NaNs - here we can fix that
threshpats = 0; % 1 = YES

% what if our data of interest don't start at "run_01" (e.g., your localizer is the 10th and 11th scan session of the day)?
% specify scan run range
% NOTE: currently used only if runs_concat = 0
runstart = 1;
runend = 2;
realrnums = runstart:1:runend;

%Subject ID/number
par.substr = ['S' Sub{1}];
S.subj_id = par.substr;

% What is your study's name / study folder name (subject folders should be
% subdirectories of this folder
S.exp_name = 'kl_valence';
study_prefix = ''; %deprecated/unused?

% What is the name of your .mat file which specificies the condition labels
% and onset times for when the events of those conditions occurred?
S.onsets_filename = [Sub{1} '_rsa_SOTs_ValMem_dur2'];

% what if we don't want to use all of the conditions specified in our condition labels filed?
% ditch unneeded indices. NOTE: this is custom for class example
if cutoffonsets == 1
conditions_range = 1:7; %specify a vector to filter names/onsets/durations; class example 1:7 means use all conditions in the .mat file from #1 to 7 -> any conditions beyond that we won't bother finding patterns for them in the time-series
end

mask = Mask;

% flags for specialized analyses to run
run_mds = 0; % run a "multimensional scaling" analysis?
run_hrclust = 1; % run a hierarchical clustering analysis?
run_simmodelfit = 0; % build "model" similarity matrices and test how well these fit observed

%% Directories
S.expt_dir = ['/Users/giova/Documents/work/My4803Folder/mvpa_sample_data/' S.exp_name '/'];%study location

par.subdir =[S.expt_dir S.subj_id];%subject location

par.funcdir =[par.subdir '/bolds/'];%subfolder for 'raw' BOLD data. Assumes BOLDs are stored in subfolders labeled 'run_01', etc)

S.workspace_dir = [par.subdir '/mvpa_workspace'];%temporary files workspace - code creates this if it doesn't exist

%model file directory (onsets.mat and betas in here) - this is still used
%when working with raw data. We must have some way to tell the classifier
%which images correspond to which classes
if strcmp(S.inputformat, 'raw')
    S.mvpa_dir = [S.expt_dir S.subj_id '/1stLevel_ValMem_s3_dur2/'];
elseif strcmp(S.inputformat, 'betas')
    S.mvpa_dir = [S.expt_dir S.subj_id '/1stLevel_ValMem_s3_dur2/'];
    if strcmp(LStype,'LSS')
        S.beta_dir = [S.expt_dir S.subj_id '/1stLevel_ValMem_s3_dur2/LSS/'];
    elseif strcmp(LStype,'LSA')
        S.beta_dir = [S.expt_dir S.subj_id '/1stLevel_ValMem_s3_dur2/LSA/'];
    elseif strcmp(LStype,'Multievent')
        S.beta_dir = [S.expt_dir S.subj_id '/1stLevel_ValMem_s3_dur2/'];
    end
end

%ROI masks (could be whole-brain mask, but the code wants a mask file
S.anat_dir = [S.expt_dir S.subj_id '/Mask'];
maskname=[S.anat_dir '/' mask '.nii'];

S.group_mvpa_dir = [S.expt_dir 'RSA_output_files'];%results .mat files are spit out in here

if ~exist(S.group_mvpa_dir)
    mkdir(S.group_mvpa_dir)
end

if ~exist([S.mvpa_dir '/RSA_data/'])
    mkdir([S.mvpa_dir '/RSA_data/'])
end


%% extract and process patterns based on flags above
if runs_concat == 1
    %load onsets
    load([S.mvpa_dir S.onsets_filename]);

    if cutoffonsets == 1
    %ditch unneeded indices. NOTE: this is custom for class example
    names = names(conditions_range);
    onsets = onsets(conditions_range);
    durations = durations(conditions_range);
    end

    rmat_condensed = [];
    onsets_TRs = [];
    names_TRs = [];%filled in if gen_onsetsTR == 1
    runsel_TRs = [];

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
            if strcmp(S.inputformat, 'raw')
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

                run_sel = cell2mat(run_sel);
                %% optional preprocessing
                figure;
                % hp filter the data - recommended
                if runhpfilt == 1
                    rmat_t = hp_filter(rmat_t,run_sel',100,2)';%2=2s TR

                    %plot for exploration
                    cm_t = corr(rmat_t);
                    subplot(2,1,1), imagesc(cm_t);
                    title('hp_filt corrmat')
                    colormap('jet'); % set the colorscheme
                    colorbar;
                    caxis([-1 1]);
                end

                % zscore within runs
                if runzscore == 1

                    pat_t = [];

                    for r = 1:length(TRsperRun) %for each run
                        activepats = rmat_t(:,logical(run_sel==r));%filter patterns to current run
                        if size(activepats,2)~=TRsperRun(r)
                            error('Your pattern count doesnt match current run length');
                        end

                        pat_t = [pat_t zscore_mvpa(activepats,2)];%2 = z-score within rows (within-voxels)
                    end

                    rmat_t = pat_t;

                    %plot for exploration
                    cm_t2 = corr(rmat_t);
                    subplot(2,1,2), imagesc(cm_t2);
                    title('hp_filt_z corrmat')
                    colormap('jet'); % set the colorscheme
                    colorbar;
                    caxis([-1 1]);
                end

                %% now compute a mean pattern for __ TRs surrounding the onset

                for n = 1:length(names)


                    time_idx = floor(onsets{n}/S.TR) + 1;%convert onsets to TRs

                    onsets_TR{n} = time_idx; %sort(horzcat(onsets_t{n}, onsets_t{n}{idxThisCond}(enoughTRs_h)));%put the onsets for cond{i} into an array,

                    theseTRWeights2 = theseTRWeights;

                    %create weighted mean pattern for that onset
                    for nvoxr = 1:size(rmat_t,1)
                        %tempvals = [];
                        for tidx = 1:length(time_idx)
                            tempvals = theseTRWeights2.*rmat_t(nvoxr,time_idx(tidx):time_idx(tidx)+(length(theseTRWeights2)-1));
                            condmat_condensed_t{nvoxr,time_idx(tidx)} = squeeze(sum(tempvals));
                        end
                    end
                    %rmat_condensed = [rmat_condensed condmat_condensed_t];

                    if gen_onsetsTR == 1 %if we need to split our names out by individual events
                        for tidx = 1:length(time_idx)
                            tempnames{time_idx(tidx)} = [names{n} '_' num2str(tidx)];
                            temprunsel{time_idx(tidx)} = run_sel(time_idx(tidx));
                        end
                    end

                    %runsel_TRs = [runsel_TRs temprunsel];
                    %names_TRs = [names_TRs tempnames];
                end

                filt = any(~cellfun('isempty', condmat_condensed_t),1);
                rmat_condensed = cell2mat(condmat_condensed_t(:,filt));
                runsel_TRs = temprunsel(:,filt);
                names_TRs = tempnames(:,filt);


            elseif strcmp(S.inputformat, 'betas')
                %note, this analysis assumes betas are grouped by
                %condition, and in the order of the conditions in your
                %onsets file. If they are in a different order, such as the
                %actual correct temporal order of events in your task
                %(e.g., conditions interleaved) then some of the code below
                %that reorders the betas will need to be updated
                allrawfilenames{1,1} = dir(fullfile(S.beta_dir, ['/*' betanames '*.nii']));%

                for idx = 1:length(allrawfilenames{1,1});
                    raw_filenames{idx,1} = [S.beta_dir allrawfilenames{1,1}(idx).name];
                end

                imgslength = length(raw_filenames);

                %% iterate through 3D frames to extract all patterns
                for i=1:imgslength
                    betamaps{i} = raw_filenames{i};%[tmp.name ',' num2str(i)];

                    [b,r] = MAP_getROI(maskname, betamaps{i}, 'vox', 0, '');
                    bmat_t(:,i) = b{1}; % returns all voxels, whether or not they have NaNs
                    rmat_t(:,i) = r; % returns voxels excluding NaNs
                    %meanbetas_t = nanmean(bmat_t(:,:));%get mean beta values from the ROI for each regressor

                end

                rmat_condensed = rmat_t;
                %% now pull out original condition names, and resort to original trial order to help index and analyze pattern data and compare with raw bold patterns. Requires input of TR numbers for study even though we are working with precomputed beta maps
                time_idx_sortvec = [];%we'll fill in the TR times for events, and use this to resort the betas into their true event order (LSS and LSA code typically create betas grouped by condition/out of temporal order)
                for n = 1:length(names)
                    time_idx = floor(onsets{n}/S.TR) + 1;%convert onsets to TRs
                    time_idx_sortvec = [time_idx_sortvec time_idx];

                    if gen_onsetsTR == 1 %if we need to split our names out by individual events
                        for tidx = 1:length(time_idx)
                            tempnames{time_idx(tidx)} = [names{n} '_' num2str(tidx)];
                            % temprunsel{time_idx(tidx)} = run_sel(time_idx(tidx));
                        end
                    end
                end

                filt = any(~cellfun('isempty', tempnames),1);
                %rmat_condensed = cell2mat(condmat_condensed_t(:,filt));
                %runsel_TRs = temprunsel(:,filt);
                names_TRs = tempnames(:,filt);

                [~,sortnames]=sort(time_idx_sortvec);%calculate sorting vector for original temporal order for the task
                rmat_condensed = rmat_condensed(:,sortnames);%let's sort the patterns into their original temporal order for the task
                %note for Kyoungeun - how do we know this "works" for your
                %data? Well - your original "names" are in the order you
                %want, and you gave "onsets" numbers asking for them to be
                %in that order. If you look a the time_idx_sortvec at line
                %353, it is successfully recomputing that order from the
                %orginal onsets -> so the code we have to handle other
                %scenarios is working ok for your scenario now too
                %(specifically, it is preserving the desired order of the
                %patterns you wanted)
            end
        end


        %% write out pattern info
        savename1 = [S.mvpa_dir '/RSA_data/' S.subj_id '_' mask '_condensedpats.mat'];
        save(savename1,'rmat_condensed');

        savename2 = [S.mvpa_dir '/RSA_data/' S.subj_id '_onsets_expanded.mat'];
        save(savename2,'names','onsets','names_TRs','durations');
    end

else %if runs are NOT concatenated %-------in debugging stage as of 1/3/2018
    tempmat = []; %initialize an empty matrix you'll append data from the runs to
    tnames = [];

    for rnum = 1:length(TRsperRun) %4 for 4 probe runs 3-6

        run = num2str(realrnums(rnum), '%02.f'); %+2 added to start at run03 instead of run01


        path = [par.funcdir '/run' run '/'];

        cd(path);

        %load onsets for the run
        load([S.onsets_filename '_' run]); %[S.mvpa_dir S.onsets_filename]);

        if ImgDims == 4 %in development

            a = [];

        elseif ImgDims == 3

            allrawfilenames{rnum,1} = dir(fullfile(['./' boldnames '*.nii']));

            %if 3D images (not recommended) check if the count matches that
            %specified for other stages of the process
            if ImgDims == 3
                if length(allrawfilenames{rnum})~=TRsperRun(rnum);
                    error('your specified run length does not match 3D file count')
                end
            end

            %             for idxf = 1:length(allrawfilenames{rnum})
            %                 allrawfilepaths{rnum,1}{rnum,1} = runfolds(rnum).name;
            %             end
            %             x = [allrawfilenames{1}(1).folder '/' allrawfilenames{1}(1).name]
            %allrawfilenames = vertcat(allrawfilenames{:});
            %             allrawfilepaths = vertcat(allrawfilepaths{:});

            raw_filenames = {};
            for idx = 1:length(allrawfilenames{rnum})%length(allrawfilenames);
                raw_filenames{idx,1} = [allrawfilenames{rnum}(idx).folder '/' allrawfilenames{rnum}(idx).name]
            end

            imgslength = length(raw_filenames);

            %% iterate through 3D frames to extract all patterns
            bmat_t = [];
            rmat_t = [];
            for i=1:imgslength
                betamaps{i} = raw_filenames{i};%[tmp.name ',' num2str(i)];

                [b,r] = MAP_getROI(maskname, betamaps{i}, 'vox', 0, '');
                bmat_t(:,i) = b{1}; % returns all voxels, whether or not they have NaNs
                rmat_t(:,i) = r; % returns voxels excluding NaNs
                %meanbetas_t = nanmean(bmat_t(:,:));%get mean beta values from the ROI for each regressor
            end

            %% optional preprocessing
            figure;
            % hp filter the data - recommended
            run_sel = ones(1,imgslength); % because we are doing one run at a time, run_sel is nonexistant at this time
            if runhpfilt == 1
                rmat_t = hp_filter(rmat_t,run_sel',100,2)';%2=2s TR

                %plot for exploration
                cm_t = corr(rmat_t);
                subplot(2,1,1), imagesc(cm_t);
                title('hp_filt corrmat')
                colormap('jet'); % set the colorscheme
                colorbar;
                caxis([-1 1]);
            end

            % zscore within runs
            if runzscore == 1

                pat_t = [];

                %    for r = 1:length(TRsperRun) %for each run
                activepats = rmat_t;%(:,logical(run_sel==r));%filter patterns to current run
                if size(activepats,2)~=TRsperRun(rnum)
                    error('Your pattern count doesnt match current run length');
                end

                pat_t = [pat_t zscore_mvpa(activepats,2)];%2 = z-score within rows (within-voxels)
                %    end

                rmat_t = pat_t;

                %plot for exploration
                cm_t2 = corr(rmat_t);
                subplot(2,1,2), imagesc(cm_t2);
                title('hp_filt_z corrmat')
                colormap('jet'); % set the colorscheme
                colorbar;
                caxis([-1 1]);
            end

            %% now compute a mean pattern for __ TRs surrounding the onset

            for n = 1:length(names)


                time_idx = floor(onsets{n}/S.TR) + 1;%convert onsets to TRs

                onsets_TR{n} = time_idx; %sort(horzcat(onsets_t{n}, onsets_t{n}{idxThisCond}(enoughTRs_h)));%put the onsets for cond{i} into an array,

                theseTRWeights2 = theseTRWeights;

                %create weighted mean pattern for that onset
                for nvoxr = 1:size(rmat_t,1)
                    %tempvals = [];
                    for tidx = 1:length(time_idx)
                        tempvals = theseTRWeights2.*rmat_t(nvoxr,time_idx(tidx):time_idx(tidx)+(length(theseTRWeights2)-1));
                        condmat_condensed_t{nvoxr,time_idx(tidx)} = squeeze(sum(tempvals));
                    end
                end
                %rmat_condensed = [rmat_condensed condmat_condensed_t];

                if gen_onsetsTR == 1 %if we need to split our names out by individual events
                    for tidx = 1:length(time_idx)
                        tempnames{time_idx(tidx)} = [names{n} '_' run '_' num2str(tidx)];
                        temprunsel{time_idx(tidx)} = run_sel(time_idx(tidx));
                    end
                end

                %runsel_TRs = [runsel_TRs temprunsel];
                %names_TRs = [names_TRs tempnames];
            end

            filt = any(~cellfun('isempty', condmat_condensed_t),1);
            rmat_condensed = cell2mat(condmat_condensed_t(:,filt));
            runsel_TRs = temprunsel(:,filt);
            names_TRs = tempnames(:,filt);
            tempnames=[];%clear for next loop
            temprunsel=[];
            condmat_condensed_t=[];
        end

        tempmat = [tempmat rmat_condensed]; % now concatenate with other runs' data
        tnames = [tnames names_TRs];
    end

    rmat_condensed = tempmat;
    names_TRs = tnames;

end

%% for visualization and statistical analyses, create indices for patterns of interest in corr matrices

if gen_onsetsTR == 1
    names = names_TRs;
end

%optional first rearrange data according to category to help with visualization
if sortmatbycat == 1
    [names_TRs y] = natsort(names);%requires natsort in path
    names = names_TRs;
    rmat_condensed = rmat_condensed(:,y);
end

% first index different valence categories (collapsing across subs memory)
for n = 1:length(names)
    if contains(names{n},{'Neg'})%
        Neg_idx(n)=1;
        Pos_idx(n)=0;
        Neu_idx(n)=0;
    elseif contains(names{n},{'Pos'})%
        Neg_idx(n)=0;
        Pos_idx(n)=1;
        Neu_idx(n)=0;
    elseif contains(names{n},{'Neu'})%
        Neg_idx(n)=0;
        Pos_idx(n)=0;
        Neu_idx(n)=1;
    end
end

% Now index different sub memory categories (collapsing across valence)
for n = 1:length(names)
    if contains(names{n},{'Hit'})%
        Hit_idx(n)=1;
        Miss_idx(n)=0;
    elseif contains(names{n},{'Miss'})%
        Hit_idx(n)=0;
        Miss_idx(n)=1;
    end
end


%example for how you could compute emotional idx (vs neutral)
for n = 1:length(names)
    if contains(names{n},{'Neg','Pos'})%strfind(names{n},'EA') | strfind(names{n}, 'AA')
        Emo_idx(n) = 1;
    else
        Emo_idx(n) = 0;
    end
end

%~~~~~~~Find intersections of the instances to examine more specific results
%scrambled faces
Pos_Hits = Pos_idx.*Hit_idx;% ".*" syntax means multiply the corresponding elements of each matrix or vector
Pos_Misses = Pos_idx.*Miss_idx;


%% create correlation matrices

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

%% global similarity measures - what is the representational match of a category to other exemplars of itself, or to those of a superordinate category (e.g., face similarity with ALL visual stimuli [vs. auditory stimuli])
res.Pos_w_Pos = cm(logical(Pos_idx),logical(Pos_idx));
res.Pos_w_Pos_mean = nanmean(res.Pos_w_Pos(:));

res.Pos_w_Neg = cm2(logical(Pos_idx),logical(Neg_idx));
res.Pos_w_Neg_mean = nanmean(res.Pos_w_Neg(:));


%% stability "item" specific effects (e.g., how similar is a particular face to itself across repetitions, or a planning period in a particular environment across repetitions?)
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

%% univariate controls - could set up analyses computing signal level ACROSS voxels in ROI, to explore how RSA results might correlate with univariate amplitudes

meanbetas = nanmean(rmat_condensed(:,:));%get mean activity values from the ROI for each regressor (only betas if patterns are actually betas...)


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
figure;
subplot(2,2,1), imagesc(cm);
title('Within-cond corrmat');
colormap('jet'); % set the colorscheme
caxis([-1 1]);
colorbar; % enable colorbar

subplot(2,2,2), imagesc(cm2);
title('Overall corrmat');
colormap('jet'); % set the colorscheme
caxis([-1 1]);
colorbar; % enable colorbar

subplot(2,2,3), imagesc(cm(logical(Pos_idx),logical(Pos_idx)));
title('Pos with Pos across runs');
colormap('jet'); % set the colorscheme
caxis([-1 1]);
colorbar; % enable colorbar

subplot(2,2,4), imagesc(cm2(logical(Pos_idx),logical(Neg_idx)));
title('Pos with Neg across runs');
colormap('jet'); % set the colorscheme
caxis([-1 1]);
colorbar; % enable colorbar

% save plot
plot_savename = [S.group_mvpa_dir '/Rcorrs_' S.subj_id '_' mask '_' weights_str '_' S.exp_name '_corrmats.png'];
saveas(gcf,plot_savename);

%% examine correlation structure between two specific classes
testinds = Pos_Hits+Pos_Misses;
custlbls = names(logical(testinds));
rmat_condensed_short = rmat_condensed(:,logical(testinds));
cm_2c = corr(rmat_condensed_short);
figure;
imagesc(cm_2c);
colormap('jet');
caxis([-1 1]);
colorbar;
set(gca, 'XTicklabel', custlbls, 'XTick', [1:length(custlbls)]);
set(gca, 'YTicklabel', custlbls, 'YTick', [1:length(custlbls)]);
res.cm_2c = cm_2c;

% is there a rhyme or reason to the similarity accorting to category?
schemaball(cm2,names_TRs);

%% clustering and multidimensional scaling
cm3 = corr(rmat_condensed);%generate a corrmat without NaNs
cm4 = 1-cm3;%generate DISsimilarity matrix - some clustering algorithms assume distance, not proximity, is the significance of the numbers in the matrix

res.cm4=cm4;%store in res structure for posterity

%visualize  a DISsimilarity matrix
figure;
subplot(2,1,1),imagesc(cm4);
colormap('jet');
colorbar;
set(gca, 'YTicklabel', names, 'YTick', [1:length(names)]);



%% MDS analysis
if run_mds == 1

    %let's see what's going on just within the Positive condition
    [Y,eigvals] = cmdscale(cm4(logical(Pos_idx),logical(Pos_idx)));
    figure;
    subplot(1,2,1), plot(1:length(eigvals),eigvals,'bo-');
    line([1,length(eigvals)],[0 0],'LineStyle',':','XLimInclude','off',...
        'Color',[.7 .7 .7])
    axis([1,length(eigvals),min(eigvals),max(eigvals)*1.1]);
    xlabel('Eigenvalue number');
    ylabel('Eigenvalue');

    labels = names_TRs(logical(Pos_idx));
    subplot(1,2,2), plot(Y(:,1),Y(:,2),'bx');
    axis(max(max(abs(Y))) * [-1.1,1.1,-1.1,1.1]); axis('square');
    text(Y(:,1),Y(:,2),labels,'HorizontalAlignment','left');
    line([-1,1],[0 0],'XLimInclude','off','Color',[.7 .7 .7])
    line([0 0],[-1,1],'YLimInclude','off','Color',[.7 .7 .7])

end


%% hierarchical clustering analysis

if run_hrclust == 1; % run a hierarchical clustering analysis?

    %create vector of distances between instances in the cm
    distm = pdist(cm3,'correlation');% tell matlab metric is pearson r
    Z1 = linkage(distm,'average');%compute dendrogram, using average distance within clusters for agglomeration

    %color code select original classes to help you evaluate clustering?
    colorcodeZ1 = 1; %1 = yes, let's color code, 0 = no; you may turn this off if working out the color coding is unnecessary
    if colorcodeZ1 == 1
        xz(1:size(names'),1) = {[0 0 0]};%create dendrogram condition colors using RGB code (default = [0 0 0], black)
        xz(logical(AA_intact)) = {[1 0 0]};
        xz(logical(EA_intact)) = {[0 1 0]};
        xz(logical(AA_scrambled)) = {[1 0.5 0]};
        xz(logical(EA_scrambled)) = {[0 1 .5]};
        xz(logical(Scene_idx)) = {[0 0 1]};
        xz(logical(Obj_idx)) = {[1 0 1]};
        %h = cell2mat(userOptions.conditionColours)
        userOptions.conditionColours = cell2mat(xz);
    end

    % compute dendrogram. This example code labels each branch termination
    % according to the name of that event in the 'names' model file you
    % provided
    subplot(2,1,2),[H_ignore T_ignore labelReordering] = dendrogram(Z1,0,'labels',names,'Orientation','left');%display dendrogram. 0 = show all items in correlation structure (default would limit to 30 clusters)

    if colorcodeZ1 == 1
        color_t = xz(labelReordering);
        hold on;
        x = xlim(gca);
        for condition = 1:size(cm3,1)%size(squareRDM(cm4), 1)
            plot(x(1), condition, 'o', 'MarkerFaceColor', color_t{condition, :}, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
        end%for:condition
        hold off;
    end

    % save plot
    plot_savename = [S.group_mvpa_dir '/Rcorrs_' S.subj_id '_' mask '_' weights_str '_' S.exp_name '_hierarchicalclustering.png'];
    saveas(gcf,plot_savename);

    % ~~~~ that level of visualization can be too hard to make sense of

    % Let's break it down with a couple of examples to simplify interpretation

    % First, let's explore which classes belong to clusters at level __ in the dendrogram
    clustlvl = 3;%define level of dendrogram to inspect. E.g., 3 means the level from the top where there are 3 clusters
    cl_content{clustlvl} = cluster(Z1,'maxclust',clustlvl); %what are the items in cluster level 'cutoff'?
    names_cl1 = names(cl_content{clustlvl}==1); % which classes are in cluster 1 at this level of the dendrogram?

    % Second, let's try clustering just on a smaller subset of the data.
    % Let's use the same example from above using cm_2c and custlbls

    %create vector of distances between instances in the cm
    distm_2 = pdist(cm_2c,'correlation');% tell matlab metric is pearson r
    Z2 = linkage(distm_2,'average');%compute dendrogram, using average distance within clusters for agglomeration
    % compute dendrogram.
    figure;
    subplot(1,1,1),[H_ignore T_ignore labelReordering] = dendrogram(Z2,0,'labels',custlbls,'Orientation','left');%display dendrogram. 0 = show all items in correlation structure (default would limit to 30 clusters)
    %color code select original classes to help you evaluate clustering?
    colorcodeZ2 = 1; %1 = yes, let's color code, 0 = no; you may turn this off if working out the color coding is unnecessary
    if colorcodeZ2 == 1
        xz2(1:size(custlbls'),1) = {[0 0 0]};%create dendrogram condition colors using RGB code (default = [0 0 0], black)
        xz2(logical(EA_intact(logical(testinds)))) = {[0 1 0]};%NOTE - before using our EA_intact idx here, we have to shorten it to just have 1s and 0s that agree with the size of our shrunken-down matrix; in this case, since our cm_2c drops values associated with AA, etc., we need to drop 0s in our EA_intact index that are associated with those same 0s. We can do this easily by filtering our EA_intact idx by the 'testinds' idx we created above for this example, which has 0s for every pattern which is dropped out of cm_2c and 1s for every pattern which is stil in cm_2c
        xz2(logical(Scene_idx(logical(testinds)))) = {[0 0 1]};
        userOptions.conditionColours2 = cell2mat(xz2);
        color_t2 = xz2(labelReordering);
        hold on;
        x = xlim(gca);
        for condition = 1:size(cm_2c,1)
            plot(x(1), condition, 'o', 'MarkerFaceColor', color_t2{condition, :}, 'MarkerEdgeColor', 'none', 'MarkerSize', 8);
        end%for:condition
        hold off;
    end

    % how do you know if your clustering fit the data well? one way is to
    % ensure the linkage heights it came up with resemble the raw pairwise
    % distances in the original matrix.
    %c = cophenet(Z2,distm_2); % is c, the cophenetic correlation very high? ideally yes. if not, your linkage algorithm may not handle the data well and you might try something else (e.g., 'average' vs 'complete')

end

%% Test models of similarity structure
if run_simmodelfit == 1; % build "model" similarity matrices and test how well these fit observed
    modfits = CMmodelcomparison(cm2,EA_intact,AA_intact,Obj_idx,othercond_idx,Scene_idx,scrambled_idx);
end

%% Save data
savename = [S.group_mvpa_dir '/Rcorrs_' S.subj_id '_' mask '_' weights_str '_' S.exp_name '.mat'];
save(savename, 'res');

end

function data = hp_filter(pat,sel,cutoff,tr)

% max(sel) amounts to the maximum number of runs there could
% be. any values in the runs selector that are <= 0 will be ignored
% by the for loop. won't mind if you're lacking a particular run in
% the middle either

nRuns = max(sel);
data  = pat';	% transposition required for spm_filter

for r = 1:nRuns
    progress(r,nRuns);
    this_run = sel==r; % select current run

    K(r).row = find(this_run == 1);
    K(r).RT  = tr;
    K(r).HParam = cutoff;

    data = spm_filter(K,data);

end
end

function modfits = CMmodelcomparison(cm,idx1,idx2,idx3,idx4,idx5,idx6)
%% construct Four sample CM models
%normally, these would be built from some theoretical or computational models of how information is organized in the
%brain's perceptual, mnemonic, etc, systems

cm2 = cm;
EA_intact = idx1;
AA_intact = idx2;
Obj_idx = idx3;
othercond_idx = idx4;
Scene_idx = idx5;
scrambled_idx = idx6;

cm2_testo1 = cm2;% model 1 - "all faces and objects are "items" to this area - scenes are scenes, and distinct from objects"
cm2_testo1(logical(EA_intact),logical(EA_intact))=0.5;
cm2_testo1(logical(AA_intact),logical(AA_intact))=0.5;
cm2_testo1(logical(Obj_idx),logical(Obj_idx))=0.5;
cm2_testo1(logical(EA_intact),logical(AA_intact))=0.5;
cm2_testo1(logical(EA_intact),logical(Obj_idx))=0.5;
cm2_testo1(logical(AA_intact),logical(Obj_idx))=0.5;
cm2_testo1(logical(AA_intact),logical(EA_intact))=0.5;
cm2_testo1(logical(Obj_idx),logical(EA_intact))=0.5;
cm2_testo1(logical(Obj_idx),logical(AA_intact))=0.5;

cm2_testo1(logical(othercond_idx),logical(othercond_idx))=0.5;

cm2_testo1(logical(othercond_idx),logical(EA_intact))=0.5;
cm2_testo1(logical(othercond_idx),logical(AA_intact))=0.5;
cm2_testo1(logical(othercond_idx),logical(Obj_idx))=0.5;
cm2_testo1(logical(EA_intact),logical(othercond_idx))=0.5;
cm2_testo1(logical(AA_intact),logical(othercond_idx))=0.5;
cm2_testo1(logical(Obj_idx),logical(othercond_idx))=0.5;

cm2_testo1(logical(Scene_idx),logical(Scene_idx))=0.5;

cm2_testo1(logical(EA_intact),logical(Scene_idx))=0.1;
cm2_testo1(logical(AA_intact),logical(Scene_idx))=0.1;
cm2_testo1(logical(Obj_idx),logical(Scene_idx))=0.1;
cm2_testo1(logical(Scene_idx),logical(EA_intact))=0.1;
cm2_testo1(logical(Scene_idx),logical(AA_intact))=0.1;
cm2_testo1(logical(Scene_idx),logical(Obj_idx))=0.1;
cm2_testo1(logical(othercond_idx),logical(Scene_idx))=0.1;
cm2_testo1(logical(Scene_idx),logical(othercond_idx))=0.1;

cm2_testo1(logical(scrambled_idx),logical(scrambled_idx))=0.1;

cm2_testo1(logical(Scene_idx),logical(scrambled_idx))=0.1;
cm2_testo1(logical(Obj_idx),logical(scrambled_idx))=0.1;
cm2_testo1(logical(EA_intact),logical(scrambled_idx))=0.1;
cm2_testo1(logical(AA_intact),logical(scrambled_idx))=0.1;
cm2_testo1(logical(othercond_idx),logical(scrambled_idx))=0.1;
cm2_testo1(logical(scrambled_idx),logical(Scene_idx))=0.1;
cm2_testo1(logical(scrambled_idx),logical(Obj_idx))=0.1;
cm2_testo1(logical(scrambled_idx),logical(EA_intact))=0.1;
cm2_testo1(logical(scrambled_idx),logical(AA_intact))=0.1;
cm2_testo1(logical(scrambled_idx),logical(othercond_idx))=0.1;
%set diagonal back to NaN
cm2_testo1(logical(eye(size(cm2_testo1)))) = NaN;

cm2_testo2 = cm2;% model 2 - "faces > bodies > inanimate objects gradient - scenes are scenes, and distinct from items"
cm2_testo2(logical(EA_intact),logical(EA_intact))=0.5;
cm2_testo2(logical(AA_intact),logical(AA_intact))=0.5;
cm2_testo2(logical(Obj_idx),logical(Obj_idx))=0.5;
cm2_testo2(logical(EA_intact),logical(AA_intact))=0.5;
cm2_testo2(logical(EA_intact),logical(Obj_idx))=0.25;
cm2_testo2(logical(AA_intact),logical(Obj_idx))=0.25;
cm2_testo2(logical(AA_intact),logical(EA_intact))=0.5;
cm2_testo2(logical(Obj_idx),logical(EA_intact))=0.25;
cm2_testo2(logical(Obj_idx),logical(AA_intact))=0.25;

cm2_testo2(logical(othercond_idx),logical(othercond_idx))=0.5;

cm2_testo2(logical(othercond_idx),logical(EA_intact))=0.35;
cm2_testo2(logical(othercond_idx),logical(AA_intact))=0.35;
cm2_testo2(logical(othercond_idx),logical(Obj_idx))=0.25;
cm2_testo2(logical(EA_intact),logical(othercond_idx))=0.35;
cm2_testo2(logical(AA_intact),logical(othercond_idx))=0.35;
cm2_testo2(logical(Obj_idx),logical(othercond_idx))=0.25;

cm2_testo2(logical(Scene_idx),logical(Scene_idx))=0.5;

cm2_testo2(logical(EA_intact),logical(Scene_idx))=0.1;
cm2_testo2(logical(AA_intact),logical(Scene_idx))=0.1;
cm2_testo2(logical(Obj_idx),logical(Scene_idx))=0.15;
cm2_testo2(logical(Scene_idx),logical(EA_intact))=0.1;
cm2_testo2(logical(Scene_idx),logical(AA_intact))=0.1;
cm2_testo2(logical(Scene_idx),logical(Obj_idx))=0.15;
cm2_testo2(logical(othercond_idx),logical(Scene_idx))=0.1;
cm2_testo2(logical(Scene_idx),logical(othercond_idx))=0.1;

cm2_testo2(logical(scrambled_idx),logical(scrambled_idx))=0.1;

cm2_testo2(logical(Scene_idx),logical(scrambled_idx))=0.1;
cm2_testo2(logical(Obj_idx),logical(scrambled_idx))=0.1;
cm2_testo2(logical(EA_intact),logical(scrambled_idx))=0.1;
cm2_testo2(logical(AA_intact),logical(scrambled_idx))=0.1;
cm2_testo2(logical(othercond_idx),logical(scrambled_idx))=0.1;
cm2_testo2(logical(scrambled_idx),logical(Scene_idx))=0.1;
cm2_testo2(logical(scrambled_idx),logical(Obj_idx))=0.1;
cm2_testo2(logical(scrambled_idx),logical(EA_intact))=0.1;
cm2_testo2(logical(scrambled_idx),logical(AA_intact))=0.1;
cm2_testo2(logical(scrambled_idx),logical(othercond_idx))=0.1;
%set diagonal back to NaN
cm2_testo2(logical(eye(size(cm2_testo2)))) = NaN;


%compare with vectorized unique correlations of cm2
truem = cm2(triu(true(size(cm2)),1));
testm1 = cm2_testo1(triu(true(size(cm2_testo1)),1));
[modfits.r_trvste1 modfits.p_trvste1] = corr(truem,testm1,'Type','Spearman');

testm2 = cm2_testo2(triu(true(size(cm2_testo2)),1));
[modfits.r_trvste2 modfits.p_trvste2] = corr(truem,testm2,'Type','Spearman');
%at the group level, you can use statistics like Wilcoxin signed-rank test
%to test for significance of specific model-CM mappings

%% create variant with two categorical predictors for "encoding" regression test (see Mur et al., 2013)
%construct two sample CM models - normally, these would be built from some
%theoretical or computational models of how information is organized in the
%brain's perceptual, mnemonic, etc, systems
cm2_testo3 = cm2;% model 3 - categorical animate-inanimate
cm2_testo3(logical(EA_intact),logical(EA_intact))=1;
cm2_testo3(logical(AA_intact),logical(AA_intact))=1;
cm2_testo3(logical(Obj_idx),logical(Obj_idx))=1;
cm2_testo3(logical(EA_intact),logical(AA_intact))=1;
cm2_testo3(logical(EA_intact),logical(Obj_idx))=-1;
cm2_testo3(logical(AA_intact),logical(Obj_idx))=-1;
cm2_testo3(logical(AA_intact),logical(EA_intact))=1;
cm2_testo3(logical(Obj_idx),logical(EA_intact))=-1;
cm2_testo3(logical(Obj_idx),logical(AA_intact))=-1;

cm2_testo3(logical(othercond_idx),logical(othercond_idx))=1;

cm2_testo3(logical(othercond_idx),logical(EA_intact))=1;
cm2_testo3(logical(othercond_idx),logical(AA_intact))=1;
cm2_testo3(logical(othercond_idx),logical(Obj_idx))=-1;
cm2_testo3(logical(EA_intact),logical(othercond_idx))=1;
cm2_testo3(logical(AA_intact),logical(othercond_idx))=1;
cm2_testo3(logical(Obj_idx),logical(othercond_idx))=-1;

cm2_testo3(logical(Scene_idx),logical(Scene_idx))=1;

cm2_testo3(logical(EA_intact),logical(Scene_idx))=-1;
cm2_testo3(logical(AA_intact),logical(Scene_idx))=-1;
cm2_testo3(logical(Obj_idx),logical(Scene_idx))=1;
cm2_testo3(logical(Scene_idx),logical(EA_intact))=-1;
cm2_testo3(logical(Scene_idx),logical(AA_intact))=-1;
cm2_testo3(logical(Scene_idx),logical(Obj_idx))=1;
cm2_testo3(logical(othercond_idx),logical(Scene_idx))=-1;
cm2_testo3(logical(Scene_idx),logical(othercond_idx))=-1;

cm2_testo3(logical(scrambled_idx),logical(scrambled_idx))=1;

cm2_testo3(logical(Scene_idx),logical(scrambled_idx))=0;
cm2_testo3(logical(Obj_idx),logical(scrambled_idx))=0;
cm2_testo3(logical(EA_intact),logical(scrambled_idx))=0;
cm2_testo3(logical(AA_intact),logical(scrambled_idx))=0;
cm2_testo3(logical(othercond_idx),logical(scrambled_idx))=0;
cm2_testo3(logical(scrambled_idx),logical(Scene_idx))=0;
cm2_testo3(logical(scrambled_idx),logical(Obj_idx))=0;
cm2_testo3(logical(scrambled_idx),logical(EA_intact))=0;
cm2_testo3(logical(scrambled_idx),logical(AA_intact))=0;
cm2_testo3(logical(scrambled_idx),logical(othercond_idx))=0;
%set diagonal back to NaN
cm2_testo3(logical(eye(size(cm2_testo3)))) = NaN;


cm2_testo4 = cm2;% model 4 - categorical scene-"item" (face,body,object)
cm2_testo4(logical(EA_intact),logical(EA_intact))=1;
cm2_testo4(logical(AA_intact),logical(AA_intact))=1;
cm2_testo4(logical(Obj_idx),logical(Obj_idx))=1;
cm2_testo4(logical(EA_intact),logical(AA_intact))=1;
cm2_testo4(logical(EA_intact),logical(Obj_idx))=1;
cm2_testo4(logical(AA_intact),logical(Obj_idx))=1;
cm2_testo4(logical(AA_intact),logical(EA_intact))=1;
cm2_testo4(logical(Obj_idx),logical(EA_intact))=1;
cm2_testo4(logical(Obj_idx),logical(AA_intact))=1;

cm2_testo4(logical(othercond_idx),logical(othercond_idx))=1;

cm2_testo4(logical(othercond_idx),logical(EA_intact))=1;
cm2_testo4(logical(othercond_idx),logical(AA_intact))=1;
cm2_testo4(logical(othercond_idx),logical(Obj_idx))=1;
cm2_testo4(logical(EA_intact),logical(othercond_idx))=1;
cm2_testo4(logical(AA_intact),logical(othercond_idx))=1;
cm2_testo4(logical(Obj_idx),logical(othercond_idx))=1;

cm2_testo4(logical(Scene_idx),logical(Scene_idx))=1;

cm2_testo4(logical(EA_intact),logical(Scene_idx))=-1;
cm2_testo4(logical(AA_intact),logical(Scene_idx))=-1;
cm2_testo4(logical(Obj_idx),logical(Scene_idx))=-1;
cm2_testo4(logical(Scene_idx),logical(EA_intact))=-1;
cm2_testo4(logical(Scene_idx),logical(AA_intact))=-1;
cm2_testo4(logical(Scene_idx),logical(Obj_idx))=-1;
cm2_testo4(logical(othercond_idx),logical(Scene_idx))=-1;
cm2_testo4(logical(Scene_idx),logical(othercond_idx))=-1;

cm2_testo4(logical(scrambled_idx),logical(scrambled_idx))=1;

cm2_testo4(logical(Scene_idx),logical(scrambled_idx))=0;
cm2_testo4(logical(Obj_idx),logical(scrambled_idx))=0;
cm2_testo4(logical(EA_intact),logical(scrambled_idx))=0;
cm2_testo4(logical(AA_intact),logical(scrambled_idx))=0;
cm2_testo4(logical(othercond_idx),logical(scrambled_idx))=0;
cm2_testo4(logical(scrambled_idx),logical(Scene_idx))=0;
cm2_testo4(logical(scrambled_idx),logical(Obj_idx))=0;
cm2_testo4(logical(scrambled_idx),logical(EA_intact))=0;
cm2_testo4(logical(scrambled_idx),logical(AA_intact))=0;
cm2_testo4(logical(scrambled_idx),logical(othercond_idx))=0;
%set diagonal back to NaN
cm2_testo4(logical(eye(size(cm2_testo3)))) = NaN;

testm3 = cm2_testo3(triu(true(size(cm2_testo3)),1));
testm4 = cm2_testo4(triu(true(size(cm2_testo4)),1));

modpredicts = horzcat(testm3,testm4);
%add constant term to modpredicts for valid F and P stats
modpredicts(:,3)=1;%matlab assumes there's a constant column of 1

[modfits.b,bint,r,rint,modfits.stats] = regress(truem,modpredicts);%for this example, run a simple multiple regression.

figure;
subplot(2,2,1), imagesc(cm2_testo1);
title('complex model1')
colormap('jet'); % set the colorscheme
caxis([-1 1])
colorbar; % enable colorbar

subplot(2,2,2), imagesc(cm2_testo2);
title('complex model2')
colormap('jet'); % set the colorscheme
caxis([-1 1])
colorbar; % enable colorbar

subplot(2,2,3), imagesc(cm2_testo3);
title('animate-inanimate')
colormap('jet'); % set the colorscheme
caxis([-1 1])
colorbar; % enable colorbar

subplot(2,2,4), imagesc(cm2_testo4);
title('item-scene')
colormap('jet'); % set the colorscheme
caxis([-1 1])
colorbar; % enable colorbar

end
