function [ output_args ] = rsa_MakeBetaMaps(Sub)
%Reads in onsets and regs from existing model and creates and runs
%   additional GLMs
%
% written TIB, 2/18/2014

%% set paths and obtain sub-specific params

Anlz = 0; %are the data nifti or analyze (img+hdr)?

% set base path
startdir = pwd;

% define subject(s)
sub = Sub;
%S.patRegDir = ['/biac4/wagner/wagner/thackery/Circmaze/subs_preprocessed/cm' sub '/model_cue_plan_goal_forbetas/'];
S.patRegDir = ['/media/sf_host/mvpa_sample_data/CM_localizer/CM' sub '/results01/LSSshort/'];
%S.moddir = ['/biac4/wagner/wagner/thackery/Circmaze/subs_preprocessed/cm' sub '/model_cue_plan_goal_forbetas/'];
if ~exist(S.patRegDir)
    mkdir(S.patRegDir);
end
%S.patRegDir = ['/Users/thackery/Work/rsa/' sub '/Results_model2/Study/RSA_test' ]
% which model? (still, want to make sure this is set in ar_params)
%model = 4; % R7M

% obtain params specific to subject
%[S par idx] = ar_ResponseParams(sub);

% load in onsets and regs of non-interest
onsetsmat = ['/media/sf_host/mvpa_sample_data/CM_localizer/CM' sub '/results01/CM' sub '_localizer_onsets_test_short.mat'];
regsmat = ['/media/sf_host/mvpa_sample_data/CM_localizer/CM' sub '/results01/CM' sub '_MvmtParamsConcat_AND_sessFX.mat'];

load(onsetsmat);
namesOrig = names;
onsetsOrig = onsets;
dursOrig = durations;
load(regsmat);

%create expanded list of onsets, durations, names
allonsets = [];
allcondnames = {};  %use this OR allcondnames
all_stimnames = []; %use this OR allcondnames
alldurations = [];
tempinx=0;

%function to get a list of the stimulus name for each onset (alternative to
%condition+trial# naming scheme for betas)
%[stimlist] = extract_stimnames(sub); %doesn't work on matlab 2012?
%stimlist_file =['/biac4/wagner/wagner/thackery/CM/Data/CM' sub '/Results_model2/safety/stimlist.mat']
%load(stimlist_file)

for i = 1:length(onsets)
    allonsets= [allonsets, onsets{i}];
    % all_stimnames=[all_stimnames, stimlist{i}];%list of the stimulus name for each onset (alternative to
    %condition+trial# naming scheme for betas)
    for j = 1:length(onsets{i})
        allcondnames{tempinx+j}= names{i};
        alldurations = [alldurations, durations{i}];
    end
    tempinx=tempinx+length(onsets{i}) ;
    
end

%% specify and run GLMs for each individual trial

for a = 1:length(allonsets);
    
    thisOnset = allonsets(a);
    thisOnsetName = [prepend(num2str(a),4)];
    %     thisStimName = all_stimnames{a}
    %thisDur = 0; % model all trials of interest as stick functions
    thisDur = alldurations(a); % model all trials of interest with durations provided above
    condName = allcondnames{a};
    
    % create new cond for thisOnset in names/onsets/durs
    names = namesOrig;
    names(2:end+1) = names(1:end);
    names{1} = [thisOnsetName '_' condName '_event' ];% 'single';
    onsets = onsetsOrig;
    onsets(2:end+1) = onsets(1:end);
    onsets{1} = thisOnset;
    durations = dursOrig;
    durations(2:end+1) = durations(1:end);
    durations{1} = thisDur;
    
    % remove thisOnset from its original location in onsets
    condNum = find(strcmp(names,condName));
    tNum = find(onsets{condNum} == thisOnset);
    onsets{condNum}(tNum) = [];
    
    %note, for a study with multiple durations, you need a similar command
    %for durations.
    
    cd(S.patRegDir);
    
    % specify model
    mod_spec_mvpa(sub,S,names,onsets,durations,R);
    
    % estimate model
    mod_est_mvpa(sub,S);
    
    % move and rename first beta, delete all other files
    cd(S.patRegDir);
    if ~exist('betas', 'dir')
        mkdir('betas');
    end
    !mv beta_0001* betas
    !rm beta*
    !rm Res*
    !rm SPM.mat
    cd('betas');
    if strcmp(condName,'id') % stupid hack to create file names of equal length
        condName = 'ide';
    end
    
    
    if Anlz == 1
        tmpimg=dir('beta*.img');
        tmphdr=dir('beta*.hdr');
        
    else
        tmpimg=dir('beta*.nii');
    end
    
 %   for n = 1:length(names);
 %       currname = names{n};
 %       for onsidx = 1:length(onsets{n})
            if Anlz == 1
                betaName = [names{1} '.img'];%[currname '_' num2str(onsidx) '.img']; %names all betas similarly based on the condition and trial number
                betaNamehd = [names{1} '.hdr'];%[currname '_' num2str(onsidx) '.hdr'];
                
                bimg = tmpimg(1,1).name;
                bhdr = tmphdr(1,1).name;
                copyfile(bhdr, betaNamehd);
                copyfile(bimg, betaName);
            else
                betaName = [names{1} '.nii'];%[currname '_' num2str(onsidx) '.nii']; %names all betas similarly based on the condition and trial number
                
                bimg = tmpimg(1,1).name;
                copyfile(bimg, betaName);
            end
 %       end
        cd(startdir);
        
        % clear vars to avoid cumulatively adding to them
        close all
        clear names onsets durations
        
    %end
    
     cd(S.patRegDir);
     cd('betas');
     !mv * ../
     cd(S.patRegDir);
     ! rm -r betas/
     cd(startdir);
    
end




