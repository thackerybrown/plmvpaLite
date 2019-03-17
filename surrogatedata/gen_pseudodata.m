%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate_pseudodata.m                      %
% generate pseudo data for testing mvpa code %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% k larocque, august 28, 2014
% dependencies: spm_hrf; shuffle (Princeton MVPA toolbox)

%% params

% output params
outputdir = '/home/brain/host/mvpa_sample_data/CM_pseudodat1/'; % where data will be saved
outprefix = 'pseudo_test'; % prefix for niftis, onsets, & plots

% imaging params
mimicimg = '/home/brain/host/mvpa_sample_data/CM_localizer/CM001/bolds/run_01/run_01_004.nii'; % generate data with same dim as this image (nifti file)
msk = '/home/brain/host/mvpa_sample_data/CM_localizer/CM001/Masks/HVisCtx_1.nii'; %[]: use no mask

% design params
ncond = 3; % how many conditions
ntrials = [25 25 25]; % how many trials per condition
duration = 2; % how long are stimuli on the screen
iti = 8; % how long in between stimuli
tr = 2; % repetition time
weighttrs = [0 .25 .5 0.25 0]; % how are TRs weighted in pattern classification (must be of length (duration + iti)/tr)

% signal params
condmu = [1 1 1]; % mean signal for each condition
condvar = [1 1 1]; % variance across voxels for individual stimuli within each condition
condcov = [0.5 0.25 0.25; 0.25 0.5 0.1; 0.25 0.1 0.5];%[0.8 0.5 0; 0.5 0.8 0.5; 0 0.5 0.8]; % covariance between conditions
noisevar = .05; % variance of noise (mean = 0)
lineardrift = .2; % total linear drift across scan
hrf = spm_hrf(tr); % hrf for convolving data

% check that all is okay
if ~exist(outputdir,'dir'), system(sprintf('mkdir %s',outputdir)); end

assert(length(ntrials) == ncond,                'trial counts mismatch number of conditions');
assert(length(weighttrs) == (duration+iti)/tr,  'tr weights mismatch trial duration');
assert(length(condmu) == ncond,                 'condition means mismatch number of conditions');
assert(length(condvar) == ncond,                'condition variances mismatch number of conditions');
assert(all(condvar >= 0),                       'condition variances must be greater than or equal to zero');
assert(size(condcov,1) == ncond,                'covariance matrix mismatches number of conditions');
assert(size(condcov,1) == size(condcov,2),      'covariance matrix must be square');

% TO DO: add in differential repetition suppression effects across classes


%% generate pseudo data

fprintf('Generating data ... ')

% get dimensions + header info
Vh = spm_vol(mimicimg);
V = spm_read_vols(Vh);
dim = size(V);
clear('V');

% load mask if one is specified
if ~isempty(msk)
    Mh = spm_vol(msk);
    M = spm_read_vols(Mh);
    assert(all(size(M)==dim),'mask dimensions do not match image dimensions!');
end

% create design structure
design=[];
for cx = 1:ncond
    design = [design repmat(cx,1,ntrials(cx))]; %#ok<AGROW>
end
design = shuffle(design);%from Princeton MVPA toolbox %Shuffle(design);
onst = (0:(length(design)-1)) .* (duration + iti);

% make big mu + sigma
bigmu = nan(1,sum(ntrials));
bigsigma = nan(sum(ntrials),sum(ntrials));

for cx = 1:ncond
    
    bigmu(design==cx) = condmu(cx);
    
    for cy = 1:ncond
        
        bigsigma(design==cx,design==cy) = condcov(cx,cy);
        
    end
    
end

for i = 1:length(bigsigma)
    bigsigma(i,i) = condvar(design(i));
end

% generate trial patterns
p = mvnrnd(bigmu,bigsigma,dim(1)*dim(2)*dim(3));

% simulate TRs
tcgt = zeros(size(p,1),(length(design) .* (duration+iti))./tr);
for i = 1:length(design)
    
    tcgt(:,((onst(i)/tr):((onst(i)+duration-1))/tr)+1) = repmat(p(:,i),1,duration./tr);
    
end

% convolve
tcconv = conv2(tcgt',hrf);
tcconv = tcconv';
tcconv = tcconv(:,1:size(tcgt,2));

% add linear drift + noise
tcfinal = tcconv + ...
    repmat(linspace(0,lineardrift,size(tcconv,2)),size(tcconv,1),1) + ...
    sqrt(noisevar) .* randn(size(tcconv,1),size(tcconv,2));

fprintf('generated\n');

fprintf('Outputting data ... ');

% output patterns + onset file
for i = 1:size(tcfinal,2)
    tmp = tcfinal(:,i);
    tmp = reshape(tmp,dim);
    Vt = Vh;
    Vt.fname = fullfile(outputdir,sprintf('%s_%04d.nii',outprefix,i));
    spm_write_vol(Vt,tmp);
end

for c = 1:ncond
    onsets{c} = onst(design==c); %#ok<SAGROW>
    durations{c} = repmat(duration,1,sum(design==c)); %#ok<SAGROW>
    names{c} = sprintf('cond%d',c); %#ok<SAGROW>
end
save(fullfile(outputdir,sprintf('%s_onsets.mat',outprefix)),'onsets','durations','names');

fprintf('outputted\n');

% sanity check: estimate based on TRs
for i = 1:length(design)
    sbst = tcfinal(:,(onst(i)/tr+1):(onst(i)/tr+(duration+iti)/tr));
    est(:,i) = sbst * weighttrs'; %#ok<SAGROW>
end

%% plot

fprintf('Plotting data ... ');

h = figure();

colorscale = [-1 1];

% apply mask to p & est
if ~isempty(msk)
    
    M = M(:);
    p = p(M~=0,:);
    est = est(M~=0,:);
    
end

% plot underlying correlations
subplot(2,2,1)
gtc = corr(p);
imagesc(gtc);
title('ground truth correlation')
%caxis(colorscale);
colorbar();
colormap(hot);
% mask out identity correlations
for i = 1:size(gtc,1)
    gtc(i,i)=NaN;
end

% get mean correlation for each pair type
gtcov = nan(ncond,ncond);
for cx=1:ncond
    for cy=1:ncond
        sbst = gtc(design==cx,design==cy);
        gtcov(cx,cy) = nanmean(sbst(:));
    end
end

subplot(2,2,2);
imagesc(gtcov);
title('ground truth correlation, collapsed')
%caxis(colorscale);
colorbar();
colormap(hot);

% plot estimated correlations
subplot(2,2,3)
estc = corr(est);
imagesc(estc);
title('estimated correlation')
%caxis(colorscale);
colorbar();
colormap(hot);
% mask out identity correlations
for i = 1:size(gtc,1)
    estc(i,i)=NaN;
end

% get mean correlation for each pair type
estcov = nan(ncond,ncond);
for cx=1:ncond
    for cy=1:ncond
        sbst = estc(design==cx,design==cy);
        estcov(cx,cy) = nanmean(sbst(:));
    end
end

subplot(2,2,4);
imagesc(estcov);
title('estimated correlation, collapsed')
%caxis(colorscale);
colorbar();
colormap(hot);

textStrings = num2str(estcov(:),'%0.3f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:3);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
    'HorizontalAlignment','center');
set(hStrings,'Color',[0,0,1])

saveas(h,fullfile(outputdir,sprintf('%s_fig.jpg',outprefix)));
%close(h);

fprintf('plotted\n');
clear;