function [] = rsa_simple()
% code for RSA analysis simplified!

% example call with 'mvpa_sample_data' - rsa_simple()

% this script is meant to demo RSA basics without the complex parameter
% lines associated with large fMRI datasets

%% data predefined
% for this simple example, bypass having separate model files and reading
% in image names and masks in a complex manner...

studydir = 'C:\Users\giova\Documents\work\prosthesis\';

%specify the "model" - simplified here, this is just condition labels, 
% 1 image per condition per subject
names = {['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact'],['pros'],['intact_post'],['pros_post'],['intact_pre'],['pros_pre'],['intact_post'],['pros_post'],['intact_pre'],['pros_pre']};
onsets = {[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[13],[14],[15],[16],[17],[18],[19],[20],[21],[22],[23],[24],[25],[26],[27],[28],[29],[30],[31],[32],[33],[34],[35],[36],[37],[38],[39],[40]};

%find the MRI files you want to load in as patterns
datadir = 'data';
raw_filenames = dir(fullfile([studydir '\' datadir],'*.img'));

%specify a mask to analyze the patterns within (must be same resolution as
%MRI files you want to analyze
maskdir = 'Masks';
maskname = 'rNativeGM_BOLDres.nii';

masknameandpath = [studydir '\' maskdir '\' maskname];

%% extract all patterns, iterating through 3D MRI frames
for i=1:length(raw_filenames)
    brainmaps{i} = [raw_filenames(i).folder '\' raw_filenames(i).name];%this is built from the filename and the folder path

    [b,r] = MAP_getROI(masknameandpath, brainmaps{i}, 'vox', 0, '');
    bmat_t(:,i) = b{1}; % returns all voxels, whether or not they have NaNs
    rmat_t(:,i) = r; % returns voxels excluding NaNs
    %meanbetas_t = nanmean(bmat_t(:,:));%get mean beta values from the ROI for each regressor

end


%% STOP!  ~Quality Assurance~
%before proceeding...
%is your Mask REALLY binary (i.e., do you get the data from the number of
%voxels you think you do?

%check # of rows in bmat_t -> if your ROI was 100 voxels this number should
%be = 100 

%% Flags and parameters
%optional - toss values below a certain number as well (e.g., maybe a zero
%is equivalent to a NaN for your study. Say you are using a mask on "raw"
%BOLD data that does nothing to account for signal drop-out and extra-brain
%voxels are zeros instead of NaNs - here we can fix that
threshpats = 1; % 1 = YES
threshold = 0; % you come up with your scheme.  example: 0.01*mean(mean(rmat_condensed')); - this example will threshold out anything <99% of the average (quite liberal thresholding)

runzscore = 0;%1=yes. Standard but controversial preprocessing.

%% preprocessing (using flags) above

if threshpats == 1
    thresh = threshold;
    x1 =[]; %vector of filtered intensity values
    for p = 1:length(rmat_t(1,:))
        x1 = [x1 (rmat_t(:,p)>thresh)];
    end
    t = mean(x1');% average across columns - anything less than 1 indicates there are some zeros
    t1 = t'>=1; %threshold once more to voxels that had signal passing our threshold for EVERY pattern

    rmat_condensed = rmat_t(t1,:);
end


% zscore within runs
if runzscore == 1
    figure;
    %plot for exploration
    subplot(2,1,1), imagesc(rmat_condensed);
    title('original pattern')
    colormap('jet'); % set the colorscheme
    colorbar;
    %caxis([-1 1]);
    hold on

    pat_t = [];
    pat_t = [pat_t zscore_mvpa(rmat_condensed,2)];%2 = z-score within rows (within-voxels)
    rmat_condensed = pat_t;

    subplot(2,1,2), imagesc(rmat_condensed);
    title('z_scored pattern')
    colormap('jet'); % set the colorscheme
    colorbar;
    caxis([-1 1]);
end



%% create correlation matrices

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



%% index different categories
%first, let's find "intact" vs "prosthetic" events
for n = 1:length(names)
    if contains(names{n},{'intact'})%
        intact_idx(n)=1;
        pros_idx(n)=0;
    elseif contains(names{n},{'pros'})%
        intact_idx(n)=0;
        pros_idx(n)=1;
    end
end

%now, let's find "pre" vs "post" events
for n = 1:length(names)
    if contains(names{n},{'pre'})%
        pre_idx(n)=1;
        post_idx(n)=0;
    elseif contains(names{n},{'post'})%
        pre_idx(n)=0;
        post_idx(n)=1;
    end
end

%finally, filter subject groups to pre, post, no training
int_notrain_idx = intact_idx-(intact_idx.*pre_idx)-(intact_idx.*post_idx); %subtract out unwanted...
pros_notrain_idx = pros_idx-(pros_idx.*pre_idx)-(pros_idx.*post_idx); %subtract out unwanted...

int_pretrain_idx = intact_idx.*pre_idx; %intersection of idx...
pros_pretrain_idx = pros_idx.*pre_idx; %intersection of idx...

int_posttrain_idx = intact_idx.*post_idx; %intersection of idx...
pros_posttrain_idx = pros_idx.*post_idx; %intersection of idx...

%% compute some RS...

% global similarity measures
res.intact_w_intact = cm(logical(int_notrain_idx),logical(int_notrain_idx));
res.intact_w_intact_mean = nanmean(res.intact_w_intact(:));

res.pros_w_pros = cm(logical(pros_notrain_idx),logical(pros_notrain_idx));
res.pros_w_pros_mean = nanmean(res.pros_w_pros(:));

res.intact_w_pros = cm2(logical(int_notrain_idx),logical(pros_notrain_idx));
res.intact_w_pros_mean = nanmean(res.intact_w_pros(:));

% how does training impact the two subjects?
% first - how similar are intacts before prosthetic training?
res.intact_w_pretrainint = cm2(logical(int_notrain_idx),logical(int_pretrain_idx));
res.intact_w_pretrainint_mean = nanmean(res.intact_w_pretrainint(:));
% does that change with prosthetic training?
res.intact_w_posttrainint = cm2(logical(int_notrain_idx),logical(int_posttrain_idx));
res.intact_w_posttrainint_mean = nanmean(res.intact_w_posttrainint(:));


%ok - do prosthetic patterns "grow away" from intact with training?
% pre-training...
res.intact_w_pretrainpros = cm2(logical(int_notrain_idx),logical(pros_pretrain_idx));
res.intact_w_pretrainpros_mean = nanmean(res.intact_w_pretrainpros(:));
% post-training...
res.intact_w_posttrainpros = cm2(logical(int_notrain_idx),logical(pros_posttrain_idx));
res.intact_w_posttrainpros_mean = nanmean(res.intact_w_posttrainpros(:));

%ok - do prosthetic patterns "grow away" from naive prosthetic with training?
% pre-training...
res.pros_w_pretrainpros = cm2(logical(pros_notrain_idx),logical(pros_pretrain_idx));
res.pros_w_pretrainpros_mean = nanmean(res.pros_w_pretrainpros(:));
% post-training...
res.pros_w_posttrainpros = cm2(logical(pros_notrain_idx),logical(pros_posttrain_idx));
res.pros_w_posttrainpros_mean = nanmean(res.pros_w_posttrainpros(:));



% %% Plots
% %plot matrices of interest
% figure;
% subplot(2,2,1), imagesc(cm);
% title('Within-cond corrmat');
% colormap('jet'); % set the colorscheme
% caxis([-1 1]);
% colorbar; % enable colorbar
% 
% subplot(2,2,2), imagesc(cm2);
% title('Overall corrmat');
% colormap('jet'); % set the colorscheme
% caxis([-1 1]);
% colorbar; % enable colorbar
% 
% subplot(2,2,3), imagesc(cm(logical(EA_intact),logical(EA_intact)));
% title('EA with EA across blocks and runs');
% colormap('jet'); % set the colorscheme
% caxis([-1 1]);
% colorbar; % enable colorbar
% 
% subplot(2,2,4), imagesc(cm2(logical(EA_intact),logical(Scene_idx)));
% title('EA with Scene across blocks and runs');
% colormap('jet'); % set the colorscheme
% caxis([-1 1]);
% colorbar; % enable colorbar
% 
% % save plot
% plot_savename = [S.group_mvpa_dir '/Rcorrs_' S.subj_id '_' mask '_' weights_str '_' S.exp_name '_corrmats.png'];
% saveas(gcf,plot_savename);




%% end of script; Save data
savename = [studydir '/Rcorrs_' maskname '.mat'];
save(savename, 'res');

end


