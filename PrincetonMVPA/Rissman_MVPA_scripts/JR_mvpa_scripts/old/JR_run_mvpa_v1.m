function [subj results] = JR_run_mvpa_v1(mvpa_workspace)

tic

if ~exist('mvpa_workspace')

    %%%%%%% specify user-defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subj_id = 's104';
    exp_name = 'PAST';
    roi_file = 'whole_brain_no_L_motor.img';
    anova_p_thresh = .005;  %p threshold for feature selection ANOVA
    generate_importance_maps = 1; % 1 = yes please, 2 = no thanks
    vol_info = spm_vol(['/Users/Jesse/fMRI/data/PAST/fMRI/' subj_id '/anatomical/s104_mean_func.img']); %get functional data resolution info for spm .img writing

    num_runs = 9;
    num_TP_per_run = 203;
    total_TP = num_runs * num_TP_per_run;

    % load user-created filename and onsets lists into workspace
    load('raw_filename_list_SRA.mat')
    load(['/Users/Jesse/fMRI/data/PAST/behavioral/FACE_expt/data/' subj_id '/onsets.mat']);

    num_conds = size(onsets,2);
    all_regs = zeros(num_conds,total_TP); % initialize regs matrix as conditions x timepoints

    for cond = 1: num_conds
        for trial = 1: length(onsets{cond})
            time_idx = onsets{cond}(trial)/2+1; % divide by 2 and add 1 to convert back from sec to TRs (first timepoint = 0 sec; first TR = 1)
            all_regs(cond,time_idx) = 1;
        end
    end
    % condnames = names; %specify names 'OLD_recollect'    'OLD_hc_old'    'OLD_lc_old'    'OLD_lc_new'    'OLD_hc_new'    'NEW_recollect' 'NEW_hc_old'    'NEW_lc_old'    'NEW_lc_new'    'NEW_hc_new'    'no_resp'

    % define conditions of interest
    objective_old = sum(all_regs(1:5,:)); 
    objective_new = sum(all_regs(7:10,:)); % changed for s104 from objective_new = sum(all_regs(6:10,:));
    
    subjective_old_all = sum(all_regs([1 2 3 6 7 8],:));
    subjective_new_all = sum(all_regs([4 5 9 10],:));
    
    subjective_old_hc_only = sum(all_regs([1 2 6 7],:));
    subjective_new_hc_only = sum(all_regs([5 10],:));
    
    correct_recollect = all_regs(1,:);
    correct_hc_old = all_regs(2,:);    
    
    no_resp = all_regs(11,:); %excluded from analysis

    %choose conditions to train/test classifier on
    regs(1,:) = objective_old;
    regs(2,:) = objective_new;
   
    condnames =  {'objective_old','objective_new'}; %{'old','new'};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % initialize subj structure
    subj = init_subj(exp_name,subj_id);

    % load mask file
    subj = load_analyze_mask(subj,roi_file,roi_file);

    
    % load functional data
    subj = load_analyze_pattern(subj,'spiral',roi_file,raw_filenames);

    % move pattern to hard drive
    % subj = move_pattern_to_hd(subj, 'spiral');

    % make runs vector
    subj = init_object(subj,'selector','runs');

    trial_idx = 1;
    for r = 1:num_runs

        runs(trial_idx:203*r) = r;
        trial_idx = trial_idx + num_TP_per_run;
    end

    subj = set_mat(subj,'selector','runs',runs);

    % detrend the timeseries data
    subj = detrend_runs(subj,'spiral','runs');  % not in mvpa tutorial, but seems important to do

    % move pattern to hard drive
    %subj = move_pattern_to_hd(subj, 'spiral_d');
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','spiral');
    
    % high-pass filter the timeseries data
    subj = hpfilter_runs(subj,'spiral_d','runs',100,2); % remove frequencies below .01 Hz

     % clean up workspace
    subj = remove_mat(subj,'pattern','spiral_d');
    
    % move pattern to hard drive
    %subj = move_pattern_to_hd(subj, 'spiral_d_hp');
    
    % zscore the data from each run
    subj = zscore_runs(subj,'spiral_d_hp','runs'); % Z-score the data
    
    % clean up workspace
    subj = remove_mat(subj,'pattern','spiral_d_hp');
    
    %save_cmd = ['save ' subj_id '_spiral_d_hp_z_' date '.mat'];
    %eval(save_cmd); 

else
    load(mvpa_workspace); %load in pre-saved workspace files and continue
end

% select TRs of interest (to correspond with peak post-stim BOLD response)

all_trials = sum(regs,1); % vector of all trial   (patterns{4} contains fully preprocessed data)

TR0 = subj.patterns{4}.mat(:,find(all_trials)-1); % -1 TR (-2-0 sec)
TR1 = subj.patterns{4}.mat(:,find(all_trials)+0); % 1st TR (0-2 sec)
TR2 = subj.patterns{4}.mat(:,find(all_trials)+1); % 2nd TR (2-4 sec)
TR3 = subj.patterns{4}.mat(:,find(all_trials)+2); % 3rd TR (4-6 sec)
TR4 = subj.patterns{4}.mat(:,find(all_trials)+3); % 4th TR (6-8 sec)
TR5 = subj.patterns{4}.mat(:,find(all_trials)+4); % 5th TR (8-10 sec)

mean_data = (TR3+TR4+TR5)/3; %take uniform average of peak signal TRs

clear TR0 TR1 TR2 TR3 TR4 TR5; %clean up matlab workspace to save memory



% BASELINE CORRECT MEAN DATA
% mean_baseline = (t0+t1)/2;
% mean_data = mean_data - mean_baseline;


% plot mean timecourses for whole ROI (useful check of data quality & overall activity pattern)

new_ind = find(objective_new);
for t = 1: length(new_ind)
     new_ts(t,:) = mean(subj.patterns{4}.mat(:,new_ind(t):new_ind(t)+6),1);
end

old_ind = find(objective_old);
for t = 1: length(old_ind)
     old_ts(t,:) = mean(subj.patterns{4}.mat(:,old_ind(t):old_ind(t)+6),1);
end

figure;
plot(mean(old_ts,1),'b');
hold on;
plot(mean(new_ts,1),'r');


% condense regs by removing zeros

trial_counter = 1;
for i = 1: size(regs,2)
    if ~isempty(find(regs(:,i)))
        condensed_regs(:,trial_counter) = regs(:,i);
        condensed_all(:,trial_counter) = all_regs(:,i); %save
        condensed_runs(1,trial_counter) = runs(i);
        trial_counter = trial_counter + 1;
    end
end

subj = init_object(subj,'regressors','conds');
subj = set_mat(subj,'regressors','conds',condensed_regs);
subj = set_objfield(subj,'regressors','conds','condnames',condnames);

% add new condensed activation pattern
subj = duplicate_object(subj,'pattern','spiral_d_hp_z','spiral_d_hp_z_condensed');
subj = set_mat(subj,'pattern','spiral_d_hp_z_condensed',mean_data);

zhist = sprintf('Pattern ''%s'' created by custom code','spiral_d_hp_z_condensed');
subj = add_history(subj,'pattern','spiral_d_hp_z_condensed',zhist,true);

 
% clean up workspace
subj = remove_mat(subj,'pattern','spiral_d_hp_z');


% update run vector to condensed format
subj.selectors{1}.mat = condensed_runs;
subj.selectors{1}.matsize = size(condensed_runs);

% create cross-validation indices
subj = create_xvalid_indices(subj,'runs');

% % remove rest timepoints
%
% temp_sel = ones(1,size(regs,2));
% temp_sel(find(sum(regs)==0)) = 0; % remove all rest timepoints
% subj = init_object(subj,'selector','no_rest');
% subj = set_mat(subj,'selector','no_rest',temp_sel);
% subj = create_xvalid_indices(subj,'runs','actives_selname','no_rest');

% run feature selection ANOVA
subj = feature_select(subj,'spiral_d_hp_z_condensed','conds','runs_xval','thresh',anova_p_thresh);

% run classifier (hidden layer netlab backprop algorithm)
% class_args.train_funct_name = 'train_bp_netlab';
% class_args.test_funct_name = 'test_bp_netlab';
% class_args.nHidden = 10;
% [subj results] = cross_validation(subj,'spiral_d_hp_z','conds','runs_xval',subj.masks{2}.group_name,class_args);

% run classifier (No hidden layer NN Toolbox backprop algorithm)
class_args.train_funct_name = 'train_bp';
class_args.test_funct_name = 'test_bp';
class_args.nHidden = 0;
for i = 1:4
[subj results{i}] = cross_validation(subj,'spiral_d_hp_z_condensed','conds','runs_xval',subj.masks{2}.group_name,class_args);
end

correct_vector = [];
desireds_vector = [];
guesses_vector = [];
acts_diff_vector = [];
for a = 1:num_runs
    correct_vector = horzcat(correct_vector,results.iterations(a).perfmet.corrects);
    desireds_vector = horzcat(desireds_vector,results.iterations(a).perfmet.desireds);
    guesses_vector = horzcat(guesses_vector,results.iterations(a).perfmet.guesses);
    acts_diff_vector = horzcat(acts_diff_vector, results.iterations(a).acts(1,:)-results.iterations(a).acts(2,:));
end

classification_accuracy_by_resp = [];
for b = 1:10
    classification_accuracy_by_resp(b) = mean(correct_vector(find(condensed_all(b,:))));
    mean_acts_diffs_by_resp(b) = mean(acts_diff_vector(find(condensed_all(b,:))));
end

classification_accuracy_by_resp
mean_acts_diffs_by_resp

[sorted_diffs ind] = sort(acts_diff_vector,2,'descend');
correct_sorted = correct_vector(ind);
desireds_sorted = desireds_vector(ind);

% for i = 1:length(sorted_diffs);
% hit_rate(i) = length(correct_sorted(intersect(find(desireds_sorted == 1),[1:i]))) / length(find(desireds_sorted == 1));
% fa_rate(i) = length(correct_sorted(intersect(find(desireds_sorted == 2),[1:i]))) / length(find(desireds_sorted == 2));
% end

% created binned ROC
counter = 1;
for i = 1:8
hit(i) = length(correct_sorted(intersect(find(desireds_sorted == 1),[1:counter+45]))) / length(find(desireds_sorted == 1));
fa(i) = length(correct_sorted(intersect(find(desireds_sorted == 2),[1:counter+45]))) / length(find(desireds_sorted == 2));
counter = counter+46;
end
figure
plot(fa,hit,'.-')
hold on
plot([0 1],[0 1],'r')
xlabel('P(Old|New)')
ylabel('P(Old|Old)')

if generate_importance_maps == 1;
           
    subj = interpret_weights(subj, results);
   
     
    immap1 = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 1
    immap2 = zeros(vol_info.dim); %initialize appropriately sized matrix for importance map 2
    
    for j = 1:num_runs
        temp1 = zeros(vol_info.dim); %initialize appropriately sized matrix
        temp2 = zeros(vol_info.dim); %initialize appropriately sized matrix
        voxel_inds{j} = find(subj.masks{end-num_runs+j}.mat); %get mask voxel indices
        temp1(voxel_inds{j})=subj.patterns{end-num_runs+j}.mat(:,1); %store immap values at appropriate voxel indices
        temp2(voxel_inds{j})=subj.patterns{end-num_runs+j}.mat(:,2); %store immap values at appropriate voxel indices
        immap1 = immap1+temp1; %add values cumulatively across iterations
        immap2 = immap2+temp2;
    end
    
    %sum across masks to get composite mask (where value of each voxel =
    %number of runs for which that voxel was included)
    composite_mask = zeros(vol_info.dim);
    for i = 2:size(subj.masks,2)  %exclude first mask (it's the starting ROI)
        composite_mask = composite_mask+subj.masks{i}.mat;
    end
    voxels_to_exclude = find(composite_mask<6);
    immap1(voxels_to_exclude)=0;
    immap2(voxels_to_exclude)=0;
    
    immap1_avg = immap1./composite_mask * 1000;  % divide by number of observations contributing to each sum (to get avg) and multiply by 1000 for scaling
    immap2_avg = immap2./composite_mask * 1000;  % divide by number of observations contributing to each sum (to get avg) and multiply by 1000 for scaling
    
    vol_info.fname = [condnames{1} '_p' num2str(anova_p_thresh) '_' mvpa_datetime '.img'];
    spm_write_vol(vol_info,immap1_avg);
    vol_info.fname = [condnames{2} '_p' num2str(anova_p_thresh) '_' mvpa_datetime '.img'];
    spm_write_vol(vol_info,immap2_avg);
end
    
time2finish = toc/60;
disp(['time2finish = ' num2str(time2finish) ' minutes']);

display('all done');