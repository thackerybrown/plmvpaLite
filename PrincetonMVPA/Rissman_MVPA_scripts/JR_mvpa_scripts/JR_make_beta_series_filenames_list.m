function [] = JR_make_beta_series_filenames_list(subj_id,total_num_of_trials)

expt_path = '/Users/Jesse/fMRI/data/PAST/fMRI';

mvpa_path = [expt_path '/' subj_id '/mvpa'];
results_path = [expt_path '/' subj_id '/results_beta_series'];
results_IBS_path = [expt_path '/' subj_id '/results_beta_series_IBS'];

IBS_canon_inds = [1:3:total_num_of_trials*3];

for i = 1:total_num_of_trials
        
    beta_series_filenames{i} = [results_path '/beta_' prepend(num2str(i),4) '.img'];
    beta_series_IBS_filenames{i} = [results_IBS_path '/beta_' prepend(num2str(IBS_canon_inds(i)),4) '.img'];
end

raw_filenames = beta_series_filenames;
save_cmd = ['save ' mvpa_path '/beta_series_filenames.mat raw_filenames'];
eval(save_cmd);

raw_filenames = beta_series_IBS_filenames;
save_cmd = ['save ' mvpa_path '/beta_series_IBS_filenames.mat raw_filenames'];
eval(save_cmd);







