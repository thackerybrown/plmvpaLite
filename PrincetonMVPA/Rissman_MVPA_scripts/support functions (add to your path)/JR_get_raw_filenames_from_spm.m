function [] = JR_get_raw_filenames_from_spm(subj)

load SPM.mat

temp = SPM.xY.P;

raw_filenames = cellstr(temp(:,1:end-2));

cd(['/Users/Jesse/fMRI/data/PAST/fMRI/' subj '/mvpa/']);

save raw_filename_list_SRA.mat raw_filenames;
