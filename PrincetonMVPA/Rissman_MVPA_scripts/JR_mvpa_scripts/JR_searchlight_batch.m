function []=JR_searchlight_batch(subjects_to_include,cond1,cond2,balanced)

topdir = pwd;

for a=1:length(subjects_to_include)
    b = subjects_to_include(a);
    if b<10
        subj_id=strcat('s10', num2str(b))
    else
        subj_id=strcat('s1', num2str(b))
    end
    
    cd([topdir '/' subj_id '/mvpa']);
    %JR_run_mvpa_v4_searchlight(subj_id, cond1,cond2,balanced)
    JR_run_mvpa_v4_searchlight(subj_id,[subj_id '_merged_AAL_ROIs_unsmoothed.mat'], cond1,cond2,balanced)
    %JR_run_mvpa_v4_searchlight(subj_id,[subj_id '_merged_AAL_ROIs_8mm_smoothing.mat'])
    
    cd(topdir)
end