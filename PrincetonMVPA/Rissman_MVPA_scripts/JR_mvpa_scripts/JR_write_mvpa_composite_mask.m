function []= JR_write_mvpa_composite_mask(subj)

temp = zeros(size(subj.masks{2}.mat));

for i = 2:size(subj.masks,2)  %exclude first mask (it's the starting ROI)
    temp = temp+subj.masks{i}.mat;
end
V = spm_vol('whole_brain.img');

thresh_str = num2str(subj.masks{2}.thresh);
V.fname = ['composite_mask_p' thresh_str(3:end) '.img']; %get rid of decimal point in thresh_str
spm_write_vol(V,temp);
