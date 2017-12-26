function [] = mask_img_with_thresholded_img(img_to_thresh, img_to_mask, thresh)

T = spm_vol(img_to_thresh);
t_data = spm_read_vols(T);
t_below_thresh = find(t_data<thresh);

V = spm_vol(img_to_mask);
v_data = spm_read_vols(V);
v_data(t_below_thresh)=0;

V.fname=[img_to_mask(1:end-4) '_t_thresh_' num2str(thresh) '.img'];

spm_write_vol(V,v_data);