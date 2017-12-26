function []=JR_exclude_motor_mask

M = spm_vol('rbilat_motor_handdrawn.nii');
mvox = find(spm_read_vols(M));
W = spm_vol('whole_brain.img');
wdata = spm_read_vols(W);
wdata(mvox) = 0;
W.fname = 'whole_brain_no_motor.img';
spm_write_vol(W,wdata)