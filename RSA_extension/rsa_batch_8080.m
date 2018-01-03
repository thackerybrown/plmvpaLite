%addpath('/biac4/wagner/wagner/thackery/scripts', '/biac4/wagner/wagner/thackery/CM/masks')
sub_grp = {'001' '003'}' 
TRsperRun = {[114,114] [114,114]};

for i = 1:length(sub_grp)
    sub = sub_grp{i};
   
   %% ROI by ROI analysis - generates similarity data for each ROI in loop
   GMmask = ['NativeGM_BOLDres'];
   rsa_CM_Localizer(sub_grp(i), GMmask, TRsperRun{i});

   %add additional ROIs here
%   vtcmask = ['bilat-parahipp_fusi_latocc_inftemp'];   
end

