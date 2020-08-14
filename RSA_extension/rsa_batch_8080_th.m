sub_grp = {'001'}' 
TRsperRun = {[800]};

for i = 1:length(sub_grp)
    sub = sub_grp{i};
   
   %% ROI by ROI analysis - generates similarity data for each ROI in loop
   %GMmask = ['NativeGM_BOLDres'];
   %rsa_CM_Localizer(sub_grp(i), GMmask, TRsperRun{i});

   vtcmask = ['latIPS'];   
   rsa_CM_Localizer_th(sub_grp(i), vtcmask, TRsperRun{i});
   
   %add additional ROIs here

end

