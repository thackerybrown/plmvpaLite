function [] = JR_replace_zeros_with_nans(wildcard)

if exist('wildcard')
    P = spm_get(inf,['*' wildcard '*.img'],'Select files to average...')
else
    P = spm_get(inf,'*.img','Select files to average...')
end




for i = 1:size(P,1)
    
    V = spm_vol(P(i,:));
    V_data = spm_read_vols(V);
    V_data_zeros = find(V_data == 0);
    V_data(V_data_zeros) = NaN;
    [path name ext versn] = fileparts(P(i,:));    
    V.fname = [path '/nan_' name ext];
    spm_write_vol(V,V_data);    
end
    
