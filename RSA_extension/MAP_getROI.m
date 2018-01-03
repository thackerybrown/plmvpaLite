function [Betas, results, mask_name]=MAP_getROI(maskname, betamaps, operation, threshold, condition)

%GETBETA usage: [Betas]=getROI(mask, betamaps)
%Returns the betas (as a cell array, one cell per image) from images listed in the string array betamaps, based on the mask.
%               [Betas, results]=getROI(mask, betamaps, operation, threshold)
%Returns also the result of a requested operation of the masked values
%
%Input
%mask       a file name (full path if not in the current directory) as a string
%betamaps   a string array of file names (full paths if not in the current directory)
%operation  (optional) a manipulation to do on the values extracted. Possible values include:
%           'mean'    Mean of all non NaN voxels within the mask
%           'over'    Count, percent, and mean of suprathreshold voxels (requires a threshold)
%           'below'   Count,  Percent, and mean of subthreshold voxels (requires a
%                       threshold)
%           'none'    (default)
%threshold  (needed for some operations) a threshold to use for count and percent 
%Output:
%Betas      A cell array, where each cell is a vector containing all the values from the betamap that were indexed in the mask 
%results    a vector with the required operation result for each input map.
%If an operation is specified and two output arguments are given, a tab delimited text
%file will also be written to disk 
%
% An error occurrs if the mask and one of the maps is not the same voxel size, origin, etc.
%
%see STRVCAT for how to make string arrays, and try spm_get as well.

% 9/2003 by Leon Deouell, The Hebrew University of Jerusalem
% msleon@mscc.huji.ac.il
%revision 9/25/03:  over and below return the mean supra/subthreshol voxels
%                   too. The 'mean' operation returns the mean of all non
%                   Nan voxels in the mask.

if nargin<3
    operation='none';
    if nargout>2
        warning('nothing to do with the second output argument. It will be an empty matrix');
        results='';
    end
end
if nargin > 2 && nargout < 2
    error ('Not enough output arguments for operations. See help');
end
% disp(['Reading mask: ' maskname])
% [mask.V, mask.data]=spmBIC_readimg(mask);

mask.V = spm_vol(maskname);
mask.data = spm_read_vols(mask.V);

mask.include=find(mask.data);
mask.size=size(mask.include, 1);
for i=1:size(betamaps,1)
%     disp(betamaps(i,:))
% 	[Beta.V, Beta.data]=spmBIC_readimg(strtok(betamaps(i,:)));
    Beta.V = spm_vol(strtok(betamaps(i,:)));
    Beta.data = spm_read_vols(Beta.V);
    
    
    % check that the images are similar
	if Beta.V.mat~=mask.V.mat
        error('not the same voxel dimensions or origin')
	end
    if Beta.V.dim~=mask.V.dim
        error('not the same file dimensions')
	end    
    % Get the betas
	Betas{i}=Beta.data(mask.include);
    current=Betas{i};
    switch lower(operation)
        case 'mean'
            if any(isnan(current))
%                 disp(['dropping ' num2str(sum(isnan(current))) ' NaN voxels'])
            end
            current=current(~isnan(current));
            results(i,1)=size(current,1);
            results(i,2)=mean(current); %this would average all voxels within the mask together; we want to separate them.           
            results(i,3)=std(current);
            
        case 'vox'   %use for splitting mask into separate voxels, will spit out t-value for each voxel within mask.
            if any(isnan(current))
%                 disp(['dropping ' num2str(sum(isnan(current))) ' NaN voxels'])
            end
            current=current(~isnan(current));
            %added this to output all voxels separately within mask.
            for a=1:size(current,1)
                results(i,a)=(current(a));  %added another dimension to vector to hold voxel data.
            end
            
        case 'over'
            if ~exist('threshold','var')
                error('Need threshold for this operation. See help')
            end
            suprathreshold=current(find(current>threshold));
            
            mask_voxels = find(mask.data);
            supra_voxels_in_mask = find(current>threshold);
            supra_voxels_in_brain_space = mask_voxels(supra_voxels_in_mask);
            
            beta_map = spm_vol(deblank(betalist(i,:)));
            beta_data = spm_read_vols(beta_map);
            all_betas_in_mask = beta_data(mask_voxels);
            suprathreshold_betas = beta_data(supra_voxels_in_brain_space);           
            
            if isempty(suprathreshold)
                disp('---------------------------------------')
                warning(['No voxels over threshold for ' betamaps(i,:)])
                disp('---------------------------------------')
				results(i,1:4)=[0 0 NaN mean(all_betas_in_mask) ];
            else
                                             
                results(i,4)=mean(all_betas_in_mask);
                %results(i,4)=mean(suprathreshold_betas)*1000;
                results(i,3)=mean(suprathreshold);
                results(i,2)=size(suprathreshold,1)*100/mask.size;
                results(i,1)=size(suprathreshold,1);
            end
            
        case 'below'
            if ~exist('threshold','var')
                error('Need threshold for this operation. See help')
            end
            subthreshold=current(find(current<threshold));
            if isempty(subthreshold)
                disp('---------------------------------------')
                warning(['No voxels below threshold for ' betamaps(i,:)])
                disp('---------------------------------------')       
                results(i,1:3)=[0 0 NaN];
            end
            results(i,3)=mean(subthreshold);   
            results(i,2)=size(subthreshold,1)*100/mask.size;
            results(i,1)=size(subthreshold,1);
        otherwise
            error(['No such operation: ' operation])
    end
        
end

if ~exist('results', 'var') %in case the mask is empty
    results = [];
end

% if exist('results', 'var')
% 	%[f,p]=uiputfile('*.txt','Save results to file')
% 	%if f~=0
%         %outfid=fopen(fullfile(p, f),'w');
%         threshstring = num2str(threshold);
%         outputfile = ['/r/d3/despo2/jeff/varian_data/MAP/ROI_stats/' mask_name(54:58) '_' condition '.txt'];  %% my addition
%         outfid=fopen(outputfile,'a');
%         fprintf(outfid, 'Mask: %s\r\n',mask.V.fname);
% 		disp('Writing text file...')
% 
%         switch lower(operation)
%             case 'mean'
%                 fprintf(outfid, 'Operation: %s\r\n', operation);;
%                 fprintf(outfid, '%s\t%s\t%s\r\n','name','count (non NaNs)', 'mean')   
%             case {'over', 'below'}
% 				fprintf(outfid, 'Operation: %s (%s)\r\n', operation, num2str(threshold));;
%                 fprintf(outfid, '%s\t%s\t%s\t%s\t%s\t%s\r\n','subject','count', 'percent','supra_t','supra_b','all_betas')
%         end
%         for i=1:size(betamaps,1)
%             fprintf(outfid, '%s\t', strtok(betalist(i,35:37)));
%             fprintf(outfid, '%.2f\t', results(i,:));
% 			fprintf(outfid, '\r\n');
%         end
%         fclose(outfid)
%     end
% end
    
