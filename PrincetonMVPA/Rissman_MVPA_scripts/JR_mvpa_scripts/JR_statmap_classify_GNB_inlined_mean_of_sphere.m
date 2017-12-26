function [acts_vector corrects_vector scratch] = JR_statmap_classify_GNB_inlined(pat,regs,scratch)

% Create a statmap based on cross-validation performance
%
% [VAL SCRATCH] = STATMAP_CLASSIFY(PAT,REGS,SCRATCH)
%
% This contains the logic for running classification on each
% voxel/sphere and creates a statmap of their
% results. Should work equally for a single voxel or bag of
% voxels. More or less reproduces the logic of
% cross-validation for this pat and this regs.
%
% This isn't designed to be called directly. Instead, use
% statmap_searchlight and specify that you want to use this
% as the objective function. You'll need to set the
% statmap_searchlight 'scratch' argument to contain the
% following:
%
% PAT3 (required). Contains the pattern timepoints that were
% marked with 3s in the SELNAME selector passed to
% FEATURE_SELECT. These are the validation timepoints
% that we're going to be testing on.
%
% REGS3 (required). As per PAT3
%
% CLASS_ARGS (required). Will be fed into the train and
% test functions, exactly as in cross_validation.
%
% PERFMET_NAME (required). Which performance metric to
% use, e.g. 'perfmet_maxclass'. You can only specify one.
%
% PERFMET_ARGS (required). Can be empty.

% License:
%=====================================================================
%
% This is part of the Princeton MVPA toolbox, released under
% the GPL. See http://www.csbmb.princeton.edu/mvpa for more
% information.
% 
% The Princeton MVPA toolbox is available free and
% unsupported to those who might find it useful. We do not
% take any responsibility whatsoever for any problems that
% you have related to the use of the MVPA toolbox.
%
% ======================================================================


pat3 = scratch.pat3;
regs3 = scratch.regs3;

if size(pat3,1)>1
    pat3 = mean(pat3,1);
end

[nConds nTimepoints] = size(regs3);
nVox = scratch.nVox;

class_args = scratch.class_args;
perfmet_funct= scratch.perfmet_funct;
perfmet_args = scratch.perfmet_args;

% train and test the classifier
%train_funct_hand = str2func(class_args.train_funct_name);
%class_scratch = train_funct_hand(pat,regs,class_args,[]);

%%% begin inlined train_gnb code %%%

trainpats = mean(pat,1);
%trainpats = pat;
traintargs = regs;
args.uniform_prior = true;

%defaults.uniform_prior = true;
%args = mergestructs(in_args, defaults);
nConds = size(traintargs,1);
[nVox nTimepoints] = size(trainpats);

% find a gaussian distribution for each voxel for each category

scratch_internal = [];
scratch_internal.mu = NaN(nVox, nConds);
scratch_internal.sigma = NaN(nVox, nConds);

for k = 1:nConds

  % grab the subset of the data with a label of category k
    k_idx = find(traintargs(k, :) == 1);

    if numel(k_idx) < 1
      error('Condition %g has no data points.', k);
    end
    
    data = trainpats(:, k_idx);

    % calculate the maximum likelihood estimators (mean and variance)
    [ mu_hat, sigma_hat] = normfit(data');

    scratch_internal.mu(:,k) = mu_hat;
    scratch_internal.sigma(:,k) = sigma_hat;
end

%calculate the priors based on occurence in the training set
scratch_internal.prior = NaN(nConds, 1);
if (args.uniform_prior)
  scratch_internal.prior = ones(nConds,1) / nConds;
else
  
  for k = 1:nConds  
    scratch_internal.prior(k) = (1 + numel( find(traintargs(k, :) == 1))) / ...
        (nConds + nTimepoints);    
  end
  
end


%%% end inlined train_gnb code %%%



%test_funct_hand = str2func(class_args.test_funct_name);
%[acts class_scratch] = test_funct_hand(pat3,regs3,class_scratch;


%%% begin inlined trest_gnb code %%%

testpats = pat3;
testtargs = regs3;

nConds = size(testtargs,1);

[nVox nTimepoints] = size(testpats);

% To make a prediction for a given test pattern, we compute the
% posterior likelihood of observing the label for the given
% pattern:
%
% Pr(Y = y | X = x) ~ Pr( X = x | Y = y) * Pr (Y = y)

warning('off');
% compute the likelihood of the data under the MLE estimated
% gaussian model for each category
for k = 1:nConds

  % we calculate the proper mu and sigma for each voxel and each
  % timepoint:
  
  % scratch.mu is a nVox x 1 vector of voxel means for that condition
  mu = repmat(scratch_internal.mu(:,k), 1, nTimepoints);

  % scratch.sigma is a nVox x 1 vector of voxel variances for a condition 
  sigma = repmat(scratch_internal.sigma(:,k), 1, nTimepoints);

  % compute the likelihood of the data under label k (this just the
  % value of the gaussian equation using the estimated means and variances)
  raw_likelihood = normpdf(testpats, mu, sigma);
  
  % GNB assumes all voxels are independent, so we multiply these probabilities
  % together to get the likelihood of each data point: we do this in
  % the log domain by summing to avoid underflow
  log_likelihood(k, :) =  nansum(log(raw_likelihood), 1);  % JR mod:  changed sum to nansum
end

% with uniform priors, the final log posterior is just the
% normalized log likelihood.  We would normally try to normalize
% within logs to avoid underflow, but that doesn't seem to be
% necessary here.

% (we don't technically need to change variable names here, but it
% makes it more correct if we want to add a prior)
log_posterior = log_likelihood + repmat(log(scratch_internal.prior), 1, nTimepoints);

% normalize the likelihoods in the real domain
acts = exp(log_posterior);
acts = acts ./ repmat(sum(acts,1), nConds, 1);

%%% end inlined test_gnb code %%%

class_scratch = scratch_internal;


% check whether we have already created a 3D MULTI_ACTS field to
% store the activations from each sphere's classifier. if we have,
% place our ACTS for the latest sphere in the right place in
% MULTI_ACTS. if we haven't, then create a MULTI_ACTS matrix of the
% right size.
%
% this complicated procedure is required because STATMAP_SEARCHLIGHT
% (our caller function) doesn't know that we need an ACTS matrix
% if ~isfield(scratch,'multi_acts')
%   % MULTI_ACTS field doesn't exist, so create it
%   scratch.multi_acts = NaN(nVox,nConds,nTimepoints);  
% end % created ACTS field
% scratch.multi_acts(scratch.v_counter,:,:) = acts;
% 
% % place our latest nConds x nTimepoints ACTS in the right part of the
% % MULTI_ACTS 3D nVox x nConds x nTimepoints matrix
% 
% perfmet_funct_hand = str2func(perfmet_funct);
% perfmet = perfmet_funct_hand(acts,regs3,class_scratch, ...
%                              perfmet_args);
% 
% % keep appending the the latest classifier scratch information from
% % the most recent time STATMAP_CLASSIFY was run
% if ~isfield(scratch,'class_scratch')
%   scratch.class_scratch = class_scratch;
% else
%   scratch.class_scratch(end+1) = class_scratch;
% end
% 
% val = perfmet.perf;

acts_vector = acts(1,:)-acts(2,:); %% JR mod:  save vector of signed activation differences
corrects_vector = zeros(1,length(acts_vector));
for cc = 1:length(acts_vector) 
    if acts_vector(cc)>0 % if activation favors classA
        if testtargs(1,cc)==1;  
            corrects_vector(cc)=1;
        end
    else % if activation favors classB
         if testtargs(2,cc)==1;  
            corrects_vector(cc)=1;
         end
    end
end

%corrects_vector = perfmet.corrects;  %% JR mod:  save vector of corrects (1=correct; 0=incorrect)


% if isnan(val)
%   warning('Val is NaN')
% end

