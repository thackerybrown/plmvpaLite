%% In Development. Do not use right now.
% WARNING - this can be extremely computationally intensive. 
% We recommend a jitter of no more than 3 units (jitter window = 2)
% More jitter could benefit enormously from parallelization

%code assumes rows = features. If your features are represented by columns,
%transpose the matrix first

%% Specify temporal parameters
seedtime = 2; % Integer referring to rows or columns in your matrix. What post-onset timepoint should we start with for every feature? This is in the units of matrix cells - so if each column is 1 BOLD TR then a 1 means we will be examining a time window from TR #1 forward

jitterwindow = 2; % Integer referring to additional rows or columns in your matrix. What time window will we allow the code to jitter over? This is in the units of matrix cells. A 1 means features can come from time window n and n+1. A 2 adds n+2. Etc.

percent_jitter = 50; % this is a very important parameter. What percent of the final feature set must include the seed time? It can DRAMATICALLY reduce classification time. Ask yourself this: if more than half of the feature times are NOT your seed time, shouldn't you use a different seed time?
%% Generate a random matrix to work with
RanM = randi(10,10,12);

inputM = RanM;
%% Generate time indices for each feature
% Generates all possible jitter window*feature # permutations with
% repetition.

% Again, this does not scale well computationally. At all. Consider the enormity of the
% permuation problem.

x = [seedtime:(seedtime+jitterwindow)];

k = size(inputM,1);

C = cell(k,1);

[C{:}] = ndgrid(x);

y = cellfun(@(x){x(:)},C);

y = [y{:}];


%% select only time indices that include your seed time
% to help with the scale of the classification problem, only use feature
% sets that include some % of the seed time

a = sum(y==seedtime,2); % how many seedtimes are in each row

targ_perm_ids = a>=(percent_jitter/100)*size(inputM,1);

targ_perms = y(targ_perm_ids,:);

%% apply the target permutations to the patterns

idx = bsxfun(@eq, cumsum(ones(size(inputM)),2), targ_perms(2,:)');

inputMadj = sum(inputM.*idx,2); %snip out a 
