%% === ANALYSIS SECTION ===

%accuracy calculator for permutation testing

featureset = 1; %this is a name for what subject or dataset you've used; 
% it is needed because the data are stored in different indices in the res. 
% structure depending on the "subnum" you ran the code with 

% 
% %first, read in performance for every test
scrambledAcc= [];
% Loop through each of the 1000 cells
for cellIdx = 1:length(res.subj{1,featureset}.penalty.nVox.weights.iter)
    % Access the current cell's content using curly braces {}
    currentIterations = res.subj{1,featureset}.penalty.nVox.weights.iter{1, cellIdx}.iterations;
    
    % Loop through the iterations within the current cell
    for i = 1:length(currentIterations)
        % Grab the 'perf' value and append it to the vector
        % Use dot notation to access the field of the structure
        scrambledAcc = [scrambledAcc currentIterations(i).perf]; 
    end
end

% Now, count how many times the total perf value is repated in the
% scrambled vector of 1000 values

% --- Count the frequency of each unique value ---
% Method 
% groupcounts returns the counts (B) and the corresponding unique grouping variables (BG)
[perfsCounter.counts, perfsCounter.uniqueValues] = groupcounts(scrambledAcc'); 

% Combine them into a single variable (a matrix) with two columns
perfsData = [perfsCounter.uniqueValues, perfsCounter.counts]; % perfsData returns all the accuracy values and how many times they are repeated

% Sort the matrix by the second column (counts) in descending order to
% quickly identify which classfier accuracy value is repeated the most
% The -1 indicates sortiuniqueValuesng the first column (unique values) in descending order
sortedPerfsData = sortrows(perfsData, -1);

% --- Define your target value ---
targetVal = 0.525;  % change this to any accuracy threshold = the "real" tr1teo totasl_perf/accuracy value without scrambling

% --- Find all values >= target ---
higherIdx = sortedPerfsData(:,1) >= targetVal;

% --- Compute p-value ---
totalHigherOrEqual = sum(sortedPerfsData(higherIdx, 2));  % total count of those ≥ target
pValue = totalHigherOrEqual / length(scrambledAcc);  % divide by total tests across permutations

fprintf('Target value: %.3f | p-value = %.4f (based on %d / number of tests)\n', ...
        targetVal, pValue, totalHigherOrEqual);


%count in scrambledacc how many values = real result or better, then / 1000

