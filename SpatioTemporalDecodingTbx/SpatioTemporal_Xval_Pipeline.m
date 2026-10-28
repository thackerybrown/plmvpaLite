%% This is pre-release, non-running code focused on modernizing the 
% classification branch to handle temporal variance across space in the 
% information contained by features.

% The logic: altough fMRI is inherently a slow measurement tool for two
% reasons - (1) TR speed is limited (typically 1s+), (2) even with fast
% TRs, the BOLD response to events itself is lagged, temporally protracted,
% and temporally smooth.

% However, there is temporal variability in the BOLD response across brain
% areas, such that a large ROI for classification forces a sub-optimal
% contstraint on a decoder if it assumes information peaks from pooling
% data at the same TR(s) across all voxels. Moreover, an increasing number
% of empirical papers show you can detect shifts in information using
% decoders and other tools (e.g., RSA) from one TR to the next, suggesting
% that changes in the BOLD response at the level of individual TRs can
% reflect meaningful shifts in processing state.

% Consequently, the purpose of this code is deploy spatiotemporal
% flattening with regularization to pattern data around a target time
% window relative to an event, rather than using an arbitrary fixed timing
% for the patterns. The code (when finalized) will generate importance maps
% for the spatial information, but also temporal importance information,
% revealing the temporal jitter across space which maximized information
% for the classification problem.

% The approach using reglarization over a higher-dimensional pattern which 
% includes both space and time tries to avoid a major mathematical hurdle: 
% combinatorial explosion. If you have 10,000 voxels and want to test 5 
% possible time shifts per voxel, there are $5^{10000}$ possible patterns. 
% We cannot reasonably expect to optimize the shift for each voxel 
% simultaneously through brute force.



% Assume the following variables from your existing pipeline:
% N = number of trials
% V = number of voxels in the ROI/mask
% S.TR_window = [3, 4, 5, 6]; % The target time window relative to event onset
% T = length(S.TR_window);
% labels = [N x 1] vector of condition labels for classification

% %% 1. Extract and Flatten Spatiotemporal Features
% X_3D = zeros(N, V, T); 
% 
% for i = 1:N
%     % EXTRACT: Get the [V x T] matrix for this specific trial
%     % (Replace this with your actual data extraction logic)
%     trial_data = get_trial_data(S, i); 
%     X_3D(i, :, :) = trial_data;
% end
% 
% % SAFETY NET: Create the traditional STATIC feature set (the mean across time)
% % By concatenating the traditional static features with the new spatiotemporal 
% % features, the L1 (LASSO) algorithm is forced to choose. If the static average 
% % is a more robust predictor for a specific voxel than any individual TR, 
% % LASSO will assign a weight to the average feature and zero out the individual TRs.
% % squeeze() ensures the result is [N x V]
% 
% X_static = squeeze(mean(X_3D, 3));
% 
% % FLATTEN: Reshape to [N x (V*T)]
% % A single trial now represents all spatial features across all allowed timepoints
% X_flat = reshape(X_3D, N, V * T);
% 
% % CONCATENATE: The new feature space is [N x (V * (T + 1))]
% X_combined = [X_flat, X_static];
% 
% %% 2. Train Regularized Classifier %****<<see comment below on adapting to a nested xvalidation to optimize lambda
% % Using fitclinear for SVM/Logistic Regression with L1 regularization
% % 'lasso' forces sparsity, pushing useless TR lags to exactly 0 weight
% Mdl = fitclinear(X_combined, labels, ...
%     'Learner', 'logistic', ... % or 'svm'
%     'Regularization', 'lasso');
% 
% % Extract the trained weights (Beta)
% % Mdl trained on X_combined...
% all_weights = Mdl.Beta;
% 
% % Split the weights based on their known indices
% weights_spatiotemporal = all_weights(1 : (V*T));
% weights_static = all_weights((V*T + 1) : end); % Size [V x 1]
% 
% %% 3. Generate the "Where" and "When" Importance Maps
% % Reshape weights back to [Voxels x Time]
% % focus on spatiotemporal weights to [V x T]
% W_3D = reshape(weights_spatiotemporal, V, T);
% 
% %note: We can now compare max(abs(W_3D), [], 2) against abs(weights_static)
% % to determine if a voxel was best classified by a specific TR or the overall mean.
% 
% % Find the maximum absolute weight across the time dimension (dim 2)
% [max_weights, max_time_indices] = max(abs(W_3D), [], 2);
% 
% % --- THE "WHERE" MAP ---
% % This shows the spatial importance of features, regardless of when they fired
% spatial_importance_map = max_weights; 
% 
% % --- THE "WHEN" MAP ---
% % Convert the matrix index back to the actual TR lag from your window
% temporal_lag_map = S.TR_window(max_time_indices)';
% 
% % Handle voxels dropped by the classifier (Weight = 0 across all TRs)
% % We don't want to assign a "When" to a voxel that wasn't used.
% dropped_voxels = (max_weights == 0);
% temporal_lag_map(dropped_voxels) = NaN; % Mark unused voxels as NaN
% 
% %% 4. Save to NIfTI / Brain Map
% % You can now map 'spatial_importance_map' and 'temporal_lag_map' 
% % back to your 3D brain space using your repository's existing affine headers.
% 
% 
% 
% % Assume X_flat [N x (V*T)] and labels [N x 1] are already created 
% % from Step 1 in the previous code block.
% 
% 
% 
% 
% %% 1. Setup Cross-Validation
% num_folds = 5; 
% % cvpartition ensures balanced class representation in each fold
% cvp = cvpartition(labels, 'KFold', num_folds); 
% 
% % Preallocate arrays to store predictions, accuracies, and weights
% predictions = zeros(N, 1);
% fold_accuracies = zeros(num_folds, 1);
% % We store weights for every fold to average them later
% fold_weights = zeros(V * T, num_folds); 
% 
% %% 2. Cross-Validation Loop
% for fold = 1:num_folds
% 
%     % --- Split Data ---
%     train_idx = training(cvp, fold);
%     test_idx = test(cvp, fold);
% 
%     X_train = X_flat(train_idx, :);
%     y_train = labels(train_idx);
% 
%     X_test = X_flat(test_idx, :);
%     y_test = labels(test_idx);
% 
%     % --- Train Classifier ---
%     % Train only on the data in the current training fold
%     Mdl = fitclinear(X_train, y_train, ...
%         'Learner', 'logistic', ... % Or 'svm'
%         'Regularization', 'lasso'); 
% 
%     % --- Test Classifier on Left-Out Data ---
%     % Predict labels for the unseen test data
%     y_pred = predict(Mdl, X_test);
%     predictions(test_idx) = y_pred;
% 
%     % Calculate and store accuracy for this specific fold
%     fold_accuracies(fold) = sum(y_pred == y_test) / length(y_test);
% 
%     % Store the beta weights learned in this fold
%     fold_weights(:, fold) = Mdl.Beta;
% 
%     fprintf('Fold %d Accuracy: %.2f%%\n', fold, fold_accuracies(fold) * 100);
% end
% 
% %% 3. Evaluate Overall Performance
% overall_accuracy = sum(predictions == labels) / N;
% fprintf('\nOverall Cross-Validated Accuracy: %.2f%%\n', overall_accuracy * 100);
% 
% %% 4. Generate Stable "Where" and "When" Maps
% % Average the weights across all folds to find the features that were 
% % consistently important, preventing fold-specific noise from dominating the maps.
% mean_weights = mean(fold_weights, 2);
% 
% % Reshape the averaged weights back to [Voxels x Time]
% W_3D = reshape(mean_weights, V, T);
% 
% % Find the maximum absolute weight across the time dimension (dim 2)
% [max_weights, max_time_indices] = max(abs(W_3D), [], 2);
% 
% % --- THE "WHERE" MAP ---
% spatial_importance_map = max_weights; 
% 
% % --- THE "WHEN" MAP ---
% temporal_lag_map = S.TR_window(max_time_indices)';
% 
% % Handle unused voxels
% % Because we averaged across folds, a voxel must have an average weight of 0 
% % (meaning it was dropped in all or most folds) to be considered unused.
% % You might want to use a small threshold (e.g., 1e-6) instead of strict 0.
% dropped_voxels = (max_weights < 1e-6); 
% temporal_lag_map(dropped_voxels) = NaN;
% 
% %Hyperparameter Tuning (Lambda): The code above uses the default regularization strength (Lambda) in fitclinear. In a rigorous pipeline, you should use nested cross-validation (passing 'OptimizeHyperparameters', 'Lambda' to fitclinear) to find the optimal penalty term on the training data without peeking at the test data.
% %Scaling/Z-Scoring: L1 regularization is highly sensitive to the scale of the features. Because BOLD signal amplitudes can vary across voxels, you should z-score the features (across trials) within your CV loop. Calculate the mean and standard deviation on X_train, apply them to X_train, and then use those same training parameters to scale X_test. Do not z-score X_flat before the CV loop, as this leaks information from the test set into the training phase.


% Assume X_3D [N x V x T] and labels [N x 1] are already extracted.
% S.TR_window contains the TR indices (e.g., [3, 4, 5, 6]).

[N, V, T] = size(X_3D);

%% 1. Build Combined Feature Space (Spatiotemporal + Static)
X_static   = squeeze(mean(X_3D, 3));        % [N x V] Traditional static average
X_flat     = reshape(X_3D, N, V * T);      % [N x (V*T)] Flattened temporal window
X_combined = [X_flat, X_static];           % [N x (V*(T+1))] Concatenated feature matrix

%% 2. Setup Outer Cross-Validation Loop
num_outer_folds = 5;
cvp_outer = cvpartition(labels, 'KFold', num_outer_folds);

predictions     = zeros(N, 1);
fold_accuracies = zeros(num_outer_folds, 1);
fold_weights    = zeros(V * (T + 1), num_outer_folds);
best_lambdas    = zeros(num_outer_folds, 1);

% Define Lambda grid for inner CV search
lambda_grid = logspace(-5, 1, 15); 

%% 3. Outer Cross-Validation Loop
for fold = 1:num_outer_folds
    
    % --- Step A: Outer Train/Test Split ---
    train_idx = training(cvp_outer, fold);
    test_idx  = test(cvp_outer, fold);
    
    X_train = X_combined(train_idx, :);
    y_train = labels(train_idx);
    
    X_test  = X_combined(test_idx, :);
    y_test  = labels(test_idx);
    
    % --- Step B: Z-Score Transformation (Fit on Train ONLY) ---
    mu_train    = mean(X_train, 1);
    sigma_train = std(X_train, [], 1);
    
    % Avoid division by zero for zero-variance voxels
    sigma_train(sigma_train == 0) = 1; 
    
    X_train_scaled = (X_train - mu_train) ./ sigma_train;
    X_test_scaled  = (X_test  - mu_train) ./ sigma_train;
    
    % --- Step C: Inner CV for Lambda Hyperparameter Tuning ---
    num_inner_folds = 4;
    cvp_inner = cvpartition(y_train, 'KFold', num_inner_folds);
    
    inner_cv_loss = zeros(length(lambda_grid), 1);
    
    % Sweep through candidate Lambdas inside training set
    for i = 1:length(lambda_grid)
        lam = lambda_grid(i);
        
        % Cross-validate model on inner folds
        Mdl_inner = fitclinear(X_train_scaled, y_train, ...
            'Learner', 'logistic', ...
            'Regularization', 'lasso', ...
            'Lambda', lam, ...
            'CVPartition', cvp_inner);
            
        inner_cv_loss(i) = kfoldLoss(Mdl_inner);
    end
    
    % Select optimal Lambda with minimum cross-validation loss
    [~, best_idx] = min(inner_cv_loss);
    best_lambda   = lambda_grid(best_idx);
    best_lambdas(fold) = best_lambda;
    
    % --- Step D: Train Final Outer Model with Best Lambda ---
    Mdl_final = fitclinear(X_train_scaled, y_train, ...
        'Learner', 'logistic', ...
        'Regularization', 'lasso', ...
        'Lambda', best_lambda);
        
    % --- Step E: Evaluate on Left-Out Test Fold ---
    y_pred = predict(Mdl_final, X_test_scaled);
    
    predictions(test_idx) = y_pred;
    fold_accuracies(fold) = mean(y_pred == y_test);
    fold_weights(:, fold) = Mdl_final.Beta;
    
    fprintf('Outer Fold %d | Best Lambda: %.5f | Accuracy: %.2f%%\n', ...
        fold, best_lambda, fold_accuracies(fold) * 100);
end

%% 4. Overall CV Metrics & Map Reconstruction
overall_accuracy = mean(predictions == labels);
fprintf('\nOverall Nested CV Accuracy: %.2f%%\n', overall_accuracy * 100);

% Average learned weights across all outer folds
mean_weights = mean(fold_weights, 2);

% Split spatiotemporal weights from static weights
spatiotemporal_weights = mean_weights(1 : (V * T));
static_weights         = mean_weights((V * T + 1) : end);

% Reshape spatiotemporal weights to [Voxels x Time]
W_3D = reshape(spatiotemporal_weights, V, T);

% Extract "Where" and "When" importance maps
[max_spatiotemporal_weights, max_time_indices] = max(abs(W_3D), [], 2);

% --- Spatial Importance ("Where") ---
% Compare highest spatiotemporal weight against static feature weight
spatial_importance_map = max(max_spatiotemporal_weights, abs(static_weights));

% --- Temporal Lag ("When") ---
temporal_lag_map = S.TR_window(max_time_indices)';

% Set lag to NaN if static feature was preferred OR voxel was zeroed out
preferred_static = abs(static_weights) >= max_spatiotemporal_weights;
uninformative    = spatial_importance_map < 1e-6;

temporal_lag_map(preferred_static | uninformative) = NaN;

%% exporting to SPM
% Read reference header
V_ref = spm_vol(reference_nifti);

% Define 'WHERE' Map Header
V_where = V_ref(1); % Take first volume structure
V_where.fname = 'spatial_importance_map.nii';
V_where.descrip = 'MVPA Spatial Importance Map';
V_where.dt = [spm_type('float32'), spm_platform('bigend')];
spm_write_vol(V_where, where_vol);

% Define 'WHEN' Map Header
V_when = V_ref(1);
V_when.fname = 'temporal_lag_map.nii';
V_when.descrip = 'MVPA Temporal Peak Lag (TR index)';
V_when.dt = [spm_type('float32'), spm_platform('bigend')];
spm_write_vol(V_when, when_vol);

%% standalone visualization of spatial and temporal importance maps (does not rely on SPM)
%% 1. File Paths & Inputs
t1_filename    = 'path/to/structural_T1.nii';
where_filename = 'spatial_importance_map.nii';
when_filename  = 'temporal_lag_map.nii';

%% 2. Load Volumes & Resample to T1 Space (Native MATLAB)
% Read Structural T1
t1_info = niftiinfo(t1_filename);
t1_vol  = double(niftiread(t1_info));

% Read Importance Maps
where_info = niftiinfo(where_filename);
where_vol  = double(niftiread(where_info));

when_info  = niftiinfo(when_filename);
when_vol   = double(niftiread(when_info));

% Resample functional maps to T1 volume space if dimensions differ
if ~isequal(size(t1_vol), size(where_vol))
    [Xq, Yq, Zq] = meshgrid(1:size(t1_vol,2), 1:size(t1_vol,1), 1:size(t1_vol,3));
    
    % Interpolate Where Map (Linear)
    where_vol = interp3(where_vol, Xq, Yq, Zq, 'linear', 0);
    
    % Interpolate When Map (Nearest-Neighbor to preserve integer TR values)
    when_vol  = interp3(when_vol, Xq, Yq, Zq, 'nearest', NaN);
end

%% 3. Set Display Slices (In Voxel Coordinates of T1)
% Choose representative axial slices across the volume
z_slices = round(linspace(size(t1_vol,3)*0.3, size(t1_vol,3)*0.7, 4));
num_slices = length(z_slices);

%% 4. Render "WHERE" Overlay (Spatial Importance)
figure('Name', 'Spatial Importance Map (Where)', 'Color', 'w', 'Position', [100, 100, 1200, 400]);
colormap_where = hot(256); % Hot colormap for weights

for i = 1:num_slices
    z = z_slices(i);
    subplot(1, num_slices, i);
    
    % Display Anatomical Background
    t1_slice = squeeze(t1_vol(:, :, z))';
    imagesc(t1_slice);
    colormap(gca, gray(256));
    axis image off; hold on;
    
    % Overlay Importance Weights
    where_slice = squeeze(where_vol(:, :, z))';
    
    % Threshold: Hide zero/negligible weights
    thresh = 0.05 * max(where_vol(:)); 
    alpha_mask = double(where_slice > thresh) * 0.7; % 70% opacity
    
    h_overlay = imagesc(where_slice);
    set(h_overlay, 'AlphaData', alpha_mask);
    
    % Apply hot colormap strictly to overlay layer using an auxiliary axes
    caxis([thresh, max(where_vol(:))]);
    title(sprintf('Slice Z = %d', z), 'FontSize', 10);
end

% Add Colorbar
cb1 = colorbar('Position', [0.92, 0.2, 0.015, 0.6]);
colormap(cb1, colormap_where);
ylabel(cb1, 'Classifier Weight (Absolute Max)');

%% 5. Render "WHEN" Overlay (Temporal Peak Lag)
figure('Name', 'Temporal Peak Lag Map (When)', 'Color', 'w', 'Position', [100, 550, 1200, 400]);

% Custom Discrete Colormap for TR Lags (e.g., Jet or Turbo)
unique_lags = unique(when_vol(~isnan(when_vol)));
num_lags    = length(unique_lags);
colormap_when = jet(max(num_lags, 1));

for i = 1:num_slices
    z = z_slices(i);
    subplot(1, num_slices, i);
    
    % Display Anatomical Background
    t1_slice = squeeze(t1_vol(:, :, z))';
    imagesc(t1_slice);
    colormap(gca, gray(256));
    axis image off; hold on;
    
    % Overlay Lag Index
    when_slice = squeeze(when_vol(:, :, z))';
    
    % Alpha mask hides NaNs (voxels where static feature won or weight = 0)
    alpha_mask = double(~isnan(when_slice)) * 0.75; 
    
    h_overlay = imagesc(when_slice);
    set(h_overlay, 'AlphaData', alpha_mask);
    
    if ~isempty(unique_lags)
        caxis([min(unique_lags), max(unique_lags)]);
    end
    title(sprintf('Slice Z = %d', z), 'FontSize', 10);
end

% Add Discrete Colorbar for Time Lags
if ~isempty(unique_lags)
    cb2 = colorbar('Position', [0.92, 0.2, 0.015, 0.6]);
    colormap(cb2, colormap_when);
    ylabel(cb2, 'Peak TR Lag (TR Index)');
    
    % Align ticks to discrete TR values
    cb2.Ticks = unique_lags;
    cb2.TickLabels = arrayfun(@(x) sprintf('TR %d', x), unique_lags, 'UniformOutput', false);
end