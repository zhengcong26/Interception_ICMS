
function [r2,r2_shf,h] = decoding_xy(x,y,y_shf,foldnum)

% Perform ridge regression decoding with 10-fold cross-validation
% Inputs:
%   x: Predictor data (n_samples x n_features x n_models)
%   y: True target data
%   y_shf: Shuffled target data
%   foldnum: Number of folds for cross-validation
% Outputs:
%   r2: R-squared values for true targets
%   r2_shf: R-squared values for shuffled targets
%   h: Statistical test results comparing true vs shuffled

% Create cross-validation partitions
kf = cvpartition(size(x, 1), 'KFold', foldnum); % kf.TrainSize % kf.TestSize
n = 1;

for k = 1:size(x,2)
    tic
    for i = 1:kf.NumTestSets
        % Training and test split
        X = squeeze(x(:,k,:));
        train_index = kf.training(i);
        test_index = kf.test(i);
        
        % Training data
        X_train = X(train_index, :);
        y_train = y(train_index,:);
        y_shf_train = y_shf(train_index,:);
        
        % Test data
        X_test = X(test_index, :);
        y_test = y(test_index,:);
        y_shf_test = y_shf(test_index,:);
        
        % Ridge regression for true target
        lambdas = logspace(-5, 5, 100); 
        model = fitrlinear(X_train, y_train, 'Learner', 'leastsquares',...
            'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
        [~, idx] = min(kfoldLoss(model)); % Find best lambda
        bestLambda = lambdas(idx); 
        model = fitrlinear(X_train, y_train, 'Learner', 'leastsquares',...
            'Regularization', 'ridge', 'Lambda', bestLambda);
        y_pred = predict(model, X_test);

        rmse(n,i) = sqrt(mean((y_test - y_pred).^2));
        r2(n,i) = 1 - sum((y_test - y_pred).^2) / sum((y_test- mean(y_test)).^2);
        
        % Ridge regression for shuffled target
        lambdas = logspace(-5, 5, 100); 
        model = fitrlinear(X_train, y_shf_train, 'Learner', 'leastsquares', ...
            'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
        [~, idx] = min(kfoldLoss(model));
        bestLambda = lambdas(idx); 
        model = fitrlinear(X_train, y_shf_train, 'Learner', 'leastsquares', ...
            'Regularization', 'ridge', 'Lambda', bestLambda);
        y_shf_pred = predict(model, X_test);
        rmse_shf(n,i) = sqrt(mean((y_shf_test - y_shf_pred).^2));
        r2_shf(n,i) = 1 - sum((y_shf_test - y_shf_pred).^2) / sum((y_shf_test - mean(y_shf_test)).^2);
        
    end
    
    % statistics
    [~,h(n,1)]=ranksum(r2(n,:),r2_shf(n,:));
    
    n = n+1;
    toc
end


