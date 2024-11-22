
function [r2_angle,r2_angle_shf,h_angle,cosine_dist, cosine_dist_shf,h_cosine, cosine_similarity, cosine_similarity_shf, h_similarity] = decoding_cosine(x,y1,y2,y_shf_1,y_shf_2,foldnum)
% Decoding with ridge regression using cosine distance and angle metrics
% Inputs:
%   x: Predictor data
%   y1, y2: True target data
%   y_shf_1, y_shf_2: Shuffled target data
%   foldnum: Number of folds for cross-validation
% Outputs:
%   r2_angle, r2_angle_shf: R-squared for angles (true vs shuffled)
%   h_angle: Statistical test results for angles
%   cosine_dist, cosine_dist_shf: Cosine distances (true vs shuffled)
%   h_cosine: Statistical test results for cosine distances
%   cosine_similarity, cosine_similarity_shf: Cosine similarities (true vs shuffled)
%   h_similarity: Statistical test results for cosine similarities

kf = cvpartition(size(x, 1), 'KFold', foldnum); % kf.TrainSize % kf.TestSize
n = 1;

for k = 1:size(x,2)
    tic
    for i = 1:kf.NumTestSets
        % Split data into training and test sets
        X = squeeze(x(:,k,:));
        train_index = kf.training(i);
        test_index = kf.test(i);
        
        X_train = X(train_index, :);
        y_train_1 = y1(train_index,:);
        y_shf_train_1 = y_shf_1(train_index,:);
        y_train_2 = y2(train_index,:);
        y_shf_train_2 = y_shf_2(train_index,:);
        
        X_test = X(test_index, :);
        y_test_1 = y1(test_index,:);
        y_shf_test_1 = y_shf_1(test_index,:);
        y_test_2 = y2(test_index,:);
        y_shf_test_2 = y_shf_2(test_index,:);
        
        %% true y1 and y2 prediction
        lambdas = logspace(-5, 5, 100); 
        model = fitrlinear(X_train, y_train_1, 'Learner', 'leastsquares',...
            'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
        [~, idx] = min(kfoldLoss(model)); 
        bestLambda = lambdas(idx); 
        model = fitrlinear(X_train, y_train_1, 'Learner', 'leastsquares',...
            'Regularization', 'ridge', 'Lambda', bestLambda);
        y_pred_1 = predict(model, X_test);
        
        
        model = fitrlinear(X_train, y_train_2, 'Learner', 'leastsquares',...
            'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
        [~, idx] = min(kfoldLoss(model)); 
        bestLambda = lambdas(idx);
        model = fitrlinear(X_train, y_train_2, 'Learner', 'leastsquares',...
            'Regularization', 'ridge', 'Lambda', bestLambda);
        y_pred_2 = predict(model, X_test);
        
        %% true angle mse r2
        theta_pred = atan2d(y_pred_2, y_pred_1);  % Predicted angle in degrees
        theta_pred(theta_pred < 0) = theta_pred(theta_pred < 0) + 360;  % Ensure angle is in [0, 360]
        theta_true = atan2d(y_test_2, y_test_1);  % True angle in degrees
        theta_true(theta_true < 0) = theta_true(theta_true < 0) + 360;  % Ensure angle is in [0, 360]
        
        mse_angle(n, i) = mean((theta_true - theta_pred).^2);  % Mean squared error for angles
        r2_angle(n, i) = 1 - sum((theta_true - theta_pred).^2) / sum((theta_true - mean(theta_true)).^2);  % R2 for angles
        
        %% cosine distance
        [cosine_dist(n,i),cosine_similarity(n,i)] = cosine_distance(y_pred_1, y_pred_2, y_test_1, y_test_2);
        
        %% shuffle y1 and y2 prediction
        lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
        model = fitrlinear(X_train, y_shf_train_1, 'Learner', 'leastsquares', ...
            'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
        [~, idx] = min(kfoldLoss(model)); % 找到最小的 MSE 所对应的索引
        bestLambda = lambdas(idx); % 最优的 lambda
        model = fitrlinear(X_train, y_shf_train_1, 'Learner', 'leastsquares',...
            'Regularization', 'ridge', 'Lambda', bestLambda);
        y_shf_pred_1 = predict(model, X_test);

        lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
        model = fitrlinear(X_train, y_shf_train_2, 'Learner', 'leastsquares',...
            'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
        [~, idx] = min(kfoldLoss(model)); % 找到最小的 MSE 所对应的索引
        bestLambda = lambdas(idx); % 最优的 lambda
        model = fitrlinear(X_train, y_shf_train_2, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);
        y_shf_pred_2 = predict(model, X_test);
        
        %% shuffled angle mse r2
        theta_shf_pred = atan2d(y_shf_pred_2, y_shf_pred_1);  % Predicted angle in degrees
        theta_shf_pred(theta_shf_pred < 0) = theta_shf_pred(theta_shf_pred < 0) + 360;  % Ensure angle is in [0, 360]
        theta_shf_true = atan2d(y_shf_test_2, y_shf_test_1);  % True angle in degrees
        theta_shf_true(theta_shf_true < 0) = theta_shf_true(theta_shf_true < 0) + 360;  % Ensure angle is in [0, 360]
        
        mse_angle_shf(n, i) = mean((theta_shf_true - theta_shf_pred).^2);  % Mean squared error for angles
        r2_angle_shf(n, i) = 1 - sum((theta_shf_true - theta_shf_pred).^2) / sum((theta_shf_true - mean(theta_shf_true)).^2);  % R2 for angles
        
        %% cosine distance
        [cosine_dist_shf(n,i),cosine_similarity_shf(n,i)] = cosine_distance(y_shf_pred_1, y_shf_pred_2, y_shf_test_1, y_shf_test_2);
        
    end
    
    %% statistics
    [~,h_angle(n,1)]=ranksum(r2_angle(n,:),r2_angle_shf(n,:));
    [~,h_cosine(n,1)]=ranksum(cosine_dist(n,:),cosine_dist_shf(n,:));
    [~,h_similarity(n,1)]=ranksum(cosine_similarity(n,:),cosine_similarity_shf(n,:));
    
    n = n+1;
    toc
end
end

function [mean_cosine_dist,mean_cosine_similarity] = cosine_distance(pred_x, pred_y, true_x, true_y)

for i = 1:size(pred_x,1)
    pred_vector = [pred_x(i), pred_y(i)];
    true_vector = [true_x(i), true_y(i)];
    pred_norm = norm(pred_vector);
    true_norm = norm(true_vector);
    dot_product = dot(pred_vector, true_vector);
    cosine_similarity(i) = dot_product / (pred_norm * true_norm);
    cosine_dist(i) = 1 - cosine_similarity(i);
end
mean_cosine_dist = mean(cosine_dist);
mean_cosine_similarity = mean(cosine_similarity);
end
