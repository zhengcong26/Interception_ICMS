
function [r2_x_TAR,r2_x_TAR_shf,h_x_TAR,r2_y_TAR,r2_y_TAR_shf,h_y_TAR,...
    cosine_dist,cosine_dist_shf,h_cosine,cosine_similarity,cosine_similarity_shf,h_similarity]...
    = decoding_target(x,initial_angle,initial_angle_shf,t,foldnum)

% Create cross-validation partition
kf = cvpartition(size(x, 1), 'KFold', foldnum); % kf.TrainSize % kf.TestSize

n = 1;
for k = 1:size(x,2)
    tic
    for i = 1:kf.NumTestSets
        X = squeeze(x(:,k,:));
        train_index = kf.training(i);
        test_index = kf.test(i);
        
        %% target location compute
        tt = t(k)/1000;
        
        y=rad2deg(-2*pi/3*tt)+initial_angle;
        y(y<0,:)=y(y<0,:)+360;
        y(y>360,:)=y(y>360,:)-360;
        y_x = cosd(y)*16;
        y_y = sind(y)*16;
        
        y_shf=rad2deg(-2*pi/3*tt)+initial_angle_shf;
        y_shf(y_shf<0,:)=y_shf(y_shf<0,:)+360;
        y_shf(y_shf>360,:)=y_shf(y_shf>360,:)-360;
        y_shf_x = cosd(y_shf)*16;
        y_shf_y = sind(y_shf)*16;
        
        %% datasest
        X_train = X(train_index, :);
        y_x_train = y_x(train_index,:);
        y_y_train = y_y(train_index,:);
        y_shf_x_train = y_shf_x(train_index,:);
        y_shf_y_train = y_shf_y(train_index,:);
        
        X_test = X(test_index, :);
        y_x_test = y_x(test_index,:);
        y_y_test = y_y(test_index,:);
        y_shf_x_test = y_shf_x(test_index,:);
        y_shf_y_test = y_shf_y(test_index,:);
        
        %% true
        % Predict x-coordinates
        y_x_pred = train_and_predict(X_train, y_x_train, X_test);
        rmse_x_TAR(n,i) = sqrt(mean((y_x_test - y_x_pred).^2));
        r2_x_TAR(n,i) = 1 - sum((y_x_test - y_x_pred).^2) / sum((y_x_test - mean(y_x_test)).^2);
        
        % Predict y-coordinates
        y_y_pred = train_and_predict(X_train, y_y_train, X_test);
        rmse_y_TAR(n,i) = sqrt(mean((y_y_test - y_y_pred).^2));
        r2_y_TAR(n,i) = 1 - sum((y_y_test - y_y_pred).^2) / sum((y_y_test - mean(y_y_test)).^2);
        
        %% true angle mse r2
        theta_pred = atan2d(y_y_pred, y_x_pred);  % Predicted angle in degrees
        theta_pred(theta_pred < 0) = theta_pred(theta_pred < 0) + 360;  % Ensure angle is in [0, 360]
        theta_true = atan2d(y_y_test, y_x_test);  % True angle in degrees
        theta_true(theta_true < 0) = theta_true(theta_true < 0) + 360;  % Ensure angle is in [0, 360]
        
        % Calculate MSE and R2 for angles
        mse_angle(n, i) = mean((theta_true - theta_pred).^2);  % Mean squared error for angles
        r2_angle(n, i) = 1 - sum((theta_true - theta_pred).^2) / sum((theta_true - mean(theta_true)).^2);  % R2 for angles
        
        [cosine_dist(n,i),cosine_similarity(n,i)] = cosine_distance(y_x_pred, y_y_pred, y_x_test, y_y_test);
        
        %% shuffle
        y_shf_x_pred = train_and_predict(X_train, y_shf_x_train, X_test);
        rmse_x_TAR_shf(n,i) = sqrt(mean((y_shf_x_test - y_shf_x_pred).^2));
        r2_x_TAR_shf(n,i) = 1 - sum((y_shf_x_test - y_shf_x_pred).^2) / sum((y_shf_x_test - mean(y_shf_x_test)).^2);
        
        y_shf_y_pred = train_and_predict(X_train, y_shf_y_train, X_test);
        rmse_y_TAR_shf(n,i) = sqrt(mean((y_shf_y_test - y_shf_y_pred).^2));
        r2_y_TAR_shf(n,i) = 1 - sum((y_shf_y_test - y_shf_y_pred).^2) / sum((y_shf_y_test - mean(y_shf_y_test)).^2);
        
        %% true angle mse r2
        theta_pred = atan2d(y_shf_y_pred, y_shf_x_pred);  % Predicted angle in degrees
        theta_pred(theta_pred < 0) = theta_pred(theta_pred < 0) + 360;  % Ensure angle is in [0, 360]
        theta_true = atan2d(y_shf_y_test, y_shf_x_test);  % True angle in degrees
        theta_true(theta_true < 0) = theta_true(theta_true < 0) + 360;  % Ensure angle is in [0, 360]
        
        mse_angle_shf(n, i) = mean((theta_true - theta_pred).^2);  % Mean squared error for angles
        r2_angle_shf(n, i) = 1 - sum((theta_true - theta_pred).^2) / sum((theta_true - mean(theta_true)).^2);  % R2 for angles
        
        [cosine_dist_shf(n,i),cosine_similarity_shf(n,i)] = ...
            cosine_distance(y_shf_x_pred, y_shf_y_pred, y_shf_x_test, y_shf_y_test);
        
    end
    
    %% statistics
    h_x_TAR(n, 1) = ranksum(r2_x_TAR(n, :), r2_x_TAR_shf(n, :));
    h_y_TAR(n, 1) = ranksum(r2_y_TAR(n, :), r2_y_TAR_shf(n, :));
    h_cosine(n, 1) = ranksum(cosine_dist(n, :), cosine_dist_shf(n, :));
    h_similarity(n, 1) = ranksum(cosine_similarity(n, :), cosine_similarity_shf(n, :));
    
    n = n+1;
    toc
end
end

%% Helper function to train and predict with ridge regression
function y_pred = train_and_predict(X_train, y_train, X_test)
lambdas = logspace(-5, 5, 100); % Range of regularization parameters
model = fitrlinear(X_train, y_train, 'Learner', 'leastsquares', ...
    'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
[~, idx] = min(kfoldLoss(model)); % 找到最小的 MSE 所对应的索引
bestLambda = lambdas(idx); % 最优的 lambda
final_model = fitrlinear(X_train, y_train, 'Learner', 'leastsquares', ...
    'Regularization', 'ridge', 'Lambda', bestLambda);
y_pred = predict(final_model, X_test);
end

%% Helper Function to compute cosine similarity and distance
function [mean_cosine_dist,mean_cosine_similarity] = cosine_distance(pred_x, pred_y, true_x, true_y)

for i = 1:size(pred_x,1)
    % Convert predicted and true coordinates into vectors
    pred_vector = [pred_x(i), pred_y(i)];
    true_vector = [true_x(i), true_y(i)];
    
    % Normalize the vectors (calculate their magnitudes)
    pred_norm = norm(pred_vector);
    true_norm = norm(true_vector);
    
    % Compute the cosine similarity
    dot_product = dot(pred_vector, true_vector);
    cosine_similarity(i) = dot_product / (pred_norm * true_norm);
    
    % Convert to cosine distance (1 - cosine similarity)
    cosine_dist(i) = 1 - cosine_similarity(i);
end

mean_cosine_dist = mean(cosine_dist);
mean_cosine_similarity = mean(cosine_similarity);

end
