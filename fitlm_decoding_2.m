
function [r2_angle,r2_angle_shf,h_angle,cosine_dist, cosine_dist_shf,h_cosine, cosine_similarity, cosine_similarity_shf, h_similarity] = fitlm_decoding_2(x,y1,y2,y_shf_1,y_shf_2,foldnum,bin)

% 10折交叉验证
kf = cvpartition(size(x, 1), 'KFold', foldnum); % kf.TrainSize % kf.TestSize
n = 1;
if bin == 2
    for k = 0:5:size(x,2)-10 % 10ms step
        tic
        for i = 1:kf.NumTestSets
            
            X = squeeze(mean(x(:,k+1:k+10,:),2)); % 2*10 = 20ms window size
            
            % datasest
            train_index = kf.training(i);
            X_train = X(train_index, :);
            y_train_1 = y1(train_index,:);
            y_shf_train_1 = y_shf_1(train_index,:);
            
            test_index = kf.test(i);
            X_test = X(test_index, :);
            y_test_1 = y1(test_index,:);
            y_shf_test_1 = y_shf_1(test_index,:);
            
            %% true
            %         layers = [
            %             featureInputLayer(size(X_train, 2))
            %             fullyConnectedLayer(50)
            %             reluLayer
            %             fullyConnectedLayer(1)
            %             regressionLayer];
            %         options = trainingOptions('adam', 'MaxEpochs', 100);
            %         net = trainNetwork(X_train, y_train, layers, options);
            %         y_pred = predict(net, X_test);
            
            %          model = fitrgp(X_train, y_train);
            %         model = fitrlinear(X_train, y_train);
            %         model = fitrsvm(X_train, y_train, 'KernelFunction', 'linear');
            %         model = fitlm(X_train, y_train);
            
            
            lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
            model = fitrlinear(X_train, y_train_1, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            
            model = fitrlinear(X_train, y_train_1, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);
            y_pred = predict(model, X_test);
            rmse_CO(n,i) = sqrt(mean((y_test_1 - y_pred).^2));
            r2_CO(n,i) = 1 - sum((y_test_1 - y_pred).^2) / sum((y_test_1- mean(y_test_1)).^2);
            
            %% shuffle
            %         layers = [
            %             featureInputLayer(size(X_train, 2))
            %             fullyConnectedLayer(50)
            %             reluLayer
            %             fullyConnectedLayer(1)
            %             regressionLayer];
            %         options = trainingOptions('adam', 'MaxEpochs', 100);
            %         net = trainNetwork(X_train, y_shf_train, layers, options);
            %         y_pred = predict(net, X_test);
            
            %          model = fitrgp(X_train, y_shf_train);
            %         model = fitrlinear(X_train, y_shf_train);
            %         model = fitrsvm(X_train, y_shf_train, 'KernelFunction', 'linear');
            %         model = fitlm(X_train, y_shf_train);
            
            lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
            model = fitrlinear(X_train, y_shf_train_1, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            
            model = fitrlinear(X_train, y_shf_train_1, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);
            y_shf_pred_1 = predict(model, X_test);
            rmse_CO_shf(n,i) = sqrt(mean((y_shf_test_1 - y_shf_pred_1).^2));
            r2_CO_shf(n,i) = 1 - sum((y_shf_test_1 - y_shf_pred_1).^2) / sum((y_shf_test_1 - mean(y_shf_test_1)).^2);
            
        end
        
        %% statistics
        A = r2_CO(n,:);
        B = r2_CO_shf(n,:);
        [~,h_angle(n,1)]=ranksum(A,B);
        
        n = n+1;
        toc
    end
elseif bin == 20
    for k = 1:size(x,2)
        tic
        for i = 1:kf.NumTestSets
            
            X = squeeze(x(:,k,:)); % 2*10 = 20ms window size
            
            % datasest
            train_index = kf.training(i);
            X_train = X(train_index, :);
            y_train_1 = y1(train_index,:);
            y_shf_train_1 = y_shf_1(train_index,:);
            y_train_2 = y2(train_index,:);
            y_shf_train_2 = y_shf_2(train_index,:);
            
            test_index = kf.test(i);
            X_test = X(test_index, :);
            y_test_1 = y1(test_index,:);
            y_shf_test_1 = y_shf_1(test_index,:);
            y_test_2 = y2(test_index,:);
            y_shf_test_2 = y_shf_2(test_index,:);
            
            %% true y1
            lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
            model = fitrlinear(X_train, y_train_1, 'Learner', 'leastsquares',...
                'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            
            model = fitrlinear(X_train, y_train_1, 'Learner', 'leastsquares',...
                'Regularization', 'ridge', 'Lambda', bestLambda);
            y_pred_1 = predict(model, X_test);
            %% true y2
            model = fitrlinear(X_train, y_train_2, 'Learner', 'leastsquares',...
                'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            
            model = fitrlinear(X_train, y_train_2, 'Learner', 'leastsquares',...
                'Regularization', 'ridge', 'Lambda', bestLambda);
            y_pred_2 = predict(model, X_test);
            
            %% true angle mse r2
            % Convert predicted coordinates to angles
            theta_pred = atan2d(y_pred_2, y_pred_1);  % Predicted angle in degrees
            theta_pred(theta_pred < 0) = theta_pred(theta_pred < 0) + 360;  % Ensure angle is in [0, 360]
            
            % Convert true test coordinates to angles
            theta_true = atan2d(y_test_2, y_test_1);  % True angle in degrees
            theta_true(theta_true < 0) = theta_true(theta_true < 0) + 360;  % Ensure angle is in [0, 360]
            
            % Calculate MSE and R2 for angles
            mse_angle(n, i) = mean((theta_true - theta_pred).^2);  % Mean squared error for angles
            r2_angle(n, i) = 1 - sum((theta_true - theta_pred).^2) / sum((theta_true - mean(theta_true)).^2);  % R2 for angles
            
            %% cosine distance
            [cosine_dist(n,i),cosine_similarity(n,i)] = cosine_distance(y_pred_1, y_pred_2, y_test_1, y_test_2);
            
            %% shuffle y1
            
            lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
            model = fitrlinear(X_train, y_shf_train_1, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            
            model = fitrlinear(X_train, y_shf_train_1, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);
            y_shf_pred_1 = predict(model, X_test);
            %% shuffle y2
            
            lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
            model = fitrlinear(X_train, y_shf_train_2, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            
            model = fitrlinear(X_train, y_shf_train_2, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);
            y_shf_pred_2 = predict(model, X_test);
            
            %% shuffled angle mse r2
            % Convert predicted coordinates to angles
            theta_shf_pred = atan2d(y_shf_pred_2, y_shf_pred_1);  % Predicted angle in degrees
            theta_shf_pred(theta_shf_pred < 0) = theta_shf_pred(theta_shf_pred < 0) + 360;  % Ensure angle is in [0, 360]
            
            % Convert true test coordinates to angles
            theta_shf_true = atan2d(y_shf_test_2, y_shf_test_1);  % True angle in degrees
            theta_shf_true(theta_shf_true < 0) = theta_shf_true(theta_shf_true < 0) + 360;  % Ensure angle is in [0, 360]
            
            % Calculate MSE and R2 for angles
            mse_angle_shf(n, i) = mean((theta_shf_true - theta_shf_pred).^2);  % Mean squared error for angles
            r2_angle_shf(n, i) = 1 - sum((theta_shf_true - theta_shf_pred).^2) / sum((theta_shf_true - mean(theta_shf_true)).^2);  % R2 for angles
            
            %% cosine distance
            [cosine_dist_shf(n,i),cosine_similarity_shf(n,i)] = cosine_distance(y_shf_pred_1, y_shf_pred_2, y_shf_test_1, y_shf_test_2);
            
        end
        
        %% statistics
        A = r2_angle(n,:);
        B = r2_angle_shf(n,:);
        [~,h_angle(n,1)]=ranksum(A,B);
        
        A = cosine_dist(n,:);
        B = cosine_dist_shf(n,:);
        [~,h_cosine(n,1)]=ranksum(A,B);
        
        A = cosine_similarity(n,:);
        B = cosine_similarity_shf(n,:);
        [~,h_similarity(n,1)]=ranksum(A,B);
        
        n = n+1;
        toc
    end
end

end

% Function to compute cosine similarity
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
