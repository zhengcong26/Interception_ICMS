
function [r2_x_TAR,r2_x_TAR_shf,h_x_TAR,r2_y_TAR,r2_y_TAR_shf,h_y_TAR,cosine_dist,cosine_dist_shf,h_cosine,cosine_similarity,cosine_similarity_shf,h_similarity] = fitlm_decoding_target(x,initial_angle,initial_angle_shf,t,foldnum,bin)

kf = cvpartition(size(x, 1), 'KFold', foldnum); % kf.TrainSize % kf.TestSize

n = 1;
if bin == 2
    for k = 0:5:size(x,2)-10 % 10ms step
        tic
        for i = 1:kf.NumTestSets
            
            X = squeeze(mean(x(:,k+1:k+10,:),2)); % 2*10 = 20ms window size
            
            %% target location compute
            if k==0
                tt = t(k+1)/1000;
            else
                tt = t(k)/1000;
            end
            
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
            train_index = kf.training(i);
            X_train = X(train_index, :);
            y_x_train = y_x(train_index,:);
            y_y_train = y_y(train_index,:);
            y_shf_x_train = y_shf_x(train_index,:);
            y_shf_y_train = y_shf_y(train_index,:);
            
            test_index = kf.test(i);
            X_test = X(test_index, :);
            y_x_test = y_x(test_index,:);
            y_y_test = y_y(test_index,:);
            y_shf_x_test = y_shf_x(test_index,:);
            y_shf_y_test = y_shf_y(test_index,:);
            
            %% true
            %         model = fitlm(X_train, y_x_train);
            lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
            model = fitrlinear(X_train, y_x_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            model = fitrlinear(X_train, y_x_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);
            y_x_pred = predict(model, X_test);
            rmse_x_TAR(n,i) = sqrt(mean((y_x_test - y_x_pred).^2));
            r2_x_TAR(n,i) = 1 - sum((y_x_test - y_x_pred).^2) / sum((y_x_test - mean(y_x_test)).^2);
            
            %         model = fitlm(X_train, y_y_train);
            lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
            model = fitrlinear(X_train, y_y_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            model = fitrlinear(X_train, y_y_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);
            y_y_pred = predict(model, X_test);
            rmse_y_TAR(n,i) = sqrt(mean((y_y_test - y_y_pred).^2));
            r2_y_TAR(n,i) = 1 - sum((y_y_test - y_y_pred).^2) / sum((y_y_test - mean(y_y_test)).^2);
            
            % shuffle
            %         model = fitlm(X_train, y_shf_x_train);
            lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
            model = fitrlinear(X_train, y_shf_x_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            model = fitrlinear(X_train, y_shf_x_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);
            y_shf_x_pred = predict(model, X_test);
            rmse_x_TAR_shf(n,i) = sqrt(mean((y_shf_x_test - y_shf_x_pred).^2));
            r2_x_TAR_shf(n,i) = 1 - sum((y_shf_x_test - y_shf_x_pred).^2) / sum((y_shf_x_test - mean(y_shf_x_test)).^2);
            
            %         model = fitlm(X_train, y_shf_y_train);
            lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
            model = fitrlinear(X_train, y_shf_y_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            model = fitrlinear(X_train, y_shf_y_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);
            y_shf_y_pred = predict(model, X_test);
            rmse_y_TAR_shf(n,i) = sqrt(mean((y_shf_y_test - y_shf_y_pred).^2));
            r2_y_TAR_shf(n,i) = 1 - sum((y_shf_y_test - y_shf_y_pred).^2) / sum((y_shf_y_test - mean(y_shf_y_test)).^2);
        end
        
        %% statistics
        A = r2_x_TAR(n,:);
        B = r2_x_TAR_shf(n,:);
        [~,h_x_TAR(n,1)]=ranksum(A,B);
        
        A = r2_y_TAR(n,:);
        B = r2_y_TAR_shf(n,:);
        [~,h_y_TAR(n,1)]=ranksum(A,B);
        
        n = n+1;
        toc
    end
    
elseif bin == 20
    for k = 1:size(x,2)
        tic
        for i = 1:kf.NumTestSets
            
            X = squeeze(x(:,k,:));
            
            %% target location compute
            if k==0
                tt = t(k+1)/1000;
            else
                tt = t(k)/1000;
            end
            
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
            train_index = kf.training(i);
            X_train = X(train_index, :);
            y_x_train = y_x(train_index,:);
            y_y_train = y_y(train_index,:);
            y_shf_x_train = y_shf_x(train_index,:);
            y_shf_y_train = y_shf_y(train_index,:);
            
            test_index = kf.test(i);
            X_test = X(test_index, :);
            y_x_test = y_x(test_index,:);
            y_y_test = y_y(test_index,:);
            y_shf_x_test = y_shf_x(test_index,:);
            y_shf_y_test = y_shf_y(test_index,:);
            
            %% true
            %         model = fitlm(X_train, y_x_train);
            lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
            model = fitrlinear(X_train, y_x_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            model = fitrlinear(X_train, y_x_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);
            y_x_pred = predict(model, X_test);
            
            rmse_x_TAR(n,i) = sqrt(mean((y_x_test - y_x_pred).^2));
            r2_x_TAR(n,i) = 1 - sum((y_x_test - y_x_pred).^2) / sum((y_x_test - mean(y_x_test)).^2);
            
            %         model = fitlm(X_train, y_y_train);
            lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
            model = fitrlinear(X_train, y_y_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            model = fitrlinear(X_train, y_y_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);
            y_y_pred = predict(model, X_test);
            
            rmse_y_TAR(n,i) = sqrt(mean((y_y_test - y_y_pred).^2));
            r2_y_TAR(n,i) = 1 - sum((y_y_test - y_y_pred).^2) / sum((y_y_test - mean(y_y_test)).^2);
            
            %% true angle mse r2
            % Convert predicted coordinates to angles
            theta_pred = atan2d(y_y_pred, y_x_pred);  % Predicted angle in degrees
            theta_pred(theta_pred < 0) = theta_pred(theta_pred < 0) + 360;  % Ensure angle is in [0, 360]
            
            % Convert true test coordinates to angles
            theta_true = atan2d(y_y_test, y_x_test);  % True angle in degrees
            theta_true(theta_true < 0) = theta_true(theta_true < 0) + 360;  % Ensure angle is in [0, 360]
            
            % Calculate MSE and R2 for angles
            mse_angle(n, i) = mean((theta_true - theta_pred).^2);  % Mean squared error for angles
            r2_angle(n, i) = 1 - sum((theta_true - theta_pred).^2) / sum((theta_true - mean(theta_true)).^2);  % R2 for angles
            
            %% cosine distance
            [cosine_dist(n,i),cosine_similarity(n,i)] = cosine_distance(y_x_pred, y_y_pred, y_x_test, y_y_test);
            
            %% shuffle
            %         model = fitlm(X_train, y_shf_x_train);
            lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
            model = fitrlinear(X_train, y_shf_x_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            model = fitrlinear(X_train, y_shf_x_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);
            y_shf_x_pred = predict(model, X_test);
            rmse_x_TAR_shf(n,i) = sqrt(mean((y_shf_x_test - y_shf_x_pred).^2));
            r2_x_TAR_shf(n,i) = 1 - sum((y_shf_x_test - y_shf_x_pred).^2) / sum((y_shf_x_test - mean(y_shf_x_test)).^2);
            
            %         model = fitlm(X_train, y_shf_y_train);
            lambdas = logspace(-5, 5, 100); % 从 10^-4 到 10^4 之间的 100 个 lambda 值
            model = fitrlinear(X_train, y_shf_y_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', lambdas, 'CrossVal', 'on');
            minMSE = kfoldLoss(model); % 计算每个 lambda 下的平均均方误差
            [~, idx] = min(minMSE); % 找到最小的 MSE 所对应的索引
            bestLambda = lambdas(idx); % 最优的 lambda
            model = fitrlinear(X_train, y_shf_y_train, 'Learner', 'leastsquares', 'Regularization', 'ridge', 'Lambda', bestLambda);
            y_shf_y_pred = predict(model, X_test);
            rmse_y_TAR_shf(n,i) = sqrt(mean((y_shf_y_test - y_shf_y_pred).^2));
            r2_y_TAR_shf(n,i) = 1 - sum((y_shf_y_test - y_shf_y_pred).^2) / sum((y_shf_y_test - mean(y_shf_y_test)).^2);
            
            %% true angle mse r2
            
            theta_pred = atan2d(y_shf_y_pred, y_shf_x_pred);  % Predicted angle in degrees
            theta_pred(theta_pred < 0) = theta_pred(theta_pred < 0) + 360;  % Ensure angle is in [0, 360]
            
            theta_true = atan2d(y_shf_y_test, y_shf_x_test);  % True angle in degrees
            theta_true(theta_true < 0) = theta_true(theta_true < 0) + 360;  % Ensure angle is in [0, 360]
            
            mse_angle_shf(n, i) = mean((theta_true - theta_pred).^2);  % Mean squared error for angles
            r2_angle_shf(n, i) = 1 - sum((theta_true - theta_pred).^2) / sum((theta_true - mean(theta_true)).^2);  % R2 for angles
            
            %% cosine distance
            [cosine_dist_shf(n,i),cosine_similarity_shf(n,i)] = cosine_distance(y_shf_x_pred, y_shf_y_pred, y_shf_x_test, y_shf_y_test);
            
        end
        
        %% statistics
        A = r2_x_TAR(n,:);
        B = r2_x_TAR_shf(n,:);
        [~,h_x_TAR(n,1)]=ranksum(A,B);
        
        A = r2_y_TAR(n,:);
        B = r2_y_TAR_shf(n,:);
        [~,h_y_TAR(n,1)]=ranksum(A,B);
        
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