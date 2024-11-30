
%%
function [c,xstars] = prms_v1(pr, W, ystar, c0, xstars0) % 
% function [c,xstars] = prms_v1(pr, W, ystar) % 

% ntrial = pr.ntrial;
n_obs = pr.NN;
n_mov = pr.ntrial;
xstars_std = 0.2;
sp = pr.x_sp;

xstars0 = randn(n_obs,n_mov)*0.1/sqrt(n_mov);
c0 = randn(2,pr.NN)*0.1/sqrt(2);

xstars0_i = xstars0';
x0 = [c0;xstars0_i];

%% 
I = eye(size(W));
Q = lyap((W-I)',2*I);
[eig_vector,eig_value] = eig(Q); % eigenvectors are arranged in columns in the matrix eig_vector
a = fliplr(eig_vector(:,end-(n_obs-1):end)); % pick up 8 eigenvectors

%% Ocaml para.
n_muscles = 2;
m1_slice = pr.NN;

lambda_traj = 1/n_mov;
lambda_reg = 1/(m1_slice*n_muscles)*0.001;
gamma = I(1:m1_slice,:);

%%
% xstars (200*8) = top_obs (200*200) * xstars (200*8);
% z (200*9) = [xstar + sp (200*8), sp (200*1)];
% xstars_motor (160*9) = gamma (160*200) * z (200*9);

%% cost function
func = @(x) lambda_reg*norm(fun1(unpack(x)),'fro')^2 ...
    + lambda_traj * norm( fun1(unpack(x))*max(0,X_after(pr, W, fun2(unpack(x), 3)))-ystar{1},'fro')^2 ...
    + lambda_traj * norm( fun1(unpack(x))*max(0,X_after(pr, W, fun2(unpack(x), 4)))-ystar{2},'fro')^2 ...
    + lambda_traj * norm( fun1(unpack(x))*max(0,X_after(pr, W, fun2(unpack(x), 5)))-ystar{3},'fro')^2 ...
    + lambda_traj * norm( fun1(unpack(x))*max(0,X_after(pr, W, fun2(unpack(x), 6)))-ystar{4},'fro')^2 ...
    + lambda_traj * norm( fun1(unpack(x))*max(0,X_after(pr, W, fun2(unpack(x), 7)))-ystar{5},'fro')^2 ...
    + lambda_traj * norm( fun1(unpack(x))*max(0,X_after(pr, W, fun2(unpack(x), 8)))-ystar{6},'fro')^2 ...
    + lambda_traj * norm( fun1(unpack(x))*max(0,X_after(pr, W, fun2(unpack(x), 9)))-ystar{7},'fro')^2 ...
    + lambda_traj * norm( fun1(unpack(x))*max(0,X_after(pr, W, fun2(unpack(x), 10)))-ystar{8},'fro')^2 ...
    + lambda_traj * norm( fun1(unpack(x))*max(0,X_after(pr, W, fun2(unpack(x), 11)))-ystar{9},'fro')^2 ...
    + lambda_traj * norm( fun1(unpack(x))*max(0,X_after(pr, W, fun2(unpack(x), 12)))-ystar{10},'fro')^2 ...
    + lambda_traj * norm( fun1(unpack(x))*max(0,X_after(pr, W, fun2(unpack(x), 13)))-ystar{11},'fro')^2 ...
    + lambda_traj * norm( fun1(unpack(x))*max(0,X_after(pr, W, fun2(unpack(x), 14)))-ystar{12},'fro')^2;

    function y = fun1(x)
        y = x(1:2,:);
    end

    function y = fun2(x, n)
        y = transpose(x(n,:));
    end

    function xstar_c = unpack(xstar_c)
        C = xstar_c(1:2,:);
        x_star = transpose(xstar_c(3:end,:));
        
        top_obs = a;  % 观测度量的前n_obs个特征向量
        z = n_obs*n_mov*xstars_std^2;  % 用于标准化xstars
        
        xstar = top_obs*x_star;
        xstar = sqrt(z / norm(xstar)^2) * xstar; 
%         xstar = sqrt(z / sum(sum(xstar.^2))).*xstar;   % 标准化xstars

%         % nullspace
        xstars_motor = [xstar+sp,sp];
        xstars_motor = gamma * xstars_motor;
        xstars_motor_t = transpose(xstars_motor);
        h = linsolve(xstars_motor_t * xstars_motor, xstars_motor_t);
%         h = linsolve(xstars_motor.')*xstars_motor\(xstars_motor.');
%         disp('1')
        C = C - C*xstars_motor*h;

        xstar_c = [C;(xstar+sp)'];
          
    end

    function X_aft = X_after(pr, W, xstar)
        
        h = pr.h;
        tspan = linspace(0, pr.t_move,pr.n_move); % 20230526 constant MT
        
        initial_cond = xstar;
        Xafter = integrate_dynamics2(pr, W, h, initial_cond);
        Xafter = Xafter'; % transpose, so that the time is in the direction of the ROW vector
        
        X_aft = Xafter;
        
        % slove ODE 2
        function X = integrate_dynamics2(pr, W, h, initial_cond)
            [~,X] = ode45(@(t,X) rate_dynamics2_ode(t, X, pr, W, h), ...
                tspan, initial_cond);
        end
        
        % define ODE 2（after MO）
        function x_dot = rate_dynamics2_ode(t, X, pr, W, h)
            h_bar = pr.x_sp - W*max(0,pr.x_sp); % phi(x)=max(x,0)
            x_dot = 1/pr.tau* ...
                ( ...
                -X ...
                +W*max(0,X) ...
                +h_bar ...
                +interp1(tspan,h,t) ... % no u(t) after MO
                );
        end
        
    end

% options = optimoptions('fmincon','Display','iter','MaxIter',10000,'MaxFunEvals',10000); % 
% x = fmincon(func,x0,[],[],[],[],[],[],[],options);

options = optimoptions('fminunc','Display','iter','MaxIter',500000,'MaxFunEvals',500000); % 
x = fminunc(func,x0,options);

c = x(1:2,:);
xstars = transpose(x(3:end,:));

end

