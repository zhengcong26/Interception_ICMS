

function [C,xstar] = Unpack(C,x_star,W,pr)
% C = xstar_c(1:2,:);
% x_star = transpose(xstar_c(3:end,:));

n_obs = pr.NN;
n_mov = pr.ntrial;
xstars_std = 0.2;
sp = pr.x_sp;
I = eye(200);
gamma = I(1:200,:);

Q = lyap((W-I)',2*I);
[eig_vector,~] = eig(Q); % eigenvectors are arranged in columns in the matrix eig_vector
a = fliplr(eig_vector(:,end-(n_obs-1):end)); % pick up 8 eigenvectors

top_obs = a;  % 观测度量的前n_obs个特征向量
z = n_obs*n_mov*xstars_std^2;  % 用于标准化xstars

xstar = top_obs*x_star;
xstar = sqrt(z / sum(sum(xstar.^2))).*xstar;   % 标准化xstars

% nullspace
xstars_motor = [xstar + sp,sp];
xstars_motor = gamma * xstars_motor;
xstars_motor_t = transpose(xstars_motor);
h = linsolve(xstars_motor_t * xstars_motor, xstars_motor_t);
%         h = linsolve(xstars_motor.')*xstars_motor\(xstars_motor.');
%         disp('1')
C = C - C*xstars_motor*h;
xstar = xstar+sp;

end