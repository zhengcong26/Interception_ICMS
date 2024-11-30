function Xafter = cal_Xafter(pr, W, xstar)
% Calculate post-movement states for multiple trials
% Inputs:
%   pr - Parameters structure
%   W - Weight matrix, size (NN, NN)
%   xstar - Initial states, size (NN, ntrial)
% Outputs:
%   Xafter - Cell array of post-movement states, size (1, ntrial)
%            Each cell contains a matrix (NN, n_timepoints)

h_bar = pr.x_sp - W*max(0,pr.x_sp); % phi(x)=max(x,0)
tspan = linspace(0, pr.t_move,pr.n_move); % 20230526 constant MT

% Preallocate output cell array
Xafter = cell(1, pr.ntrial);
    
for i = 1:pr.ntrial
    initial_cond = xstar(:,i);
    Xafter{i} = integrate_dynamics2(pr, W, h_bar, pr.h, initial_cond, tspan);
    Xafter{i} = Xafter{i}'; % transpose, so that the time is in the direction of the ROW vector
end

end

% Function to solve ODE for a single trial
function X = integrate_dynamics2(pr, W, h_bar, h, initial_cond, tspan)
[~,X] = ode45(@(t,X) rate_dynamics2_ode(t, X, pr, W, h_bar, h, tspan), ...
    tspan, initial_cond);
end

% Define ODE for post-movement dynamics
function x_dot = rate_dynamics2_ode(t, X, pr, W, h_bar, h, tspan)
x_dot = 1/pr.tau* ...
    ( ...
    -X ...
    +W*max(0,X) ...
    +h_bar ...
    +interp1(tspan,h,t) ...
    );
end

