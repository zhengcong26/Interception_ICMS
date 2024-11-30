
function [Xbefore, u_opt, du_opt, u_optstar, motor_cost, energy_cost] = cal_Xbefore(pr, W, C, xstar, delay, initial_cond_0, init_final, init_t, count, Xbefore, h_pt, motor_cost_prev)
    % Main function to calculate system dynamics and cost optimization
    % Inputs:
    %   pr - Parameters structure
    %   W, C, xstar - Matrices for system dynamics and steady-state targets
    %   delay - Time delay in the system
    %   initial_cond_0 - Initial condition for the simulation
    %   init_final - Initial state indices at MO
    %   init_t - Initial state indices every 100 ms
    %   count - Iteration count
    %   Xbefore - Storage for previous dynamics results
    %   h_pt - External perturbations
    %   motor_cost_prev - Costs from previous iterations
    % Outputs:
    %   Xbefore - Updated system states
    %   u_opt, du_opt - Optimized control inputs
    %   u_optstar - Optimal steady-state input
    %   motor_cost, energy_cost - Calculated cost metrics

    % Constants
    I = eye(size(W));             % Identity matrix
    A = W - I;                    % Adjusted system matrix
    Q = lyap(A', C' * C);         % Lyapunov solution for Q
    B = ones(size(A)) / pr.lambda; % Control matrix (adjusted with lambda)
    P = are(A, B, Q);             % Algebraic Riccati Equation solution

    % Calculate Y (state transition matrix)
    Acl = A - P / pr.lambda;
    BB = (1 / pr.lambda)^2 * P * P;
    Y = lyap(Acl', BB);

    % Time span for integration
    delta_t = 100;                 % Target speed = 250 degrees/sec
    x_inter = pr.tfinal / delta_t; % Number of time intervals for integration
    ahead = pr.t_move / delta_t + 1; % predicted state for interception

    % Generate time spans for integration
    tspan = arrayfun(@(i) linspace(delta_t * (i - 1), delta_t * i, delta_t), 1:x_inter, 'UniformOutput', false);
    tspan{x_inter} = linspace(pr.tfinal - delta_t, delay, delay - (pr.tfinal - delta_t));

    % Set initial condition
    if count == 1
        initial_cond = initial_cond_0;   
    else
        initial_cond = Xbefore(:,(count-1)*100);
%         xbf = squeeze(Xbefore(count-1,:,:));
%         initial_cond = xbf(:,end);
    end
    
    % Calculate baseline inputs
    h_bar = pr.x_sp - W * max(0, pr.x_sp);

    %% run dynamic model
    naive = false;
    if init_final == init_t % static
        xstar_k = xstar(:,init_t);
        if count == 1 || motor_cost_prev(end) > 200
            xbf = integrate_dynamics(pr, W, xstar_k, P, initial_cond, tspan{count}, h_pt(:, count), h_bar);
        else
            xbf = integrate_dynamics(pr, W, xstar_k, [], initial_cond, tspan{count}, h_pt(:, count), h_bar);
            naive = true;
        end
    else % moving
        if count == 1 || motor_cost_prev(end) > 200
            xstar_k = xstar(:,init_t);
        else
            init_ahead = mod(init_t+ahead,12);
            init_ahead = init_ahead + (init_ahead == 0) * 12;
            xstar_k = xstar(:,init_ahead);
        end
        xbf = integrate_dynamics(pr, W, xstar_k, P, initial_cond, tspan{count}, h_pt(:, count), h_bar);
    end
    Xbefore(:,(1+100*(count-1)):(1+100*(count-1)+size(xbf,2)-1)) = xbf;
    
    %% Calculate motor and energy costs
    
    % Preallocate cost arrays
    n_steps = size(xbf, 2);
    motor_cost = zeros(1, n_steps);
    energy_cost = zeros(1, n_steps);
    
    dx = xbf-xstar_k;
    for k = 1:n_steps
        motor_cost(k) = dx(:,k)'* (P - pr.lambda * Y) * dx(:, k);
        energy_cost(k) = dx(:,k)' * Y * dx(:, k);
    end

    % Calculate input optimization
    u_optstar = xstar_k - W * max(0, xstar_k);
    if naive == false
        u_opt = u_optstar - h_bar - (1 / pr.lambda) * P * (xbf - xstar_k);
        du_opt = -(1 / pr.lambda) * P * (xbf - xstar_k);
    else
        u_opt = repmat(u_optstar - h_bar, 1, n_steps);
        du_opt = zeros(size(u_opt));
    end

end

% Integrate dynamics using ODE solver
function X = integrate_dynamics(pr, W, xstar_k, P, initial_cond, tspan, h_pt, h_bar)
    if isempty(P)
        P = 0;
    end
    [~, X] = ode45(@(t, X) rate_dynamics(t, X, pr, W, xstar_k, P, h_pt, h_bar), tspan, initial_cond);
    X = X';
end

% Define the system's rate dynamics
function x_dot = rate_dynamics(t, X, pr, W, xstar_k, P, h_pt, h_bar)
    x_dot = 1/pr.tau* ...
        ( ...
        -X ...
        +W*max(0,X) ...
        +h_bar ...
        +h_pt ...
        +xstar_k - W*max(0,xstar_k) - h_bar ... 
        -1/pr.lambda*P*(X-xstar_k) ...
        );
end
