% define parameters

function pr = def_Params()

% (1) Time parameters
pr.tfinal = 700; % delay
pr.n_timepoints = 700; % delay

pr.t_move = 300; % constant MT
pr.n_move = 300; 

pr.ntrial = 12;

% (2) Parameters of the M1 circuit model
pr.NN = 200;
pr.NE = 100;
pr.tau = 150;
pr.tau_rise = 50;
pr.tau_decay = 500;

pr.x_sp = normrnd(20,2,pr.NN,1);

% h
pr.tmove = 0;
pr.h = h_func(linspace(0,pr.t_move,pr.n_move), pr); 
    function output = h_func(t,pr)
        A = 1;
        output = A * (exp(-(t-pr.tmove)/pr.tau_decay)...
            - exp(-(t-pr.tmove)/pr.tau_rise));
        A = 5/max(output); % make sure h(t) has a maximum of 5, cong
        output = A * (exp(-(t-pr.tmove)/pr.tau_decay)...
            - exp(-(t-pr.tmove)/pr.tau_rise));
    end

% (3) Parameters of the arm mechanics and hand trajectories
pr.L1 = 30; % unit£ºcm
pr.L2 = 33;
pr.M1 = 1.4; % unit£ºkg
pr.M2 = 1.0;
pr.D2 = 16;
pr.I1 = 0.025E4;
pr.I2 = 0.045E4;
reach_angle = [15,-15:-30:-165,165:-30:45];
pr.theta_reach = reach_angle/180*pi;
pr.d_reach = 20;
pr.tau_reach = 120;

% (4) Parameters of the LQR algorithm
pr.lambda = 0.1;

% (5) Parameters of penalization of Schimel 2023
pr.a_effort = 5e-7;
pr.a_null = 1;

end
