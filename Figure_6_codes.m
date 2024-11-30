
%% Run the model
% model_1_optimization.m
% model_2_simulation.m

%% Figure 6c

load('model_trj.mat') 
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
plot_hand_all(Hand, Y_CO, colo(1,:))
plot_hand_all(Hand, Y_INT, colo(2,:))

%% Figure 6d

load('model_PSTH.mat') 
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
k = 197;
plot_PSTH_oneunit(X_CO,k,colo(1,:))
plot_PSTH_oneunit(X_INT,k,colo(2,:))

%% Figure 6e

load('model_motor_cost.mat') % for Figure 6e and 6f

motor_cost_CO_n = vertcat(motor_cost_CO{:});
motor_cost_INT_n = vertcat(motor_cost_INT{:});

% Define colors and line styles
colors = [0.06275, 0.30588, 0.5451; 0.80392, 0.33333, 0.33333]; % [CO; INT]
linestyle = {'-'};

figure; hold on;
plotMotorCost(motor_cost_INT_n, colors(2,:), linestyle{1}, 20);
plotMotorCost(motor_cost_CO_n, colors(1,:), linestyle{1}, 20);

% Customize figure
set(gca, 'yscale', 'log');
xlabel('Time (ms)');
ylabel('Prospective motor error');
box off;
set(gca, 'TickDir', 'out', 'FontSize', 11, 'color', 'none', 'TickLength', [0.015, 0.025], 'LineWidth', 1);
set(gcf, 'position', [30, 300, 1000, 500]);
xlim([-50 700]); ylim([-10 1e6]);
xticks(0:100:1000); yticks([1, 1e4]);

%% Figure 6f

load('model_motor_cost.mat') % for Figure 6e and 6f

motor_cost_CO_n = vertcat(motor_cost_CO{:});
motor_cost_INT_n = vertcat(motor_cost_INT{:});

%% Figure 6h

load('model_pert_r2_rmse.mat')

for k = 1:size(r2_co,1)
    for i = 1:size(r2_co,2)
        for j = 1:size(r2_co,3)
            try
                rt_1(k,i,j)=find(squeeze(r2_co(k,i,j,:))>0.9,1);
            catch
                rt_1(k,i,j)=size(r2_co,4); %
            end
        end
    end
end

for k = 1:size(r2_int,1)
    for i = 1:size(r2_int,2)
        for j = 1:size(r2_int,3)
            try
                rt_2(k,i,j)=find(squeeze(r2_int(k,i,j,:))>0.9,1);
            catch
                rt_2(k,i,j)=size(r2_int,4); %
            end
        end
    end
end

rt_co=(rt_1-1)*20;
rt_int=(rt_2-1)*20;

rt_co_1 = (reshape(rt_co, 4, 120))';
rt_int_1 = (reshape(rt_int, 4, 120))';

%% Figure 6i

load('model_u.mat')
 
it = [102 120 167];

u = u_INT{1,1};
plot_u(u,it)

u = u_CO{1,1};
plot_u(u,it)
