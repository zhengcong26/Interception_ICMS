
%% Figure 6c

load('model_trj_PSTH.mat') % for Figure 6c and 6d

colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];

plot_hand_all(hand, Y_check_CO, colo(1,:))
plot_hand_all(hand, Y_check_INT, colo(2,:))


%% Figure 6d

colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];

k = 197;
plot_PSTH_oneunit(X_check_CO,k,colo(1,:))
plot_PSTH_oneunit(X_check_INT,k,colo(2,:))

%% Figure 6e

load('model_u.mat')
 
it = [102 120 167];

u = u_INT{1,1};
plot_u(u,it)

u = u_CO{1,1};
plot_u(u,it)

%% Figure 6f



%% Figure 6h



%% Figure 6i