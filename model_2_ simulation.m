
%% load
load('model_workspace.mat')

%% simulating unperturb condition

delay = 700;
repetition = 1;

X_CO = nan(repetition,pr.ntrial,pr.NN,pr.tfinal+pr.t_move);
X_INT = nan(repetition,pr.ntrial,pr.NN,pr.tfinal+pr.t_move);

Y_CO = nan(repetition,pr.ntrial,2,pr.t_move);
Y_INT = nan(repetition,pr.ntrial,2,pr.t_move);

u_CO = nan(repetition,pr.ntrial,pr.NN,pr.tfinal);
u_INT = nan(repetition,pr.ntrial,pr.NN,pr.tfinal);

motor_cost_CO = nan(repetition,pr.ntrial,pr.tfinal);
motor_cost_INT = nan(repetition,pr.ntrial,pr.tfinal);

r2_co = nan(repetition,pr.ntrial);
r2_int = nan(repetition,pr.ntrial);

rmse_co = nan(repetition,pr.ntrial);
rmse_int = nan(repetition,pr.ntrial);

for rep = 1:repetition
    parfor r_dir = 1:pr.ntrial
        
        
        initial_cond_0 = normrnd(20,0.15,pr.NN,1);

        h_pt = zeros(200,pr.tfinal/100);

        Xbefore_opt_CO = nan(pr.NN, pr.tfinal);
        u_opt_CO = nan(pr.NN, pr.tfinal);
        motor_cost_opt_CO = nan(1, pr.tfinal);
        
        Xbefore_opt_INT = nan(pr.NN, pr.tfinal);
        u_opt_INT = nan(pr.NN, pr.tfinal);
        motor_cost_opt_INT = nan(1, pr.tfinal);
        
        tic
        %% static
        u_opt = nan(pr.NN, 100);
        motor_cost_opt = nan(1, 100);
        r2 = 0;
        rmse = 0;
        
        state_start = r_dir;
        for i = 1:pr.tfinal/100
            
            state_t = r_dir;

            [Xbefore_opt_CO, u_opt, ~, ~, motor_cost_opt, ~] = ...
                cal_Xbefore(pr, Wsoc, CC1, Xstar, delay, initial_cond_0, state_start, state_t, i, Xbefore_opt_CO, h_pt, motor_cost_opt);

            u_opt_CO(:,(1+100*(i-1)):(1+100*(i-1)+99)) = u_opt;
            motor_cost_opt_CO(:,(1+100*(i-1)):(1+100*(i-1)+99)) = motor_cost_opt;
            
        end
        
        [X_opt_CO,Y_opt_CO,r2,rmse] = cal_XY(pr, Wsoc, CC1, Hand, Xbefore_opt_CO, Xafter, r_dir, state_start, 0);
        X_CO(rep,r_dir,:,:)=X_opt_CO;
        Y_CO(rep,r_dir,:,:)=Y_opt_CO;
        u_CO(rep,r_dir,:,:)=u_opt_CO;
        motor_cost_CO(rep,r_dir,:)=motor_cost_opt_CO;
        r2_co(rep,r_dir) = r2;
        rmse_co(rep,r_dir) = rmse;
        
        %% moving
        u_opt = nan(pr.NN, 100);
        motor_cost_opt = nan(1, 100);
        r2 = 0;
        rmse = 0;
        
        state_start = mod(r_dir+2,12);
        state_start = state_start + (state_start == 0) * 12;
        for i = 1:pr.tfinal/100

            state_t = mod(state_start+i-1,12);
            state_t = state_t + (state_t == 0) * 12;
            
            [Xbefore_opt_INT, u_opt, ~, ~, motor_cost_opt, ~] = ...
                cal_Xbefore(pr, Wsoc, CC1, Xstar, delay, initial_cond_0, state_start, state_t, i, Xbefore_opt_INT, h_pt, motor_cost_opt);
            
            u_opt_INT(:,(1+100*(i-1)):(1+100*(i-1)+99)) = u_opt;
            motor_cost_opt_INT(:,(1+100*(i-1)):(1+100*(i-1)+99)) = motor_cost_opt;
            
        end

        [X_opt_INT,Y_opt_INT,r2,rmse] = cal_XY(pr, Wsoc, CC1, Hand, Xbefore_opt_INT, Xafter, r_dir, state_start, 0);
        X_INT(rep,r_dir,:,:)=X_opt_INT;
        Y_INT(rep,r_dir,:,:)=Y_opt_INT;
        u_INT(rep,r_dir,:,:)=u_opt_INT;
        motor_cost_INT(rep,r_dir,:)=motor_cost_opt_INT;
        r2_int(rep,r_dir) = r2;
        rmse_int(rep,r_dir) = rmse;
        
        toc
        
    end
end

%% simulating perturbed condition
 
pr.lambda = 0.1;
repetition = 1;
delay = 700:20:1000;
pert_amp = [0 20 30 40];

mdelay = pr.t_move+delay(end);

X_CO = nan(length(pert_amp),pr.ntrial,repetition,length(delay),pr.NN,mdelay);
X_INT = nan(length(pert_amp),pr.ntrial,repetition,length(delay),pr.NN,mdelay);

Y_CO = nan(length(pert_amp),pr.ntrial,repetition,length(delay),2,pr.t_move);
Y_INT = nan(length(pert_amp),pr.ntrial,repetition,length(delay),2,pr.t_move);

u_CO = nan(length(pert_amp),pr.ntrial,repetition,length(delay),pr.NN,delay(end));
u_INT = nan(length(pert_amp),pr.ntrial,repetition,length(delay),pr.NN,delay(end));

motor_cost_CO = nan(length(pert_amp),pr.ntrial,repetition,length(delay),1,delay(end));
motor_cost_INT = nan(length(pert_amp),pr.ntrial,repetition,length(delay),1,delay(end));

r2_co = nan(length(pert_amp),pr.ntrial,repetition);
rmse_co = nan(length(pert_amp),pr.ntrial,repetition);

r2_int = nan(length(pert_amp),pr.ntrial,repetition);
rmse_int = nan(length(pert_amp),pr.ntrial,repetition);

for n = 1:length(pert_amp)
    pert_t = 500;
    pert_index = randsample(200,100);
    pt = ones(200,1)*pert_amp(n);
    pt(pert_index)=0;
    h_pt = zeros(200,pr.tfinal/100);
    h_pt(:, pert_t/100+1) = pt;
            
    for r_dir = 1:pr.ntrial
        for rep = 1:repetition
            
            initial_cond_0 = normrnd(20,1,pr.NN,1);

            parfor d = 1:length(delay)
                
                Xbefore_opt_CO = nan(pr.NN, delay(d));
                u_opt_CO = nan(pr.NN, delay(d));
                motor_cost_opt_CO = nan(1, delay(d));
                
                
                Xbefore_opt_INT = nan(pr.NN, delay(d));
                u_opt_INT = nan(pr.NN, delay(d));
                motor_cost_opt_INT = nan(1, delay(d));
                
                tic
              %% static
                motor_cost_opt = nan(1, 100);
                
                state_start = r_dir;
                for i = 1:pr.tfinal/100
                    
                    state_t = r_dir;
                    
                    [Xbefore_opt_CO, u_opt, ~, ~, motor_cost_opt, ~] = ...
                        cal_Xbefore(pr, Wsoc, CC1, Xstar, delay(d), initial_cond_0, state_start, state_t, i, Xbefore_opt_CO, h_pt, motor_cost_opt);
                    
                    if i < pr.tfinal/100
                        u_opt_CO(:,(1+100*(i-1)):(1+100*(i-1)+99)) = u_opt;
                        motor_cost_opt_CO(:,(1+100*(i-1)):(1+100*(i-1)+99)) = motor_cost_opt;
                    else
                        u_opt_CO(:,(1+100*(i-1)):(1+100*(i-1)+size(u_opt,2)-1)) = u_opt;
                        motor_cost_opt_CO(:,(1+100*(i-1)):(1+100*(i-1)+size(u_opt,2)-1)) = motor_cost_opt;
                    end
                    
                end 
                
                [X_opt_CO,Y_opt_CO,r2,rmse] = cal_XY(pr, Wsoc, CC1, Hand, Xbefore_opt_CO, Xafter, r_dir, state_start, 1);
                
                extra_cols = mdelay - size(X_opt_CO, 2);
                
                if extra_cols > 0
                    X_opt_CO = padarray(X_opt_CO, [0, extra_cols], NaN, 'post');
                    u_opt_CO = padarray(u_opt_CO, [0, extra_cols], NaN, 'post');
                    motor_cost_opt_CO = padarray(motor_cost_opt_CO, [0, extra_cols], NaN, 'post');
                end

                X_CO(n,r_dir,rep,d,:,:)=X_opt_CO;
                Y_CO(n,r_dir,rep,d,:,:)=Y_opt_CO;
                u_CO(n,r_dir,rep,d,:,:)=u_opt_CO;
                motor_cost_CO(n,r_dir,rep,d,:)=motor_cost_opt_CO;
                r2_co(n,r_dir,rep,d) = r2;
                rmse_co(n,r_dir,rep,d) = rmse;
                
                %% moving
                motor_cost_opt = nan(1, 100);
                
                state_start = mod(r_dir+2,12);
                state_start = state_start + (state_start == 0) * 12;
                for i = 1:pr.tfinal/100
                    
                    state_t = mod(state_start+i-1,12);
                    state_t = state_t + (state_t == 0) * 12;
                    
                    [Xbefore_opt_INT, u_opt, ~, ~, motor_cost_opt, ~] = ...
                        cal_Xbefore(pr, Wsoc, CC1, Xstar, delay(d), initial_cond_0, state_start, state_t, i, Xbefore_opt_INT, h_pt, motor_cost_opt);
                    
                    if i < pr.tfinal/100
                        u_opt_INT(:,(1+100*(i-1)):(1+100*(i-1)+99)) = u_opt;
                        motor_cost_opt_INT(:,(1+100*(i-1)):(1+100*(i-1)+99)) = motor_cost_opt;
                    else
                        u_opt_INT(:,(1+100*(i-1)):(1+100*(i-1)+size(u_opt,2)-1)) = u_opt;
                        motor_cost_opt_INT(:,(1+100*(i-1)):(1+100*(i-1)+size(u_opt,2)-1)) = motor_cost_opt;
                    end

                end

                [X_opt_INT,Y_opt_INT,r2,rmse] = cal_XY(pr, Wsoc, CC1, Hand, Xbefore_opt_INT, Xafter, r_dir, state_start, 1);
                
                extra_cols = mdelay - size(X_opt_INT, 2);
                
                if extra_cols > 0
                    X_opt_INT = padarray(X_opt_INT, [0, extra_cols], NaN, 'post');
                    u_opt_INT = padarray(u_opt_INT, [0, extra_cols], NaN, 'post');
                    motor_cost_opt_INT = padarray(motor_cost_opt_INT, [0, extra_cols], NaN, 'post');
                end
                
                X_INT(n,r_dir,rep,d,:,:)=X_opt_INT;
                Y_INT(n,r_dir,rep,d,:,:)=Y_opt_INT;
                u_INT(n,r_dir,rep,d,:,:)=u_opt_INT;
                motor_cost_INT(n,r_dir,rep,d,:)=motor_cost_opt_INT;
                r2_int(n,r_dir,rep,d) = r2;
                rmse_int(n,r_dir,rep,d) = rmse;

                toc
            end
        end
    end
end
 