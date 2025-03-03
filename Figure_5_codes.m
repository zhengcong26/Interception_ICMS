

%% Figure 5b TDR target location space

% % Load data
load('G_-120_tdr.mat')
% load('L_-120_tdr.mat')

speed = -120;
C = squeeze(median(J_TO_targ_t,1));

t = -2000:100:2000-100;
t0 = find(t==0);
t_start = 0;
t_end = 1200;
t1 = t0+t_start/100;
t2 = t0+t_end/100;
tt = t1:1:t2;
T = (t_start+100)/1000:0.1:(t_end+100)/1000;

target_loc_mov = zeros(length(T), size(C, 2)); % Initialize target locations
for k = 1:size(C,2)
    initial_angle = C(2,k);
    for i = 1:length(T)
        loc=initial_angle+speed*T(i);
        loc(loc<0,:)=loc(loc<0,:)+360; % Wrap angles below 0
        loc(loc>360,:)=loc(loc>360,:)-360; % Wrap angles above 360
        target_loc_mov(i,k) = loc;
    end
end
target_loc_static = C(1,:).*ones(size(target_loc_mov,1),size(target_loc_mov,2));

% Extract neural data before learning
beforeLearningN = cell(length(tt), 1);
for n = 1:length(tt)
    MO_N = J_TO_N_t(:, 1, :, tt(n)); % Neural data for condition 1
    beforeLearningN{n} = squeeze(permute(MO_N, [3, 1, 2])); % Reformat data
end

% Convert target angles to Cartesian coordinates (before learning)
A = C(1,:);
for i = 1:length(A)
    beforeLearningX(i,1) = cosd(A(i))*10.2;
    beforeLearningX(i,2) = sind(A(i))*10.2;
end

% Extract neural data after learning
afterLearningN = cell(length(tt), 1);
for n = 1:length(tt)
    MO_N = J_TO_N_t(:, 2, :, tt(n)); % Neural data for condition 2
    afterLearningN{n} = squeeze(permute(MO_N, [3, 1, 2])); % Reformat data
end

% Convert target angles to Cartesian coordinates (after learning)
A = target_loc_mov;
afterLearningX = zeros(size(A, 1), size(A, 2), 2);
for k = 1:size(A,1)
    for i = 1:size(A,2)
        afterLearningX(k,i,1) = cosd(A(k,i))*10.2;
        afterLearningX(k,i,2) = sind(A(k,i))*10.2;
    end
end

% static
N_mean = mean(cell2mat(beforeLearningN), 1);
beforeLearningN_center = beforeLearningN{1,1} - N_mean;
afterLearningN_center = cellfun(@(x) (x-N_mean), beforeLearningN, 'UniformOutput', false);

TDRoption = 1;
projectedStates_1 = [];
for i = 1:size(afterLearningN_center,1)
    [betaNeural2Behav,betaNeural2BehavOrth,projectedStates_1(i,:,:)] = buildTDRSubspace(beforeLearningX,beforeLearningX,beforeLearningN_center,afterLearningN_center{i,1},TDRoption);
end
plot_TDR_target(projectedStates_1, target_loc_static, 1)

% moving
N_mean = mean(cell2mat(afterLearningN), 1);
beforeLearningN_center = afterLearningN{1,1} - N_mean;
afterLearningN_center = cellfun(@(x) (x-N_mean), afterLearningN, 'UniformOutput', false);

TDRoption = 1;
projectedStates_2 = [];
for i = 1:size(afterLearningN_center,1)
    [betaNeural2Behav,betaNeural2BehavOrth,projectedStates_2(i,:,:)] = buildTDRSubspace(squeeze(afterLearningX(3,:,:)),squeeze(afterLearningX(i,:,:)),beforeLearningN_center,afterLearningN_center{i,1},TDRoption);
end
plot_TDR_target(projectedStates_2, target_loc_mov, 2)

%% Figure 5cd TDR reach direction space 

B = squeeze(mean(J_TO_X_t,1));

t = -2000:100:2000-100;
t0 = find(t==0);
t_start = -100;
t1 = t0+t_start/100;

MO_N = J_MO_N_t(:,1,:,t1);
beforeLearningN_1 = squeeze(permute(MO_N,[3,1,2]));
A = B(1,:);
beforeLearningX_1 = [cosd(A)' * 16, sind(A)' * 16]; % Convert angles to Cartesian coordinates

MO_N = J_MO_N_t(:,2,:,t1);
beforeLearningN_2 = squeeze(permute(MO_N,[3,1,2]));
A = B(2,:);
beforeLearningX_2 = [cosd(A)' * 16, sind(A)' * 16]; % Convert angles to Cartesian coordinates

t_start = 0;
t_end = 1200;
t1 = t0+t_start/100;
t2 = t0+t_end/100;
tt = t1:1:t2;

afterLearningN_1 = cell(length(tt), 1); % Initialize storage
for i = 1:length(tt)
    MO_N = J_TO_N_t(:,1,:,tt(i));
    afterLearningN_1{i,1} = squeeze(permute(MO_N,[3,1,2]));
end
afterLearningX_1 = beforeLearningX_1;

afterLearningN_2 = cell(length(tt), 1); % Initialize storage
for i = 1:length(tt)
    MO_N = J_TO_N_t(:,2,:,tt(i));
    afterLearningN_2{i,1} = squeeze(permute(MO_N,[3,1,2]));
end
afterLearningX_2 = beforeLearningX_2;

% static
N_mean = mean([beforeLearningN_1;cell2mat(afterLearningN_1)], 1);
beforeLearningN_center = beforeLearningN_1 - N_mean;
afterLearningN_center = cellfun(@(x) (x-N_mean), afterLearningN_1, 'UniformOutput', false);

TDRoption = 1;
projectedStates_1 = [];
for i = 1:size(afterLearningN_center)
    [betaNeural2Behav,betaNeural2BehavOrth,projectedStates_1(i,:,:)] = buildTDRSubspace(beforeLearningX_1,afterLearningX_1,beforeLearningN_center,afterLearningN_center{i,1},TDRoption);
end
plot_TDR_reach(projectedStates_1, target_loc_static, 1)

% moving
N_mean = mean([beforeLearningN_2;cell2mat(afterLearningN_2)], 1);
beforeLearningN_center = beforeLearningN_2 - N_mean;
afterLearningN_center = cellfun(@(x) (x-N_mean), afterLearningN_2, 'UniformOutput', false);

TDRoption = 1;
projectedStates_2 = [];
for i = 1:size(afterLearningN_center)
    [betaNeural2Behav,betaNeural2BehavOrth,projectedStates_2(i,:,:)] = buildTDRSubspace(beforeLearningX_2,afterLearningX_2,beforeLearningN_center,afterLearningN_center{i,1},TDRoption);
end
plot_TDR_reach(projectedStates_2, target_loc_mov, 2)


