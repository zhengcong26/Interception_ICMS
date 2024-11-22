

%% Figure 5b

% Detour index for one example session 
load('G_0914_2msBin_25msGaussian_4s.mat')

detour_l = [];
num=1;
trialconstraint = 675;

% Select relevant trials based on criteria
stim_list_select = sortrows(stim_list(stim_list(:,4)>0&stim_list(:,5)==0&stim_list(:,1)<=trialconstraint,:),1);

CO = stim_list_select(stim_list_select(:,3)==0,:);
INT = cell(size(CO,1),1);

for s = 1:size(CO,1)
    condition_non = stim_list_select(stim_list_select(:,17)>CO(s,17)-20 & stim_list_select(:,17)<CO(s,17)+20 &...
        stim_list_select(:,3)~=0 & stim_list_select(:,4)>CO(s,4)-100 & stim_list_select(:,4)<CO(s,4)+200,:);
    INT{s,1} = [INT{s,1};condition_non];
end

% Remove duplicate rows in INT and sort
INT_mat = cell2mat(INT);
[~,ia,~] = unique(INT_mat(:,1),'rows');
INT_unique = INT_mat(ia,:);
INT_unique = sortrows(INT_unique,1);

tran = [CO;INT_unique];

% Assign hand locations based on angular information
hand_loc=zeros(size(tran,1),1);
hand_loc(0<tran(:,17) & tran(:,17)<90)=1;
hand_loc(90<tran(:,17) & tran(:,17)<180)=4;
hand_loc(180<tran(:,17) & tran(:,17)<270)=3;
hand_loc(270<tran(:,17) & tran(:,17)<360)=2;
hand_loc(360<tran(:,17))=1;

% Extract unique movement speeds
sp=unique(tran(:,3));
movingspeed = setdiff(sp,0);
sp=[0,movingspeed];

TargetBins = J_GO_bin;
Target_l=[];
for i = 1:2
    for j = 1:4
        tran_group_l{j,i}=tran(tran(:,4)>500&tran(:,3)==sp(i)&hand_loc==j,:);% long
        for k = 1:size(tran_group_l{j,i},1)
            Target_l{j,i}(k,:,:) = TargetBins{stim_list(:,1)==tran_group_l{j,i}(k,1),1};
            RT_l{j,i}(k,1) = stim_list(stim_list(:,1)==tran_group_l{j,i}(k,1),13);
        end
    end
end

% Define time windows for processing
t = -2000:2:2000;
t0 = find(t==0);
t1_l = t0-1000/2;
tt2_l = 300;
t2_l = t0+tt2_l/2; 
dt_l = t2_l-t1_l+1;

Target_l_t = cellfun(@(x) (x(:,t1_l:t2_l,:)),Target_l,'UniformOutput',0);

% Calculate detour index
t_l = -1000:2:tt2_l;
[detour_l(num,1,:),detour_l(num,2,:)] = detour_index(Target_l_t, t_l, 1000);

%% Figure 5cd TDR target location space

% Load data
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

target_loc = zeros(length(T), size(C, 2)); % Initialize target locations
for k = 1:size(C,2)
    initial_angle = C(2,k);
    for i = 1:length(T)
        loc=initial_angle+speed*T(i);
        loc(loc<0,:)=loc(loc<0,:)+360; % Wrap angles below 0
        loc(loc>360,:)=loc(loc>360,:)-360; % Wrap angles above 360
        target_loc(i,k) = loc;
    end
end

% Extract neural data before learning
beforeLearningN = cell(length(tt), 1);
for n = 1:length(tt)
    MO_N = J_TO_N_t(:, 1, :, tt(n)); % Neural data for condition 1
    beforeLearningN{n} = squeeze(permute(MO_N, [3, 1, 2])); % Reformat data
end

% Convert target angles to Cartesian coordinates (before learning)
A = C(1,:);
for i = 1:length(A)
    beforeLearningX(i,1) = cosd(A(i))*16;
    beforeLearningX(i,2) = sind(A(i))*16;
end

% Extract neural data after learning
afterLearningN = cell(length(tt), 1);
for n = 1:length(tt)
    MO_N = J_TO_N_t(:, 2, :, tt(n)); % Neural data for condition 2
    afterLearningN{n} = squeeze(permute(MO_N, [3, 1, 2])); % Reformat data
end

% Convert target angles to Cartesian coordinates (after learning)
A = target_loc;
afterLearningX = zeros(size(A, 1), size(A, 2), 2);
for k = 1:size(A,1)
    for i = 1:size(A,2)
        afterLearningX(k,i,1) = cosd(A(k,i))*16;
        afterLearningX(k,i,2) = sind(A(k,i))*16;
    end
end

% static
N_mean = mean(cell2mat(beforeLearningN), 1);
beforeLearningN_center = beforeLearningN{1,1} - N_mean;
afterLearningN_center = cellfun(@(x) (x-N_mean), beforeLearningN, 'UniformOutput', false);

TDRoption = 1;
projectedStates_1 = [];
for i = 1:size(afterLearningN_center)
    [betaNeural2Behav,betaNeural2BehavOrth,projectedStates_1(i,:,:)] = buildTDRSubspace(beforeLearningX,beforeLearningX,beforeLearningN_center,afterLearningN_center{i,1},TDRoption);
end
plot_TDR_target(projectedStates_1)

% moving
N_mean = mean(cell2mat(afterLearningN), 1);
beforeLearningN_center = afterLearningN{1,1} - N_mean;
afterLearningN_center = cellfun(@(x) (x-N_mean), afterLearningN, 'UniformOutput', false);

TDRoption = 1;
projectedStates_2 = [];
for i = 1:size(afterLearningN_center)
    [betaNeural2Behav,betaNeural2BehavOrth,projectedStates_2(i,:,:)] = buildTDRSubspace(squeeze(afterLearningX(3,:,:)),squeeze(afterLearningX(i,:,:)),beforeLearningN_center,afterLearningN_center{i,1},TDRoption);
end
plot_TDR_target(projectedStates_2)

%% Figure 5e TDR reach direction space 

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
plot_TDR_reach(projectedStates_1)

% moving
N_mean = mean([beforeLearningN_2;cell2mat(afterLearningN_2)], 1);
beforeLearningN_center = beforeLearningN_2 - N_mean;
afterLearningN_center = cellfun(@(x) (x-N_mean), afterLearningN_2, 'UniformOutput', false);

TDRoption = 1;
projectedStates_2 = [];
for i = 1:size(afterLearningN_center)
    [betaNeural2Behav,betaNeural2BehavOrth,projectedStates_2(i,:,:)] = buildTDRSubspace(beforeLearningX_2,afterLearningX_2,beforeLearningN_center,afterLearningN_center{i,1},TDRoption);
end
plot_TDR_reach(projectedStates_2)

%% Figure 5f decode
% detailed decoeding process refer to decoding_xy_cosine.m

load('G_0914_decoding_20ms.mat')
load('L_0815_decoding_20ms.mat')

t = -1000:20:500;
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;1	0.64706	0];
figure
set(gcf,'unit','centimeters','position',[10 5 5.8 3.6]);
hold on
plot_decoding(r2_CO_x,r2_CO_shf_x,h_CO_x,t,colo(1,:),0.01)
plot_decoding(r2_INT_x,r2_INT_shf_x,h_INT_x,t,colo(2,:),0.06)
plot_decoding(r2_TAR_x,r2_TAR_shf_x,h_TAR_x,t,colo(3,:),0.11)

figure
set(gcf,'unit','centimeters','position',[10 5 5.8 3.6]);
hold on
plot_decoding(r2_CO_y,r2_CO_shf_y,h_CO_y,t,colo(1,:),0.01)
plot_decoding(r2_INT_y,r2_INT_shf_y,h_INT_y,t,colo(2,:),0.06)
plot_decoding(r2_TAR_y,r2_TAR_shf_y,h_TAR_y,t,colo(3,:),0.11)

figure
set(gcf,'unit','centimeters','position',[10 5 5.8 3.6]);
hold on
plot_decoding(cosine_similarity_CO, cosine_similarity_CO_shf, h_similarity_CO,t,colo(1,:),0.01)
plot_decoding(cosine_similarity_INT, cosine_similarity_INT_shf, h_similarity_INT,t,colo(2,:),0.06)
plot_decoding(cosine_similarity_TAR, cosine_similarity_TAR_shf, h_similarity_TAR,t,colo(3,:),0.11)
