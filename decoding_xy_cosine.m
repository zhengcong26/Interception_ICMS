
clc;clear;close all;

load('G_0914_GO_decoding.mat')
stim_list_select = stim_list(stim_list(:,1)<=580,:);

% load('L_0815_GO_decoding.mat')
% stim_list_select = stim_list(stim_list(:,1)<=675,:);

t = -2000:20:2000;
t0 = find(t==0);
t1 = t0-1000/20;
t2 = t0+500/20;

GO_t = GO(:,t1:t2,:);

CO = stim_list_select(stim_list_select(:,3)==0,:);
INT = cell(size(CO,1),1);

for s = 1:size(CO,1)
    condition_non = stim_list_select(stim_list_select(:,17)>CO(s,17)-20 & stim_list_select(:,17)<CO(s,17)+20 &...
        stim_list_select(:,3)~=0 & stim_list_select(:,4)>CO(s,4)-100 & stim_list_select(:,4)<CO(s,4)+200,:);
    INT{s,1} = [INT{s,1};condition_non];
end

INT_mat = cell2mat(INT);

[~,ia,~] = unique(INT_mat(:,1),'rows');
INT_unique = INT_mat(ia,:);
INT_unique = sortrows(INT_unique,1);

tran = [CO;INT_unique];

%% CO regress
tran_all = tran(tran(:,3)==0&tran(:,4)>500,:); 

x = [];
for k = 1:size(tran_all,1)
    x(k,:,:) = GO_t(stim_list(:,1)==tran_all(k,1),:,:);
end

y1 = tran_all(:,20);
y2 = tran_all(:,21);
p = randperm(size(y1,1));
y_shf_1 = y1(p,:);
y_shf_2 = y2(p,:);

[r2_CO_x,r2_CO_shf_x,h_CO_x] = decoding_xy(x,y1,y_shf_1,10);
[r2_CO_y,r2_CO_shf_y,h_CO_y] = decoding_xy(x,y2,y_shf_2,10);
[r2_CO_angle,r2_CO_angle_shf,h_CO_angle,cosine_dist_CO, cosine_dist_CO_shf,h_cosine_CO, cosine_similarity_CO, cosine_similarity_CO_shf, h_similarity_CO]...
    = decoding_cosine(x,y1,y2,y_shf_1,y_shf_2,10);

t = -1000:20:500;
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;1	0.64706	0];
figure
set(gcf,'position',[100,100,360,250]);
hold on
plot_decoding(r2_CO_x,r2_CO_shf_x,h_CO_x,t,colo(1,:),0.01)
plot_decoding(r2_CO_y,r2_CO_shf_y,h_CO_y,t,colo(1,:),0.01)
plot_decoding(r2_CO_angle,r2_CO_angle_shf,h_CO_angle,t,colo(1,:),0.01)
plot_decoding(cosine_dist_CO,cosine_dist_CO_shf,h_cosine_CO,t,colo(1,:),0.01)
plot_decoding(cosine_similarity_CO, cosine_similarity_CO_shf, h_similarity_CO,t,colo(1,:),0.01)

%% INT regress
tran_all = tran(tran(:,3)~=0&tran(:,4)>500,:); 

x = [];
for k = 1:size(tran_all,1)
    x(k,:,:) = GO_t(stim_list(:,1)==tran_all(k,1),:,:);
end

y1 = tran_all(:,20);
y2 = tran_all(:,21);
p = randperm(size(y1,1));
y_shf_1 = y1(p,:);
y_shf_2 = y2(p,:);

[r2_INT_x,r2_INT_shf_x,h_INT_x] = decoding_xy(x,y1,y_shf_1,10);
[r2_INT_y,r2_INT_shf_y,h_INT_y] = decoding_xy(x,y2,y_shf_2,10);
[r2_INT_angle,r2_INT_angle_shf,h_INT_angle,cosine_dist_INT, cosine_dist_INT_shf,h_cosine_INT, cosine_similarity_INT, cosine_similarity_INT_shf, h_similarity_INT]...
    = decoding_cosine(x,y1,y2,y_shf_1,y_shf_2,10);

t = -1000:20:500;
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;1	0.64706	0];
figure
set(gcf,'position',[100,100,360,250]);
hold on
plot_decoding(r2_INT_x,r2_INT_shf_x,h_INT_x,t,colo(2,:),0.01)
plot_decoding(r2_INT_y,r2_INT_shf_y,h_INT_y,t,colo(2,:),0.01)
plot_decoding(r2_INT_angle,r2_INT_angle_shf,h_INT_angle,t,colo(1,:),0.01)
plot_decoding(cosine_dist_INT,cosine_dist_INT_shf,h_cosine_INT,t,colo(1,:),0.01)
plot_decoding(cosine_similarity_INT, cosine_similarity_INT_shf, h_similarity_INT,t,colo(1,:),0.01)

%% INT target regress
tran_all = tran(tran(:,3)~=0&tran(:,4)>500,:); % CO vs. INT delay

x = [];
for k = 1:size(tran_all,1)
    x(k,:,:) = GO_t(stim_list(:,1)==tran_all(k,1),:,:);
end

initial_angle = tran_all(:,15);
p = randperm(size(initial_angle,1));
initial_angle_shf = initial_angle(p,:);

[r2_TAR_x,r2_TAR_shf_x,h_TAR_x,r2_TAR_y,r2_TAR_shf_y,h_TAR_y,cosine_dist_TAR, cosine_dist_TAR_shf,h_cosine_TAR, cosine_similarity_TAR, cosine_similarity_TAR_shf, h_similarity_TAR]...
    = decoding_target(x,initial_angle,initial_angle_shf,t,10);

colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;1	0.64706	0];
figure
set(gcf,'position',[100,100,360,250]);
hold on
plot_decoding(r2_TAR_x,r2_TAR_shf_x,h_TAR_x,t,colo(3,:),0.01)
plot_decoding(r2_TAR_y,r2_TAR_shf_y,h_TAR_y,t,colo(3,:),0.01)
plot_decoding(cosine_dist_TAR,cosine_dist_TAR_shf,h_cosine_TAR,t,colo(1,:),0.01)
plot_decoding(cosine_similarity_TAR, cosine_similarity_TAR_shf, h_similarity_TAR,t,colo(1,:),0.01)
