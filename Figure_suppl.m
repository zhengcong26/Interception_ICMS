
%% Suppl Fig 9

% panel 3b
% load('G_0914_PC_NSST_PreGO.mat')

% panel 3c
load('G_0914_PC_NSST_PreMO.mat')

t_start = -200;
t_end = 400;
bin = 20;

ses = 1;
para = 4;
q = 1;

% PC
CO_NS_PC = squeeze(J_C_m(ses,para,:,:));
CO_ST_PC = squeeze(J_C_m(ses,para+4,:,:));
INT_NS_PC = squeeze(J_C_m(ses,para+8,:,:));
INT_ST_PC = squeeze(J_C_m(ses,para+12,:,:));

% CI
CO_NS_ciL = squeeze(J_ciLower(ses,para,:,:));
CO_ST_ciL = squeeze(J_ciLower(ses,para+4,:,:));
INT_NS_ciL = squeeze(J_ciLower(ses,para+8,:,:));
INT_ST_ciL = squeeze(J_ciLower(ses,para+12,:,:));

CO_NS_ciU = squeeze(J_ciUpper(ses,para,:,:));
CO_ST_ciU = squeeze(J_ciUpper(ses,para+4,:,:));
INT_NS_ciU = squeeze(J_ciUpper(ses,para+8,:,:));
INT_ST_ciU = squeeze(J_ciUpper(ses,para+12,:,:)); 

% ICMS
CO_icms_1 = squeeze(J_ICMSp(ses,1,para,6));
CO_icms_2 = squeeze(J_ICMSp(ses,2,para,6));
INT_icms_1 = squeeze(J_ICMSp(ses,1,para,8));
INT_icms_2 = squeeze(J_ICMSp(ses,2,para,8));

% RT MT
CO_NS_RT = squeeze(J_RTMT_mean(ses,1,para,5));
CO_NS_MT = squeeze(J_RTMT_mean(ses,2,para,5));
CO_ST_RT = squeeze(J_RTMT_mean(ses,1,para,6));
CO_ST_MT = squeeze(J_RTMT_mean(ses,2,para,6));

INT_NS_RT = squeeze(J_RTMT_mean(ses,1,para,7));
INT_NS_MT = squeeze(J_RTMT_mean(ses,2,para,7));
INT_ST_RT = squeeze(J_RTMT_mean(ses,1,para,8));
INT_ST_MT = squeeze(J_RTMT_mean(ses,2,para,8));

% significance
p_value_CO=[];p_value_INT=[];
for t = 1:size(CO_NS_PC,1)
    for i = 1:3
        [~, p_value_CO(t,i)] = ttest_from_ci(CO_NS_PC(t,i), [CO_NS_ciL(t,i), CO_NS_ciU(t,i)], 100, CO_ST_PC(t,i), [CO_ST_ciL(t,i), CO_ST_ciU(t,i)], 100, 0.05);
        [~, p_value_INT(t,i)] = ttest_from_ci(INT_NS_PC(t,i), [INT_NS_ciL(t,i), INT_NS_ciU(t,i)], 100, INT_ST_PC(t,i), [INT_ST_ciL(t,i), INT_ST_ciU(t,i)], 100, 0.05);
    end
end

p_value_CO(p_value_CO>=0.05)=0;
p_value_CO(p_value_CO<0.05&p_value_CO>0)=1;

p_value_INT(p_value_INT>=0.05)=0;
p_value_INT(p_value_INT<0.05&p_value_INT>0)=1;

t = t_start:bin:t_end;
plot_PC_4_2(J_C_m,J_ciLower,J_ciUpper,J_ICMSp,J_RTMT_mean,ses,t)

%% Suppl Fig 11 PSTH

load('unit_G09175003.mat')
% load('unit_G09171401.mat')
% load('unit_L09092401.mat')

t = -2000:50:2000-50;
t0 = find(t==0);
m = 1;

t_TO = [-0.05 0.75];
t_GO = [-0.15 0.15];
t_MO = [-0.05 0.25];
plot_psth_without_raster(spike_time_TO,t_TO,spike_time_GO,t_GO,spike_time_MO,t_MO,...
    cor_tri_sp_sl,hand_Sur_loc_4_sl)

%% Suppl Fig 11 PD

load('G_0917_PD.mat')
n = 13; % G09175003
% n = 1; % G09171401

% load('L_0909_PD.mat')
% n = 3; % L09092401

pd_dif = [];
for i = 1:size(pd_b,1)
    a = pd_b(i,:,1);
    b = pd_b(i,:,2);
    for k = 1:length(a)-1
        if a(k)>270 && a(k+1)<90
            aa = a(k+1)+360;
            pd_dif(i,k,1) = aa-a(k);
        elseif a(k+1)>270 && a(k)<90
            aa = a(k)+360;
            pd_dif(i,k,1) = a(k+1)-aa;
        else
            pd_dif(i,k,1) = a(k+1)-a(k);
        end
        
        if b(k)>270 && b(k+1)<90
            bb = b(k+1)+360;
            pd_dif(i,k,2) = bb-b(k);
        elseif b(k+1)>270 && b(k)<90
            bb = b(k)+360;
            pd_dif(i,k,2) = b(k+1)-bb;
        else
            pd_dif(i,k,2) = b(k+1)-b(k);
        end
    end
end

A = squeeze(pd_b(n,:,1));
B = squeeze(pd_b(n,:,2));
a = abs(squeeze(pd_dif(n,:,1)));
b = abs(squeeze(pd_dif(n,:,2)));

colo=[0.09412	0.4549	0.80392;0.80392	0.2	0.2];
figure
set(gcf,'Position',[100,300,700,200])

h1 = subplot(2,1,1);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1), pos1(2), pos1(3), pos1(4)]); % [left bottom width height]
hold on
t1 = 1:1:length(tt);
plot_pd_arrow(A,colo(1,:))
plot_pd_arrow(B,colo(2,:))
xlim([-1,33])
x0 = xlim;
plot([x0(1) x0(2)-1],[0,0],'k:','linewidth',1)
ylim([-1.5 1]);
y0 = ylim;
plot([17.5 17.5],[y0(1) y0(2)],'k:','linewidth',1.2)
plot([24.5 24.5],[y0(1) y0(2)],'k:','linewidth',1.2)
scatter(2,y0(1),25,'ko','filled')
scatter(21,y0(1),25,'ko','filled')
scatter(26,y0(1),25,'ko','filled')

plot_pd_arrow(A,colo(1,:))
plot_pd_arrow(B,colo(2,:))

ax = gca;
ax.TickDir='out';
axis off

subplot(2,1,2)
hold on
t1 = 2:1:length(tt);
plot(t1,smoothdata(a,'gaussian',5),'color',colo(1,:),'linewidth',2.3)
plot(t1,smoothdata(b,'gaussian',5),'color',colo(2,:),'linewidth',2.3)

y0 = ylim;
dy = y0(2)/5;
ylim([-dy,y0(2)])
y0 = ylim;
%             ylim([y0(1)-50 y0(2)])
%             y0 = ylim;
plot([17.5 17.5],[y0(1) y0(2)],'k:','linewidth',1.2)
plot([24.5 24.5],[y0(1) y0(2)],'k:','linewidth',1.2)
xlim([-1,33])
x0 = xlim;
plot([x0(2)-1 x0(2)-1],[0 20],'k','linewidth',2)
scatter(2,y0(1),25,'ko','filled')
scatter(21,y0(1),25,'ko','filled')
scatter(26,y0(1),25,'ko','filled')
plot([x0(1)+1 x0(2)-1],[0,0],'k:','linewidth',1)
%             plot([6,8],[y0(1),y0(1)],'k','linewidth',1)

plot(t1,smoothdata(a,'gaussian',5),'color',colo(1,:),'linewidth',2.3)
plot(t1,smoothdata(b,'gaussian',5),'color',colo(2,:),'linewidth',2.3)

ax = gca;
ax.TickDir='out';
axis off

%% Suppl Fig 13

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

%% Suppl Fig 15abc
% detailed decoeding process refer to decoding_xy_cosine.m

% load('G_0914_decoding_20ms.mat')
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

%% Suppl Fig 18i

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

%% Suppl Fig 18j

load('model_u.mat')
 
it = [102 120 167];

u = u_INT{1,1};
plot_u(u,it)

u = u_CO{1,1};
plot_u(u,it)

%% Suppl Fig 19a 
load('G_eye_distance.mat')

figure
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
set(gcf,'Position',[100,100,1100,200]);

for p = 1:4
    subplot(1,4,p)
    hold on
    
    % Set time range and x-axis limits
    if p == 1
        t = 0:1:200;
        xlim([0,200])
    else
        t = -100:1:100;
        xlim([-100,100])
    end

    % Extract data for the current subplot
    a = squeeze(dist_tg_CO(:,p,:));
    b = squeeze(dist_ct_CO(:,p,:));
    A = squeeze(dist_tg_INT(:,p,:));
    B = squeeze(dist_ct_INT(:,p,:));
    
    % Identify and remove outliers where values exceed the threshold
    a_r = any(a > 30, 2);
    b_r = any(b > 30, 2);
    A_r = any(A > 30, 2);
    B_r = any(B > 30, 2);
    
    remov_1 = a_r | b_r;
    remov_2 = A_r | B_r; 
    
    a = a(~remov_1, :);
    b = b(~remov_1, :);
    A = A(~remov_2, :);
    B = B(~remov_2, :);
    
     % Statistical testing (Rank-Sum Test) for each time point
    h_a = [];h_b = [];
    for tt = 1:size(a,2)
        [~,h_a(tt,1)]=ranksum(a(:,tt),A(:,tt)); % tg
        [~,h_b(tt,1)]=ranksum(b(:,tt),B(:,tt)); % ct
    end
    
    plot_eye_distance(a,t,colo(1,:),'--')
    plot_eye_distance(b,t,colo(1,:),'-')
    plot_eye_distance(A,t,colo(2,:),'--')
    plot_eye_distance(B,t,colo(2,:),'-')   
    
    ylim([0 22])
    y0 = ylim;
    
    % Mark significant differences with horizontal lines
    for tt = 1:size(h_a,1)
        if h_a(tt,1)==1
            plot([t(tt)-2 t(tt)+2],[y0(2)-1 y0(2)-1],'color',[0.62745	0.12549	0.94118],'LineWidth',1.6)
        end
        if h_b(tt,1)==1
            plot([t(tt)-2 t(tt)+2],[y0(1)+1 y0(1)+1],'color',[0.62745	0.12549	0.94118],'LineWidth',1.6)
        end
    end
    
    plot([0 0],[y0(1) y0(2)],'k:')
    
    % Customize axis appearance
    ax = gca;
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    ax.TickLength = [0.03 0.025];
    yticklabels([]) % Remove y-axis tick labels
    set(gca,'XColor', 'none')
    if p >=2
        set(gca,'YColor', 'none')
    end

end

h1 =  subplot(1,4,2);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1)*0.9, pos1(2), pos1(3), pos1(4)]); % [left bottom width height]

h1 = subplot(1,4,3);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1)*0.88, pos1(2), pos1(3), pos1(4)]); % [left bottom width height]

h1 = subplot(1,4,4);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1)*0.87, pos1(2), pos1(3), pos1(4)]); % [left bottom width height]

%% Suppl Fig 19b

load('saccade_probability.mat') 

% % GO -0.9  -0.5
x_1 = cos_thetas_CO_NS(locs_CO_NS>-900&locs_CO_NS<-500);
x_2 = cos_thetas_INT_NS(locs_INT_NS>-900&locs_INT_NS<-500);
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];

a = -1;
b = 1;
binEdges = linspace(a, b, 20); 
x_1 = x_1(x_1>=a&x_1<=b);
x_2 = x_2(x_2>=a&x_2<=b);

h1 = histogram(x_1, 'BinEdges', binEdges);
hold on
h2 = histogram(x_2, 'BinEdges', binEdges);
dbin = h1.BinWidth;
data1 = h1.Values;
data2 = h2.Values;
close all

counts1 = data1; 
counts2 = data2;

figure
set(gcf,'Position',[100,100,150,100]);
hold on
h1 = histogram(x_1, 'BinEdges', binEdges, 'Normalization', 'probability');
h1.FaceColor = colo(1,:);
h1.EdgeColor = [0.5 0.5 0.5];
h2 = histogram(x_2, 'BinEdges', binEdges, 'Normalization', 'probability');
h2.FaceColor = colo(2,:);
h2.EdgeColor = [0.5 0.5 0.5];
y0 = ylim;
ylim([y0(1) y0(2)+0.01]);
y0 = ylim;
for i = 1:length(counts1)
    table = [counts1(i), sum(counts1) - counts1(i);
        counts2(i), sum(counts2) - counts2(i)];
    
    [h, ~] = fishertest(table,'Tail','right','Alpha',0.05);
    if h == 1
        plot([binEdges(i)+dbin/2-0.02 binEdges(i)+dbin/2+0.02],[y0(2) y0(2)],'color',colo(1,:),'linewidth',2)
    end
    [h, ~] = fishertest(table,'Tail','left','Alpha',0.05);
    if h == 1
        plot([binEdges(i)+dbin/2-0.02 binEdges(i)+dbin/2+0.02],[y0(2) y0(2)],'color',colo(2,:),'linewidth',2)
    end
end

ax = gca;
ax.LineWidth = 0.5;
ax.TickLength = [0.02 0.025];
ax.TickDir = 'out';

%% Suppl Fig 19c

% % б└240бу/s
load('G_eye_suppl_c_1.mat') 

% % ICMS < GO-500 ms
load('G_eye_suppl_c_2.mat') 

% % ICMS > GO
load('G_eye_suppl_c_3.mat') 

colo=[0.4	0.80392	0.66667;0.4902	0.14902	0.80392;1	0.44706	0.33725];

a = -1;
b = 1;
binEdges = linspace(a, b, 20); 
x_1 = x_1(x_1>=a&x_1<=b);
x_2 = x_2(x_2>=a&x_2<=b);

h1 = histogram(x_1, 'BinEdges', binEdges);
hold on
h2 = histogram(x_2, 'BinEdges', binEdges);
dbin = h1.BinWidth;
data1 = h1.Values;
data2 = h2.Values;
close all

counts1 = data1; 
counts2 = data2;

figure
set(gcf,'Position',[100,100,150,100]);
hold on
h1 = histogram(x_1, 'BinEdges', binEdges, 'Normalization', 'probability');
h1.FaceColor = colo(1,:);
h1.EdgeColor = [0.5 0.5 0.5];
h2 = histogram(x_2, 'BinEdges', binEdges, 'Normalization', 'probability');
h2.FaceColor = colo(2,:);
h2.EdgeColor = [0.5 0.5 0.5];
y0 = ylim;
ylim([y0(1) y0(2)+0.01]);
y0 = ylim;
for i = 1:length(counts1)
    table = [counts1(i), sum(counts1) - counts1(i);
        counts2(i), sum(counts2) - counts2(i)];
    
    [h, ~] = fishertest(table,'Tail','right','Alpha',0.05);
    if h == 1
        plot([binEdges(i)+dbin/2-0.02 binEdges(i)+dbin/2+0.02],[y0(2) y0(2)],'color',colo(1,:),'linewidth',2)
    end
    [h, ~] = fishertest(table,'Tail','left','Alpha',0.05);
    if h == 1
        plot([binEdges(i)+dbin/2-0.02 binEdges(i)+dbin/2+0.02],[y0(2) y0(2)],'color',colo(2,:),'linewidth',2)
    end
end

ax = gca;
ax.LineWidth = 0.5;
ax.TickLength = [0.02 0.025];
ax.TickDir = 'out';

