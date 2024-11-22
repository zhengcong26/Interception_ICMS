
%% Figure 4a
% plot_units_heatmap.py

%% Figure 4b
load('unit_G09146402.mat')
% load('unit_G10125601.mat')
% load('unit_L08190801.mat')
% load('unit_L09096001.mat')

plot_psth_with_raster(cor_tri_sp_sl,spike_time_sl,hand_Sur_loc_4_sl,cor_tri_Times_sl)

%% Figure 4c histogram
% load('G_psth_p_value.mat')
load('L_psth_p_value.mat')

d_step = 50;  % Step size (50ms)
x_1=[];x_2=[];
for i = 1:size(J_p_1,1)
    % Extract p-values for two datasets and find where p < 0.05
    d_dir_1 = find((J_p_1(i,:)<0.05)==1)*d_step;
    d_dir_1(d_dir_1>1500)=[];
    x_1 = [x_1,d_dir_1/1000-1];
    
    d_dir_2 = find((J_p_2(i,:)<0.05)==1)*d_step;
    d_dir_2(d_dir_2>1500)=[];
    x_2 = [x_2,d_dir_2/1000-1];
end

binEdges = linspace(-1, 0.6, 30); 
h1 = histogram(x_1, 'BinEdges', binEdges);
hold on
h2 = histogram(x_2, 'BinEdges', binEdges);
dbin = h1.BinWidth;
counts1 = h1.Values;
counts2 = h2.Values;
close all

figure
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
set(gcf,'Position',[100,100,400,120]);
hold on

% Plot histograms for both datasets
h1 = histogram(x_1, 'BinEdges', binEdges, 'Normalization', 'probability');
h1.FaceColor = colo(1,:);
h1.EdgeColor = [0.5 0.5 0.5];
h2 = histogram(x_2, 'BinEdges', binEdges, 'Normalization', 'probability');
h2.FaceColor = colo(2,:);
h2.EdgeColor = [0.5 0.5 0.5];

% Add statistical significance markers (Fisher's exact tes
y0 = ylim;
for i = 1:length(counts1)
    table = [counts1(i), sum(counts1) - counts1(i);
        counts2(i), sum(counts2) - counts2(i)];
    
    % Fisher's test for right-tail significance
    [h, ~] = fishertest(table,'Tail','right','Alpha',0.05);
    if h == 1
        plot([binEdges(i)+dbin/2-0.02 binEdges(i)+dbin/2+0.02],[y0(2) y0(2)],'color',colo(1,:),'linewidth',2)
    end
    
    % Fisher's test for left-tail significance
    [h, ~] = fishertest(table,'Tail','left','Alpha',0.05);
    if h == 1
        plot([binEdges(i)+dbin/2-0.02 binEdges(i)+dbin/2+0.02],[y0(2) y0(2)],'color',colo(2,:),'linewidth',2)
    end
end

ax = gca;
ax.LineWidth = 1;
ax.TickDir = 'out';
xlim([-1.05 0.6])

%% Figure 4d PSTH

% load('unit_G09175003.mat')
% load('unit_G09171401.mat')
load('unit_L09092401.mat')

t = -2000:50:2000-50;
t0 = find(t==0);
m = 1;

t_TO = [-0.05 0.75];
t_GO = [-0.15 0.15];
t_MO = [-0.05 0.25];
plot_psth_without_raster(spike_time_TO,t_TO,spike_time_GO,t_GO,spike_time_MO,t_MO,...
    cor_tri_sp_sl,hand_Sur_loc_4_sl)
%% Figure 4d PD

% load('G_0917_PD.mat')
% n = 13; % G09175003
% n = 1; % G09171401

load('L_0909_PD.mat')
n = 3; % L09092401

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
