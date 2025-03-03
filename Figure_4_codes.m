
%% Figure 4a
% plot_units_heatmap.py

%% Figure 4b
% load('unit_G09146402.mat')
load('unit_G10125601.mat')
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

