
function plot_psth_with_raster(cor_tri_sp,spike_time_MO,hand_Sur_loc,event_time)

speed_con=sort(unique(cor_tri_sp),'descend');%正240逆时针，CO，负240顺时针

d_step = 50;  % Step size (50ms)
d_step1 = 1:d_step:3951;  % Start points for each step
d_step2 = d_step1 + d_step - 1;  % End points for each step

% Prepare spike time data
spike_time_MO_zero = spike_time_MO;
spike_time_MO_zero(isnan(spike_time_MO_zero))=0; % Replace NaN values with 0
index_bin_MO = ceil(spike_time_MO*1000+2000); % Convert time to bin index
spike_bin_MO  = zeros(size(cor_tri_sp,1),4000); % Initialize spike bin array

% Fill spike bins
for i = 1:size(index_bin_MO,1)
    spike_bin_MO(i,index_bin_MO(i,~isnan(index_bin_MO(i,:))) ) = 1;
end

for i_t = 1:length(d_step1)
    % Calculate firing rate for each step
    FR_MO_50(:,i_t) = sum(spike_bin_MO(:,d_step1(i_t):d_step2(i_t)),2);
    % Speed comparison (ranksum test)
    x = FR_MO_50 (cor_tri_sp==speed_con(1),i_t);
    y = FR_MO_50 (cor_tri_sp==speed_con(2),i_t);
    p_2sp(:,i_t) = ranksum(x,y);
    % Direction comparison (kruskalwallis test)
    p_dir_1(:,i_t) = kruskalwallis(FR_MO_50(cor_tri_sp==speed_con(1),i_t), hand_Sur_loc(cor_tri_sp==speed_con(1),:),'off');%每个d_step内比较八方向有没有显著差异
    p_dir_2(:,i_t) = kruskalwallis(FR_MO_50(cor_tri_sp==speed_con(2),i_t), hand_Sur_loc(cor_tri_sp==speed_con(2),:),'off');%每个d_step内比较八方向有没有显著差异
end

% Event time correction for GO
tri_time=event_time(:,3:7)-event_time(:,5)*ones(1,5); % GO

colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
c1 = gradient_ramp(colo(1,:),4,1);
c2 = gradient_ramp(colo(2,:),4,1);

figure
set(gcf,'Position',[400,100,250,600]);

%% static
subplot(3,1,1)
for i_dir=1:4
    set(gca,'XTick',[-1,-0.5,0,0.5], 'XLim',[-1 0.6], 'FontSize',10);
    spk_dir_MO{i_dir,1}=spike_time_MO_zero(hand_Sur_loc==i_dir&cor_tri_sp==speed_con(1) ,:);
    psth_NeurAn(spk_dir_MO{i_dir,1},0.05,c1(i_dir,:));
    hold on
end
axis off

% Find significant points for direction comparison
tt = -2:0.05:2-0.05;
d_dir_1 = find(p_dir_1 < 0.05);
x_q_1 = tt(d_dir_1)-d_step/1000/2; 
x_q2_1 = x_q_1+d_step/1000; 
y_q_1 = ylim;

%% moving
subplot(3,1,2)
for i_dir=1:4
    set(gca,'XTick',[-1,-0.5,0,0.5], 'XLim',[-1 0.6],'FontSize',10);
    spk_dir_MO{i_dir,1}=spike_time_MO_zero(hand_Sur_loc==i_dir&cor_tri_sp==speed_con(2) ,:);
    psth_NeurAn(spk_dir_MO{i_dir,1},0.05,c2(i_dir,:));
    hold on
end
axis off

d_dir_2 = find(p_dir_2 < 0.05);
x_q_2 = tt(d_dir_2)-d_step/1000/2; 
x_q2_2 = x_q_2+d_step/1000; 
y_q_2 = ylim;

%% compare
subplot(3,1,3)
for k=1:length(speed_con)
    set(gca,'XTick',[-1,-0.5,0,0.5], 'XLim',[-1 0.6],'FontSize',10);
    if speed_con(k) == 0
        spk_speed_MO{k,1}=spike_time_MO_zero(cor_tri_sp==speed_con(k) & any(spike_time_MO_zero,2),:);
        psth_NeurAn(spk_speed_MO{k,1},0.05,[0.06275	0.30588	0.5451]);
        hold on
    else
        spk_speed_MO{k,1}=spike_time_MO_zero(cor_tri_sp==speed_con(k) & any(spike_time_MO_zero,2),:);
        psth_NeurAn(spk_speed_MO{k,1},0.05,[0.80392	0.33333	0.33333]);
        hold on
    end
end
axis off

tt = -2:0.05:2-0.05;
d_5sp = find(p_2sp < 0.05);
x_p = tt(d_5sp)-d_step/1000/2; 
x_p2 = x_p+d_step/1000; 
y_p = ylim;

%% Adjust axis limits
y_max = max ([y_p, y_q_1, y_q_2]);
y_min = min ([y_p, y_q_1, y_q_2]);

% static plot markers
subplot(3,1,1)
set(gca,'YLim', [y_min y_max*1.2], 'FontSize',10);
for i_dir=4:-1:1
    rasters_speed_timemarker(spike_time_MO(hand_Sur_loc==i_dir&cor_tri_sp==speed_con(1),:),tri_time(hand_Sur_loc==i_dir&cor_tri_sp==speed_con(1),:),colo(1,:),'GO');
end
for i_p= 1:length(x_q_1)
    plot([x_q_1(i_p) x_q2_1(i_p)],[y_max y_max],'k','LineWidth',2);
end
plot([0.6 0.6],[y_min y_min+20],'k','LineWidth',2);
scatter(0,y_min,15,'k','filled')
scatter(-0.9,y_min,15,'k','filled')
scatter(0.25,y_min,15,'k','filled')

% moving plot markers
subplot(3,1,2)
set(gca,'YLim', [y_min y_max*1.2], 'FontSize',10);
for i_dir=4:-1:1
    rasters_speed_timemarker(spike_time_MO(hand_Sur_loc==i_dir&cor_tri_sp==speed_con(2),:),tri_time(hand_Sur_loc==i_dir&cor_tri_sp==speed_con(2),:),colo(2,:),'GO');
end
for i_p= 1:length(x_q_2)
    plot([x_q_2(i_p) x_q2_2(i_p)],[y_max*1.1 y_max*1.1],'k','LineWidth',2);
end
plot([0.6 0.6],[y_min y_min+20],'k','LineWidth',2);
scatter(0,y_min,15,'k','filled')
scatter(-0.9,y_min,15,'k','filled')
scatter(0.25,y_min,15,'k','filled')

% Compare plot markers
subplot(3,1,3)
for i_p= 1:length(x_p)
    plot([x_p(i_p) x_p2(i_p)],[y_max y_max],'k','LineWidth',2);
end
plot([0.6 0.6],[y_min y_min+20],'k','LineWidth',2);
scatter(0,y_min,15,'k','filled')
scatter(-0.9,y_min,15,'k','filled')
scatter(0.25,y_min,15,'k','filled')

%% Adjust subplot positions for better layout
h1 = subplot(3,1,1);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1), pos1(2)*0.9, pos1(3), pos1(4)*1.3]); % [left bottom width height]

h1 = subplot(3,1,2);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1), pos1(2)*0.8, pos1(3), pos1(4)*1.3]); % [left bottom width height]

h1 = subplot(3,1,3);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1), pos1(2), pos1(3), pos1(4)*0.6]); % [left bottom width height]

