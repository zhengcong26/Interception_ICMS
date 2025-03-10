
%% Figure 7a eye 
load('G_eye_example.mat')

dtheta = atand(10.2/30)/10.2;
scale = 10.2/16*dtheta;

% for n = 1:num 

% for n = 42 % static
for n = 2 % moving
    %% Initialize figure
    figure;
    set(gcf, 'unit', 'centimeters', 'position', [10 5 11 8]);
    
    %% Plot Eye X
    subplot('Position', [0.13, 0.56, 0.6, 0.37]); % Upper plot
    hold on;
    
    % Plot velocity and Eye X position
    plot(eye{n,5}(:,1), smoothdata(eye{n,5}(:,2), 'gaussian', 30) / 22, 'color', [0.7, 0.7, 0.7], 'Linewidth', 1.4); % Velocity
    plot(eye{n,5}(:,1), eye{n,1}(:,1) * scale, 'b', 'Linewidth', 1.4); % Eye X position
    
    ylim([-26, 26]);
    xlabel("Time from center show (ms)");
    ylabel("Eye X (deg)");
    
    %% Calculate and plot target trajectory
    x_t = eye{n,5}(eye{n,1}(:,1) == eye{n,2}(3,1) & eye{n,1}(:,2) == eye{n,2}(3,2), 1);
    T_TO = atan2d(eye{n,4}(3,2), eye{n,4}(3,1));
    target_time = x_t(1):eye{n,5}(end,1);
    
    for i = 1:numel(target_time)
        eye{n,6}(i,1) = target_time(i);
        [eye{n,6}(i,2), eye{n,6}(i,3)] = pol2cart(deg2rad(T_TO + stim_list(n,3) * i / 1000), 16);
        eye{n,6}(i,4) = T_TO + stim_list(n,3) * i / 1000;
    end
    
    plot(eye{n,6}(:,1), eye{n,6}(:,2) * scale, 'color', [0.13333, 0.5451, 0.13333], 'Linewidth', 1.4);
    
    %% Plot significant events (TC, TO, GO, Stim, TA)
    ylimit = ylim;
    event_labels = {'TC', 'TO', 'GO', 'TA'};
    event_indices = [2, 3, 4, 6]; % Event indices in eye{n,2}
    
    for idx = 1:numel(event_indices)
        x_event = eye{n,5}(eye{n,1}(:,1) == eye{n,2}(event_indices(idx),1) & ...
            eye{n,1}(:,2) == eye{n,2}(event_indices(idx),2), 1);
        plot([x_event(1), x_event(1)], ylimit, ':k', 'linewidth', 1.2);
        text(x_event(1), ylimit(2) + 1, event_labels{idx});
        if idx == 4 % TA event
            scatter(x_event(1), eye{n,3}(6,1) * scale, 'r', "filled");
        end
    end
    
    % MO
    x = eye{n,7}(1,6)-eye{n,7}(1,1);
    plot([x(1) x(1)],ylimit,'k:','linewidth',1.2)
    text(x(1),ylimit(2)+1,'MO');
    
    % Add stimulation markers if present
    if stim_list(n,27) ~= 0
        stim_start = x_event(1) + stim_list(n,27);
        for j = stim_start:10:stim_start + 60
            plot([j, j], ylimit, 'r', 'linewidth', 1);
        end
        text(stim_start, ylimit(2) + 1, 'Stim');
    end
    
    xlim([-inf, x_event(1) + 5]);
    
    %% formatting
    ax = gca;
    box off
    ax.TickDir='out';
    ax.FontName='Arial';
    ax.FontSize=10;
    set(gca,'XColor', 'none')
    
    %% Plot Eye Y
    subplot('Position', [0.13, 0.11, 0.6, 0.37]); % Lower plot
    hold on;
    
    % Plot velocity and Eye Y position
    plot(eye{n,5}(:,1), smoothdata(eye{n,5}(:,3), 'gaussian', 30) / 22, 'color', [0.7, 0.7, 0.7], 'Linewidth', 1.4); % Velocity
    plot(eye{n,5}(:,1), eye{n,1}(:,2) * scale, 'b', 'Linewidth', 1.4); % Eye Y position
    plot(eye{n,6}(:,1), eye{n,6}(:,3) * scale, 'color', [0.13333, 0.5451, 0.13333], 'Linewidth', 1.4); % Target trajectory
    
    ylim([-26, 26]);
    xlabel("Time from center show (ms)");
    ylabel("Eye Y (deg)");
    
    %% Plot significant events (TC, TO, GO, Stim, TA)
    ylimit = ylim;
    event_indices = [2, 3, 4, 6]; % Event indices in eye{n,2}
    
    for idx = 1:numel(event_indices)
        x_event = eye{n,5}(eye{n,1}(:,1) == eye{n,2}(event_indices(idx),1) & ...
            eye{n,1}(:,2) == eye{n,2}(event_indices(idx),2), 1);
        plot([x_event(1), x_event(1)], ylimit, ':k', 'linewidth', 1.2);
        if idx == 4 % TA event
            scatter(x_event(1), eye{n,3}(6,2) * scale, 'r', "filled");
        end
    end
    
    % MO
    x = eye{n,7}(1,6)-eye{n,7}(1,1);
    plot([x(1) x(1)],ylimit,'k:','linewidth',1.2)
    
    % Add stimulation markers if present
    if stim_list(n,27) ~= 0
        stim_start = x_event(1) + stim_list(n,27);
        for j = stim_start:10:stim_start + 60
            plot([j, j], ylimit, 'r', 'linewidth', 1);
        end
        text(stim_start, ylimit(2) + 1, 'Stim');
    end
    
    xlim([-inf, x_event(1) + 5]);
    
    %% formatting
    ax = gca;
    box off
    ax.TickDir='out';
    ax.FontName='Arial';
    ax.FontSize=10;
    set(gca,'XColor', 'none')
    
    %% save
%     saveas(1,['Trial ' num2str(n)],'jpg'); % individual file
%     close all
    
end

%% Figure 7b
load('G_eye_RT.mat')

colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
plot_eye_RT_pdf(eye_RT_CO,colo(1,:))
plot_eye_RT_pdf(eye_RT_INT,colo(2,:))

%% Figure 7c

load('saccade_probability.mat') % for Figure c, f, e, g

x_1 = locs_CO_NS;
x_2 = locs_INT_NS;
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
a = -900;
b = 100;
binEdges = linspace(a, b, 50); 

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
set(gcf,'Position',[100,100,500,100]);
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
        plot([binEdges(i)+dbin/2-5 binEdges(i)+dbin/2+5],[y0(2) y0(2)],'color',colo(1,:),'linewidth',2)
    end
    [h, ~] = fishertest(table,'Tail','left','Alpha',0.05);
    if h == 1
        plot([binEdges(i)+dbin/2-5 binEdges(i)+dbin/2+5],[y0(2) y0(2)],'color',colo(2,:),'linewidth',2)
    end
end

ax = gca;
ax.LineWidth = 0.5;
ax.TickLength = [0.008 0.025];
ax.TickDir = 'out';

%% Figure 7f

% % static
% x_1 = locs_CO_NS;
% x_2 = locs_CO_ST;

% moving
x_1 = locs_INT_NS;
x_2 = locs_INT_ST;

colo=[0.4	0.80392	0.66667;0.4902	0.14902	0.80392;1	0.44706	0.33725];

a = -100;
b = 100;
binEdges = linspace(a, b, 30); 

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
% set(gcf,'Position',[100,100,400,120]);
set(gcf,'Position',[100,100,220,100]);
hold on
h1 = histogram(x_1, 'BinEdges', binEdges, 'Normalization', 'probability');
h1.FaceColor = colo(1,:);
h1.EdgeColor = [0.5 0.5 0.5];
h2 = histogram(x_2, 'BinEdges', binEdges, 'Normalization', 'probability');
h2.FaceColor = colo(2,:);
h2.EdgeColor = [0.5 0.5 0.5];
y0 = ylim;
ylim([y0(1) y0(2)+0.05]);
y0 = ylim;
for i = 1:length(counts1)
    table = [counts1(i), sum(counts1) - counts1(i);
        counts2(i), sum(counts2) - counts2(i)];
    
    [h, ~] = fishertest(table,'Tail','right','Alpha',0.05);
    if h == 1
        plot([binEdges(i)+dbin/2-2 binEdges(i)+dbin/2+2],[y0(2) y0(2)],'color',colo(1,:),'linewidth',2)
    end
    [h, ~] = fishertest(table,'Tail','left','Alpha',0.05);
    if h == 1
        plot([binEdges(i)+dbin/2-2 binEdges(i)+dbin/2+2],[y0(2) y0(2)],'color',colo(2,:),'linewidth',2)
    end
end

ax = gca;
ax.LineWidth = 0.5;
ax.TickDir = 'out';

%% Figure 7e, g

% %  Figure 7e GO -0.5  +0.1
x_1 = cos_thetas_CO_NS(locs_CO_NS>-500&locs_CO_NS<100);
x_2 = cos_thetas_INT_NS(locs_INT_NS>-500&locs_INT_NS<100);
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];

% % Figure 7g
% % static
% x_1 = cos_thetas_CO_NS(locs_CO_NS>-100&locs_CO_NS<100);
% x_2 = cos_thetas_CO_ST(locs_CO_ST>-100&locs_CO_ST<100);
% colo=[0.4	0.80392	0.66667;0.4902	0.14902	0.80392;1	0.44706	0.33725];

% % Figure 7g
% % moving
% x_1 = cos_thetas_INT_NS(locs_INT_NS>-100&locs_INT_NS<100);
% x_2 = cos_thetas_INT_ST(locs_INT_ST>-100&locs_INT_ST<100);
% colo=[0.4	0.80392	0.66667;0.4902	0.14902	0.80392;1	0.44706	0.33725];

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

