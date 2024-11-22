
function plot_TDR_reach(projectedStates)

t_start = 0;
t_end = 1200;
T = (t_start+100)/1000:0.1:(t_end+100)/1000;

%% plot 1
n_cond = size(projectedStates,2) / 2; % Assume there're equal numbers of conditions before and after learning.
c1 = [0.80392	0.16078	0.56471;0.80392	0.52157	0;0.18039	0.5451	0.34118;0.09412	0.4549	0.80392];
subp = [2,4,3,1];

figure
set(gcf,'position',[50,200,800,800]);
hold on
for i=1:n_cond
    %     subplot(2,2,subp(i))
    hold on
    for k = 1:size(projectedStates,1)
        c2 = gradient_ramp(c1(i,:),size(projectedStates,1),2);
        plot(projectedStates(k,i,1),projectedStates(:,i,2),'o','MarkerSize',25,'MarkerFaceColor',c1(i,:),'MarkerEdgeColor','k');
        plot(projectedStates(k,i+n_cond,1),projectedStates(k,i+n_cond,2),"d",'MarkerSize',15,'MarkerFaceColor',c2(k,:),'MarkerEdgeColor','k');
    end
    plot(projectedStates(:,i,1),projectedStates(:,i,2),'color',c1(i,:));
    plot(projectedStates(:,i+n_cond,1),projectedStates(:,i+n_cond,2),'color',c2(end,:));
    axis square
    box on
    set(gca, 'XTick', [], 'YTick', [])
%     axis off
end
x1 = xlim;
y1 = ylim;

%% Plot 4

figure
set(gcf,'position',[50,200,800,800]);
% hold on
for i=1:n_cond
    subplot(2,2,subp(i))
    hold on
    for k = 1:size(projectedStates,1)
        c2 = gradient_ramp(c1(i,:),size(projectedStates,1),2);
        plot(projectedStates(k,i,1),projectedStates(:,i,2),'o','MarkerSize',20,'MarkerFaceColor',c1(i,:),'MarkerEdgeColor','k');
        plot(projectedStates(k,i+n_cond,1),projectedStates(k,i+n_cond,2),"s",'MarkerSize',12,'MarkerFaceColor',c2(k,:),'MarkerEdgeColor','k');
    end
    plot(projectedStates(:,i,1),projectedStates(:,i,2),'color',c1(i,:));
    plot(projectedStates(:,i+n_cond,1),projectedStates(:,i+n_cond,2),'color',c2(end,:));
    xlim(x1);
    ylim(y1);
    axis square
    box on
    set(gca, 'XTick', [], 'YTick', [])
%     axis off
    
end

h1 =  subplot(2,2,1);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1), pos1(2)*0.78, pos1(3), pos1(4)]); % [left bottom width height]

h1 = subplot(2,2,2);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1)*0.83, pos1(2)*0.78, pos1(3), pos1(4)]); % [left bottom width height]

h1 = subplot(2,2,4);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1)*0.83, pos1(2), pos1(3), pos1(4)]); % [left bottom width height]

%% Plot degree 

theta = [];
for k = 1:size(projectedStates,1)
    for i=1:size(projectedStates,2)
        
        x = projectedStates(k,i,1);
        y = projectedStates(k,i,2);
        
        the = atan2(y, x);
        
        theta(k,i) = rad2deg(the);
        
        if theta(k,i) < 0
            theta(k,i) = theta(k,i) + 360;
        end
        
    end
end

figure
set(gcf,'position',[50,200,350,300]);
hold on

for i=1:n_cond
    plot(T(2:end-1),theta(2:end-1,i+n_cond),'Color',c1(i,:),'linewidth',3);
end

xlim([0.1,T(end)])
ylim([0 360])
ax = gca;
ax.TickDir = 'out';
set(gca, 'XColor', 'none'); 
xticks([0.2 0.6 1])
yticks([0 90 180 270 360])
plot([0.9 0.9],[0 360],':','color',[0.18039	0.5451	0.34118],'linewidth',1.2)
plot([1.1 1.1],[0 360],':','color',[0.80392	0.2	0.2],'linewidth',1.2)

hx = offsetaxis(ax, 'y', 0.02); % Default, y-offset = 0.1 and yloc = 'l'
hy = offsetaxis(ax, 'x', 0.02); % Default, y-offset = 0.1 and yloc = 'l'
hx.TickLength = [0.03 0.025];
hy.TickLength = [0.03 0.025];
hx.LineWidth = 1;
hy.LineWidth = 1;
