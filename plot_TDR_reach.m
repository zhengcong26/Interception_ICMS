
function plot_TDR_reach(projectedStates, target_loc, condition)

t_start = 0;
t_end = 1200;
T = (t_start+100)/1000:0.1:(t_end+100)/1000;

%% plot 1
n_cond = size(projectedStates,2) / 2; % Assume there're equal numbers of conditions before and after learning.
c1 = [0.80392	0.16078	0.56471;0.80392	0.52157	0;0.18039	0.5451	0.34118;0.09412	0.4549	0.80392];
subp = [2,4,3,1];

n = 1;temp_1=[];temp_2=[];
figure
set(gcf,'position',[50,200,300,300]);
hold on
for i=1:n_cond
    %     subplot(2,2,subp(i))
    hold on
    for k = 1:size(projectedStates,1)
        c2 = gradient_ramp(c1(i,:),size(projectedStates,1),2);
        plot(projectedStates(k,i,1),projectedStates(:,i,2),'o','MarkerSize',10,'MarkerFaceColor',c1(i,:),'MarkerEdgeColor','k');
        plot(projectedStates(k,i+n_cond,1),projectedStates(k,i+n_cond,2),"d",'MarkerSize',6,'MarkerFaceColor',c2(k,:),'MarkerEdgeColor','k');
        temp_1(k,n) = projectedStates(k,i+n_cond,1);
        temp_2(k,n) = projectedStates(k,i+n_cond,2);
    end
    n = n+1;
    plot(projectedStates(:,i,1),projectedStates(:,i,2),'color',c1(i,:));
    plot(projectedStates(:,i+n_cond,1),projectedStates(:,i+n_cond,2),'color',c2(end,:),'linewidth',1.2);
    axis square
    box on
    set(gca, 'XTick', [], 'YTick', [])
    %     axis off
end

x1 = xlim;
y1 = ylim;

% x1(1) = x1(1)-3;
% y1(2) = y1(2)+3;
% xlim(x1)
% ylim(y1)

%% Plot 4
figure
set(gcf,'position',[50,200,400,400]);
% hold on
for i=1:n_cond
    subplot(2,2,subp(i))
    hold on
    for k = 1:size(projectedStates,1)
        c2 = gradient_ramp(c1(i,:),size(projectedStates,1),2);
        plot(projectedStates(k,i,1),projectedStates(:,i,2),'o','MarkerSize',12,'MarkerFaceColor',c1(i,:),'MarkerEdgeColor','k');
        plot(projectedStates(k,i+n_cond,1),projectedStates(k,i+n_cond,2),"s",'MarkerSize',8,'MarkerFaceColor',c2(k,:),'MarkerEdgeColor','k');
    end
    plot(projectedStates(:,i,1),projectedStates(:,i,2),'color',c1(i,:));
    plot(projectedStates(:,i+n_cond,1),projectedStates(:,i+n_cond,2),'color',c2(end,:),'linewidth',1.2);
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

if condition == 1
    
    for i=1:n_cond
        for t = 1:size(theta,1)-1
            if theta(t,i+n_cond)<90 && theta(t+1,i+n_cond)>270
                theta(t+1,i+n_cond) = theta(t+1,i+n_cond)-360;
                continue
            end
        end
    end
    
    figure
    set(gcf,'position',[50,200,400,300]);
    hold on
    
    for i=1:n_cond
        plot(T(2:end-1),theta(2:end-1,i+n_cond),'Color',c1(i,:),'linewidth',4);
        tgt = plot(T(2:end-1),target_loc(2:end-1,i),'Color',c1(i,:),'linewidth',4);
        tgt.Color(4) = 0.2;
    end
    
    xlim([0.1,T(end)])
    ylim([-45 360])
    ax = gca;
    ax.TickDir = 'out';
    set(gca, 'XColor', 'none');
    xticks([0.2 0.6 1])
    yticks([0 90 180 270 360])
    plot([0.9 0.9],[-45 360],':','color',[0.18039	0.5451	0.34118],'linewidth',1.2)
    plot([1.1 1.1],[-45 360],':','color',[0.80392	0.2	0.2],'linewidth',1.2)
    
    hx = offsetaxis(ax, 'y', 0.02); % Default, y-offset = 0.1 and yloc = 'l'
    hy = offsetaxis(ax, 'x', 0.02); % Default, y-offset = 0.1 and yloc = 'l'
    hx.TickLength = [0.03 0.025];
    hy.TickLength = [0.03 0.025];
    hx.LineWidth = 1.2;
    hy.LineWidth = 1.2;
    
else
    for i=1:n_cond
%         % G
%         if i == 3
%             theta(:,i+n_cond) = theta(:,i+n_cond)-360;
%         end
%         
        for t = 1:size(theta,1)-1
            if theta(t,i+n_cond)<90 && theta(t+1,i+n_cond)>270
                theta(t+1:end,i+n_cond) = theta(t+1:end,i+n_cond)-360;
                continue
            end
        end
        
        % L
        if i == 4
            theta(:,i+n_cond) = theta(:,i+n_cond)+360;
        end
    end
    
    for i = 1:n_cond
        for t = 1:size(target_loc,1)-1
            if target_loc(t,i)<90 && target_loc(t+1,i)>270
                target_loc(t+1:end,i) = target_loc(t+1:end,i)-360;
                continue
            end
        end
    end
    
    figure
    set(gcf,'position',[50,200,400,300]);
    hold on
    
    for i=1:n_cond
        plot(T(2:end-1),theta(2:end-1,i+n_cond),'Color',c1(i,:),'linewidth',4);
        tgt = plot(T(2:end-1),target_loc(2:end-1,i),'Color',c1(i,:),'linewidth',4);
        tgt.Color(4) = 0.2;
    end
    
    xlim([0.1,T(end)])
    ylim([-180 360])
    ax = gca;
    ax.TickDir = 'out';
    set(gca, 'XColor', 'none');
    xticks([0.2 0.6 1])
    yticks([-180 -90 0 90 180 270 360])
    plot([0.9 0.9],[-180 360],':','color',[0.18039	0.5451	0.34118],'linewidth',1.2)
    plot([1.1 1.1],[-180 360],':','color',[0.80392	0.2	0.2],'linewidth',1.2)
    
    hx = offsetaxis(ax, 'y', 0.02); % Default, y-offset = 0.1 and yloc = 'l'
    hy = offsetaxis(ax, 'x', 0.02); % Default, y-offset = 0.1 and yloc = 'l'
    hx.TickLength = [0.03 0.025];
    hy.TickLength = [0.03 0.025];
    hx.LineWidth = 1.2;
    hy.LineWidth = 1.2;
    
end
