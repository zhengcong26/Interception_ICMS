%% Figure 1c
% depends on plot_session_structure.m

% load('E:\manuscript\20230323 manuscript\version NC revision-1\example codes\G_20220914_J_bin_20msBin_4s.mat')
% trialconstraint = 675; 

load('E:\Laplace task collection\J\L_20220818_J_bin_20msBin_4s.mat')
trialconstraint = 644; 

CO = stim_list(stim_list(:,2) == 0&stim_list(:,3) == 0&stim_list(:,19) ~= 1&stim_list(:,1) < trialconstraint,1);
INT = stim_list(stim_list(:,2) == 0&stim_list(:,3) ~= 0&stim_list(:,19) ~= 1&stim_list(:,1) < trialconstraint,1);
RAND = stim_list(stim_list(:,2) == 0&stim_list(:,19) == 1&stim_list(:,1) < trialconstraint,1);

figure
set(gcf,'position',[30,300,500,500]);
colo = [0.06275	0.30588	0.5451;0.5451	0.1451	0;0.7 0.7 0.7];
hold on
plot_session_structure(RAND,velocity_list,colo(3,:))
plot_session_structure(CO,velocity_list,colo(1,:))
plot_session_structure(INT,velocity_list,colo(2,:))

axis equal
axis off

CO(:,2)=1;
INT(:,2)=2;
RAND(:,2)=3;
A = sortrows([CO;INT;RAND],1);
length(CO)/length(A)
length(INT)/length(A)
length(RAND)/length(A)

%% Figure 1d
% depends on ErrorEllipse.m

% monkey G
load('E:\manuscript\20230323 manuscript\2023-07-17_G_J_COINT_behaviors.mat')

% monkey L
load('E:\manuscript\20230323 manuscript\2023-07-17_L_J_COINT_behaviors.mat')

figure
set(gcf,'unit','centimeters','position',[10 5 15 15]);
pos = [2,1,3,4];
colo = [0	0.69804	0.93333;1	0.49804	0;0.06275	0.30588	0.5451;0.5451	0.1451	0];
linestyle = {'-','-','-','-'};
for q = 1:4
    subplot(2,2,pos(q))
    hold on
    if q == 1
        M=[cos(pi/4),sin(pi/4);-sin(pi/4),cos(pi/4)];
    elseif q == 2
        M=[cos(pi/4),-sin(pi/4);sin(pi/4),cos(pi/4)];
    elseif q == 3
        M=[cos(3*pi/4),-sin(3*pi/4);sin(3*pi/4),cos(3*pi/4)];
    elseif q == 4
        M=[cos(3*pi/4),sin(3*pi/4);-sin(3*pi/4),cos(3*pi/4)];
    end
    %% arc
    R1=[[-5 5];[0,0]];
    R2=M*R1;%旋转后坐标
    plot(R2(1,:),R2(2,:),':k','LineWidth',0.8)
    
    R1=[[0,0];[-5 5]];
    R2=M*R1;%旋转后坐标
    plot(R2(1,:),R2(2,:),':k','LineWidth',0.8)
    
    % plot arc
    theta1 = pi/4;theta2 = 3*pi/4;
    the = theta1:pi/180:theta2;
    r = 16;
    x0 = 0;y0 = -16;
    x = x0 + r*cos(the);
    y = y0 + r*sin(the);
    R1=[x;y];R2=M*R1;%旋转后坐标
    plot(R2(1,:),R2(2,:),':k','LineWidth',1)
    
    % plot error tolerance
    theta1 = 0;theta2 = 2*pi;
    the = theta1:pi/180:theta2;
    r = 4.5; % tolerance
    x0 = 0;y0 = 0;
    x = x0 + r*cos(the);
    y = y0 + r*sin(the);
    R1=[x;y];R2=M*R1;%旋转后坐标
    fill(R2(1,:), R2(2,:), 'b', 'FaceAlpha', 0.08, 'EdgeColor', 'none'); % 填充半透明蓝色
    
    %%
    for c = 1:4
        X = J_tran_group{q,c}(:,24);
        Y = J_tran_group{q,c}(:,25);
        
        [~,~,~,X0,Y0,r_ellipse] = ErrorEllipse([X Y],0.95,colo,'-','k',1,0);
        
        R1=[r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0];
        R2=M*R1';%旋转后坐标
        plot(R2(1,:),R2(2,:),'Color',colo(c,:),'LineStyle',linestyle{c},'LineWidth',2.2)
        
        R1=[X0;Y0];
        R2=M*R1;%旋转后坐标
        scatter(R2(1),R2(2),70,'b','d','MarkerEdgeColor',colo(c,:),'LineWidth',2.2)
        
        %%
        axis([-7 7 -7 7])
        axis square
        box off
        axis off;
        set(gca,'LineWidth',2)
        
        if q == 3
            plot([-6 -6],[-4 -4+2*16/10.2],'k','LineWidth',1)
            %             text(-6.5,-8,'2 cm','FontSize',10)
        end
    end
end

%% Figure 1e
% depends on rotatetrj.m

% monkey G
cd 'E:\manuscript\20230323 manuscript\2024-04-18_G_COINT_rotrj'
load('E:\manuscript\20230323 manuscript\2023-07-17_G_J_COINT_behaviors.mat')

% mongkey L
cd 'E:\manuscript\20230323 manuscript\2024-04-18_L_COINT_rotrj'
load('E:\manuscript\20230323 manuscript\2023-07-17_L_J_COINT_behaviors.mat')

currentFolder = pwd;
namelist = dir('**/*.mat');

figure
set(gcf,'position',[30,300,400,300])
hold on

% plot direction
r = 102;
theta_rad = deg2rad(45);  % 将角度转换为弧度
x0 = r .* cos(theta_rad);  % 计算 x 坐标
y0 = r .* sin(theta_rad);  % 计算 y 坐标
plot([0,x0],[0,y0],':k','linewidth',1)

theta_rad = deg2rad(90+45);
x0 = r .* cos(theta_rad);
y0 = r .* sin(theta_rad);
plot([0,x0],[0,y0],':k','linewidth',1)

theta_rad = deg2rad(180+45);
x0 = r .* cos(theta_rad);
y0 = r .* sin(theta_rad);
plot([0,x0],[0,y0],':k','linewidth',1)

theta_rad = deg2rad(360-45);
x0 = r .* cos(theta_rad);
y0 = r .* sin(theta_rad);
plot([0,x0],[0,y0],':k','linewidth',1)

% plot CI
colo = [0	0.69804	0.93333;1	0.49804	0;0.06275	0.30588	0.5451;0.5451	0.1451	0;...
    0	0.69804	0.93333;1	0.49804	0;0.06275	0.30588	0.5451;0.5451	0.1451	0;...
    0	0.69804	0.93333;1	0.49804	0;0.06275	0.30588	0.5451;0.5451	0.1451	0;...
    0	0.69804	0.93333;1	0.49804	0;0.06275	0.30588	0.5451;0.5451	0.1451	0];
linestyle = {'-','-','-','-','-','-','-','-','-','-','-','-','-','-','-','-'};
n = 1;s_name=[];
for s = 1:length(namelist)
    if ~isempty(strfind(namelist(s).name,'rotrj_ci')) %&& ~isempty(strfind(namelist(ses).name,'0815'))
        tic
        load(namelist(s).name,'x','y','z')
        s_name{n,1} = namelist(s).name;
        
        z = squeeze(z);
        [C,h]=contour(x,y,z);
        h.LineColor = colo(n,:);
        h.LineStyle = linestyle{n};
        h.Visible = 'off';
        npoints = C(2,1);
        fill(C(1,2:npoints+1),C(2,2:npoints+1),colo(n,:),'EdgeColor','none','FaceAlpha',0.1) % 填充
        
        n = n+1;
        toc
    end
end

axis off
axis square

% get mean
XYr = [];XYr_trial = [];
for q = 1:4
    for c = 1:4
        
        n = 1;tjC = [];
        for k = 1:size(J_tran_group{q,c},1)
            for i = 1:length(J_velocity_list)
                if ~isempty(J_velocity_list{i,1})
                    if J_velocity_list{i,1}(1,3) == J_tran_group{q,c}(k,1) && J_velocity_list{i,1}(1,4) == J_tran_group{q,c}(k,32)
                        try
                            x0 = -J_velocity_list{i,3}(1,1);
                            y0 = J_velocity_list{i,3}(1,2);
                            tjC(1:10,n) = -J_velocity_list{i,3}(1:10,1)-x0;
                            tjC(11:20,n) = J_velocity_list{i,3}(1:10,2)-y0;
                            n = n+1;
                            %                                         plot(-velocity_list{i,3}(1:10,1)-x0, velocity_list{i,3}(1:10,2)-y0,'k','linewidth',1.5) ;
                        catch
                        end
                    end
                end
            end
        end
        
        xC = mean(tjC(1:10,:),2,'omitnan'); % mean
        yC = mean(tjC(11:end,:),2,'omitnan');
        %% rotate
        [XYr{q,c},~] = rotatetrj(xC,yC);
        XYr{q,c} = [smoothdata(XYr{q,c}(1,:)','gaussian',4),smoothdata(XYr{q,c}(2,:)','gaussian',4)];
        
        for i = 1:size(tjC,2)
            [XYr_trial_temp,~] = rotatetrj(tjC(1:10,i),tjC(11:20,i));
            XYr_trial{q,c}{1,1}(:,i) = smoothdata(XYr_trial_temp(1,:)','gaussian',4);
            XYr_trial{q,c}{1,2}(:,i) = smoothdata(XYr_trial_temp(2,:)','gaussian',4);
        end
        
    end
end

% plot mean
colo = [0	0.69804	0.93333;1	0.49804	0;0.06275	0.30588	0.5451;0.5451	0.1451	0];
linestyle = {'-','-','-','-'};
for i = 1:4
    for j = 1:4
        if i == 1
            M=[cos(pi/4),-sin(pi/4);sin(pi/4),cos(pi/4)];
        elseif i == 2
            M=[cos(3*pi/4),-sin(3*pi/4);sin(3*pi/4),cos(3*pi/4)];
        elseif i == 3
            M=[cos(3*pi/4),sin(3*pi/4);-sin(3*pi/4),cos(3*pi/4)];
        elseif i == 4
            M=[cos(pi/4),sin(pi/4);-sin(pi/4),cos(pi/4)];
        end
        
        R1=[XYr{i,j}(:,1),XYr{i,j}(:,2)];
        R2=M*R1';%旋转后坐标
        plot(R2(1,:),R2(2,:),'color',colo(j,:),'LineStyle',linestyle{j},'linewidth',2)
        
    end
end

plot([-50 -30],[-70 -70],'k','LineWidth',1)
% text(-50,-80,'2 cm','FontSize',14)

%% Figure 1f

% monkey G
load('E:\manuscript\20230323 manuscript\2023-07-17_G_J_COINT_behaviors.mat')

% monkey L
load('E:\manuscript\20230323 manuscript\2023-07-17_L_J_COINT_behaviors.mat')

tran_v = J_tran_group;
trj = [];
t = -200:10:200;

figure
set(gcf,'unit','centimeters','position',[10 5 8 10]);
pos = [2,1,3,4];
colo = [0	0.69804	0.93333;1	0.49804	0;0.06275	0.30588	0.5451;0.5451	0.1451	0];
linestyle = {'-','-','-','-'};
for a = 1:4
    subplot(2,2,pos(a))
    hold on
    for b = 1:4
        n = 1;
        tj = [];
        for k = 1:size(tran_v{a,b})
            for i = 1:length(J_velocity_list)
                if ~isempty(J_velocity_list{i,1})
                    try
                        if J_velocity_list{i,1}(1,3) == tran_v{a,b}(k,1) && J_velocity_list{i,1}(1,4) == tran_v{a,b}(k,32)
                            tj(:,n) = J_velocity_list{i,2}(:,1);
                            n = n+1;
                        end
                    catch
                    end
                end
            end
        end
        
        exceed_indices = any(tj > 1000, 1);
        tj = tj(:, ~exceed_indices);
        
        if a == 1 || a == 2
            tj(21,:) = NaN;
            nanIndices = isnan(tj);
            vq2 = tj;
            vq2(nanIndices) = interp1(find(~nanIndices), tj(~nanIndices), find(nanIndices), 'spline');
            tj = vq2;
        end
        
        trj{a,b}=tj';
        
        y0 = mean(tj,2,'omitnan');
        y1 = y0';
        plot(t,y1,'color',colo(b,:),'LineStyle',linestyle{b},'LineWidth',1.5); % black
        
        statFunc = @(x) mean(x,'omitNaN');
        ci = bootci(1000, {statFunc, tj'}, 'alpha', 0.05, 'type', 'percentile');
        g=fill([t,t(end:-1:1)],[ci(1,:),ci(2,end:-1:1)],'b');
        set(g,'FaceColor',colo(b,:),'FaceAlpha',0.2,'EdgeColor','none');
        
        % %         axis([-200 200 0 120]) % L
        axis([-200 200 0 140]) % G
        axis off
        if a == 3
            %             plot([-200 -200],[50 75],'k','LineWidth',0.8)
            %             plot([-200 -140],[50 50],'k','LineWidth',0.8)
            plot([-200 -200],[20 40],'k','LineWidth',2)
            plot([0 100],[0 0],'k','LineWidth',2)
            
            xlabel('Time to Peak velocity (ms)')
            ylabel('Hand speed (cm/s)')
        end
        
    end
end

h1 = subplot(2,2,1);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1), pos1(2)*0.95, pos1(3)*1.1, pos1(4)]); % [left bottom width height]

h1 = subplot(2,2,2);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1), pos1(2)*0.95, pos1(3)*1.1, pos1(4)]); % [left bottom width height]

h1 = subplot(2,2,3);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1), pos1(2), pos1(3)*1.1, pos1(4)]); % [left bottom width height]

h1 = subplot(2,2,4);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1), pos1(2), pos1(3)*1.1, pos1(4)]); % [left bottom width height]


% %% difference between static and moving
% for i = 1:size(trj,1)
%     for k = 1:size(trj{1,1},2)
%         [~,h_s(i,k)] = ranksum(trj{i,1}(:,k),trj{i,2}(:,k));
%         [~,h_l(i,k)] = ranksum(trj{i,3}(:,k),trj{i,4}(:,k));
%         [dif_s(i,k), ciLower_s(i,k), ciUpper_s(i,k)] = bootstrapDifference(trj{i,2}(:,k), trj{i,1}(:,k), 1000);
%         [dif_l(i,k), ciLower_l(i,k), ciUpper_l(i,k)] = bootstrapDifference(trj{i,4}(:,k), trj{i,3}(:,k), 1000);
%     end
% end
%
% figure
% set(gcf,'unit','centimeters','position',[10 5 7 10]);
% colo = [0.62745	0.32157	0.17647; 0.6	0.19608	0.8];
% for i = 1:4
%     subplot(2,2,pos(i))
%     hold on
%
%     plot(t,dif_s(i,:),'color',colo(1,:),'LineStyle',linestyle{b},'LineWidth',1.5); % black
%     plot(t,dif_l(i,:),'color',colo(2,:),'LineStyle',linestyle{b},'LineWidth',1.5); % black
%
%     g=fill([t,t(end:-1:1)],[ciLower_s(i,:),ciUpper_s(i,end:-1:1)],'b');
%     set(g,'FaceColor',colo(1,:),'FaceAlpha',0.2,'EdgeColor','none');
%
%     g=fill([t,t(end:-1:1)],[ciLower_l(i,:),ciUpper_l(i,end:-1:1)],'b');
%     set(g,'FaceColor',colo(2,:),'FaceAlpha',0.2,'EdgeColor','none');
%
%     y0 = ylim;
%     for k = 1:size(h_s,2)
%         if h_s(i,k)==1
%             plot([t(k)-5 t(k)+5],[y0(2) y0(2)],'color',colo(1,:),'LineWidth',2)
%         end
%     end
%     for k = 1:size(h_l,2)
%         if h_l(i,k)==1
%             plot([t(k)-5 t(k)+5],[y0(2)-1 y0(2)-1],'color',colo(2,:),'LineWidth',2)
%         end
%     end
%
%     plot([-200 200],[0 0],'color',[0.5 0.5 0.5],'LineWidth',0.8)
%     axis off
% %     if a == 3
% %         plot([-200 -200],[50 75],'k','LineWidth',0.8)
% %         plot([-200 -140],[50 50],'k','LineWidth',0.8)
% %
% %         xlabel('Time to Peak velocity (ms)')
% %         ylabel('Hand speed (cm/s)')
% %     end
% end

%% Figure 1g

% monkey G
load('E:\manuscript\20230323 manuscript\2023-07-17_G_J_COINT_behaviors.mat')

% % monkey L
% load('E:\manuscript\20230323 manuscript\2023-07-17_L_J_COINT_behaviors.mat')

para = 13; % RT
tran_para = cellfun(@(x) (x(:,para)),J_tran_4,'UniformOutput',0);
tran_group_para = cellfun(@(x) (x(:,para)),J_tran_group,'UniformOutput',0);

tran_para4 = nan(2000,4);
for i=1:4
    tran_para4(1:size(tran_para{1,i},1),i) = tran_para{1,i};
end

n = 1;
tran_para4_group = nan(500,16);
for i=1:4
    for j = 1:4
        tran_para4_group(1:size(tran_group_para{i,j},1),n) = tran_group_para{i,j};
        n = n+1;
    end
end

% cdf
colo = [0	0.69804	0.93333;1	0.49804	0;0.06275	0.30588	0.5451;0.5451	0.1451	0];
figure
ax = axes;
for i = 1:4
    x = tran_para{1,i};
    h1=cdfplot_my(x);
    set(h1,'LineStyle', '-', 'Color', colo(i,:),'LineWidth',1) % CO short
    hold on
    scatter(median(x),0.05,30,'o','MarkerEdgeColor',colo(i,:))
end

xlim([100 400])
box off; grid off
ylabel('cdf','FontName','Arial');
xlabel('reaction time (ms)','FontName','Arial');
ax.TickDir='out';
set(gcf,'unit','centimeters','position',[5 10 15 4]);
set(gca,'Position',[.2 .2 .7 .7]);
offsetaxis(ax, 'y', 0.02); % Default, y-offset = 0.1 and yloc = 'l'
offsetaxis(ax, 'x', 0.02); % Default, y-offset = 0.1 and yloc = 'l'
xticks([100,200,300,400,500]);

A = [];
proportion_150ms = [];
for i = 1:4
    A{1,i} = J_tran_4{1,i}(J_tran_4{1,i}(:,13)<150,:);
    proportion_150ms(i)=length(find(A{1,i}(:,19)~=1))*100/size(A{1,i},1);
end

%% Figure 2f
% G
[num,~,~] = xlsread('E:\manuscript\20230323 manuscript\manuscript table 202301203 correct session median.xlsx',1);

map1 = 1:1:18;
map2 = 18:-1:1;

num2 = [];
for i = 1:length(find(~isnan(num(:,6))))
    num2(i) = map2(map1==num(i,6));
end

radius_static = [12.34;-19.78;-4.85;11.73;25.99;27.36;18.64;35.59;0.79;19.88;22.23;20.01;15.60;6.75];
radius_moving = [-14.69;-34.48;-3.16;14.71;-19.94;5.22;-44.54;20.73;-31.39;-33.92;-29.86;13.32;-15.89;-36.19];
% radius = radius_static;
radius = radius_moving;

colo = [0	0.60392	0.80392;0.81961	0.37255	0.93333];
figure
set(gcf,'Position',[100,300,500,500])
hold on
for i = 1:size(radius,1)
    if radius(i)<0
        scatter3(num2(i),num(i,7),-num(i,10),abs(radius(i))*10,'filled','markerfacecolor',colo(1,:),'markeredgecolor','k')
    else
        scatter3(num2(i),num(i,7),-num(i,10),abs(radius(i))*10,'filled','markerfacecolor',colo(2,:),'markeredgecolor','k')
    end
    plot3([num2(i),num2(i)],[num(i,7),num(i,7)],[0,-num(i,10)],'k-','linewidth',1.5)
    plot3([num2(i),num2(i)],[num(i,7),num(i,7)],[-num(i,10),-1.5],'k:','linewidth',1)
end

plot3([17.5,18.5],[19 19],[-1.5,-1.5],'Color',[0 0 0],'linewidth',2)

% img2 = imread('E:\manuscript\20230323 manuscript\Golden chamber MRI flip.jpg');
img2 = imread('E:\manuscript\20230323 manuscript\version NC revision-1\example codes\G brain chamber.jpg');
xImage = [2 20; 2 20];  % upperleft, upperright; bottomleft, bottomright
yImage = [20 20; 2 2];
zImage = [-1.5 -1.5; -1.5 -1.5];
surf(xImage,yImage,zImage,'CData',img2, 'FaceColor','texturemap','facealpha',0.5);% Plot the bottom

axis off
axis square

% L
[num,txt,raw] = xlsread('E:\manuscript\20230323 manuscript\manuscript table 202301203 correct session median.xlsx',2);

radius_static = ([9.37	16.41	15.29	28.87	16.89	40.30	1.64	-7.30	-2.40	29.85	53.52	11.64 26.98 7.45 10.05])';
radius_moving = ([10.36	9.33	-9.57	-16.86	14.95	-14.21	-48.37	-30.09	-4.06	-44.49	-3.89	-13.90 -29.47 -17.21 -30.54])';
% radius = radius_static;
radius = radius_moving;

colo = [0	0.60392	0.80392;0.81961	0.37255	0.93333];
figure
set(gcf,'Position',[100,300,500,500])
hold on
for i = 1:size(radius,1)
    if radius(i)<0
        scatter3(num(i,6),num(i,7),-num(i,10),abs(radius(i))*10,'filled','markerfacecolor',colo(1,:),'markeredgecolor','k')
    else
        scatter3(num(i,6),num(i,7),-num(i,10),abs(radius(i))*10,'filled','markerfacecolor',colo(2,:),'markeredgecolor','k')
    end
    plot3([num(i,6),num(i,6)],[num(i,7),num(i,7)],[0,-num(i,10)],'k-','linewidth',1.5)
    plot3([num(i,6),num(i,6)],[num(i,7),num(i,7)],[-num(i,10),-2],'k:','linewidth',1)
end

plot3([15,16],[17,17],[-2,-2],'Color',[0 0 0],'linewidth',2)

img2 = imread('E:\manuscript\20230323 manuscript\version NC revision-1\example codes\L brain chamber.jpg');
xImage = [0 18; 0 18];   % The x data for the image corners
yImage = [18 18; 0 0];   % The y data for the image corners
zImage = [-2 -2; -2 -2];   % The z data for the image corners
surf(xImage,yImage,zImage,'CData',img2, 'FaceColor','texturemap','facealpha',0.5);% Plot the bottom

axis off
axis square

%% Figure 2g
% location_ICMS_20231212.pzfx

