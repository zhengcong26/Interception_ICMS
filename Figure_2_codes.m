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