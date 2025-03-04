
%% Figure 2c

% Uncomment the following line to load the data file
% load('G_J_NSST_behavior.mat')
load('L_J_NSST_behavior.mat')

tran1 = cellfun(@(x,y) ([x;y]),J_GO_raw(:,1:4),J_GO_raw(:,5:8),'UniformOutput',0);
tran2 = cellfun(@(x,y) ([x;y]),J_MO_raw(:,1:4),J_MO_raw(:,5:8),'UniformOutput',0);
J_tran_group = tran1;
J_tran_group(:,2) = cellfun(@(x,y) ([x;y]),J_tran_group(:,2),tran2(:,2),'UniformOutput',0);
J_tran_group(:,4) = cellfun(@(x,y) ([x;y]),J_tran_group(:,4),tran2(:,4),'UniformOutput',0);

figure;
set(gcf, 'unit', 'centimeters', 'position', [10, 5, 15, 15]);
pos = [2, 1, 3, 4]; % Subplot positions
colo = [0.4	0.80392	0.66667;0.4	0.80392	0.66667;0.4902	0.14902	0.80392;0.4902	0.14902	0.80392];
linestyle = {':','-',':','-'};

for q = 1:4
    subplot(2, 2, pos(q)); hold on;
    
    % Define rotation matrix M based on the quadrant
    if q == 1
        M = [cos(pi/4), sin(pi/4); -sin(pi/4), cos(pi/4)];
    elseif q == 2
        M = [cos(pi/4), -sin(pi/4); sin(pi/4), cos(pi/4)];
    elseif q == 3
        M = [cos(3*pi/4), -sin(3*pi/4); sin(3*pi/4), cos(3*pi/4)];
    else
        M = [cos(3*pi/4), sin(3*pi/4); -sin(3*pi/4), cos(3*pi/4)];
    end
    
    % Plot reference arc lines with rotation
    R1 = [[-5, 5]; [0, 0]]; R2 = M * R1;
    plot(R2(1, :), R2(2, :), ':k', 'LineWidth', 0.8);
    R1 = [[0, 0]; [-5, 5]]; R2 = M * R1;
    plot(R2(1, :), R2(2, :), ':k', 'LineWidth', 0.8);
    
    % Plot arc for trajectory
    the = pi/4:pi/180:3*pi/4; r = 16;
    x = r * cos(the); y = -16 + r * sin(the);
    R1 = [x; y]; R2 = M * R1;
    plot(R2(1, :), R2(2, :), ':k', 'LineWidth', 1);
    
    % Plot error tolerance region
    the = 0:pi/180:2*pi; r = 4.5;
    x = r * cos(the); y = r * sin(the);
    R1 = [x; y]; R2 = M * R1;
    fill(R2(1, :), R2(2, :), 'b', 'FaceAlpha', 0.08, 'EdgeColor', 'none');
    
    % Plot data for each condition
    for c = 1:4
        X = J_tran_group{q, c}(:, 24);
        Y = J_tran_group{q, c}(:, 25);
        
        % Plot error ellipse
        [~,~,~,X0,Y0,r_ellipse] = ErrorEllipse([X, Y], 0.95, colo, '-', 'k', 1, 0);
        R1 = [r_ellipse(:, 1) + X0, r_ellipse(:, 2) + Y0]; R2 = M * R1';
        plot(R2(1, :), R2(2, :), 'Color', colo(c, :), 'LineStyle', linestyle{c}, 'LineWidth', 2.2);
        
        % Plot center point
        R1 = [X0; Y0]; R2 = M * R1;
        scatter(R2(1), R2(2), 70, 'b', 'd', 'MarkerEdgeColor', colo(c, :), 'LineWidth', 2.2);
        
        % Set axis properties
        axis([-7 7 -7 7]); axis square; box off; axis off;
        set(gca, 'LineWidth', 2);
        
        % Scale bar in bottom-left plot
        if q == 3
            plot([-6, -6], [-4, -4 + 2 * 16 / 10.2], 'k', 'LineWidth', 1);
        end
    end
end

C = J_tran_group;
col_indices = [24, 25]; 
C = C(:);  
max_rows = max(cellfun(@(x) size(x,1), C));
num_cells = numel(C);
result_matrix = NaN(max_rows, num_cells * length(col_indices)); 
for i = 1:num_cells
    data = C{i}(:, col_indices); 
    rows = size(data, 1); 
    result_matrix(1:rows, (i-1)*2 + (1:2)) = data;
end

%% Figure 2e

% monkey G
load('sesbehaviors_G_14ses.mat')

% monkey L
% load('sesbehaviors_L_15ses.mat')

% Main scatter plot
figure;
set(gcf, 'unit', 'centimeters', 'position', [10 5 4.5 4.5]);
hold on;

% Plot reference lines for x=0 and y=0
plot([0 0], [-150 150], 'k:', 'linewidth', 1);
plot([-150 150], [0 0], 'k:', 'linewidth', 1);

% Scatter and error ellipse for GO condition
scatter(x_GO, y_GO, 25, 'o', 'MarkerEdgeColor', [1 0.64706 0.3098]);
ErrorEllipse([x_GO, y_GO], 0.95, [1 0.64706 0.3098], '-', 'k', 1, 1);

% Scatter and error ellipse for MO condition
scatter(x_MO, y_MO, 25, 'o', 'MarkerEdgeColor', [0.5451 0.27059 0.07451]);
ErrorEllipse([x_MO, y_MO], 0.95, [0.5451 0.27059 0.07451], '-', 'k', 1, 1);

% Set axis properties
axis square;
box on;
xlim([-70, 70]);
ylim([-70, 70]);

% Beautify the axis
h = gca;
h.FontSize = 8;
h.TickLength = [0.02, 0.25];
xticks([-60 -30 0 30 60]);
yticks([-60 -30 0 30 60]);

% Scatter plot with histogram (scatterhist)
figure;
set(gcf, 'unit', 'centimeters', 'position', [10 5 6.35 6.35]);
x = [x_GO; x_MO];
y = [y_GO; y_MO];

% Create grouping labels
group = cell(length(x), 1);
group(1:length(x_GO)) = {'1'};  % GO group
group(length(x_GO)+1:end) = {'2'};  % MO group

% Scatter plot with histogram and kernel density estimate
s = scatterhist(x, y, 'Group', group, 'Kernel', 'on', 'Location', 'NorthEast', ...
    'Direction', 'out', 'Color', [1 0.64706 0.3098; 0.5451 0.27059 0.07451], ...
    'LineStyle', {'-', '-'}, 'LineWidth', [1, 1]);

% Set limits for each subplot
for i = 1:3
    s(i).XLim = [-70, 70];
end
s(1).YLim = [-70, 70];
s(2).YLim = [0, 0.06];
s(3).YLim = [0, 0.06];

axis square;

%% Figure 2f
% G
% Load session information data
[num,~,~] = xlsread('session information table.xlsx',1);

% Define mapping from session indices to positions in the 3D plot
map1 = 1:1:18;
map2 = 18:-1:1;

% Map session indices based on table entries
num2 = []; n = 1;
for i = 1:size(num,1)
    if ~isnan(num(i,6))
        num2(i) = map2(map1 == num(i,6));
    end
end

% deltaRT
radius_static = [12.34; -19.78; -4.85; 11.73; 25.99; 27.36; 18.64; 35.59; 0.79; 19.88; 22.23; 20.01; 15.60; 6.75];
radius_moving = [-14.69; -34.48; -3.16; 14.71; -19.94; 5.22; -44.54; 20.73; -31.39; -33.92; -29.86; 13.32; -15.89; -36.19];

% Choose which data to use (static or moving)
radius = radius_static;
% radius = radius_moving;

colo = [0 0.60392 0.80392; 0.81961 0.37255 0.93333];
figure;
set(gcf, 'Position', [100, 300, 500, 500]);
hold on;

% Loop through each radius and plot it in 3D space
for i = 1:size(radius, 1)
    % Determine color based on radius sign (positive or negative)
    if radius(i) < 0
        scatter3(num2(i), num(i,7), -num(i,10), abs(radius(i)) * 10, 'filled', ...
                 'MarkerFaceColor', colo(1,:), 'MarkerEdgeColor', 'k');
    else
        scatter3(num2(i), num(i,7), -num(i,10), abs(radius(i)) * 10, 'filled', ...
                 'MarkerFaceColor', colo(2,:), 'MarkerEdgeColor', 'k');
    end
    
    % Draw a line from each point to the z-axis plane (z = 0)
    plot3([num2(i), num2(i)], [num(i,7), num(i,7)], [0, -num(i,10)], 'k-', 'LineWidth', 1.5);
    % Dashed line to extend downward from each point
    plot3([num2(i), num2(i)], [num(i,7), num(i,7)], [-num(i,10), -1.5], 'k:', 'LineWidth', 1);
end

% Add scale line as a visual reference
plot3([17.5, 18.5], [19, 19], [-1.5, -1.5], 'Color', [0 0 0], 'LineWidth', 2);

% Load and display an image of the brain chamber on a specific plane
img2 = imread('G brain chamber.jpg');
xImage = [2 20; 2 20];   % Upper left, upper right; bottom left, bottom right coordinates
yImage = [20 20; 2 2];
zImage = [-1.5 -1.5; -1.5 -1.5];  % Set z coordinates for the image plane
surf(xImage, yImage, zImage, 'CData', img2, 'FaceColor', 'texturemap', 'FaceAlpha', 0.5);
axis off;  
axis square;  

% L
[num, txt, raw] = xlsread('session information table.xlsx', 2);

% deltaRT
radius_static = ([9.37, 16.41, 15.29, 28.87, 16.89, 40.30, 1.64, -7.30, -2.40, 29.85, 53.52, 11.64, 26.98, 7.45, 10.05])';
radius_moving = ([10.36, 9.33, -9.57, -16.86, 14.95, -14.21, -48.37, -30.09, -4.06, -44.49, -3.89, -13.90, -29.47, -17.21, -30.54])';

% Choose which radius data to use (static or moving)
% radius = radius_static;
radius = radius_moving;

colo = [0 0.60392 0.80392; 0.81961 0.37255 0.93333];
figure;
set(gcf, 'Position', [100, 300, 500, 500]);
hold on;

% Loop through each data point and plot in 3D space
for i = 1:size(radius, 1)
    % Determine color based on sign of the radius (positive or negative)
    if radius(i) < 0
        % Plot negative radius values in color specified by colo(1,:)
        scatter3(num(i,6), num(i,7), -num(i,10), abs(radius(i)) * 10, 'filled', ...
                 'MarkerFaceColor', colo(1,:), 'MarkerEdgeColor', 'k');
    else
        % Plot positive radius values in color specified by colo(2,:)
        scatter3(num(i,6), num(i,7), -num(i,10), abs(radius(i)) * 10, 'filled', ...
                 'MarkerFaceColor', colo(2,:), 'MarkerEdgeColor', 'k');
    end
    
    % Plot vertical line from each point to z = 0 plane
    plot3([num(i,6), num(i,6)], [num(i,7), num(i,7)], [0, -num(i,10)], 'k-', 'LineWidth', 1.5);
    % Dashed line extending downward from each point to z = -2
    plot3([num(i,6), num(i,6)], [num(i,7), num(i,7)], [-num(i,10), -2], 'k:', 'LineWidth', 1);
end

% Add a scale bar for reference
plot3([15, 16], [17, 17], [-2, -2], 'Color', [0 0 0], 'LineWidth', 2);

% Load and display an image (e.g., MRI scan) on the bottom of the 3D plot
img2 = imread('L brain chamber.jpg');
xImage = [0 18; 0 18];   % X data for the image corners
yImage = [18 18; 0 0];   % Y data for the image corners
zImage = [-2 -2; -2 -2]; % Z data for the image corners (places image on z = -2 plane)
surf(xImage, yImage, zImage, 'CData', img2, 'FaceColor', 'texturemap', 'FaceAlpha', 0.5);
axis off;
axis square;
