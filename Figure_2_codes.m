
%% Figure 2e

% monkey G
load('sesbehaviors_G_14ses.mat')

% monkey L
% load('sesbehaviors_L_12ses')

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
% radius = radius_static;
radius = radius_moving;

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
