
%% Figure 1c
% Uncomment the following line to load the data file and set the trial constraint value
% load('G_0914_J_bin_20msBin_4s.mat'); trialconstraint = 675; % monkey G
% load('L_0818_J_bin_20msBin_4s.mat'); trialconstraint = 644; % monkey L

% designed static condition
CO = stim_list(stim_list(:,2) == 0 & stim_list(:,3) == 0 & ...
               stim_list(:,19) ~= 1 & stim_list(:,1) < trialconstraint, 1);

% designed moving condition
INT = stim_list(stim_list(:,2) == 0 & stim_list(:,3) ~= 0 & ...
                stim_list(:,19) ~= 1 & stim_list(:,1) < trialconstraint, 1);

% random condition
RAND = stim_list(stim_list(:,2) == 0 & stim_list(:,19) == 1 & ...
                 stim_list(:,1) < trialconstraint, 1);

% Set up the figure window with specified position and size
figure;
set(gcf, 'Position', [30, 300, 500, 500]);
colo = [0.06275, 0.30588, 0.5451;   
        0.5451,  0.1451,  0;        
        0.7,     0.7,     0.7];

hold on; 
plot_session_structure(RAND, velocity_list, colo(3,:)); % Plot for RAND
plot_session_structure(CO, velocity_list, colo(1,:));   % Plot for CO
plot_session_structure(INT, velocity_list, colo(2,:));  % Plot for INT

axis equal;
axis off;

% Add category labels to the trial data and concatenate for sorting
CO(:,2) = 1;   
INT(:,2) = 2; 
RAND(:,2) = 3; 
A = sortrows([CO; INT; RAND], 1); % Concatenate and sort by trial number

% Calculate and display the proportions of each condition
disp(['Proportion of CO trials: ', num2str(length(CO) / length(A))]);
disp(['Proportion of INT trials: ', num2str(length(INT) / length(A))]);
disp(['Proportion of RAND trials: ', num2str(length(RAND) / length(A))]);

%% Figure 1d

% Uncomment the following line to load the data file
% load('G_J_COINT_behaviors.mat')
% load('L_J_COINT_behaviors.mat')

figure;
set(gcf, 'unit', 'centimeters', 'position', [10, 5, 15, 15]);
pos = [2, 1, 3, 4]; % Subplot positions
colo = [0, 0.69804, 0.93333; 1, 0.49804, 0; 0.06275, 0.30588, 0.5451; 0.5451, 0.1451, 0];
linestyle = {'-', '-', '-', '-'};

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

%% Figure 1e

% monkey G
load('G_J_COINT_behaviors.mat')
load('G_COINT_rotrj.mat')

% % mongkey L
load('L_J_COINT_behaviors.mat')
load('L_COINT_rotrj.mat')

figure
set(gcf, 'position', [30, 300, 400, 300])
hold on

% Plot directional lines at 45-degree intervals
r = 102;
angles_deg = [45, 135, 225, 315];  % Define angles in degrees
for theta_deg = angles_deg
    theta_rad = deg2rad(theta_deg);  % Convert angle to radians
    x0 = r * cos(theta_rad);
    y0 = r * sin(theta_rad);
    plot([0, x0], [0, y0], ':k', 'linewidth', 1);  % Plot each directional line
end

% Define colors and line styles for confidence intervals (CI)
colo = repmat([0 0.698 0.933; 1 0.498 0; 0.063 0.306 0.545; 0.545 0.145 0], 4, 1);
linestyle = repmat({'-'}, 1, 16);

% Plot confidence intervals for each dataset
for n = 1:size(X, 1)
    x = squeeze(X(n, :, :));
    y = squeeze(Y(n, :, :));
    z = squeeze(Z(n, :, :));
    [C, h] = contour(x, y, z);  % Calculate contour
    h.LineColor = colo(n, :);
    h.LineStyle = linestyle{n};
    h.Visible = 'off';
    npoints = C(2, 1);
    fill(C(1, 2:npoints+1), C(2, 2:npoints+1), colo(n, :), 'EdgeColor', 'none', 'FaceAlpha', 0.1); % Fill CI area
end

axis off
axis square

% Initialize containers for rotated trajectories
XYr = [];
XYr_trial = [];

% Loop over groups and conditions to compute rotated mean trajectories
for q = 1:4
    for c = 1:4
        n = 1; tjC = [];
        for k = 1:size(J_tran_group{q, c}, 1)
            for i = 1:length(J_velocity_list)
                if ~isempty(J_velocity_list{i, 1}) && ...
                   J_velocity_list{i, 1}(1, 3) == J_tran_group{q, c}(k, 1) && ...
                   J_velocity_list{i, 1}(1, 4) == J_tran_group{q, c}(k, 32)
                    try
                        % Calculate relative trajectory with origin adjusted
                        x0 = -J_velocity_list{i, 3}(1, 1);
                        y0 = J_velocity_list{i, 3}(1, 2);
                        tjC(1:10, n) = -J_velocity_list{i, 3}(1:10, 1) - x0;
                        tjC(11:20, n) = J_velocity_list{i, 3}(1:10, 2) - y0;
                        n = n + 1;
                    catch
                        % Skip if error
                    end
                end
            end
        end
        
        % Calculate and smooth mean trajectory
        xC = mean(tjC(1:10, :), 2, 'omitnan');
        yC = mean(tjC(11:end, :), 2, 'omitnan');
        [XYr{q, c}, ~] = rotatetrj(xC, yC);
        XYr{q, c} = [smoothdata(XYr{q, c}(1, :)', 'gaussian', 4), smoothdata(XYr{q, c}(2, :)', 'gaussian', 4)];
        
        % Rotate and smooth individual trial trajectories
        for i = 1:size(tjC, 2)
            [XYr_trial_temp, ~] = rotatetrj(tjC(1:10, i), tjC(11:20, i));
            XYr_trial{q, c}{1, 1}(:, i) = smoothdata(XYr_trial_temp(1, :)', 'gaussian', 4);
            XYr_trial{q, c}{1, 2}(:, i) = smoothdata(XYr_trial_temp(2, :)', 'gaussian', 4);
        end
    end
end

% Plot rotated mean trajectories for each group and condition
colo = [0 0.698 0.933; 1 0.498 0; 0.063 0.306 0.545; 0.545 0.145 0];
linestyle = {'-', '-', '-', '-'};
for i = 1:4
    for j = 1:4
        % Define rotation matrix based on quadrant
        if i == 1
            M = [cos(pi/4), -sin(pi/4); sin(pi/4), cos(pi/4)];
        elseif i == 2
            M = [cos(3*pi/4), -sin(3*pi/4); sin(3*pi/4), cos(3*pi/4)];
        elseif i == 3
            M = [cos(3*pi/4), sin(3*pi/4); -sin(3*pi/4), cos(3*pi/4)];
        else
            M = [cos(pi/4), sin(pi/4); -sin(pi/4), cos(pi/4)];
        end
        
        % Apply rotation and plot mean trajectory
        R1 = [XYr{i, j}(:, 1), XYr{i, j}(:, 2)];
        R2 = M * R1'; % Rotate coordinates
        plot(R2(1, :), R2(2, :), 'color', colo(j, :), 'LineStyle', linestyle{j}, 'linewidth', 2);
    end
end

% Plot reference line for scale
plot([-50, -30], [-70, -70], 'k', 'LineWidth', 1);

%% Figure 1f

% monkey G
% load('G_J_COINT_behaviors.mat')

% % monkey L
load('L_J_COINT_behaviors.mat')

% Define the velocity grouping and time vector
tran_v = J_tran_group;
trj = [];
t = -200:10:200;  % Time axis

% Initialize figure with custom settings
figure
set(gcf, 'unit', 'centimeters', 'position', [10 5 8 10]);
pos = [2, 1, 3, 4];
colo = [0 0.69804 0.93333; 1 0.49804 0; 0.06275 0.30588 0.5451; 0.5451 0.1451 0];
linestyle = {'-', '-', '-', '-'};

% Loop over each subplot
for a = 1:4
    subplot(2, 2, pos(a))
    hold on
    for b = 1:4
        n = 1;
        tj = [];
        
        % Collect trajectories for each condition
        for k = 1:size(tran_v{a, b}, 1)
            for i = 1:length(J_velocity_list)
                if ~isempty(J_velocity_list{i, 1})
                    try
                        % Check conditions and gather trajectory data
                        if J_velocity_list{i, 1}(1, 3) == tran_v{a, b}(k, 1) && ...
                           J_velocity_list{i, 1}(1, 4) == tran_v{a, b}(k, 32)
                            tj(:, n) = J_velocity_list{i, 2}(:, 1);
                            n = n + 1;
                        end
                    catch
                        % Skip if an error occurs
                    end
                end
            end
        end
        
        % Remove trajectories with outlier values (>1000)
        exceed_indices = any(tj > 1000, 1);
        tj = tj(:, ~exceed_indices);
        
        % Fill missing data in middle segment for groups 1 and 2 using spline interpolation
        if a == 1 || a == 2
            tj(21, :) = NaN;  % Introduce NaN at center
            nanIndices = isnan(tj);
            tj(nanIndices) = interp1(find(~nanIndices), tj(~nanIndices), find(nanIndices), 'spline');
        end
        
        % Store the cleaned trajectory data
        trj{a, b} = tj';
        
        % Compute mean trajectory for plotting
        y0 = mean(tj, 2, 'omitnan');
        y1 = y0';
        plot(t, y1, 'color', colo(b, :), 'LineStyle', linestyle{b}, 'LineWidth', 1.5);
        
        % Bootstrap for confidence interval (CI) calculation
        statFunc = @(x) mean(x, 'omitNaN');
        ci = bootci(1000, {statFunc, tj'}, 'alpha', 0.05, 'type', 'percentile');
        g = fill([t, t(end:-1:1)], [ci(1, :), ci(2, end:-1:1)], 'b');
        set(g, 'FaceColor', colo(b, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        
        % Set plot limits and format
        axis([-200 200 0 140])  % Adjusted for monkey G
        axis off
        if a == 3
            % Add scale indicators
            plot([-200 -200], [20 40], 'k', 'LineWidth', 2);
            plot([0 100], [0 0], 'k', 'LineWidth', 2);
            xlabel('Time to Peak velocity (ms)')
            ylabel('Hand speed (cm/s)')
        end
    end
end

% Adjust subplot positions for visual alignment
for idx = 1:4
    h = subplot(2, 2, idx);
    pos = get(h, 'Position');
    set(h, 'Position', [pos(1), pos(2) * 0.95, pos(3) * 1.1, pos(4)]);  % Adjust each subplot
end

%% Figure 1g

% monkey G
% load('G_J_COINT_behaviors.mat')

% monkey L
load('L_J_COINT_behaviors.mat')

% Set parameter index to Reaction Time (RT)
para = 13;  % Reaction Time (RT) is in column 13

% Extract RT data for each group in J_tran_4 and J_tran_group
tran_para = cellfun(@(x) (x(:, para)), J_tran_4, 'UniformOutput', 0);
tran_group_para = cellfun(@(x) (x(:, para)), J_tran_group, 'UniformOutput', 0);

% Combine RT data into a matrix (4 groups for J_tran_4)
tran_para4 = nan(2000, 4);  % Initialize with NaNs for padding
for i = 1:4
    tran_para4(1:size(tran_para{1, i}, 1), i) = tran_para{1, i};
end

% Combine RT data into a matrix for each sub-group in J_tran_group
n = 1;
tran_para4_group = nan(500, 16);  % 16 sub-groups total (4 x 4)
for i = 1:4
    for j = 1:4
        tran_para4_group(1:size(tran_group_para{i, j}, 1), n) = tran_group_para{i, j};
        n = n + 1;
    end
end

% Plot CDF of Reaction Time (RT) for each main group
colo = [0 0.69804 0.93333; 1 0.49804 0; 0.06275 0.30588 0.5451; 0.5451 0.1451 0]; % Colors for each group
figure
ax = axes;
for i = 1:4
    x = tran_para{1, i};  % Extract RT data for current group
    h1 = cdfplot_my(x);  % Custom CDF plot function
    set(h1, 'LineStyle', '-', 'Color', colo(i, :), 'LineWidth', 1);  % Customize line style and color
    hold on
    
    % Add a scatter plot for the median RT value
    scatter(median(x), 0.05, 30, 'o', 'MarkerEdgeColor', colo(i, :))
end

% Customize plot limits and appearance
xlim([100 400]);
box off; grid off;
ylabel('CDF', 'FontName', 'Arial');
xlabel('Reaction time (ms)', 'FontName', 'Arial');
ax.TickDir = 'out';
set(gcf, 'unit', 'centimeters', 'position', [5 10 15 4]);
set(gca, 'Position', [.2 .2 .7 .7]);
offsetaxis(ax, 'y', 0.02);  % Offset y-axis for better readability
offsetaxis(ax, 'x', 0.02);  % Offset x-axis for better readability
xticks([100, 200, 300, 400, 500]);


