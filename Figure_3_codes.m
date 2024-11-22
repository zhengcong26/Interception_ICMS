%% Figure 3a
% plot_units_heatmap.py

%% Figure 3bc 

% panel 3b
% load('G_0914_PC_NSST_PreGO.mat')

% panel 3c
load('G_0914_PC_NSST_PreMO.mat')

t_start = -200;
t_end = 400;
bin = 20;

ses = 1;
para = 4;
q = 1;

% PC
CO_NS_PC = squeeze(J_C_m(ses,para,:,:));
CO_ST_PC = squeeze(J_C_m(ses,para+4,:,:));
INT_NS_PC = squeeze(J_C_m(ses,para+8,:,:));
INT_ST_PC = squeeze(J_C_m(ses,para+12,:,:));

% CI
CO_NS_ciL = squeeze(J_ciLower(ses,para,:,:));
CO_ST_ciL = squeeze(J_ciLower(ses,para+4,:,:));
INT_NS_ciL = squeeze(J_ciLower(ses,para+8,:,:));
INT_ST_ciL = squeeze(J_ciLower(ses,para+12,:,:));

CO_NS_ciU = squeeze(J_ciUpper(ses,para,:,:));
CO_ST_ciU = squeeze(J_ciUpper(ses,para+4,:,:));
INT_NS_ciU = squeeze(J_ciUpper(ses,para+8,:,:));
INT_ST_ciU = squeeze(J_ciUpper(ses,para+12,:,:)); 

% ICMS
CO_icms_1 = squeeze(J_ICMSp(ses,1,para,6));
CO_icms_2 = squeeze(J_ICMSp(ses,2,para,6));
INT_icms_1 = squeeze(J_ICMSp(ses,1,para,8));
INT_icms_2 = squeeze(J_ICMSp(ses,2,para,8));

% RT MT
CO_NS_RT = squeeze(J_RTMT_mean(ses,1,para,5));
CO_NS_MT = squeeze(J_RTMT_mean(ses,2,para,5));
CO_ST_RT = squeeze(J_RTMT_mean(ses,1,para,6));
CO_ST_MT = squeeze(J_RTMT_mean(ses,2,para,6));

INT_NS_RT = squeeze(J_RTMT_mean(ses,1,para,7));
INT_NS_MT = squeeze(J_RTMT_mean(ses,2,para,7));
INT_ST_RT = squeeze(J_RTMT_mean(ses,1,para,8));
INT_ST_MT = squeeze(J_RTMT_mean(ses,2,para,8));

% significance
p_value_CO=[];p_value_INT=[];
for t = 1:size(CO_NS_PC,1)
    for i = 1:3
        [~, p_value_CO(t,i)] = ttest_from_ci(CO_NS_PC(t,i), [CO_NS_ciL(t,i), CO_NS_ciU(t,i)], 100, CO_ST_PC(t,i), [CO_ST_ciL(t,i), CO_ST_ciU(t,i)], 100, 0.05);
        [~, p_value_INT(t,i)] = ttest_from_ci(INT_NS_PC(t,i), [INT_NS_ciL(t,i), INT_NS_ciU(t,i)], 100, INT_ST_PC(t,i), [INT_ST_ciL(t,i), INT_ST_ciU(t,i)], 100, 0.05);
    end
end

p_value_CO(p_value_CO>=0.05)=0;
p_value_CO(p_value_CO<0.05&p_value_CO>0)=1;

p_value_INT(p_value_INT>=0.05)=0;
p_value_INT(p_value_INT<0.05&p_value_INT>0)=1;

t = t_start:bin:t_end;
plot_PC_4_2(J_C_m,J_ciLower,J_ciUpper,J_ICMSp,J_RTMT_mean,ses,t)

%% Figure 3de 3D 

load('G_0914_PC_NSST_PreGO_3D.mat')

T = -2000:2:2000; % Time vector
t_start = -1000;
t0 = find(T == 0);
t1 = find(T == 0) + t_start / 2;

para = 5; ym = ys_m_CO; % long CO
% para = 7; ym = ys_m_INT; % long INT

colo = [0.4 0.80392 0.66667; 0.4902 0.14902 0.80392; 1 0.44706 0.33725]; % Colors for different conditions
figure
set(gcf, 'Position', [30, 50, 900, 900]); % Set figure size
hold on

% Quadrant 1
q = 1;

% NS condition
ys_m_NS = ym{q, 1};
t = size(ys_m_NS, 1); % Number of time steps
GOtime = t0 - t1 + 1; % Time index for Go

% Scatter plot for NS condition
for i = 1:2:GOtime - 200
    scatter3(ys_m_NS(i, 1), ys_m_NS(i, 2), ys_m_NS(i, 3), (i / 80)^2, colo(1, :), 'filled', 'MarkerFaceAlpha', 0.3);
end
for i = GOtime - 200 + 1:2:t
    scatter3(ys_m_NS(i, 1), ys_m_NS(i, 2), ys_m_NS(i, 3), (i / 80)^2, colo(1, :), 'filled', 'MarkerFaceAlpha', 0.8);
end

% Key events (TO, GO, MO, Tch)
scatter3(ys_m_NS(1 + 100 / 2, 1), ys_m_NS(1 + 100 / 2, 2), ys_m_NS(1 + 100 / 2, 3), 200, 'o', 'c', 'filled', 'LineWidth', 10); % TO
scatter3(ys_m_NS(t0 - t1 + 1, 1), ys_m_NS(t0 - t1 + 1, 2), ys_m_NS(t0 - t1 + 1, 3), 200, 'o', 'filled', 'MarkerFaceColor', [0.23529, 0.70196, 0.44314], 'LineWidth', 10); % GO
scatter3(ys_m_NS(t0 - t1 + 1 + RT_mean(q, para), 1), ys_m_NS(t0 - t1 + 1 + RT_mean(q, para), 2), ys_m_NS(t0 - t1 + 1 + RT_mean(q, para), 3), 200, 'o', 'r', 'filled', 'LineWidth', 10); % MO
scatter3(ys_m_NS(t0 - t1 + 1 + RT_mean(q, para) + MT_mean(q, para), 1), ys_m_NS(t0 - t1 + 1 + RT_mean(q, para) + MT_mean(q, para), 2), ys_m_NS(t0 - t1 + 1 + RT_mean(q, para) + MT_mean(q, para), 3), 300, 'o', 'k', 'filled', 'LineWidth', 10); % Tch

% ST condition
icms_1 = ICMSp_1(q, para + 1);
icms_2 = ICMSp_2(q, para + 1);

ys_m_ST = ym{q, 2};
t = size(ys_m_ST, 1); % Number of time steps
GOtime = t0 - t1 + 1; % Time index for Go

% Scatter plot for ST condition
for i = 1:2:GOtime - 200
    scatter3(ys_m_ST(i, 1), ys_m_ST(i, 2), ys_m_ST(i, 3), (i / 80)^2, colo(2, :), 'filled', 'MarkerFaceAlpha', 0.3);
end
for i = GOtime - 200 + 1:2:t
    if i >= icms_1 && i <= icms_2
        continue % Skip ICMS period
    end
    scatter3(ys_m_ST(i, 1), ys_m_ST(i, 2), ys_m_ST(i, 3), (i / 80)^2, colo(2, :), 'filled', 'MarkerFaceAlpha', 0.8);
end

% Key events (TO, MO, Tch)
scatter3(ys_m_ST(1 + 100 / 2, 1), ys_m_ST(1 + 100 / 2, 2), ys_m_ST(1 + 100 / 2, 3), 200, 'o', 'c', 'filled', 'LineWidth', 10); % TO
scatter3(ys_m_ST(t0 - t1 + 1 + RT_mean(q, para + 1), 1), ys_m_ST(t0 - t1 + 1 + RT_mean(q, para + 1), 2), ys_m_ST(t0 - t1 + 1 + RT_mean(q, para + 1), 3), 200, 'o', 'r', 'filled', 'LineWidth', 10); % MO
scatter3(ys_m_ST(t0 - t1 + 1 + RT_mean(q, para + 1) + MT_mean(q, para + 1), 1), ys_m_ST(t0 - t1 + 1 + RT_mean(q, para + 1) + MT_mean(q, para + 1), 2), ys_m_ST(t0 - t1 + 1 + RT_mean(q, para + 1) + MT_mean(q, para + 1), 3), 300, 'o', 'k', 'filled', 'LineWidth', 10); % Tch

% ICMS event markers
scatter3(ys_m_ST(icms_1 - 1, 1), ys_m_ST(icms_1 - 1, 2), ys_m_ST(icms_1 - 1, 3), 200, '^', 'filled', 'MarkerFaceColor', colo(3, :), 'LineWidth', 10); % ICMS start
scatter3(ys_m_ST(icms_2 + 1, 1), ys_m_ST(icms_2 + 1, 2), ys_m_ST(icms_2 + 1, 3), 200, '^', 'filled', 'MarkerFaceColor', colo(3, :), 'LineWidth', 10); % ICMS end

% NS-ST connected lines visualization
for i = GOtime - 200 + 1:2:t - 1
    if i >= icms_1 && i <= icms_2
        continue % Skip ICMS period
    end
    
    v1 = ys_m_NS(i + 1, :) - ys_m_NS(i, :); % Velocity in NS
    v2 = ys_m_ST(i, :) - ys_m_NS(i, :); % Velocity in ST
    angle_deg = VectorAngle(v1, v2); % Angle between velocities
    Colo = (angle_deg < 90) * colo(2, :) + (angle_deg >= 90) * colo(1, :); % Select color based on angle
    plot3([ys_m_NS(i, 1), ys_m_ST(i, 1)], [ys_m_NS(i, 2), ys_m_ST(i, 2)], [ys_m_NS(i, 3), ys_m_ST(i, 3)], '-', 'Color', Colo, 'LineWidth', 0.5);
end

% Axis labels
grid off
xlabel(['PC1 v=' num2str(round(explained(1))) '%']);
ylabel(['PC2 v=' num2str(round(explained(2))) '%']);
zlabel(['PC3 v=' num2str(round(explained(3))) '%']);

%% Figure 3de 3PCs 

load('G_0914_PC_NSST_PreGO_3PCs.mat')

% 1 5 9  13
% 2 6 10 14
% 3 7 11 15
% 4 8 12 16

targetsp = -120;
t_start = -1000;
t_end = 500;
bin = 20;

ses = 1;
para = 1;

% PC
CO_NS_PC = squeeze(J_C_m(ses,para,:,:));
CO_ST_PC = squeeze(J_C_m(ses,para+4,:,:));
INT_NS_PC = squeeze(J_C_m(ses,para+8,:,:));
INT_ST_PC = squeeze(J_C_m(ses,para+12,:,:));

% CI
CO_NS_ciL = squeeze(J_ciLower(ses,para,:,:));
CO_ST_ciL = squeeze(J_ciLower(ses,para+4,:,:));
INT_NS_ciL = squeeze(J_ciLower(ses,para+8,:,:));
INT_ST_ciL = squeeze(J_ciLower(ses,para+12,:,:));

CO_NS_ciU = squeeze(J_ciUpper(ses,para,:,:));
CO_ST_ciU = squeeze(J_ciUpper(ses,para+4,:,:));
INT_NS_ciU = squeeze(J_ciUpper(ses,para+8,:,:));
INT_ST_ciU = squeeze(J_ciUpper(ses,para+12,:,:)); 

% ICMS
CO_icms_1 = squeeze(J_ICMSp(ses,1,para,6));
CO_icms_2 = squeeze(J_ICMSp(ses,2,para,6));
INT_icms_1 = squeeze(J_ICMSp(ses,1,para,8));
INT_icms_2 = squeeze(J_ICMSp(ses,2,para,8));

% RT MT
CO_NS_RT = squeeze(J_RTMT_mean(ses,1,para,5));
CO_NS_MT = squeeze(J_RTMT_mean(ses,2,para,5));
CO_ST_RT = squeeze(J_RTMT_mean(ses,1,para,6));
CO_ST_MT = squeeze(J_RTMT_mean(ses,2,para,6));

INT_NS_RT = squeeze(J_RTMT_mean(ses,1,para,7));
INT_NS_MT = squeeze(J_RTMT_mean(ses,2,para,7));
INT_ST_RT = squeeze(J_RTMT_mean(ses,1,para,8));
INT_ST_MT = squeeze(J_RTMT_mean(ses,2,para,8));

% significance
p_value_CO=[];p_value_INT=[];
for t = 1:size(CO_NS_PC,1)
    for i = 1:3
        [~, p_value_CO(t,i)] = ttest_from_ci(CO_NS_PC(t,i), [CO_NS_ciL(t,i), CO_NS_ciU(t,i)], 100, CO_ST_PC(t,i), [CO_ST_ciL(t,i), CO_ST_ciU(t,i)], 100, 0.05);
        [~, p_value_INT(t,i)] = ttest_from_ci(INT_NS_PC(t,i), [INT_NS_ciL(t,i), INT_NS_ciU(t,i)], 100, INT_ST_PC(t,i), [INT_ST_ciL(t,i), INT_ST_ciU(t,i)], 100, 0.05);
    end
end

p_value_CO(p_value_CO>=0.05)=0;
p_value_CO(p_value_CO<0.05&p_value_CO>0)=1;

p_value_INT(p_value_INT>=0.05)=0;
p_value_INT(p_value_INT<0.05&p_value_INT>0)=1;

% plot
t_start = -1000;t_end = 500;
t = t_start:bin:t_end;
plot_PC(CO_NS_PC,CO_ST_PC,CO_NS_ciL,CO_ST_ciL,CO_NS_ciU,CO_ST_ciU,CO_icms_1,CO_icms_2,CO_NS_RT,CO_NS_MT,...
    CO_ST_RT,CO_ST_MT,p_value_CO,t,bin)

plot_PC(INT_NS_PC,INT_ST_PC,INT_NS_ciL,INT_ST_ciL,INT_NS_ciU,INT_ST_ciU,INT_icms_1,INT_icms_2,INT_NS_RT,INT_NS_MT,...
    INT_ST_RT,INT_ST_MT,p_value_INT,t,bin)

% plot_PC_4(J_C_m,J_ciLower,J_ciUpper,J_ICMSp,J_RTMT_mean,ses,t)

%% Figure 3fg CCF
% G
load('G_-120_GO_CCF.mat')
load('G_-240_GO_CCF.mat')

% L
load('L_-120_GO_CCF.mat')
load('L_180_GO_CCF.mat')

for i = 1:size(J_C_m,1)
        
        score = squeeze(J_C_m(i,:,end-200:end,1));
        
        for k = 1:4
            x1 = score(k,:);
            y1 = score(k+4,:);
            x2 = score(k+8,:);
            y2 = score(k+4+8,:);
            
            [lagTime_1(i,k),ccf_1(i,k,:),lags_1(i,k,:)] = CCF(x1,y1,0.002,1);
            [lagTime_2(i,k),ccf_2(i,k,:),lags_2(i,k,:)] = CCF(x2,y2,0.002,1);
        end

end

temp = squeeze(lagTime_1(1,1,:));

lagTime(:,1) = reshape(lagTime_1,size(lagTime_1,1)*size(lagTime_1,2),1);
lagTime(:,2) = reshape(lagTime_2,size(lagTime_2,1)*size(lagTime_2,2),1);

%% Figure 3h neural distance

% Time vector and load data
t_start = -600;
t_end = 300;
t = t_start:2:t_end;

% Load data
% load('_neural_distance.mat')
% icms_2 = find(t == 0) + 8; % G
% icms_1 = find(t == -410); % G

load('L_neural_distance.mat')
icms_2 = find(t == 0)+2; % L
icms_1 = find(t==-420); % L

JJJ_EDsingle_s = [];
% Loop through the data and extract relevant time windows for each condition
for i = 1:size(JJ_EDsingle_s, 1)
    for j = 1:size(JJ_EDsingle_s, 3)
        try
            % Extract ICMS times
            icms_co = J_ICMSp(i, 2, j, 6);
            icms_int = J_ICMSp(i, 2, j, 8);

            % Extract data for each condition, using ICMS times and t_start/t_end
            JJJ_EDsingle_s(i, 1, j, :) = JJ_EDsingle_s(i, 1, j, icms_co + t_start / 2 : icms_co + t_end / 2);
            JJJ_EDsingle_s(i, 2, j, :) = JJ_EDsingle_s(i, 2, j, icms_int + t_start / 2 : icms_int + t_end / 2);
        catch
            % Handle any errors (e.g., if indices are out of bounds)
            JJJ_EDsingle_s(i, 1, j, :) = nan(1, length(icms_co + t_start / 2 : icms_co + t_end / 2));
            JJJ_EDsingle_s(i, 2, j, :) = nan(1, length(icms_co + t_start / 2 : icms_co + t_end / 2));
            disp(['Error at i=' num2str(i) ', j=' num2str(j)]);
        end
    end
end

% Convert data from mV to V for further analysis
A = JJJ_EDsingle_s / 1000;

% Reshape the data for statistical analysis
X = permute(A, [1, 3, 2, 4]);
XX = reshape(X, [size(X, 1) * size(X, 2), 2, size(X, 4)]);

h_m = [];p_m = [];
% Perform the Wilcoxon signed-rank test for each time point
for i = 1:size(XX, 3)
    [p_m(i, 1), h_m(i, 1)] = signrank(squeeze(XX(:, 1, i)), squeeze(XX(:, 2, i)));
end

% Calculate mean values across all trials
C = squeeze(mean(A, [1, 3], 'omitnan'));

% Bootstrap Confidence Intervals (CI)
statFunc = @(x) mean(x, 'omitnan');
ci_1 = bootci(1000, {statFunc, squeeze(XX(:, 1, :))}, 'alpha', 0.05, 'type', 'percentile');
ci_2 = bootci(1000, {statFunc, squeeze(XX(:, 2, :))}, 'alpha', 0.05, 'type', 'percentile');

figure
set(gcf, 'Unit', 'centimeters', 'Position', [5 10 4.5 4]); % Set figure size
colo = [0.06275, 0.30588, 0.5451; 0.80392, 0.33333, 0.33333; 1, 0.44706, 0.33725]; % Color scheme
hold on

% Plot the mean and confidence intervals for each condition
plot(t(1:icms_1), C(1, 1:icms_1), 'Color', colo(1, :), 'LineWidth', 1.5);
g = fill([t(1:icms_1), t(icms_1:-1:1)], [ci_1(1, 1:icms_1), ci_1(2, icms_1:-1:1)], 'b');
set(g, 'FaceColor', colo(1, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot the second part of the curve with a shift
shift = 300;
plot(t(icms_2:end) - shift, C(1, icms_2:end), 'Color', colo(1, :), 'LineWidth', 1.5);
g = fill([t(icms_2:end) - shift, t(end:-1:icms_2) - shift], [ci_1(1, icms_2:end), ci_1(2, end:-1:icms_2)], 'b');
set(g, 'FaceColor', colo(1, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Repeat for the second condition (C2)
plot(t(1:icms_1), C(2, 1:icms_1), 'Color', colo(2, :), 'LineWidth', 1.5);
g = fill([t(1:icms_1), t(icms_1:-1:1)], [ci_2(1, 1:icms_1), ci_2(2, icms_1:-1:1)], 'b');
set(g, 'FaceColor', colo(2, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Plot the second part of the second curve with a shift
plot(t(icms_2:end) - shift, C(2, icms_2:end), 'Color', colo(2, :), 'LineWidth', 1.5);
g = fill([t(icms_2:end) - shift, t(end:-1:icms_2) - shift], [ci_2(1, icms_2:end), ci_2(2, end:-1:icms_2)], 'b');
set(g, 'FaceColor', colo(2, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');

% Mark significant time points based on the p-value from the Wilcoxon test
y0 = ylim;
xlim([t_start t_end+50-shift])
x0 = xlim;
for k = 1:length(h_m)
    if h_m(k) == 1 && (k < icms_1 || k > icms_2) % Significant difference
        plot([t(k) - 5 - shift, t(k) + 5 - shift], [y0(2) - 0.001, y0(2) - 0.001], 'k', 'LineWidth', 3);
    end
end

% Highlight the ICMS period with a shaded patch
x1 = t(icms_1); x2 = t(icms_2);
y1 = y0(1); y2 = y0(2);
patch([x1, x1, x2 - shift, x2 - shift], [y1, y2, y2, y1], 'y', 'FaceColor', colo(3, :), 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% Add vertical and horizontal lines for reference
plot([t(icms_2) - shift, t(icms_2) - shift], [y0(1), y0(2)], ':', 'Color', colo(3, :), 'LineWidth', 1);
plot([t_start, t_end - shift], [0, 0], 'k:', 'LineWidth', 1);
plot([x0(2), x0(2)], [0, 0.01], 'k', 'LineWidth', 1.5);

% Add vertical line at the start time
plot([t_start + 10, t_start + 10 + 100], [y0(1), y0(1)], 'k', 'LineWidth', 1.5);

axis off

%% Figure 3i neural speed

% Time vector and window size
t_start = -600;
t_end = 300;
windowSize = 5; % Define window size for median filter
t = t_start:2:t_end;

% Load data
load('G_neural_speed.mat')
icms_2 = find(t == 0)+8; % G
icms_1 = find(t==-410); % G

% load('L_neural_speed.mat')
% icms_2 = find(t == 0)+2;
% icms_1 = find(t==-420);

vel = [];
for i = 1:length(J_xyz_v)
    for j = 1:4
        % Compute median velocity and convert to m/s (divide by 1000)
        vel(i,j,:,:) = cell2mat((cellfun(@(x) ((median(x,1)/1000)),J_xyz_v{i,1}(j,:),'UniformOutput',0))');
    end
end

% Process data with median filtering
JJ_xyz = [];
for i = 1:size(J_xyz_v,1)
    for j = 1:4
        try
            % Extract data around ICMS times
            icms_co = J_ICMSp(i,2,j,6);
            data = vel(i,j,1:2,icms_co+t_start/2:icms_co+t_end/2);
            for k = 1:2
                JJ_xyz(i,j,k,:) = medfilt1(data(:,:,k,:), windowSize,'omitnan');
            end
            
            % Process second part (internal stimulation)
            icms_int = J_ICMSp(i,2,j,8);
            data = vel(i,j,3:4,icms_int+t_start/2:icms_int+t_end/2);
            for k = 1:2
                JJ_xyz(i,j,k+2,:) = medfilt1(data(:,:,k,:), windowSize,'omitnan');
            end
        catch
            JJ_xyz(i,j,1:2,:) = nan(2,length(icms_co+t_start/2:icms_co+t_end/2));
            JJ_xyz(i,j,3:4,:) = nan(2,length(icms_co+t_start/2:icms_co+t_end/2));
        end
    end
end

% Reshape data for statistical testing
A = JJ_xyz;
XX = reshape(A, [size(A,1)*size(A,2), size(A,3), size(A,4)]);

% Perform Wilcoxon signed-rank tests
h_m_ns =[];p_m_ns=[];h_m_st =[];p_m_st=[];p_m_co=[];h_m_co=[];p_m_int=[];h_m_int=[];
for i = 1:size(XX,3)
    [p_m_ns(i,1),h_m_ns(i,1)] = signrank(squeeze(XX(:,1,i)),squeeze(XX(:,3,i)));
    [p_m_st(i,1),h_m_st(i,1)] = signrank(squeeze(XX(:,2,i)),squeeze(XX(:,4,i)));
    [p_m_co(i,1),h_m_co(i,1)] = signrank(squeeze(XX(:,1,i)),squeeze(XX(:,2,i)));
    [p_m_int(i,1),h_m_int(i,1)] = signrank(squeeze(XX(:,3,i)),squeeze(XX(:,4,i)));
end

% Calculate the mean across trials
C = squeeze(mean(A,[1,2],'omitnan'));

% Bootstrap Confidence Intervals (CI)
statFunc = @(x) mean(x,'omitnan');
ci = [];
for i = 1:4
    An = squeeze(XX(:,i,:));
    ci(i,:,:) = bootci(1000, {statFunc, An}, 'alpha', 0.05, 'type', 'percentile');
end

figure
set(gcf,'unit','centimeters','position',[5 10 4.5 4]);
colo=[0.06275	0.30588	0.5451;0.06275	0.30588	0.5451;0.80392	0.33333	0.33333;0.80392	0.33333	0.33333;1	0.44706	0.33725];
L = {':','-',':','-'};
hold on

for i = 1:4
    % Plot the first part (before ICMS)
    plot(t(1:icms_1),C(i,1:icms_1),'color',colo(i,:),'linestyle',L{i},'LineWidth',1.5)
    g=fill([t(1:icms_1),t(icms_1:-1:1)],[squeeze(ci(i,1,1:icms_1))',squeeze(ci(i,2,icms_1:-1:1))'],'b');
    set(g,'FaceColor',colo(i,:),'FaceAlpha',0.1,'EdgeColor','none');
    
    % Plot the second part (after ICMS with shift)
    shift = 300;
    plot(t(icms_2:end)-shift,C(i,icms_2:end),'color',colo(i,:),'linestyle',L{i},'LineWidth',1.5)
    g=fill([t(icms_2:end)-shift,t(end:-1:icms_2)-shift],[squeeze(ci(i,1,icms_2:end))',squeeze(ci(i,2,end:-1:icms_2))'],'b');
    set(g,'FaceColor',colo(i,:),'FaceAlpha',0.1,'EdgeColor','none');
end

xlim([t_start t_end+50-shift])
x0 = xlim;y0 = ylim;
dy = (y0(2)-y0(1))*0.04;

% Plot significance markers for each comparison
for k = 1:length(h_m_ns)
    if h_m_ns(k)==1 && (k < icms_1 || k > icms_2)
        plot([t(k)-5-shift t(k)+5-shift],[y0(2)-2*dy y0(2)-2*dy],'-','color',[0.5 0.5 0.5],'LineWidth',3)
    end
    
    if h_m_st(k)==1 && (k < icms_1 || k > icms_2)
        plot([t(k)-5-shift t(k)+5-shift],[y0(2)-dy y0(2)-dy],'k-','LineWidth',3)
    end
    
    if h_m_co(k)==1 && (k < icms_1 || k > icms_2)
        plot([t(k)-5-shift t(k)+5-shift],[y0(2)-3*dy y0(2)-3*dy],'-','LineWidth',3,'color',[0.06275	0.30588	0.5451])
    end
    
    if h_m_int(k)==1 && (k < icms_1 || k > icms_2)
        plot([t(k)-5-shift t(k)+5-shift],[y0(2)-4*dy y0(2)-4*dy],'-','LineWidth',3,'color',[0.80392	0.33333	0.33333])
    end
end

% Highlight ICMS period with a patch
x1 = t(icms_1); x2 = t(icms_2);   % 方块的x范围
y1 = y0(1); y2 = y0(2);   % 方块的y范围
patch([x1 x1 x2-shift x2-shift], [y1 y2 y2 y1],'y','FaceColor',colo(5,:),'EdgeColor','none','FaceAlpha',0.2)
plot([t(icms_2)-shift t(icms_2)-shift],[y0(1) y0(2)],':','color',colo(5,:),'LineWidth',1)

% Add reference lines
plot([t_start+10,t_start+10+100],[y0(1),y0(1)],'k','LineWidth',1.5)
plot([x0(2),x0(2)],[y0(1),y0(1)+0.2],'k','LineWidth',1.5)

axis off

