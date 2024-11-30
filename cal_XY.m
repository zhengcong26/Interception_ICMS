
function [X, Y, r2, rmse] = cal_XY(pr, Wsoc, C, hand, Xbefore, Xafter_true, r_dir, state_start, figon)
    % Simulate and compare system states and outputs
    % Inputs:
    %   pr - Parameters structure
    %   Wsoc - Weight matrix for social dynamics
    %   C - Output mapping matrix
    %   hand - Observed hand trajectories
    %   Xbefore - System states before movement onset (MO)
    %   Xafter - System states after MO
    %   init, init_ahead - Indices for initial and ahead states
    % Outputs:
    %   X - Simulated and concatenated states (before + after MO)
    %   Y - Simulated output after MO
    %   r2 - R-squared metric for output comparison
    %   rmse - Root mean square error metric for output comparison

    % cal true X
    X_true = [Xbefore, Xafter_true{state_start}];

    % Simulate states after MO using last state from Xbefore_INT
    Xafter = cal_example_Xafter(pr, Wsoc, Xbefore(:, end)); % Simulated after-MO states

    % cal X of Xbefore and Xafter
    X = [Xbefore, Xafter];
    
    % cal y after MO
    Y = C * max(0, Xafter);
    
    % Calculate R-squared and RMSE metrics for the output
    [r2, rmse] = cal_r(hand{state_start}, Y);

    %% Create figure for visualization
    if figon == 1
        figure
        set(gcf, 'position', [100, 150, 200, 650]);
        
        % Plot state comparison
        subplot(3, 1, 1);
        hold on
        X_true = [repmat(X_true(:,1),1,100), X_true]; % Extend the drawing timeline, before MO, cong
        X_p = [repmat(X(:,1),1,100), X];
        plot(X_true(1:3,:)','k','LineWidth',1) % Draw the activities of 3 neurons
        plot(X_p(1:3,:)','b','LineWidth',1)
        ylim([10 40]);
        box on; axis off
        title('State Comparison');
        
        % Plot output comparison
        subplot(3, 1, 2);
        hold on
        y_t = hand{r_dir};
        plot(y_t','k')
        plot(Y','r','LineWidth',1)
        box on
        xticks([]);
        title('Output Comparison');
        
        % Visualize observed hand trajectory and simulated outputs
        subplot(3, 1, 3);
        colo = [0.06275, 0.30588, 0.5451; 0.80392, 0.33333, 0.33333];
        plot_example_y(hand, Y, r_dir, state_start, colo(2, :), 0);
        title('Hand Trajectory Comparison');
    end
end
