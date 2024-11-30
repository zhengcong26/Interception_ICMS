
function plotMotorCost(data, color, linestyle, maxSamples)
    % Transpose data for easier manipulation
    y = data';
    t = 1:size(y, 1);

    % Calculate mean and standard error
    y_mean = mean(y, 2, 'omitnan');
    y_std = std(y, 0, 2, 'omitnan') / sqrt(size(y, 2) - 1);

    % Plot individual trajectories
    for i = 1:min(maxSamples, size(y, 2))
        plot(t, y(:, i), 'color', [color, 0.2], 'LineStyle', linestyle, 'LineWidth', 1);
    end

    % Plot mean trajectory
    plot(t, y_mean, 'color', color, 'LineStyle', linestyle, 'LineWidth', 2);

    % Plot shaded error area if standard deviation exists
    if any(y_std ~= 0)
        fill([t, t(end:-1:1)], [y_mean - y_std; flipud(y_mean + y_std)]', ...
            color, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
end