
function plot_session_structure(RAND, velocity_list, colo)
    
%     n = 1;
    for k = 1:length(RAND)
        for i = 1:length(velocity_list)
            % Check if the third element in the first row of velocity_list matches RAND(k)
            if velocity_list{i,1}(1,3) == RAND(k)
                
                % Initialize x and y offsets based on the initial values of the velocity data
                x0 = -velocity_list{i,3}(1,1);
                y0 = velocity_list{i,3}(1,2);
                
                % Calculate relative trajectory positions for x and y components
                tjC(1:10,1) = -velocity_list{i,3}(1:10,1) - x0;
                tjC(11:20,1) = velocity_list{i,3}(1:10,2) - y0;
                
                % Skip this iteration if any values in tjC are less than -100
                if any(tjC < -100)
                    continue
                end
                
                % Smooth the x and y trajectory data using a Gaussian window of size 5
                x = smoothdata(tjC(1:10,1), 'gaussian', 5);
                y = smoothdata(tjC(11:20,1), 'gaussian', 5);
                
                plot(x, y, 'color', colo, 'linewidth', 1.5);
                hold on
                
%                 temp(1:10,n) = x;
%                 temp(11,n) = NaN;
%                 temp(12:21,n) = y;
%                 n = n+1;
                
                % Mark the end of the trajectory with a circle marker
                scatter(x(end), y(end), 40, 'o', 'MarkerEdgeColor', colo, 'LineWidth', 1.5);
            end
        end
    end
    
end

