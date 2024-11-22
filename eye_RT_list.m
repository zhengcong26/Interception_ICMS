
function eye_RT_CO = eye_RT_list(eye_CO)

% correct degree
dtheta = atand(10.2/30)/10.2;

eye_RT_CO = zeros(size(eye_CO,1),6);
for n = 1:size(eye_CO,1)
    eye = eye_CO(n,:);
    t_TO = eye{1,5}(eye{1,1}(:,1) == eye{1,2}(3,1) & eye{1,1}(:,2) == eye{1,2}(3,2),1);
    v = eye{1,5}(t_TO:end,4); % TO-Tch
    v_f = low_pass(v)*10.2/16*dtheta;
    
    t = 1:1:length(v_f);
    [pks,locs] = findpeaks(v_f,t);
    [Y_1,I_1] = max(pks);
    pk_1 = locs(I_1);
    pks(pks == Y_1) = -Inf;
    [Y_2,I_2] = max(pks);
    pk_2 = locs(I_2);
    
    t_delay = eye{1,8}(1,4);
    eye_RT_CO(n,1) = t_delay; % 1 delay
    eye_RT_CO(n,2) = eye{1,8}(1,13); % 2 RT
    try
        eye_RT_CO(n,3) = eye{1,8}(1,48); % 3 MT
    catch
        eye{1,8}(1,48) = eye{1,8}(1,9) - eye{1,8}(1,11);
        eye_RT_CO(n,3) = eye{1,8}(1,48); % 3 MT
    end
    eye_RT_CO(n,5) = pk_1-t_delay; % #1 peak
    eye_RT_CO(n,6) = pk_2-t_delay; % #2 peak
    
    eye_RT_CO(n,7) = Y_1;
    eye_RT_CO(n,8) = Y_2;
end

for n = 1:size(eye_RT_CO,1)
    
    if eye_RT_CO(n,5) > 0 && eye_RT_CO(n,6) < 0
        eye_RT_CO(n,10) = eye_RT_CO(n,5);
        if eye_RT_CO(n,6) > -200
            eye_RT_CO(n,11) = eye_RT_CO(n,6);
        end
        
    elseif eye_RT_CO(n,5) < 0 && eye_RT_CO(n,6) > 0
        eye_RT_CO(n,10) = eye_RT_CO(n,6);
        
    elseif eye_RT_CO(n,5) < 0 && eye_RT_CO(n,6) < 0
        if eye_RT_CO(n,5) < eye_RT_CO(n,6)
            eye_RT_CO(n,10) = eye_RT_CO(n,6);
        else
            eye_RT_CO(n,10) = eye_RT_CO(n,5);
        end
        
    elseif eye_RT_CO(n,5) > 0 && eye_RT_CO(n,6) > 0
        if eye_RT_CO(n,5) < eye_RT_CO(n,6)
            eye_RT_CO(n,10) = eye_RT_CO(n,5);
        else
            eye_RT_CO(n,10) = eye_RT_CO(n,6);
        end
    end
    
    if eye_RT_CO(n,10) > eye_RT_CO(n,2)+eye_RT_CO(n,3)
        eye_RT_CO(n,10) = NaN;
    end

end

for n = 1:size(eye_RT_CO,1)
    
    if eye_RT_CO(n,11) == 0
        eye_RT_CO(n,12) = eye_RT_CO(n,10);
    else
        eye_RT_CO(n,12) = eye_RT_CO(n,11);
    end
           
end
