
function hand = cal_Hand(pr, v0)
% v0：integer
% t：row vector，(1,n_timepoints)，time after MO (movement onset)
% v：row vector，(1,n_timepoints)，conforms to the bell-shaped curve, reach 20cm within 500ms
% d：row vector，(1,n_timepoints)，absolute distance of the reach
% alpha：row vector，(1,ntrial)，the angle between the eight directions of reaching and the positive direction of the x-axis
% hand：matrix，(2*ntrial,n_timepoints)，x&y coordinates of the hand

tfinal = pr.t_move; 
n_timepoints = pr.n_move; 
alpha = pr.theta_reach;
tau_reach = pr.tau_reach;

% calculate distance: d(t)
t = linspace(0,tfinal,n_timepoints);
fun_v = @(t) v0*(t/tau_reach).^2.*exp(-0.5*(t/tau_reach).^2);
v = fun_v(t);
d = cumtrapz(t,v); % 累积梯形数值积分

% calculate hand trajectory: hand(t)
for i = 1:pr.ntrial
    hand_x = d * cos( alpha(i) ) + 0;
    hand_y = d * sin( alpha(i) ) + 0;
    hand{i} = [hand_x;hand_y];
end

end