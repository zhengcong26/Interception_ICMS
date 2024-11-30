
%% step 1
% create Wsoc
run soc_W.m

%% step 2
pr = def_Params(); % define parameters

v0 = cal_v0(pr); % velocity parameter of the hand
Hand = cal_Hand(pr, v0); % hand trajectory (ground truth, cong)
load Wsoc21

%% step 3   --- optimize xstar and C together

% step 1
[c,xstars] = prms_v1(pr, Wsoc, Hand);

% step 2-n
[c,xstars] = prms_v2(pr, Wsoc, Hand, c, xstars); 

%% step 4
[C,Xstar] = Unpack(c,xstars,Wsoc,pr);
Xafter = cal_Xafter(pr, Wsoc, Xstar);

%% step 5  --- optimize C further

% step 1
CC1 = cal_C(pr, Xafter, Xstar, C, Hand); 

% step 2-n
CC1 = cal_C(pr, Xafter, Xstar, CC1, Hand); 