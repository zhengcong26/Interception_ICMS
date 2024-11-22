
function plot_psth_without_raster(spike_time_TO,t_TO,spike_time_GO,t_GO,spike_time_MO,t_MO,cor_tri_sp,hand_Sur_loc)

speed_con=sort(unique(cor_tri_sp),'descend');
d_step = 50;

[spike_time_MO_zero,p_dir_1_MO,p_dir_2_MO] = process_spike_data(spike_time_MO,cor_tri_sp,hand_Sur_loc,speed_con);
[spike_time_GO_zero,p_dir_1_GO,p_dir_2_GO] = process_spike_data(spike_time_GO,cor_tri_sp,hand_Sur_loc,speed_con);
[spike_time_TO_zero,p_dir_1_TO,p_dir_2_TO] = process_spike_data(spike_time_TO,cor_tri_sp,hand_Sur_loc,speed_con);

%%
colo=[0.06275	0.30588	0.5451;0.80392	0.33333	0.33333];
c1 = gradient_ramp(colo(1,:),4,1);
c2 = gradient_ramp(colo(2,:),4,1);

figure
set(gcf,'Position',[100,300,700,300])

subplot(2,1,1)
hold on
plot_psth(spike_time_TO_zero,spike_time_GO_zero,spike_time_MO_zero,t_TO,t_GO,t_MO,p_dir_1_TO,p_dir_1_GO,p_dir_1_MO,hand_Sur_loc,cor_tri_sp,speed_con,d_step,c1)

subplot(2,1,2)
hold on
plot_psth(spike_time_TO_zero,spike_time_GO_zero,spike_time_MO_zero,t_TO,t_GO,t_MO,p_dir_2_TO,p_dir_2_GO,p_dir_2_MO,hand_Sur_loc,cor_tri_sp,speed_con,d_step,c2)

h1 = subplot(2,1,1);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1), pos1(2)*0.85, pos1(3), pos1(4)]); % [left bottom width height]

h1 = subplot(2,1,2);
pos1 = get(h1, 'Position');
set(h1, 'Position', [pos1(1), pos1(2), pos1(3), pos1(4)]); % [left bottom width height]

end

function [spike_time_MO_zero,p_dir_1_MO,p_dir_2_MO] = process_spike_data(spike_time_MO,cor_tri_sp,hand_Sur_loc,speed_con)

d_step = 50;
d_step1 = 1:d_step:3951;
d_step2 = d_step1+d_step-1;

spike_time_MO_zero = spike_time_MO;
spike_time_MO_zero(isnan(spike_time_MO_zero))=0;
index_bin_MO = ceil(spike_time_MO*1000+2000);
spike_bin_MO  = zeros(size(cor_tri_sp,1),4000);
for i = 1:size(index_bin_MO,1)
    spike_bin_MO(i,index_bin_MO(i,~isnan(index_bin_MO(i,:))) ) = 1;
end

for i_t = 1:length(d_step1)
    FR_MO_50(:,i_t) = sum(spike_bin_MO(:,d_step1(i_t):d_step2(i_t)),2);
    x = FR_MO_50 (cor_tri_sp==speed_con(1),i_t);
    y = FR_MO_50 (cor_tri_sp==speed_con(2),i_t);
    p_2sp(:,i_t) = ranksum(x,y);
    
    p_dir_1_MO(:,i_t) = kruskalwallis(FR_MO_50(cor_tri_sp==speed_con(1),i_t), hand_Sur_loc(cor_tri_sp==speed_con(1),:),'off');%每个d_step内比较八方向有没有显著差异
    p_dir_2_MO(:,i_t) = kruskalwallis(FR_MO_50(cor_tri_sp==speed_con(2),i_t), hand_Sur_loc(cor_tri_sp==speed_con(2),:),'off');%每个d_step内比较八方向有没有显著差异
end

end

function plot_psth(spike_time_TO_zero,spike_time_GO_zero,spike_time_MO_zero,t_TO,t_GO,t_MO,p_dir_2_TO,p_dir_2_GO,p_dir_2_MO,hand_Sur_loc,cor_tri_sp,speed_con,d_step,c2)

for i_dir=1:4
    T =[-2,2];
    t = T(1):0.05:T(2)-0.05;
    err = 2;
    
    spk_dir_TO{i_dir,1}=spike_time_TO_zero(hand_Sur_loc==i_dir&cor_tri_sp==speed_con(2) ,:);
    [R_TO,t1,E_TO,~] = psth_NeurAn(spk_dir_TO{i_dir,1},0.05,'n',T,err,t);
    
    spk_dir_GO{i_dir,1}=spike_time_GO_zero(hand_Sur_loc==i_dir&cor_tri_sp==speed_con(2) ,:);
    [R_GO,t2,E_GO,~] = psth_NeurAn(spk_dir_GO{i_dir,1},0.05,'n',T,err,t);
    
    spk_dir_MO{i_dir,1}=spike_time_MO_zero(hand_Sur_loc==i_dir&cor_tri_sp==speed_con(2) ,:);
    [R_MO,t3,E_MO,~] = psth_NeurAn(spk_dir_MO{i_dir,1},0.05,'n',T,err,t);
    
    ta = find(t1>= t_TO(1)-0.005 & t1<= t_TO(1)+0.005,1,'last');
    tb = find(t1>= t_TO(2)-0.005 & t1<= t_TO(2)+0.005,1,'first');
    t = t1(ta:tb);
    R = R_TO(ta:tb);
    E = E_TO(ta:tb);
    plot(t,R,'color',c2(i_dir,:),'LineWidth',2);
    h1=fill([t fliplr(t)],[R-2*E fliplr(R+2*E)],c2(i_dir,:));
    alpha(0.2);
    set(h1,'LineStyle','none');
    
    ta = find(t2>= t_GO(1)-0.005 & t2<= t_GO(1)+0.005,1,'last');
    tb = find(t2>= t_GO(2)-0.005 & t2<= t_GO(2)+0.005,1,'first');
    dt_1 = t(end)-t2(ta)+0.05;
    t = t2(ta:tb)+dt_1;
    R = R_GO(ta:tb);
    E = E_GO(ta:tb);
    plot(t,R,'color',c2(i_dir,:),'LineWidth',2);
    h1=fill([t fliplr(t)],[R-2*E fliplr(R+2*E)],c2(i_dir,:));
    alpha(0.2);
    set(h1,'LineStyle','none');
    
    ta = find(t3>= t_MO(1)-0.005 & t3<= t_MO(1)+0.005,1,'last');
    tb = find(t3>= t_MO(2)-0.005 & t3<= t_MO(2)+0.005,1,'first');
    dt_2 = t(end)-t3(ta)+0.05;
    t = t3(ta:tb)+dt_2;
    R = R_MO(ta:tb);
    E = E_MO(ta:tb);
    plot(t,R,'color',c2(i_dir,:),'LineWidth',2);
    h1=fill([t fliplr(t)],[R-2*E fliplr(R+2*E)],c2(i_dir,:));
    alpha(0.2);
    set(h1,'LineStyle','none');
    
    box off
end

set(gca, 'XLim',[-0.2 1.8],'FontSize',1);

tt = -2:0.05:2-0.05;
d_dir_2 = find(p_dir_2_TO<0.05);
x_TO_2 = tt(d_dir_2)-d_step/1000/2;
x_TO_2(x_TO_2>t_TO(2))=[];
x_TO_2(x_TO_2<t_TO(1))=[];
x_TO2_2 = x_TO_2+d_step/1000;

d_dir_2 = find(p_dir_2_GO<0.05);
x_GO_2 = tt(d_dir_2)-d_step/1000/2;
x_GO_2(x_GO_2>t_GO(2))=[];
x_GO_2(x_GO_2<t_GO(1))=[];
x_GO_2 = x_GO_2+dt_1;
x_GO2_2 = x_GO_2+d_step/1000;

d_dir_2 = find(p_dir_2_MO<0.05);
x_MO_2 = tt(d_dir_2)-d_step/1000/2;
x_MO_2(x_MO_2>t_MO(2))=[];
x_MO_2(x_MO_2<t_MO(1))=[];
x_MO_2 = x_MO_2++dt_1+dt_2;
x_MO2_2 = x_MO_2+d_step/1000;

axis off

y0 = ylim;

for i_p=1:length(x_TO_2)
    plot([x_TO_2(i_p) x_TO2_2(i_p)],[y0(2) y0(2)],'color',c2(i_dir,:),'LineWidth',2);
end

for i_p=1:length(x_GO_2)
    plot([x_GO_2(i_p) x_GO2_2(i_p)],[y0(2) y0(2)],'color',c2(i_dir,:),'LineWidth',2);
end

for i_p=1:length(x_MO_2)
    plot([x_MO_2(i_p) x_MO2_2(i_p)],[y0(2) y0(2)],'color',c2(i_dir,:),'LineWidth',2);
end

x0 = xlim;
plot([x0(2)-0.3 x0(2)-0.3],[y0(1) y0(1)+20],'k','LineWidth',2);
plot([t_TO(2)+0.025 t_TO(2)+0.025],ylim,'k:','LineWidth',1);
plot([t_TO(2)+t_GO(2)-t_GO(1)+0.075 t_TO(2)+t_GO(2)-t_GO(1)+0.075],ylim,'k:','LineWidth',1);

scatter(0,y0(1),18,'k','filled')
scatter(t_TO(2)+abs(t_GO(1))+0.05,y0(1),18,'k','filled')
scatter(t_TO(2)+t_GO(2)-t_GO(1)+abs(t_MO(1))+0.1,y0(1),18,'k','filled')

end

