
function plot_decoding(r2_CO_x,r2_CO_shf_x,h_CO_x,t,colo,sig)

sm = 15;
x = smoothdata(r2_CO_x',2,'gaussian',sm);
x_shf = smoothdata(r2_CO_shf_x',2,'gaussian',sm);

X = mean(x);
X_shf = mean(x_shf);
statFunc = @(x) mean(x,'omitNaN');
ci_x = bootci(1000, {statFunc, x}, 'alpha', 0.05, 'type', 'percentile');
ci_x_shf = bootci(1000, {statFunc, x_shf}, 'alpha', 0.05, 'type', 'percentile');

% t = -990:10:490;
plot(t,X,'color',colo,'LineWidth',1.2)
plot(t,X_shf,'--','color',colo,'LineWidth',0.8)
% g=fill([t,t(end:-1:1)],[X(1,:),ci_x(2,end:-1:1)],'b');
% set(g,'FaceColor',colo,'FaceAlpha',0.2,'EdgeColor','none');
g=fill([t,t(end:-1:1)],[X_shf(1,:),ci_x_shf(2,end:-1:1)],'b');
set(g,'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none');
g=fill([t,t(end:-1:1)],[ci_x(1,:),ci_x(2,end:-1:1)],'b');
set(g,'FaceColor',colo,'FaceAlpha',0.2,'EdgeColor','none');
% g=fill([t,t(end:-1:1)],[ci_x_shf(1,:),ci_x_shf(2,end:-1:1)],'b');
% set(g,'FaceColor','k','FaceAlpha',0.2,'EdgeColor','none');

diff_X = X-X_shf;

xlim([t(1) t(end)])
% ylim([-0.3 1])
ylim([-0.4 0.8])
y0 = ylim;
for i = 1:size(h_CO_x,1)
    if h_CO_x(i,1)==1 && diff_X(i)>0
        plot([t(i)-10 t(i)+10],[y0(2)-sig y0(2)-sig],'color',colo,'LineWidth',1.6)
    end
end

% scatter(-900,y0(1),5,'ko','filled')
% scatter(0,y0(1),5,'ko','filled')
% scatter(200,y0(1),5,'ko','filled')
ax = gca;
ax.TickDir='out';
ax.FontName='Arial';
ax.FontSize=8;
ax.TickLength=[0.02,0.025];
set(gca,'XColor', 'none')
% ylabel('r2')
