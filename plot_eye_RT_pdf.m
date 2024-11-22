
function plot_eye_RT_pdf(eye_RT_CO,colo)

figure
ax = axes;
x = eye_RT_CO(:,2);
h1=cdfplot_my(x);
set(h1,'LineStyle', '-', 'Color', 'k','LineWidth',1.5) % CO short
hold on
scatter(median(x,'omitnan'),0.05,30,'o','MarkerEdgeColor','k')

x = eye_RT_CO(:,10);
h1=cdfplot_my(x);
set(h1,'LineStyle', '-', 'Color', colo,'LineWidth',1.5) % CO short
scatter(median(x,'omitnan'),0.05,30,'o','MarkerEdgeColor',colo(1,:))

xlim([-200 600])
box off; grid off
ylabel('cdf','FontName','Arial');
xlabel('reaction time (ms)','FontName','Arial');
ax.TickDir='out';
set(gcf,'unit','centimeters','position',[5 10 10 4]);
set(gca,'Position',[.2 .2 .7 .7]);
offsetaxis(ax, 'y', 0.02); % Default, y-offset = 0.1 and yloc = 'l'
offsetaxis(ax, 'x', 0.02); % Default, y-offset = 0.1 and yloc = 'l'
% xticks([100,200,300,400,500]);