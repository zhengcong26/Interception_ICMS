
function plot_u(u,it)

c1 = copper(4);

figure
hold on
for i = 1:length(it)
    
    U = [repmat(u(:,1),1,100), u];
    t=linspace(-100,size(U,2)-100,size(U,2));
    plot(t',U(it(i),:)','color',c1(i,:),'LineWidth',1)
    
end

set(gcf,'position',[30,100,230,80]);
set(gca,'yscale','log')
ylim([-100 6000]);
box off
set(gca,'TickDir','out');
h = gca;
h.FontSize = 8.5;
h.TickLength = [0.025,0.03];
xticks([0 100 200 300 400 500 600 700])
    
