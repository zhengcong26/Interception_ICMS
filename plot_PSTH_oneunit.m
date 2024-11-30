

function plot_PSTH_oneunit(X_check,num,colo)

endColor = colo;
c1 = gradient_ramp(endColor,size(X_check,2),3);

figure
set(gcf,'position',[30,100,350,200]);
hold on
for i = 1:size(X_check,2)
    A = squeeze(X_check(:,i,num,:));
    B = mean(A,1);
    B = [repmat(B(:,1),1,100), B];
    t=linspace(-100,size(B,2)-100,size(B,2));
    plot(t,B','color',c1(i,:),'LineWidth',2)
    
%     statFunc = @(x) mean(x,'omitNaN');
%     ci = bootci(1000, {statFunc, A}, 'alpha', 0.05, 'type', 'percentile');
%     ci = [repmat(ci(:,1),1,100), ci];
%     
%     g=fill([t,t(end:-1:1)],[ci(1,:),ci(2,end:-1:1)],'b');
%     set(g,'FaceColor',c1(i,:),'FaceAlpha',0.2,'EdgeColor','none');
end

box off
ylabel('firing rates (Hz)');
xlim([-150 inf]);
ylim([17 25]);
y1 =  ylim;
xticks([0 100 200 300 400 500 600 700 800 900 1000])
set(gca,'TickDir','out');

h = gca;
h.FontSize = 15;
h.LineWidth = 1;
h.TickLength=[0.02,0.025];

