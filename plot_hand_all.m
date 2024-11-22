
function plot_hand_all(hand, y_check, colo)

ntrial = size(hand,2);

endColor = colo;
c1 = gradient_ramp(endColor,ntrial,5);

figure
set(gcf,'position',[30,300,300,300]);
hold on
for trial = 1:ntrial
    plot(hand{trial}(1,:),hand{trial}(2,:),'k:','linewidth',1)
    
    for i = 1:size(y_check,1)
        plot(squeeze(y_check(i,trial,1,:)), squeeze(y_check(i,trial,2,:)),'color',c1(trial,:),'linewidth',2)
    end

    box off
    axis equal
    axis off
    xlim([-20 20]);
    ylim([-20 20]);
end
    
end