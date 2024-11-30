
function plot_example_y(hand, y_check, r_dir, state_start, colo, onefig)

ntrial = size(hand,2);
for i = 1:ntrial
    plot(hand{i}(1,:),hand{i}(2,:),'k--','linewidth',1.5)
    hold on
end

if onefig == 0
    scatter(hand{state_start}(1,end),hand{state_start}(2,end), 100,'go')
    scatter(hand{r_dir}(1,end),hand{r_dir}(2,end), 100,'go','filled')
    
    plot(y_check(1,:), y_check(2,:),'color','r','linewidth',3)
else
    endColor = colo;
    c1 = gradient_ramp(endColor,ntrial,5);
    for n = 1:length(y_check)
        plot(y_check{n}(1,:), y_check{n}(2,:),'color',c1(n,:),'linewidth',4)
    end
end

box off
axis equal
axis off
xlim([-20 20]);
ylim([-20 20]);

end