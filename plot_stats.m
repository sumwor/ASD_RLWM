function plot_stats(stats, savedatafolder, plotLabel)

setup_figprop;

nPlot = 1:length(stats.c);
figure;
sgtitle(plotLabel, 'Interpreter', 'none')
subplot(4,1,1)

value2Plot = (-1).^(stats.s + 1) .* ceil(stats.s / 2);

bar(value2Plot.*(value2Plot<0),1,'FaceColor','r','EdgeColor','none');
hold on;
bar(value2Plot.*(value2Plot>0),1,'FaceColor','b','EdgeColor','none');
hold on;
for ii = 1:length(stats.r)
    if stats.r(ii,1) > 0
        if stats.c(ii) == 1
            hold on; scatter(ii,max(value2Plot)+0.25,'.','red');  %reward right
        elseif stats.c(ii) == 0
            hold on; scatter(ii,-max(value2Plot)-0.25,'.','red');
        end
    end
end

ylim([-max(value2Plot)-0.5 max(value2Plot)+0.5]);
set(gca, 'Box', 'off');
title('Behavior')

% plot estimated QA and QB
subplot(4,1,2)
plot(nPlot, stats.qA(:,1),'r-');
hold on;
plot(nPlot, stats.qB(:,1),'b-');
plot(nPlot, stats.qB(:,2), 'b--')
plot(nPlot, stats.qA(:,2), 'r--')

ylim([-0.05 1.05]);
set(gca, 'Box', 'off');
title('Action values')

% estimated PL and fitted pL
binSize = 20;
estimated_pCorrect = nan(1, length(nPlot));
simulated_pCorrect = nan(1, length(nPlot));
for ii = 1:length(estimated_pCorrect)
    if ii<=binSize/2
        estimated_pCorrect(ii) = sum(stats.r(1:ii)>0)/ii;
        simulated_pCorrect(ii) = nanmean(stats.pCorrect(1:ii));
    elseif ii >= length(nPlot)-binSize/2
        estimated_pCorrect(ii) = sum(stats.r(ii:end)>0)/(length(nPlot)-ii+1);
         simulated_pCorrect(ii) = nanmean(stats.pCorrect(ii:end));
    else
        estimated_pCorrect(ii) = sum(stats.r(ii-binSize/2:ii+binSize/2)>0)/(binSize+1);
         simulated_pCorrect(ii) = nanmean(stats.pCorrect(ii-binSize/2:ii+binSize/2));
    end
end

%% third subplot
subplot(4,1,3)
if isfield(stats,'ck')
    plot(nPlot, stats.ck(:,1),'r-');
    hold on;
    plot(nPlot, stats.ck(:,2),'b-');


    ylim([-0.05 1.05]);
    set(gca, 'Box', 'off');
    title('Choice kernels')
    
else


    plot(nPlot, estimated_pCorrect, 'k')
    hold on;
    plot(nPlot, simulated_pCorrect, 'color',[0.7, 0.7, 0.7,0.7])

    ylim([-0.05 1.05]);
    set(gca, 'Box', 'off');
    title('Probability correct')
end

if isfield(stats,'ck')
    subplot(4,1,4)
        plot(nPlot, estimated_pCorrect, 'k')
    hold on;
    plot(nPlot, simulated_pCorrect, 'color',[0.7, 0.7, 0.7,0.7])

    ylim([-0.05 1.05]);
    set(gca, 'Box', 'off');
    title('Probability correct')
end

print(gcf,'-dpng',plotLabel);    %png format
saveas(gcf, plotLabel, 'fig');