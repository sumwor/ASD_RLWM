function plot_session(resultdf,n_plot,tlabel, BehPath)
% % plot_session %
%PURPOSE:   Plot performance of a single session of two-player game
%AUTHORS:   AC Kwan 180424
%
%INPUT ARGUMENTS
%   stats:  stats of the game
%   n_plot: plot up to this number of trials
%   tlabel: Text string that will be put on top of the figure

figpath = fullfile(BehPath, ['session-beh_', tlabel,'.png']);
%if ~exist(figpath)
setup_figprop;

%%
availOdors = unique(resultdf.schedule);
figure;

value2Plot = (-1).^(resultdf.schedule + 1) .* ceil(resultdf.schedule / 2);

subplot(3,1,1);
bar(value2Plot.*(value2Plot<0),1,'FaceColor','r','EdgeColor','none');
hold on;
bar(value2Plot.*(value2Plot>0),1,'FaceColor','b','EdgeColor','none');
hold on;
for ii = 1:length(resultdf.reward)
    if resultdf.reward(ii,1) > 0
        if resultdf.actions(ii) == 1
            hold on; scatter(ii,max(value2Plot)+0.25,'.','red');  %reward right
        elseif resultdf.actions(ii) == 0
            hold on; scatter(ii,-max(value2Plot)-0.25,'.','red');
        end
    end
end

% odor code needs to be changed
% available schedule: 1,2,3,4,5,6

% look for odor switch
switchTrial = [nan,nan];
if max(value2Plot) >= 2
    switchTrial(1) = find(value2Plot == 2 | value2Plot == -2, 1);
    if max(value2Plot) == 3
        switchTrial(2) = find(value2Plot == 3 | value2Plot == -3, 1);
    end
end

leftRewardRate = zeros(1,3);
rightRewardRate = zeros(1,3);

if isnan(switchTrial(1))
    leftRewardRate(1) = sum(resultdf.schedule==1 & resultdf.reward>0)/sum(resultdf.schedule==1);
    rightRewardRate(1) = sum(resultdf.schedule==2 & resultdf.reward>0)/sum(resultdf.schedule==2);
else
    leftRewardRate(1) = sum(resultdf.schedule(1:switchTrial(1))==1 & resultdf.reward(1:switchTrial(1))>0)/sum(resultdf.schedule(1:switchTrial(1))==1);
    rightRewardRate(1) = sum(resultdf.schedule(1:switchTrial(1))==2 & resultdf.reward(1:switchTrial(1))>0)/sum(resultdf.schedule(1:switchTrial(1))==2);

    if isnan(switchTrial(2))
        leftRewardRate(2) = sum(resultdf.schedule(switchTrial(1):end)==3 & resultdf.reward(switchTrial(1):end)>0)/sum(resultdf.schedule(switchTrial(1):end)==3);
        rightRewardRate(2) = sum(resultdf.schedule(switchTrial(1):end)==4 & resultdf.reward(switchTrial(1):end)>0)/sum(resultdf.schedule(switchTrial(1):end)==4);
    else
        leftRewardRate(2) = sum(resultdf.schedule(switchTrial(1):switchTrial(2))==3 & resultdf.reward(switchTrial(1):switchTrial(2))>0)/sum(resultdf.schedule(switchTrial(1):switchTrial(2))==3);
        rightRewardRate(2) = sum(resultdf.schedule(switchTrial(1):switchTrial(2))==4 & resultdf.reward(switchTrial(1):switchTrial(2))>0)/sum(resultdf.schedule(switchTrial(1):switchTrial(2))==4);
        leftRewardRate(3) = sum(resultdf.schedule(switchTrial(2):end)==5 & resultdf.reward(switchTrial(2):end)>0)/sum(resultdf.schedule(switchTrial(1):end)==5);
        rightRewardRate(3) = sum(resultdf.schedule(switchTrial(2):end)==6 & resultdf.reward(switchTrial(2):end)>0)/sum(resultdf.schedule(switchTrial(1):end)==6);
    end
end

if isnan(switchTrial(1))
    text(n_plot-100, 2, [num2str(leftRewardRate(1))], 'FontSize', 12, 'Color', 'red');
    text(n_plot-100, -2, [ num2str(rightRewardRate(1))], 'FontSize', 12, 'Color', 'red');
elseif isnan(switchTrial(2))
    text(switchTrial(1)-100, 3, [num2str(leftRewardRate(1))], 'FontSize', 12, 'Color', 'red');
    text(switchTrial(1)-100, -3, [num2str(rightRewardRate(1))], 'FontSize', 12, 'Color', 'red');
    text(n_plot-200, 3, [ num2str(leftRewardRate(2))], 'FontSize', 12, 'Color', 'red');
    text(n_plot-200, -3, [ num2str(rightRewardRate(2))], 'FontSize', 12, 'Color', 'red');
else
    text(switchTrial(1)-100, 4, [num2str(leftRewardRate(1))], 'FontSize', 12, 'Color', 'red');
    text(switchTrial(1)-100, -4, [num2str(rightRewardRate(1))], 'FontSize', 12, 'Color', 'red');
    text(switchTrial(2)-100, 4, [num2str(leftRewardRate(2))], 'FontSize', 12, 'Color', 'red');
    text(switchTrial(2)-100, -4, [num2str(rightRewardRate(2))], 'FontSize', 12, 'Color', 'red');
    text(n_plot-100, 4, [num2str(leftRewardRate(3))], 'FontSize', 12, 'Color', 'red');
    text(n_plot-100, -4, [num2str(rightRewardRate(3))], 'FontSize', 12, 'Color', 'red');
end

xlim([0 n_plot]);

ylim([-max(value2Plot)-0.5 max(value2Plot)+0.5]);
if max(value2Plot) == 1
    set(gca,'ytick',[-1.75 -1 1 1.75]);
    set(gca,'yticklabel',{'Reward','1','2','Reward'});
elseif max(value2Plot) == 2
    plot([switchTrial(1),switchTrial(1)], [-3,3], 'black', 'LineWidth', 3);
    set(gca,'ytick',[-2.75 -2 -1 1 2 2.75]);
    set(gca,'yticklabel',{'Reward','1', '3','2', '4','Reward'});
elseif max(value2Plot) == 3
    plot([switchTrial(1),switchTrial(1)], [-4,4], 'black', 'LineWidth', 3);
    plot([switchTrial(2),switchTrial(2)], [-4,4], 'black', 'LineWidth', 3);
    set(gca,'ytick',[-3.75 -3 -2 -1 1 2 3 3.75]);
    set(gca,'yticklabel',{'Reward','1', '3', '5','2', '4','6','Reward'});
end
set(gca, 'Box', 'off');
title(tlabel);

print(gcf,'-dpng',fullfile(BehPath, ['session-beh_', tlabel]));    %png format
saveas(gcf, fullfile(BehPath, ['session-beh_', tlabel]), 'fig');
saveas(gcf, fullfile(BehPath, ['session-beh_', tlabel]),'svg');
end
%end

