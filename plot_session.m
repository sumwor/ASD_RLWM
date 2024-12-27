function plot_session(resultdf,n_plot,tlabel, BehPath)
% % plot_session %
%PURPOSE:   Plot performance of a single session of two-player game
%AUTHORS:   AC Kwan 180424
%
%INPUT ARGUMENTS
%   stats:  stats of the game
%   n_plot: plot up to this number of trials
%   tlabel: Text string that will be put on top of the figure

setup_figprop;

%%
figure;
subplot(3,1,1);
bar(-1*(resultdf.actions==0),1,'FaceColor','r','EdgeColor','none');
hold on;
bar((resultdf.actions==1),1,'FaceColor','b','EdgeColor','none');
hold on;
for ii = 1:length(resultdf.reward)
    if resultdf.reward(ii,1) > 0
        if resultdf.actions(ii) == 1
            hold on; plot([ii,ii],[1.35,1.75],'black');  %reward right
        elseif resultdf.actions(ii) == 0
            hold on; plot([ii,ii],[-1.35,-1.75],'black');
        end
    end
end

leftRewardRate = sum(resultdf.odors>6.5 & resultdf.reward>0)/sum(resultdf.odors>6.5);
rightRewardRate = sum(resultdf.odors<6.5 & resultdf.reward>0)/sum(resultdf.odors<6.5);

text(n_plot-100, 2.1, ['Reward rate: ', num2str(leftRewardRate)], 'FontSize', 12, 'Color', 'red');
text(n_plot-100, -2.1, ['Reward rate: ', num2str(rightRewardRate)], 'FontSize', 12, 'Color', 'red');

xlim([0 n_plot]);
ylim([-2.5 2.5]);
set(gca,'ytick',[-1.75 -1 1 1.75]);
set(gca,'yticklabel',{'Reward','L','R','Reward'});
set(gca, 'Box', 'off');
title(tlabel);

print(gcf,'-dpng',fullfile(BehPath, ['session-beh_', tlabel]));    %png format
saveas(gcf, fullfile(BehPath, ['session-beh_', tlabel]), 'fig');
saveas(gcf, fullfile(BehPath, ['session-beh_', tlabel]),'svg');
end

