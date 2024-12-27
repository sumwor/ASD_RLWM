function responseTime = plot_responseTimes(resultdf,edges,tlabel, BehPath)
% % plot_session %
%PURPOSE:   Plot performance of a single session of two-player game
%AUTHORS:   AC Kwan 180424
%
%INPUT ARGUMENTS
%   stats:  stats of the game
%   n_plot: plot up to this number of trials
%   tlabel: Text string that will be put on top of the figure

% get response time 
rt = resultdf.side_in-resultdf.center_in;

% left correct, right correct, left incorrect, right incorrect
rtHist = zeros(4, length(edges)-1);


for ii = 1:length(rt)
    rtIdx =ceil((rt(ii)-0.5)/0.05);
    if rtIdx > 0
        if resultdf.odors(ii) > 6.5 & resultdf.reward(ii)>0
            rtHist(1, rtIdx) = rtHist(1, rtIdx)+1;
        elseif resultdf.odors(ii) < 6.5 & isnan(resultdf.reward(ii))
            rtHist(2, rtIdx) = rtHist(2, rtIdx)+1;
        elseif resultdf.odors(ii) < 6.5 & resultdf.reward(ii)>0
            rtHist(3, rtIdx) = rtHist(3, rtIdx)+1;
        elseif resultdf.odors(ii) > 6.5 & isnan(resultdf.reward(ii))
            rtHist(4, rtIdx) = rtHist(4, rtIdx)+1;
        end
    end
end

responseTime.leftCorrect = rtHist(1,:);
responseTime.leftIncorrect = rtHist(2,:);
responseTime.rightCorrect = rtHist(3,:);
responseTime.rightIncorrect = rtHist(4,:);

setup_figprop;


%%
figure;
sgtitile('Response time (s)')
subplot(2,2,1);
bar(edges(1:end-1),rtHist(1,:),'FaceColor','r','EdgeColor','r');  
title('Left correct');
set(gca,'box','off')
subplot(2,2,2);
bar(edges(1:end-1),rtHist(2,:),'FaceColor','k', 'EdgeColor','k');
title('Left incorrect');
set(gca,'box','off')
subplot(2,2,3);
bar(edges(1:end-1),rtHist(3,:),'FaceColor','b'); 
title('Right correct');
set(gca,'box','off')
subplot(2,2,4);
bar(edges(1:end-1),rtHist(4,:),'FaceColor','k');
title('Right incorrect');
set(gca,'box','off')


print(gcf,'-dpng',fullfile(BehPath, ['response-time_', tlabel]));    %png format
%saveas(gcf, fullfile(BehPath, ['response-time_', tlabel]), 'fig');
%saveas(gcf, fullfile(BehPath, ['response-time_', tlabel]),'svg');
end

