function responseTime = plot_responseTimes(resultdf,protocol, edges,tlabel, BehPath)
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

%% AB odor always exists in every session
rtHist_AB = zeros(4, length(edges)-1);

% n_respon


for ii = 1:length(rt)
    rtIdx =ceil((rt(ii)+0.5)/0.05);
    if rtIdx > 0 & rtIdx < length(edges)
        if resultdf.schedule(ii) == 1 & resultdf.reward(ii)>0 % left rewarded
            rtHist_AB(1, rtIdx) = rtHist_AB(1, rtIdx)+1;
        elseif resultdf.schedule(ii) == 1 & isnan(resultdf.reward(ii)) % left unrewarded
            rtHist_AB(2, rtIdx) = rtHist_AB(2, rtIdx)+1;
        elseif resultdf.schedule(ii) == 2 & resultdf.reward(ii)>0 % right rewarded
            rtHist_AB(3, rtIdx) = rtHist_AB(3, rtIdx)+1;
        elseif resultdf.schedule(ii) == 2 & isnan(resultdf.reward(ii)) % right rewarded
            rtHist_AB(4, rtIdx) = rtHist_AB(4, rtIdx)+1;
        end
    end
end

responseTime.AB.leftCorrect = rtHist_AB(1,:);
responseTime.AB.leftIncorrect = rtHist_AB(2,:);
responseTime.AB.rightCorrect = rtHist_AB(3,:);
responseTime.AB.rightIncorrect = rtHist_AB(4,:);


%% CD odor
if strcmp(protocol, 'AB-CD') | strcmp(protocol, 'AB-CD-DC')
    rtHist_CD = zeros(4, length(edges)-1);

% n_respon


for ii = 1:length(rt)
    rtIdx =ceil((rt(ii)+0.5)/0.05);
    if rtIdx > 0 & rtIdx < length(edges)
        if resultdf.schedule(ii) == 3 & resultdf.reward(ii)>0 % left rewarded
            rtHist_CD(1, rtIdx) = rtHist_CD(1, rtIdx)+1;
        elseif resultdf.schedule(ii) == 3 & isnan(resultdf.reward(ii)) % left unrewarded
            rtHist_CD(2, rtIdx) = rtHist_CD(2, rtIdx)+1;
        elseif resultdf.schedule(ii) == 4 & resultdf.reward(ii)>0 % right rewarded
            rtHist_CD(3, rtIdx) = rtHist_CD(3, rtIdx)+1;
        elseif resultdf.schedule(ii) == 4 & isnan(resultdf.reward(ii)) % right rewarded
            rtHist_CD(4, rtIdx) = rtHist_CD(4, rtIdx)+1;
        end
    end
end

responseTime.CD.leftCorrect = rtHist_CD(1,:);
responseTime.CD.leftIncorrect = rtHist_CD(2,:);
responseTime.CD.rightCorrect = rtHist_CD(3,:);
responseTime.CD.rightIncorrect = rtHist_CD(4,:);

end

%% DC
if strcmp(protocol, 'AB-DC') | strcmp(protocol, 'AB-CD-DC')
    rtHist_DC = zeros(4, length(edges)-1);

% n_respon


for ii = 1:length(rt)
    rtIdx =ceil((rt(ii)+0.5)/0.05);
    if rtIdx > 0 & rtIdx < length(edges)
        if resultdf.schedule(ii) == 5 & resultdf.reward(ii)>0 % left rewarded
            rtHist_DC(1, rtIdx) = rtHist_DC(1, rtIdx)+1;
        elseif resultdf.schedule(ii) == 5 & isnan(resultdf.reward(ii)) % left unrewarded
            rtHist_DC(2, rtIdx) = rtHist_DC(2, rtIdx)+1;
        elseif resultdf.schedule(ii) == 6 & resultdf.reward(ii)>0 % right rewarded
            rtHist_DC(3, rtIdx) = rtHist_DC(3, rtIdx)+1;
        elseif resultdf.schedule(ii) == 6 & isnan(resultdf.reward(ii)) % right rewarded
            rtHist_DC(4, rtIdx) = rtHist_DC(4, rtIdx)+1;
        end
    end
end

responseTime.DC.leftCorrect = rtHist_DC(1,:);
responseTime.DC.leftIncorrect = rtHist_DC(2,:);
responseTime.DC.rightCorrect = rtHist_DC(3,:);
responseTime.DC.rightIncorrect = rtHist_DC(4,:);

end

%%
figure;
sgtitle('Response time AB (s)')
subplot(2,2,1);
bar(edges(1:end-1),rtHist_AB(1,:),'FaceColor','r','EdgeColor','r');  
title('Left correct');
set(gca,'box','off')
subplot(2,2,2);
bar(edges(1:end-1),rtHist_AB(2,:),'FaceColor','k', 'EdgeColor','k');
title('Left incorrect');
set(gca,'box','off')
subplot(2,2,3);
bar(edges(1:end-1),rtHist_AB(3,:),'FaceColor','b'); 
title('Right correct');
set(gca,'box','off')
subplot(2,2,4);
bar(edges(1:end-1),rtHist_AB(4,:),'FaceColor','k');
title('Right incorrect');
set(gca,'box','off')


print(gcf,'-dpng',fullfile(BehPath, ['response-time-AB_', tlabel]));    %png format
%saveas(gcf, fullfile(BehPath, ['response-time_', tlabel]), 'fig');
%saveas(gcf, fullfile(BehPath, ['response-time_', tlabel]),'svg');

if strcmp(protocol, 'AB-CD') | strcmp(protocol, 'AB-CD-DC')
    figure;
sgtitle('Response time CD (s)')
subplot(2,2,1);
bar(edges(1:end-1),rtHist_CD(1,:),'FaceColor','r','EdgeColor','r');  
title('Left correct');
set(gca,'box','off')
subplot(2,2,2);
bar(edges(1:end-1),rtHist_CD(2,:),'FaceColor','k', 'EdgeColor','k');
title('Left incorrect');
set(gca,'box','off')
subplot(2,2,3);
bar(edges(1:end-1),rtHist_CD(3,:),'FaceColor','b'); 
title('Right correct');
set(gca,'box','off')
subplot(2,2,4);
bar(edges(1:end-1),rtHist_CD(4,:),'FaceColor','k');
title('Right incorrect');
set(gca,'box','off')


print(gcf,'-dpng',fullfile(BehPath, ['response-time-CD_', tlabel]));   
end

if strcmp(protocol, 'AB-DC') | strcmp(protocol, 'AB-CD-DC')
    figure;
sgtitle('Response time DC (s)')
subplot(2,2,1);
bar(edges(1:end-1),rtHist_DC(1,:),'FaceColor','r','EdgeColor','r');  
title('Left correct');
set(gca,'box','off')
subplot(2,2,2);
bar(edges(1:end-1),rtHist_DC(2,:),'FaceColor','k', 'EdgeColor','k');
title('Left incorrect');
set(gca,'box','off')
subplot(2,2,3);
bar(edges(1:end-1),rtHist_DC(3,:),'FaceColor','b'); 
title('Right correct');
set(gca,'box','off')
subplot(2,2,4);
bar(edges(1:end-1),rtHist_DC(4,:),'FaceColor','k');
title('Right incorrect');
set(gca,'box','off')


print(gcf,'-dpng',fullfile(BehPath, ['response-time-DC_', tlabel]));   
end

end

