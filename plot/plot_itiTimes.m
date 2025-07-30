function itiTime = plot_itiTimes(resultdf, protocol, edges,tlabel, BehPath)
% % plot_session %
%PURPOSE:   Plot performance of a single session of two-player game
%AUTHORS:   AC Kwan 180424
%
%INPUT ARGUMENTS
%   stats:  stats of the game
%   n_plot: plot up to this number of trials
%   tlabel: Text string that will be put on top of the figure

% get response time 
iti = resultdf.center_in(2:end)-resultdf.center_in(1:end-1);

% left correct, right correct, left incorrect, right incorrect
itiHist_AB = zeros(4, length(edges)-1);


for ii = 1:length(iti)
    itiIdx =ceil((iti(ii)-0)/1);

    if itiIdx > 0 & itiIdx < length(edges)
        if resultdf.schedule(ii) == 1 & resultdf.reward(ii)>0 % left correct
            itiHist_AB(1, itiIdx) = itiHist_AB(1, itiIdx)+1;
        elseif resultdf.schedule(ii) ==1 & isnan(resultdf.reward(ii)) % left incorrect
            itiHist_AB(2, itiIdx) = itiHist_AB(2, itiIdx)+1;
        elseif resultdf.schedule(ii) ==2 & resultdf.reward(ii)>0 % right correct
            itiHist_AB(3, itiIdx) = itiHist_AB(3, itiIdx)+1;
        elseif resultdf.schedule(ii) ==2 & isnan(resultdf.reward(ii)) % right incorrect
            itiHist_AB(4, itiIdx) = itiHist_AB(4, itiIdx)+1;
        end
    end
end


itiTime.AB.leftCorrect = itiHist_AB(1,:);
itiTime.AB.leftIncorrect = itiHist_AB(2,:);
itiTime.AB.rightCorrect = itiHist_AB(3,:);
itiTime.AB.rightIncorrect = itiHist_AB(4,:);

%% CD performance
if strcmp(protocol, 'AB-CD') | strcmp(protocol, 'AB-CD-DC')
itiHist_CD = zeros(4, length(edges)-1);


for ii = 1:length(iti)
    itiIdx =ceil((iti(ii)-0)/1);

    if itiIdx > 0 & itiIdx < length(edges)
        if resultdf.schedule(ii) == 3 & resultdf.reward(ii)>0 % left correct
            itiHist_CD(1, itiIdx) = itiHist_CD(1, itiIdx)+1;
        elseif resultdf.schedule(ii) ==3 & isnan(resultdf.reward(ii)) % left incorrect
            itiHist_CD(2, itiIdx) = itiHist_CD(2, itiIdx)+1;
        elseif resultdf.schedule(ii) ==4 & resultdf.reward(ii)>0 % right correct
            itiHist_CD(3, itiIdx) = itiHist_CD(3, itiIdx)+1;
        elseif resultdf.schedule(ii) ==4 & isnan(resultdf.reward(ii)) % right incorrect
            itiHist_CD(4, itiIdx) = itiHist_CD(4, itiIdx)+1;
        end
    end
end


itiTime.CD.leftCorrect = itiHist_CD(1,:);
itiTime.CD.leftIncorrect = itiHist_CD(2,:);
itiTime.CD.rightCorrect = itiHist_CD(3,:);
itiTime.CD.rightIncorrect = itiHist_CD(4,:);
end

%% DC performance
if strcmp(protocol, 'AB-DC') | strcmp(protocol, 'AB-CD-DC')
itiHist_DC = zeros(4, length(edges)-1);


for ii = 1:length(iti)
    itiIdx =ceil((iti(ii)-0)/1);

    if itiIdx > 0 & itiIdx < length(edges)
        if resultdf.schedule(ii) == 5 & resultdf.reward(ii)>0 % left correct
            itiHist_DC(1, itiIdx) = itiHist_DC(1, itiIdx)+1;
        elseif resultdf.schedule(ii) ==5 & isnan(resultdf.reward(ii)) % left incorrect
            itiHist_DC(2, itiIdx) = itiHist_DC(2, itiIdx)+1;
        elseif resultdf.schedule(ii) ==6 & resultdf.reward(ii)>0 % right correct
            itiHist_DC(3, itiIdx) = itiHist_DC(3, itiIdx)+1;
        elseif resultdf.schedule(ii) ==6 & isnan(resultdf.reward(ii)) % right incorrect
            itiHist_DC(4, itiIdx) = itiHist_DC(4, itiIdx)+1;
        end
    end
end


itiTime.DC.leftCorrect = itiHist_DC(1,:);
itiTime.DC.leftIncorrect = itiHist_DC(2,:);
itiTime.DC.rightCorrect = itiHist_DC(3,:);
itiTime.DC.rightIncorrect = itiHist_DC(4,:);
end

%%
figure;
sgtitle('Intertrial interval AB (s)')
subplot(2,2,1);
bar(edges(1:end-1),itiHist_AB(1,:),'FaceColor','r','EdgeColor','r');  
title('Left correct');
set(gca,'box','off')
subplot(2,2,2);
bar(edges(1:end-1),itiHist_AB(2,:),'FaceColor','k', 'EdgeColor','k');
title('Left incorrect');
set(gca,'box','off')
subplot(2,2,3);
bar(edges(1:end-1),itiHist_AB(3,:),'FaceColor','b'); 
title('Right correct');
set(gca,'box','off')
subplot(2,2,4);
bar(edges(1:end-1),itiHist_AB(4,:),'FaceColor','k');
title('Right incorrect');
set(gca,'box','off')


print(gcf,'-dpng',fullfile(BehPath, ['iti-time_AB', tlabel]));    %png format
%saveas(gcf, fullfile(BehPath, ['iti-time_', tlabel]), 'fig');
%saveas(gcf, fullfile(BehPath, ['iti-time_', tlabel]),'svg');

if strcmp(protocol, 'AB-CD') | strcmp(protocol, 'AB-CD-DC')
   figure;
sgtitle('Intertrial interval CD (s)')
subplot(2,2,1);
bar(edges(1:end-1),itiHist_CD(1,:),'FaceColor','r','EdgeColor','r');  
title('Left correct');
set(gca,'box','off')
subplot(2,2,2);
bar(edges(1:end-1),itiHist_CD(2,:),'FaceColor','k', 'EdgeColor','k');
title('Left incorrect');
set(gca,'box','off')
subplot(2,2,3);
bar(edges(1:end-1),itiHist_CD(3,:),'FaceColor','b'); 
title('Right correct');
set(gca,'box','off')
subplot(2,2,4);
bar(edges(1:end-1),itiHist_CD(4,:),'FaceColor','k');
title('Right incorrect');
set(gca,'box','off')


print(gcf,'-dpng',fullfile(BehPath, ['iti-time_CD', tlabel]));    %png format
%saveas(gcf, fullfile(BehPath, ['iti-time_', tlabel]), 'fig');
%saveas(gcf, fullfile(BehPath, ['iti-time_', tlabel]),'svg');
end

if strcmp(protocol, 'AB-DC') | strcmp(protocol, 'AB-CD-DC')
   figure;
sgtitle('Intertrial interval DC (s)')
subplot(2,2,1);
bar(edges(1:end-1),itiHist_DC(1,:),'FaceColor','r','EdgeColor','r');  
title('Left correct');
set(gca,'box','off')
subplot(2,2,2);
bar(edges(1:end-1),itiHist_DC(2,:),'FaceColor','k', 'EdgeColor','k');
title('Left incorrect');
set(gca,'box','off')
subplot(2,2,3);
bar(edges(1:end-1),itiHist_DC(3,:),'FaceColor','b'); 
title('Right correct');
set(gca,'box','off')
subplot(2,2,4);
bar(edges(1:end-1),itiHist_DC(4,:),'FaceColor','k');
title('Right incorrect');
set(gca,'box','off')


print(gcf,'-dpng',fullfile(BehPath, ['iti-time_DC', tlabel]));    %png format
%saveas(gcf, fullfile(BehPath, ['iti-time_', tlabel]), 'fig');
%saveas(gcf, fullfile(BehPath, ['iti-time_', tlabel]),'svg');
end

end

