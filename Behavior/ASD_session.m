function ASD_session(resultdf, protocol, animal, session, BehPath)
% % bandit_session %
%PURPOSE:   Preparing to analyze a single session of mouse behavior
%AUTHORS:   H Atilgan and AC Kwan 191203
%
%INPUT ARGUMENTS
%   BehPath:        path for the location of the analysis folder containing
%                   the behavioral .mat file
%   LogFileName:    name of the logfile
%
%OUTPUT ARGUMENTS
%
setup_figprop;
%% load the behavioral data
disp('-----------------------------------------------------------');
disp('--- Analyzing a single behavioral session ');
disp('-----------------------------------------------------------');

% Get trial information

%% action 0 -> odor 6; action 1 -> odor 7
% response time = side_in - center_in
% intertrial interval = center_in - center_in

%% session summary

% What to put as title for some of the figures generated
tlabel=strcat('Subject=',animal,', Time=',session);
% Create a subfolder to save the images for this session
% folder named using year/month/day of file
savebehfigpath = fullfile(BehPath,'Figures',session);
if ~exist(savebehfigpath)
    mkdir(savebehfigpath)
end

 %% analysis of behavioral performance
    
    % plot behavior in raw format
    
    plot_session(resultdf,size(resultdf,1),tlabel, savebehfigpath);



    %% plot response times
    valLabel='Response time (s)';    
    edges=[-0.5:0.05:10];
    responseTime = plot_responseTimes(resultdf, protocol, edges, tlabel, savebehfigpath);


    %% plot ITI
    edges=[0:0.5:30]; %ill-defined ITI for last trial
    itiTime = plot_itiTimes(resultdf, protocol, edges, tlabel, savebehfigpath);
    
    %% check error rate
    error_rate = plot_errors(resultdf, protocol, edges, tlabel, savebehfigpath);

    %% logistic regression
    % 1) p(correct) = rewarded choice + unreward choice (stimulus is redundent here)
    % 2) p(correct) = choice + choicexreward
    tlabel=strcat('Subject=',animal,', Time=',session);

    step_back = 15;
    [lregRCUC_output, negloglike, bic, nlike]=logreg_RCUC(resultdf,step_back); 
    plot_logreg(lregRCUC_output,tlabel);
    print(gcf,'-dpng','logregRCUC');    %png format
    saveas(gcf, 'logregCRInt', 'fig');
    
    [lregCRInt_output, negloglike, bic, nlike]=logreg_CRInt(resultdf,step_back); 
    plot_logreg(lregCRInt_output,tlabel);
    print(gcf,'-dpng','logregCRInt');    %png format
    saveas(gcf, 'logregCRInt', 'fig');


    
    %[~, ~, bic.logregCRInt2, nlike.logregCRInt2]=logreg_CRInt(stats,1,2);
    %[~, ~, bic.logregCRInt5, nlike.logregCRInt5]=logreg_CRInt(stats,1,5);
    %[~, ~, bic.logregCRInt10, nlike.logregCRInt10]=logreg_CRInt(stats,1,10);
    

    saveanalysispath = fullfile(BehPath, 'behAnalysis');
    if ~exist(saveanalysispath)
        mkdir(saveanalysispath);
    end
    save(fullfile(saveanalysispath,'beh_analysis.mat'),...
            'responseTime', 'itiTime');
    close all;
    clearvars -except i dirs expData;
end