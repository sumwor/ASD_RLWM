function ASD_odor_summary(dataIndex, strain, savefigpath, savedatapath)

%% go over every session to plot performance in blocks for
% 1. first two sessions of AB
% 2. first two sessions of CD
% 3. DC reverse

%% number of AB/CD/DC trials experienced based on training session
    % training session: 3 AB + 3 AB-CD + 6 AB-(CD) - DC sessions
    % ignore the retrain sessions for now

nFiles = size(dataIndex,1);
Subjects = unique(dataIndex.Animal);
nSubjects = length(unique(dataIndex.Animal));

genotype = cell(nSubjects,1);

% number of trials
nAB_trials = nan(nSubjects, 12);
nCD_trials = nan(nSubjects, 4);
nDC_trials = nan(nSubjects, 6);

% trial timing
timing = 1:8*60; % in minutes (8 hour maximum)
perf_timing = 1:30:8*60;
Trials_time = nan(nFiles, length(timing));
performance_time = nan(nFiles, length(perf_timing));

% performance
blockLength = 100; % check performance in blocks of 100 trials
perf_AB_block = nan(nSubjects, 50);
perf_CD_block = nan(nSubjects, 50);
perf_DC_block = nan(nSubjects, 50);

perf_AB_quantile = nan(nSubjects,4, 3); % dim3: first session, second session
perf_CD_quantile = nan(nSubjects,4, 3); % dim3: session
perf_DC_quantile = nan(nSubjects,4, 6); % dim4: session

perf_AB_block_session = nan(nSubjects,20, 3); % dim3: first session, second session
perf_CD_block_session = nan(nSubjects,20, 3); % dim3: session
perf_DC_block_session = nan(nSubjects,20, 6); % dim4: session

perf_DC_session = nan(nSubjects,6);

%% trial number and performance
for ii = 1:nSubjects

    genotype(ii) = unique(dataIndex.Genotype(strcmp(dataIndex.Animal, Subjects{ii})));
    %% look for first 3 AB sessions
    subDataIndex = dataIndex(strcmp(dataIndex.Animal, Subjects{ii}),:);
    % load the first three sessions
    protocol = 'AB';

    nSessions = 3;
    results = cell(nSessions,1);
    for ss=1:nSessions
        if ss<= size(subDataIndex,1)
            csvFile = fullfile(subDataIndex.BehPath{ss}, ...
                subDataIndex.BehCSV{ss});
            if exist(csvFile)
                results{ss} = readtable(csvFile);
                nAB_trials(ii, ss) = size(results{ss},1);
            else
                results{ss} = NaN;
                nAB_trials(ii,ss) = NaN;
            end
        else
            results{ss} = NaN;
            nAB_trials(ii,ss) = NaN;
        end
    end
    % calculate performance in quantile
    perf_AB_quantile(ii,:,:) = perf_in_quantile(results, protocol);
    
    % calculate performance in quantile on a session basis
    perf_AB_block_session(ii,:,:) = perf_in_block_session(results, protocol, blockLength);

    % calculate performance in blocks
    %perf_AB_block(ii,:) = perf_in_block(results, protocol, blockLength);

    %% for CD odors
    animalMask = strcmp(dataIndex.Animal, Subjects{ii});
    stageMask = cellfun(@(x) isequal(x, [1;2;3;4]), dataIndex.OdorPresented);
    subDataIndex = dataIndex(animalMask & stageMask,:);

    protocol = 'CD';
    nSessions = 3;
    results = cell(nSessions,1);
    for ss=1:nSessions
        if ss<= size(subDataIndex,1)
            csvFile = fullfile(subDataIndex.BehPath{ss}, ...
                subDataIndex.BehCSV{ss});

            results{ss} = readtable(csvFile);       

            % find the switch point
            switchTrial = find(results{ss}.schedule == 3 | results{ss}.schedule == 4, 1);
            results{ss} = results{ss}(switchTrial:end,:);
            
            nAB_trials(ii, ss+3) = switchTrial-1;
            nCD_trials(ii, ss) = size(results{ss},1);
        else
            results{ss} = NaN;
        end
    end
    % calculate performance in quantile
    perf_CD_quantile(ii,:,:) = perf_in_quantile(results, protocol);
    
     % calculate performance in quantile on a session basis
    perf_CD_block_session(ii,:,:) = perf_in_block_session(results, protocol, blockLength);

 
    % calculate performance in blocks
    %perf_CD_block(ii,:) = perf_in_block(results, protocol,blockLength);

    %% for DC performance

    animalMask = strcmp(dataIndex.Animal, Subjects{ii});
    stageMask = cellfun(@(x) isequal(x, [1;2;5;6]), dataIndex.OdorPresented) | cellfun(@(x) isequal(x, [1;2;3;4;5;6]), dataIndex.OdorPresented);
    subDataIndex = dataIndex(animalMask & stageMask,:);

    if ~isempty(subDataIndex)
        protocol = 'DC';
        nSessions = 6;
        results = cell(nSessions,1);
        for ss =1 :nSessions
            if ss<= size(subDataIndex,1)
                csvFile = fullfile(subDataIndex.BehPath{ss}, ...
                    subDataIndex.BehCSV{ss});

                results{ss} = readtable(csvFile);

                % find the switch point
                switchCD =  find(results{ss}.schedule == 3 | results{ss}.schedule == 4, 1);
                switchTrial = find(results{ss}.schedule == 5 | results{ss}.schedule == 6, 1);
                results{ss} = results{ss}(switchTrial:end,:);

                if switchCD > 0
                    nAB_trials(ii, ss+6) = switchCD-1;
                    nCD_trials(ii,ss+3) = switchTrial-switchCD;
                else
                    nAB_trials(ii,ss+6) = switchTrial-1;
                end
                nDC_trials(ii, ss) = size(results{ss},1);
            else
                results{ss} = NaN;
            end
        end
        % calculate performance in quantile
        perf_DC_quantile(ii,:,:) = perf_in_quantile(results, protocol);
        
         % calculate performance in quantile on a session basis
        perf_DC_block_session(ii,:,:) = perf_in_block_session(results, protocol, blockLength);

        % calculate performance in blocks
        %perf_DC_block(ii,:) = perf_in_block(results, protocol,blockLength);
    end
end

%% check DC performance based on how many AB/CD trials they have experienced
% check every session before starting DC reversal, including retrain
% session

nAB_trials_total = {};
nCD_trials_total = {};
nABinDC_trials_total = {};
for ii = 1:nSubjects
    %% look for first 3 AB sessions
    subDataIndex = dataIndex(strcmp(dataIndex.Animal, Subjects{ii}),:);
    tempTrial_AB = [];
    tempTrial_CD = [];
    tempTrial_ABinDC = [];
    reversed = 0; % check if the animal had experienced reversal already
    % used to identify AB retrain sessions in during reversal
    for ss = 1:length(subDataIndex.Protocol)
        csvFile = fullfile(subDataIndex.BehPath{ss}, ...
            subDataIndex.BehCSV{ss});
        if exist(csvFile)
            result = readtable(csvFile);
        if ~ strcmp('AB-CD-DC', subDataIndex.Protocol{ss})
                if reversed==0
                    Index_AB = find(result.schedule==3 | result.schedule==4,1);
                    if isempty(Index_AB)
                        tempTrial_AB = [tempTrial_AB, size(result,1)];
                    else
                        tempTrial_AB = [tempTrial_AB, Index_AB-1];
                        Index_CD = find(result.schedule==3|result.schedule==4, 1, 'last');
                        tempTrial_CD = [tempTrial_CD, Index_CD-Index_AB+1];
                    end
                else
                    Index_AB = find(result.schedule==5 | result.schedule==6,1);
                    if isempty(Index_AB)
                        tempTrial_ABinDC = [tempTrial_ABinDC, size(result,1)];
                    else
                        tempTrial_ABinDC = [tempTrial_ABinDC, Index_AB-1];
                    end
                end
                elseif strcmp('AB-CD-DC', subDataIndex.Protocol{ss})
                    reversed =1;
                    Index_AB = find(result.schedule==3 | result.schedule==4,1);
                    if isempty(Index_AB)
                    tempTrial_AB = [tempTrial_AB, size(result,1)];
                else
                    tempTrial_AB = [tempTrial_AB, Index_AB-1];
                    Index_CD = find(result.schedule==3|result.schedule==4, 1, 'last');
                    tempTrial_CD = [tempTrial_CD, Index_CD-Index_AB+1];
                end
            end
        else
            tempTrial_AB=NaN;
            tempTrial_CD = NaN;
            tempTrial_DC = NaN;
        end
    end
    nAB_trials_total{ii}=tempTrial_AB;
    nCD_trials_total{ii} = tempTrial_CD;
    nABinDC_trials_total{ii} = tempTrial_ABinDC;
end

%% plot DC learning curve in blocks based on number of AB trials experienced
% AB_trials = zeros(nSubjects,1);
% CD_trials = zeros(nSubjects,1);
% 
% figure;
% hold on;
% 
% % Define colors for different animals
% colors = lines(nSubjects); 
% perf_DC_session = squeeze(nanmean(perf_DC_quantile,2));
% 
% subplot(1,2,1)
% hold on;
% % Plot each animal's DC learning curve
% for i = 1:nSubjects
%     AB_trials(i) = sum(nAB_trials_total{i});
%     CD_trials(i) = sum(nCD_trials_total{i});
%     if strcmp(genotype{i}, 'WT')
%         valid_idx = find(~isnan(perf_DC_session(i, :))); 
%         plot(1:6, perf_DC_session(i, :), 'Color', colors(i, :), 'LineWidth', 1.5); % DC curve
%         if ~isempty(valid_idx)
%             text(valid_idx(end), perf_DC_session(i, valid_idx(end)), sprintf('%d', sum(nAB_trials_total{i})), 'Color', colors(i, :), 'FontSize', 10); % AB trial numbers
%         end
%     end
% end
% xlabel('DC Learning Trials');
% ylabel('Performance');
% title('WT')
% 
% subplot(1,2,2)
% hold on;
% 
% if ismember('HEM', genotype)
%     hetGeno = 'HEM';
% elseif length(unique(genotype)) == 3
%     hetGeno = 'KO';
% else
%     hetGeno = 'HET';
% end
% % Plot each animal's DC learning curve
% for i = 1:nSubjects
%     if strcmp(genotype{i}, hetGeno)
%         valid_idx = find(~isnan(perf_DC_session(i, :))); 
%         plot(1:6, perf_DC_session(i, :), 'Color', colors(i, :), 'LineWidth', 1.5); % DC curve
%         if ~isempty(valid_idx)
%             text(valid_idx(end), perf_DC_session(i, valid_idx(end)), sprintf('%d', sum(nAB_trials_total{i})), 'Color', colors(i, :), 'FontSize', 10); % AB trial numbers
%         end
%     end
% end
% xlabel('DC Learning Trials');
% ylabel('Performance');
% title(hetGeno)
% % Labels and title
% 
% sgtitle('DC Learning Curve for Each Animal (Based on AB Trials Shown)');
% 
% % Adjust figure
% xlim([0 7]); % Extra space for text annotation
% hold off;
% 
% % save figure
% print(gcf,'-dpng',fullfile(savefigpath,['Reversal based on AB trials experienced']));    %png format
% saveas(gcf, fullfile(savefigpath, ['Reversal based on AB trials experienced']), 'fig');
% %savefig(fullfile(savefigpath, 'Number of trials performed.fig'));
% saveas(gcf, fullfile(savefigpath, ['Reversal based on AB trials experienced']),'svg');
% close();
% 
% % plot DC learning curve in blocks based on number of CD trials experienced
% figure;
% hold on;
% 
% % Define colors for different animals
% colors = lines(nSubjects); 
% perf_DC_session = squeeze(nanmean(perf_DC_quantile,2));
% 
% subplot(1,2,1)
% hold on;
% % Plot each animal's DC learning curve
% for i = 1:nSubjects
%     if strcmp(genotype{i}, 'WT')
%         valid_idx = find(~isnan(perf_DC_session(i, :))); 
%         plot(1:6, perf_DC_session(i, :), 'Color', colors(i, :), 'LineWidth', 1.5); % DC curve
%         if ~isempty(valid_idx)
%             text(valid_idx(end), perf_DC_session(i, valid_idx(end)), sprintf('%d', sum(nCD_trials_total{i})), 'Color', colors(i, :), 'FontSize', 10); % AB trial numbers
%         end
%     end
% end
% xlabel('DC Learning Trials');
% ylabel('Performance');
% title('WT')
% 
% subplot(1,2,2)
% hold on;
% % Plot each animal's DC learning curve
% 
% if ismember('HEM', genotype)
%     hetGeno = 'HEM';
% elseif length(unique(genotype)) == 3
%     hetGeno = 'KO';
% else
%     hetGeno = 'HET';
% end
% for i = 1:nSubjects
%     if strcmp(genotype{i}, hetGeno)
%         valid_idx = find(~isnan(perf_DC_session(i, :))); 
%         plot(1:6, perf_DC_session(i, :), 'Color', colors(i, :), 'LineWidth', 1.5); % DC curve
%         if ~isempty(valid_idx)
%             text(valid_idx(end), perf_DC_session(i, valid_idx(end)), sprintf('%d', sum(nCD_trials_total{i})), 'Color', colors(i, :), 'FontSize', 10); % AB trial numbers
%         end
%     end
% end
% xlabel('DC Learning Trials');
% ylabel('Performance');
% title(hetGeno)
% % Labels and title
% 
% sgtitle('DC Learning Curve for Each Animal (Based on CD Trials Shown)');
% 
% % Adjust figure
% xlim([0 7]); % Extra space for text annotation
% hold off;
% 
% % save figure
% print(gcf,'-dpng',fullfile(savefigpath,['Reversal based on CD trials experienced']));    %png format
% saveas(gcf, fullfile(savefigpath, ['Reversal based on CD trials experienced']), 'fig');
% %savefig(fullfile(savefigpath, 'Number of trials performed.fig'));
% saveas(gcf, fullfile(savefigpath, ['Reversal based on CD trials experienced']),'svg');
% close()
% 
% % scatter plot of number of AB and CD trials for WT and Het animals
% figure;
% subplot(1,2,1)
% scatter(AB_trials(find(strcmp(genotype, 'WT'))), ...
%     nanmean(perf_DC_session(find(strcmp(genotype, 'WT')), 5:6),2),150, 'b', 'filled');
% hold on;
% scatter(AB_trials(find(strcmp(genotype, hetGeno))), ...
%     nanmean(perf_DC_session(find(strcmp(genotype, hetGeno)), 5:6),2),150, 'r', 'filled');
% title('AB trials experienced')
% 
% WT_AB = AB_trials(find(strcmp(genotype, 'WT')));
% DC_perf_WT = nanmean(perf_DC_session(find(strcmp(genotype, 'WT')), 5:6),2);
% WT_AB_valid = WT_AB(~isnan(DC_perf_WT));
% DC_perf_WT_valid = DC_perf_WT(~isnan(DC_perf_WT));
% [r_WTAB, p_WTAB] = corr(WT_AB_valid, ...
%     DC_perf_WT_valid, 'Type', 'Pearson');
% 
% HET_AB = AB_trials(find(strcmp(genotype, hetGeno)));
% DC_perf_HET = nanmean(perf_DC_session(find(strcmp(genotype, hetGeno)), 5:6),2);
% HET_AB_valid = HET_AB(~isnan(DC_perf_HET));
% DC_perf_HET_valid = DC_perf_HET(~isnan(DC_perf_HET));
% [r_HETAB, p_HETAB] = corr(HET_AB_valid, ...
%     DC_perf_HET_valid, 'Type', 'Pearson');
% 
% text(5000, 0.8, ['p(HET)', num2str(p_HETAB)], 'FontSize', 10)
% text(5000, 0.7, ['p(WT)', num2str(p_WTAB)], 'FontSize', 10)
% 
% subplot(1,2,2)
% scatter(CD_trials(find(strcmp(genotype, 'WT'))), ...
%     nanmean(perf_DC_session(find(strcmp(genotype, 'WT')), 5:6),2),150, 'b', 'filled');
% hold on;
% scatter(CD_trials(find(strcmp(genotype, hetGeno))), ...
%     nanmean(perf_DC_session(find(strcmp(genotype, hetGeno)), 5:6),2),150, 'r', 'filled');
% title('CD trials experienced')
% 
% WT_CD = CD_trials(find(strcmp(genotype, 'WT')));
% WT_CD_valid = WT_CD(~isnan(DC_perf_WT));
% [r_WTCD, p_WTCD] = corr(WT_CD_valid, ...
%     DC_perf_WT_valid, 'Type', 'Pearson');
% HET_CD = CD_trials(find(strcmp(genotype, hetGeno)));
% HET_CD_valid = HET_CD(~isnan(DC_perf_HET));
% [r_HETCD, p_HETCD] = corr(HET_CD_valid, ...
%     DC_perf_HET_valid, 'Type', 'Pearson');
% 
% text(5000, 0.9, ['p(HET)', num2str(p_HETCD)], 'FontSize', 10)
% text(5000, 0.8, ['p(WT)', num2str(p_WTCD)], 'FontSize', 10)
% 
% print(gcf,'-dpng',fullfile(savefigpath,['Scatter plot Reversal based on CD trials experienced']));    %png format
% saveas(gcf, fullfile(savefigpath, ['Scatter plot Reversal based on CD trials experienced']), 'fig');
% %savefig(fullfile(savefigpath, 'Number of trials performed.fig'));
% saveas(gcf, fullfile(savefigpath, ['Scatter plot Reversal based on CD trials experienced']),'svg');
% close()

%% trial timing and performance
for ss=1:nFiles
   csvFile = fullfile(dataIndex.BehPath{ss}, ...
                dataIndex.BehCSV{ss});
   if exist(csvFile)
   results = readtable(csvFile);
   for tt = timing
       startSeconds = (tt-1)*60;
       endSeconds = tt*60-1;
       sessionTiming = results.center_in-results.center_in(1);
       %if endSeconds <= sessionTiming(end)
       trialMask=(sessionTiming>=startSeconds & sessionTiming<=endSeconds);
       %elseif endSeconds > sessionTiming(end) && startSeconds <=sessionTiming(end)
       %     trialMask=(sessionTiming>=startSeconds);
       %else

       Trials_time(ss,tt) = sum(trialMask);
   end

   for pp = 1:length(perf_timing)
       startSeconds = (pp-1)*30*60;
       endSeconds = pp*30*60-1;
       sessionTiming = results.center_in-results.center_in(1);
       trialMask=(sessionTiming>=startSeconds & sessionTiming<=endSeconds);
       performance_time(ss,pp) = sum(~isnan(results.reward(trialMask)))/sum(trialMask);
   end
   else
       Trials_time(ss,:) = NaN;
       performance_time(ss,:) = NaN;
   end

end

%% make the plot

setup_figprop;

%% percentage of total trials performed by time
Trials_time_percentage = cumsum(Trials_time,2)./sum(Trials_time,2);
%plot_Trial_distrubution(Trials_time, performance_time, savefigpath);
% cap at 60 minues
%% number of trials performed

plot_nTrials(nAB_trials, nCD_trials, nDC_trials,genotype, savefigpath);

%% how trials build up over time
%plot_trials_time()
%%  performance in quantile
plot_performance_quantile(perf_AB_quantile, strain, genotype, savefigpath, 'AB');
plot_performance_quantile(perf_CD_quantile, strain, genotype, savefigpath, 'CD');
plot_performance_quantile(perf_DC_quantile, strain, genotype, savefigpath, 'DC');

%% performance in block on session-basis
plot_performance_block_session(perf_AB_block_session, strain, genotype, savefigpath, 'AB');
plot_performance_block_session(perf_CD_block_session, strain, genotype, savefigpath, 'CD');
plot_performance_block_session(perf_DC_block_session, strain, genotype, savefigpath, 'DC');

%% performance in block
% plot_performance_block(perf_AB_block, blockLength, genotype, savefigpath, 'AB');
% plot_performance_block(perf_CD_block, blockLength,  genotype, savefigpath, 'CD');
% plot_performance_block(perf_DC_block, blockLength, genotype, savefigpath, 'DC');

%% for DC reversal , separately examine the animals that reversed and not reversed
% criteria: reached more than 60% in the last 2 sessions

DC_perf_session= squeeze(nanmean(perf_DC_quantile(:,:,5:6),2));
animal_reversed = find(any(DC_perf_session > 0.6, 2));
animal_notreversed = find(~any(DC_perf_session > 0.6, 2));
% plot
plot_performance_quantile(perf_AB_quantile(animal_reversed ,:,:), strain,...
    genotype(animal_reversed ,:,:), savefigpath, 'AB-reversed');
plot_performance_quantile(perf_CD_quantile(animal_reversed ,:,:), strain,...
    genotype(animal_reversed ,:,:), savefigpath, 'CD-reversed');
plot_performance_quantile(perf_DC_quantile(animal_reversed ,:,:), strain, ...
    genotype(animal_reversed ,:,:), savefigpath, 'DC-reversed');

plot_performance_block_session(perf_AB_block_session(animal_reversed ,:,:), strain,...
    genotype(animal_reversed ,:,:), savefigpath, 'AB-reversed');
plot_performance_block_session(perf_CD_block_session(animal_reversed ,:,:), strain,...
    genotype(animal_reversed ,:,:), savefigpath, 'CD-reversed');
plot_performance_block_session(perf_DC_block_session(animal_reversed ,:,:),strain, ...
    genotype(animal_reversed ,:,:), savefigpath, 'DC-reversed');

% not
% reversed---------------------------------------------------------------
plot_performance_quantile(perf_AB_quantile(animal_notreversed ,:,:),strain, ...
    genotype(animal_notreversed ,:,:), savefigpath, 'AB-notreversed');
plot_performance_quantile(perf_CD_quantile(animal_notreversed ,:,:), strain,...
    genotype(animal_notreversed ,:,:), savefigpath, 'CD-notreversed');
plot_performance_quantile(perf_DC_quantile(animal_notreversed ,:,:), strain,...
    genotype(animal_notreversed ,:,:), savefigpath, 'DC-notreversed');

plot_performance_block_session(perf_AB_block_session(animal_notreversed ,:,:), strain,...
    genotype(animal_notreversed ,:,:), savefigpath, 'AB-notreversed');
plot_performance_block_session(perf_CD_block_session(animal_notreversed ,:,:), strain,...
    genotype(animal_notreversed ,:,:), savefigpath, 'CD-notreversed');
plot_performance_block_session(perf_DC_block_session(animal_notreversed ,:,:), strain,...
    genotype(animal_notreversed ,:,:), savefigpath, 'DC-notreversed');

% compare group performance based on AUC score
plot_performance_block_session_reverse(perf_AB_block_session, strain,...
    genotype, savefigpath, 'AB-reversal', animal_reversed, animal_notreversed);
plot_performance_block_session_reverse(perf_CD_block_session, strain,...
    genotype, savefigpath, 'CD-reversal', animal_reversed, animal_notreversed);
plot_performance_block_session_reverse(perf_DC_block_session, strain,...
    genotype, savefigpath, 'DC-reversal', animal_reversed, animal_notreversed);


% calculate and save AUC here
perf_AUC_AB = squeeze(nanmean((perf_AB_block_session(:,1:10,:)-0.5), 2));
perf_AUC_CD = squeeze(nanmean((perf_CD_block_session(:,1:10,:)-0.5), 2));
perf_AUC_DC = squeeze(nanmean((perf_DC_block_session(:,1:10,:)-0.5), 2));


% correlating reversal to AB performance (separate in sessions?)
% scatter plot
perf_correlation(perf_AUC_AB, perf_AUC_CD, perf_AUC_DC, genotype, savefigpath);

%% save the data
save(fullfile(savedatapath, 'perf_quantile.mat'), 'perf_AB_quantile', ...
    'perf_CD_quantile', 'perf_DC_quantile', 'genotype', 'Subjects');
save(fullfile(savedatapath, 'perf_block.mat'), 'perf_AB_block', ...
    'perf_CD_block', 'perf_DC_block', 'genotype', 'Subjects','blockLength');
save(fullfile(savedatapath, 'perf_block_session.mat'), 'perf_AB_block_session', ...
    'perf_CD_block_session', 'perf_DC_block_session', 'genotype', 'Subjects','blockLength');

% save all data in csv
savecsvpath = fullfile(savedatapath, 'perf_AUC_AB.csv');
save_data_to_csv(perf_AUC_AB, genotype, Subjects,  {'Session1', 'Session2', 'Session3'}, savecsvpath);
savecsvpath = fullfile(savedatapath, 'perf_AUC_CD.csv');
save_data_to_csv(perf_AUC_CD, genotype, Subjects,  {'Session1', 'Session2', 'Session3'}, savecsvpath);
savecsvpath = fullfile(savedatapath, 'perf_AUC_DC.csv');
save_data_to_csv(perf_AUC_DC, genotype, Subjects,  {'Session1', 'Session2', 'Session3', 'Session4', 'Session5', 'Session6'}, savecsvpath);

% save AUC?
