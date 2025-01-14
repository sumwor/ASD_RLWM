function ASD_odor_summary(dataIndex, savefigpath)

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

% performance
blockLength = 100; % check performance in blocks of 100 trials
perf_AB_block = nan(nSubjects, 50);
perf_CD_block = nan(nSubjects, 50);
perf_DC_block = nan(nSubjects, 50);

perf_AB_quantile = zeros(nSubjects,4, 3); % dim3: first session, second session
perf_CD_quantile = zeros(nSubjects,4, 3); % dim3: session
perf_DC_quantile = zeros(nSubjects,4, 6); % dim4: session
perf_DC_session = zeros(nSubjects,6);

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

            results{ss} = readtable(csvFile);
            nAB_trials(ii, ss) = size(results{ss},1);
        else
            results{ss} = NaN;
        end
    end
    % calculate performance in quantile
    perf_AB_quantile(ii,:,:) = perf_in_quantile(results, protocol);

    % calculate performance in blocks
    perf_AB_block(ii,:) = perf_in_block(results, protocol, blockLength);

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

    % calculate performance in blocks
    perf_CD_block(ii,:) = perf_in_block(results, protocol,blockLength);

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

        % calculate performance in blocks
        perf_DC_block(ii,:) = perf_in_block(results, protocol,blockLength);
    end
end

% make the plot

setup_figprop;

%% number of trials performed

plot_nTrials(nAB_trials, nCD_trials, nDC_trials,genotype, savefigpath);

%%  performance in quantile
plot_performance_quantile(perf_AB_quantile, genotype, savefigpath, 'AB');
plot_performance_quantile(perf_CD_quantile, genotype, savefigpath, 'CD');
plot_performance_quantile(perf_DC_quantile, genotype, savefigpath, 'DC');

%% performance in block
plot_performance_block(perf_AB_block, blockLength, genotype, savefigpath, 'AB');
plot_performance_block(perf_CD_block, blockLength,  genotype, savefigpath, 'CD');
plot_performance_block(perf_DC_block, blockLength, genotype, savefigpath, 'DC');


x=1;