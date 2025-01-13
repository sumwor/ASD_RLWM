function ASD_odor_summary(dataIndex, savefigpath)

%% go over every session to plot performance in blocks for
% 1. first two sessions of AB
% 2. first two sessions of CD
% 3. DC reverse

nFiles = size(dataIndex,1);
Subjects = unique(dataIndex.Animal);
nSubjects = length(unique(dataIndex.Animal));

genotype = cell(nSubjects,1);

% performance
blockLength = 100; % check performance in blocks of 100 trials
perf_AB_block = nan(nSubjects, 50); 
perf_CD_block = nan(nSubjects, 50);
perf_DC_block = nan(nSubjects, 50);

perf_AB_quantile = zeros(nSubjects,4, 2); % dim3: first session, second session
perf_CD_quantile = zeros(nSubjects,4, 2); % dim3: session
perf_DC_quantile = zeros(nSubjects,4, 2); % dim4: session

for ii = 1:nSubjects
    
    genotype(ii) = unique(dataIndex.Genotype(strcmp(dataIndex.Animal, Subjects{ii})));
    %% look for first 2 AB sessions
    subDataIndex = dataIndex(strcmp(dataIndex.Animal, Subjects{ii}),:);
    % load the first two sessions
    result1 = readtable(fullfile(subDataIndex.BehPath{1}, ...
        subDataIndex.BehCSV{1}));
    result2 = readtable(fullfile(subDataIndex.BehPath{2}, ...
        subDataIndex.BehCSV{2}));
    
    % calculate performance in quantile
    perf_AB_quantile(ii,:,:) = perf_in_quantile(result1, result2);
    
    % calculate performance in blocks
    perf_AB_block(ii,:) = perf_in_block(result1, result2, blockLength);
    
    %% for CD odors
    animalMask = strcmp(dataIndex.Animal, Subjects{ii});
    stageMask = cellfun(@(x) isequal(x, [1;2;3;4]), dataIndex.OdorPresented);
    subDataIndex = dataIndex(animalMask & stageMask,:);

    result1 = readtable(fullfile(subDataIndex.BehPath{1}, ...
        subDataIndex.BehCSV{1}));
    result2 = readtable(fullfile(subDataIndex.BehPath{2}, ...
        subDataIndex.BehCSV{2}));
    
    % find the switch point
    switchTrial1 = find(result1.schedule == 3 | result1.schedule == 4, 1);
    switchTrial2 = find(result2.schedule == 3 | result2.schedule == 4, 1);
    
        % calculate performance in quantile
    perf_CD_quantile(ii,:,:) = perf_in_quantile(result1(switchTrial1:end,:), result2(switchTrial2:end,:));
    
    % calculate performance in blocks
    perf_CD_block(ii,:) = perf_in_block(result1(switchTrial1:end,:), result2(switchTrial2:end,:),blockLength);
    
    %% for DC performance

    animalMask = strcmp(dataIndex.Animal, Subjects{ii});
    stageMask = cellfun(@(x) isequal(x, [1;2;5;6]), dataIndex.OdorPresented);
    subDataIndex = dataIndex(animalMask & stageMask,:);
    
    if ~isempty(subDataIndex)
    result1 = readtable(fullfile(subDataIndex.BehPath{1}, ...
        subDataIndex.BehCSV{1}));
    result2 = readtable(fullfile(subDataIndex.BehPath{2}, ...
        subDataIndex.BehCSV{2}));
    
    % find the switch point
    switchTrial1 = find(result1.schedule == 5 | result1.schedule == 6, 1);
    switchTrial2 = find(result2.schedule == 5 | result2.schedule == 6, 1);
    
        % calculate performance in quantile
    perf_DC_quantile(ii,:,:) = perf_in_quantile(result1(switchTrial1:end,:), result2(switchTrial2:end,:));
    
    % calculate performance in blocks
    perf_DC_block(ii,:) = perf_in_block(result1(switchTrial1:end,:), result2(switchTrial2:end,:),blockLength);
    end
end

% make the plot

setup_figprop;



%%  performance in quantile
plot_performance_quantile(perf_AB_quantile, genotype, savefigpath, 'AB');
plot_performance_quantile(perf_CD_quantile, genotype, savefigpath, 'CD');
plot_performance_quantile(perf_DC_quantile, genotype, savefigpath, 'DC');

%% performance in block

x=1;