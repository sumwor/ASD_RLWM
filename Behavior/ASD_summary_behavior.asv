function ASD_summary_behavior(dataIndex, savefigpath)

% plot behavior summary 

nFiles = size(dataIndex,1);
animalList = unique(dataIndex.Animal);
nAnimals = length(animalList);
Genotypes = unique(dataIndex.Genotype);

%% plot number of retrain sessions in different phases
% 1. after three AB sessions, how many AB retrain sessions before reach
% first AB-CD session
% 1.1 how many retrain sessions between AB-CD 1 and AB-CD3
% 2. after three AB-CD sessions, how many retrain sessions before reach
% AB-CD-DC sessions
%.3. how many retrain sessions between AB-CD-DC and AB-DC 6
retrain = struct;
retrain_phase = {'AB', 'ABCD','ABCDDC' 'ABDC'};
for ii =1:length(retrain_phase)
    retrain.(retrain_phase{ii}) = zeros(nAnimals,1);
end
phaseMask = {}; % assign phase Mask to every session

genoList = {};
for ff = 1:nAnimals
    tempAB = nan; tempABCD = nan; tempABDC = nan;
    tempMask = {};
    subIndex = dataIndex(strcmp(dataIndex.Animal, animalList{ff}),:);
    genoList = [genoList; subIndex.Genotype(1)];
    nSubFiles = size(subIndex,1);
    % the first three sessions is always AB
    ABCDStart = find(strcmp(subIndex.Protocol,'AB-CD'),1);
    ABCDDC = find(strcmp(subIndex.Protocol,'AB-CD-DC'),1);
    if ~isempty(ABCDStart)
        tempAB = ABCDStart-4; % number of extra AB sessions before reach AB-CD
        if tempAB < 0
            tempAB = 0;
        end
    end
    if ~isempty(ABCDDC)
        tempABDC = 0;
        tempABCDDC = sum(strcmp(subIndex.Protocol(ABCDStart:ABCDDC), 'AB')); % number of AB retrain sessions before reaching ABCDDC
        tempABCD = ABCDDC - ABCDStart -1-tempABCDDC; % number of extra CD learning sessions before ABCDDC
    end
    for ss = ABCDDC+1:nSubFiles
        if ~strcmp(subIndex.Protocol{ss}, 'AB-DC')
            tempABDC = tempABDC + 1;
        end
    end
    retrain.AB(ff) = tempAB;
    retrain.ABCD(ff) = tempABCD;
    retrain.ABCDDC(ff) = tempABCDDC;
    retrain.ABDC(ff) = tempABDC;
end

% group the data for plot
data = [];
group_geno = {};
group_retrain = {};
for aa = 1:nAnimals
    for m = 1:length(retrain_phase)
        data = [data; retrain.(retrain_phase{m})(aa)];
        group_geno = [group_geno; genoList{aa}];
        group_retrain = [group_retrain; retrain_phase{m}];
    end
end

g1 = categorical(group_geno, Genotypes);
g2 = categorical(group_retrain, retrain_phase);

figure;
boxplot(data, {g2, g1}, 'notch', 'on' ,'FactorSeparator', 1, 'Colors', 'k', 'Widths', 0.6);

% Customize colors
colors = lines(3); % Define colors for WT, HET, KO
hBoxes = flip(findobj(gca, 'Tag', 'Box'));  % flip to match plotting order
nGroups = length(retrain_phase);
nGenotypes = length(Genotypes);

% Store one patch handle per genotype for legend
legendHandles = gobjects(nGenotypes,1);

for i = 1:length(hBoxes)
    genotypeIdx = mod(i-1, nGenotypes) + 1;  % 1-based index: 1,2,3
    patchHandle = patch(get(hBoxes(i), 'XData'), get(hBoxes(i), 'YData'), ...
                        colors(genotypeIdx,:), 'FaceAlpha', 0.5, ...
                        'EdgeColor', colors(genotypeIdx,:));
    
    % Save the first occurrence of each genotype's patch
    if i <= nGenotypes
        legendHandles(genotypeIdx) = patchHandle;
    end
end

% Set x-axis ticks at center of triplets
xticks = (1.5:3:12);  % every 3 positions, centered
xticklabels(retrain_phase);
set(gca, 'XTick', xticks, 'XTickLabel',retrain_phase, 'FontSize', 12);

xlabel('Training phase');
ylabel('Number of sessions');
title('Number of retrain sessions in different training phases');

box off

% ✅ Correct legend order
legendHandle = legend(legendHandles, Genotypes, 'Location', 'Best');
set(legendHandle, 'Box', 'off', 'Color', 'none');

hold on;
nGenotypes = length(Genotypes);
nMeasurements = length(retrain_phase);
groupSpacing = 1;  % spacing between each genotype within a measurement group
baseX = [];

for m = 1:nMeasurements
    for g = 1:nGenotypes
        % Compute the x-position for the current group
        groupIndex = (m - 1) * nGenotypes + g;
        x = groupIndex * groupSpacing;

        % Extract relevant data points
        mask = g1 == Genotypes{g} & g2 == retrain_phase{m};
        yVals = data(mask);

        % Jitter around x to avoid overlap
        jitteredX = x + 0.15 * (rand(size(yVals)) - 0.5);

        % Plot
        scatter(jitteredX, yVals, 20, colors(g, :), 'filled', ...
                'MarkerFaceAlpha', 0.6, 'MarkerEdgeAlpha', 0.3,'HandleVisibility','off');
    end
end

print(gcf,'-dpng',fullfile(savefigpath, 'Number of retrain sessions in different training phases'));    %png format
saveas(gcf, fullfile(savefigpath, 'Number of retrain sessions in different training phases'));

%% load response time, lick rate, and error rate from single-session file

protocol = {'AB', 'CD', 'DC'};
trialTypes = {'leftCorrect', 'leftIncorrect', 'rightCorrect', 'rightIncorrect'};
for ii = 1:nFiles
    BehPath = dataIndex.BehPath{ii};
    session = dataIndex.Session{ii};
    savedatapath = fullfile(BehPath, 'behAnalysis',num2str(session),'beh_analysis.mat');
    data = load(savedatapath);

    % load response time, intertrial interval, and error rate
    if ii == 1
        % initiate the variable
        rt = struct; iti = struct; error_rate = struct; error_reinforce = struct;
        rtLength = length(data.responseTime.AB.leftCorrect);
        itiLength = length(data.itiTime.AB.leftCorrect);
        error_size = size(data.error_rate.AB);
        for pp = 1:length(protocol)
            rt.(protocol{pp})=struct;
            iti.(protocol{pp}) = struct;
            rt.(protocol{pp}) = struct;
            error_rate.(protocol{pp}) = nan(error_size(1), error_size(2), nFiles);
            error_reinforce.(protocol{pp}) = nan(error_size(1), error_size(2), nFiles);
            for tt = 1:length(trialTypes)
                rt.(protocol{pp}).(trialTypes{tt}) = nan(rtLength, nFiles);
                iti.(protocol{pp}).(trialTypes{tt}) = nan(itiLength, nFiles);
            end
        end
    end
    
    % load the data
    for pp = 1:length(protocol)
        if isfield(data.responseTime, protocol{pp})
            for tt = 1:length(trialTypes)
                rt.(protocol{pp}).(trialTypes{tt})(:, ii) = data.responseTime.(protocol{pp}).(trialTypes{tt});
                iti.(protocol{pp}).(trialTypes{tt})(:, ii) = data.itiTime.(protocol{pp}).(trialTypes{tt});
            end
            error_rate.(protocol{pp})(:,:,ii) = data.error_rate.(protocol{pp});
            error_reinforce.(protocol{pp})(:,:,ii) = data.error_rate_reinforce.(protocol{pp});
        end
    end
    logreg_RCUC_array{ii} = data.lregRCUC_output;
    logreg_CRInt_array{ii} = data.lregCRInt_output;
end

% plot average response curve by genotype
rt_edge = [-0.5:0.05:10];
iti_edge = [0:0.5:30];




%% response time distribution

plot_responseTimes_summary(rt,dataIndex, rt_edge, 'AB', savefigpath);
plot_responseTimes_summary(rt,dataIndex, rt_edge, 'AB-CD-AB', savefigpath);
plot_responseTimes_summary(rt,dataIndex, rt_edge, 'AB-CD', savefigpath);
plot_responseTimes_summary(rt,dataIndex, rt_edge, 'AB-DC-AB', savefigpath);
plot_responseTimes_summary(rt,dataIndex, rt_edge, 'AB-DC', savefigpath);

close('all')

%% intertrial interval distributino
plot_itiTimes_summary(iti,dataIndex, iti_edge, 'AB', savefigpath);
plot_itiTimes_summary(iti,dataIndex, iti_edge, 'AB-CD-AB', savefigpath);
plot_itiTimes_summary(iti,dataIndex, iti_edge, 'AB-CD', savefigpath);
plot_itiTimes_summary(iti,dataIndex, iti_edge, 'AB-DC-AB', savefigpath);
plot_itiTimes_summary(iti,dataIndex, iti_edge, 'AB-DC', savefigpath);

%% error rate based on previous choice and reward history
err_edge = 0:6;
plot_error_summary(error_rate, 'choice', dataIndex, err_edge, 'AB', savefigpath)
plot_error_summary(error_rate, 'choice', dataIndex, err_edge, 'AB-CD-AB', savefigpath)
plot_error_summary(error_rate, 'choice', dataIndex, err_edge, 'AB-CD', savefigpath)
plot_error_summary(error_rate, 'choice', dataIndex, err_edge, 'AB-DC-AB', savefigpath)
plot_error_summary(error_rate, 'choice', dataIndex, err_edge, 'AB-DC', savefigpath)

plot_error_summary(error_reinforce, 'rewarded choice', dataIndex, err_edge, 'AB', savefigpath)
plot_error_summary(error_reinforce, 'rewarded choice', dataIndex, err_edge, 'AB-CD-AB', savefigpath)
plot_error_summary(error_reinforce, 'rewarded choice', dataIndex, err_edge, 'AB-CD', savefigpath)
plot_error_summary(error_reinforce, 'rewarded choice', dataIndex, err_edge, 'AB-DC-AB', savefigpath)
plot_error_summary(error_reinforce, 'rewarded choice', dataIndex, err_edge, 'AB-DC', savefigpath)

%% plot logistic regression result
plot_logreg_summary(logreg_RCUC_array, 'RCUC', dataIndex, 'AB', savefigpath)
plot_logreg_summary(logreg_CRInt_array, 'CRInt', dataIndex, 'AB', savefigpath)
plot_logreg_summary(logreg_RCUC_array, 'RCUC', dataIndex, 'AB-CD-AB', savefigpath)
plot_logreg_summary(logreg_CRInt_array, 'CRInt', dataIndex, 'AB-CD-AB', savefigpath)
plot_logreg_summary(logreg_RCUC_array, 'RCUC', dataIndex, 'AB-CD', savefigpath)
plot_logreg_summary(logreg_CRInt_array, 'CRInt', dataIndex, 'AB-CD', savefigpath)
plot_logreg_summary(logreg_RCUC_array, 'RCUC', dataIndex, 'AB-DC-AB', savefigpath)
plot_logreg_summary(logreg_CRInt_array, 'CRInt', dataIndex, 'AB-DC-AB', savefigpath)
plot_logreg_summary(logreg_RCUC_array, 'RCUC', dataIndex, 'AB-DC', savefigpath)
plot_logreg_summary(logreg_CRInt_array, 'CRInt', dataIndex, 'AB-DC', savefigpath)

%% mixed effect logistic regression with session ID being a factor


end
