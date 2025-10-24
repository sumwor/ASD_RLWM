function estimate_Latent_fit_psy(psymodelName, dataIndex, label, savefigpath)
setup_figprop;
% used fitted latent variables to estimate the result

psy_fit = readtable(psymodelName);
if strcmp(label,'AB-AB1')
    protocol = 'AB';
    protocolDay = 1;
elseif strcmp(label, 'AB-AB2')
    protocol = 'AB';
    protocolDay = 2;
elseif strcmp(label, 'AB-AB3')
    protocol = 'AB';
    protocolDay = 3;
elseif strcmp(label, 'AB-CD-CD1')
    protocol = 'AB-CD';
    protocolDay = 1;
    elseif strcmp(label, 'AB-CD-CD2')
    protocol = 'AB-CD';
    protocolDay = 2;
    elseif strcmp(label, 'AB-CD-CD3')
    protocol = 'AB-CD';
    protocolDay = 3;
end

subjects = psy_fit.Animal;
nSubject = length(subjects);
weighted_sum_step = -9.95:0.1:9.95;
weighted_sum_count = zeros(length(weighted_sum_step),nSubject);
n_right_psychometric = zeros(length(weighted_sum_step),nSubject);

genotypes = cell(nSubject,1);
rt_all = cell(nSubject,1);
decision_time_all = cell(nSubject,1);
sample_time_all = cell(nSubject,1);
iti_all = cell(nSubject,1);
relative_weight = cell(nSubject,1); 
abs_weight = cell(nSubject,1);

for k = 1:length(subjects)
    % load csv file for response time and intertrial interval
    tempgeno = dataIndex.Genotype(strcmp(dataIndex.Animal,num2str(subjects(k))));
    genotypes{k} = tempgeno{1};
    csvfilemask = strcmp(dataIndex.Animal, num2str(subjects(k))) & strcmp(dataIndex.Protocol,protocol) & cell2mat(dataIndex.ProtocolDay) == protocolDay;
    csvfile = dataIndex.BehCSV{csvfilemask};
    csv = readtable(fullfile(dataIndex.BehPath{csvfilemask}, csvfile));
    
    if contains(label, 'CD')
        csv.schedule = csv.schedule-2;
        csv = csv(csv.schedule>0,:);
    end
    analysis = dataIndex.BehPath(strcmp(dataIndex.Animal,num2str(subjects(k))));
    modelpath = fullfile(analysis{1},'latent',['psy_fit_',label,'.json']);

    txt = fileread(modelpath);
    psy_latent = jsondecode(txt);

    % estimate p_correct based on wMode

    reward = psy_latent.args.dat.correct;
        nTrials = length(reward);
    choice = psy_latent.args.dat.y;
    sti = psy_latent.args.dat.inputs.cBoth;
    stick = psy_latent.args.dat.inputs.stick;
    fitted_weight = psy_latent.wMode;
    x = [ones(length(sti),1), sti, stick];
    
    weighted_sum = sum(x'.*fitted_weight);
    pR = 1 ./ (1 + exp(-(sum(x' .* fitted_weight))));

    % convert to p_correct
    rightMask = sti<0;
    leftMask = sti>0;
    pCorrect = zeros(nTrials,1);
    pCorrect(rightMask) = pR(rightMask);
    pCorrect(leftMask) = 1-pR(leftMask);
    pCorrect_smooth = movmean(pCorrect, 60);
    pR_smooth = movmean(pR,60);
    % calculate average reward rate
    nTrials = length(sti);
    ave_reward_rate = zeros(nTrials,1);
    ave_pR = zeros(nTrials,1);
    windowSize = 60;
    for tt = 1:nTrials
        if tt<=windowSize/2
            ave_reward_rate(tt,1) = nanmean(reward(1:tt));
            ave_pR(tt,1) = sum(choice(1:tt)==1)/tt;
        elseif tt+windowSize/2>=nTrials
            ave_reward_rate(tt,1) = nanmean(reward(tt:end));
            ave_pR(tt,1) = sum(choice(tt:end)==1)/(nTrials-tt+1);
        else
            ave_reward_rate(tt,1) = nanmean(reward(tt-windowSize/2:tt+windowSize/2));
            ave_pR(tt,1) = sum(choice(tt-windowSize/2:tt+windowSize/2)==1)/windowSize;
        end
    end


    nPlot = 1:nTrials;
    weight_var = {'Bias', 'Stimulus', 'Stick'};
    figure;
    sgtitle('Psy model fit')
    subplot(3,1,1)
    nPlot = 1:nTrials;
    plot(nPlot,ave_reward_rate,'k')
    hold on;
    %plot(nPlot, fitData.policy, '-','Color', [0.7 0.7 0.7]);
    plot(nPlot, pCorrect_smooth, '--','Color', [0.7 0.7 0.7]);
    %legend('PCorrect', 'fitted PCorrect')
    set(gca, 'Box', 'off');
    ylabel('P_{correct}')
    ylim([0 1])

    subplot(3,1,2)
    colors = lines(length(weight_var));
    for ww = 1:length(weight_var)
        plot(nPlot, fitted_weight(ww,:),'Color', colors(ww,:))
        hold on;
    end
    set(gca, 'Box', 'off');
    lgd=legend(weight_var)
    set(lgd, 'Box', 'off', 'Color', 'none')
    ylabel('Weight');
    %ylim([-0.5 0.5])


    print(gcf,'-dpng',fullfile(analysis{1},'latent', ['psy_session-fit_', label]));    %png format
    saveas(gcf, fullfile(analysis{1},'latent', ['psy_session-fit_', label]), 'fig');
    saveas(gcf, fullfile(analysis{1}, 'latent',['psy_session-fit_', label]),'svg');
    
    close;
    %% psychometric curve
    for ww = 1:length(weighted_sum_step)
        weighted_sum_Mask = weighted_sum>=weighted_sum_step(ww)-0.05 & weighted_sum<weighted_sum_step(ww)+0.05;

        weighted_sum_count(ww,k) = sum(weighted_sum_Mask);

        n_right_psychometric(ww,k) = sum(choice(weighted_sum_Mask)==1);  % number of right choice trials
    end

    %% response time and iti
    trialsMask = ~isnan(csv.actions) & ~isnan(csv.schedule);
    sample_time= csv.center_out-csv.center_in;
    rt = csv.side_in-csv.center_in;
    decision_time = csv.side_in - csv.center_out;
    iti = [nan;csv.center_in(2:end) - csv.center_in(1:end-1)];
    rt_all{k} = rt(trialsMask);
    iti_all{k} = iti(trialsMask);
    sample_time_all{k} = sample_time(trialsMask);
    decision_time_all{k} = decision_time(trialsMask);

    relative_weight{k} = abs(fitted_weight)./sum(abs(fitted_weight),1);
    abs_weight{k} = abs(fitted_weight);
    % iti of trial n is the iti before trial n odor on 
end

%% make summary plot
All_geno = unique(genotypes);

if sum(contains(All_geno, 'KO'))>0
    mutGene = 'KO';
elseif sum(contains(All_geno,'HET'))>0
    mutGene = 'HET';
elseif sum(contains(All_geno,'HEM'))>0
    mutGene = 'HEM';
end

WTMask = strcmp(genotypes,'WT');

HETMask = strcmp(genotypes, mutGene);

predicted_curve = 1./(1+exp(-weighted_sum_step));
%1. psychometric curve
figure;

subplot(2,2,1)

yyaxis left

bar(weighted_sum_step,nanmean(weighted_sum_count(:,WTMask)./sum(weighted_sum_count(:,WTMask),1),2));
ylabel('Precentage of trials')

hold on;
yyaxis right
%ylabel('Prob(R)')
pR = n_right_psychometric(:,WTMask)./weighted_sum_count(:,WTMask);
meanPR = nanmean(pR,2);
stePR = nanstd(pR,0,2)/sqrt(sum(WTMask));
errorbar(weighted_sum_step, meanPR, stePR, 'o', ...
    'MarkerFaceColor','k','MarkerEdgeColor','k', ...
    'Color','k','LineStyle','none','CapSize',6);
plot(weighted_sum_step, predicted_curve)
xlabel(['\bf\it x \cdot w'], 'Interpreter', 'tex')
set(gca,'box', 'off')
title('WT')

subplot(2,2,2)

yyaxis left
%ylabel('Precentage of trials')
bar(weighted_sum_step,nanmean(weighted_sum_count(:,HETMask)./sum(weighted_sum_count(:,HETMask),1),2));
hold on;
yyaxis right
ylabel('Prob(R)')
pR = n_right_psychometric(:,HETMask)./weighted_sum_count(:,HETMask);
meanPR = nanmean(pR,2);
stePR = nanstd(pR,0,2)/sqrt(sum(HETMask));
errorbar(weighted_sum_step, meanPR, stePR, 'o', ...
    'MarkerFaceColor','k','MarkerEdgeColor','k', ...
    'Color','k','LineStyle','none','CapSize',6);
plot(weighted_sum_step, predicted_curve)
set(gca,'box', 'off')
title(mutGene)
xlabel(['\bf\it x \cdot w'], 'Interpreter', 'tex')
    print(gcf,'-dpng',fullfile(savefigpath,'latent', ['psy_psychometric_', label]));    %png format
    saveas(gcf, fullfile(savefigpath,'latent', ['psy_psychometric_', label]), 'fig');
    saveas(gcf, fullfile(savefigpath, 'latent',['psy_psychometric_', label]),'svg');

%% plot relative weight
% average weight from the last 100 trials from previous session
% and average weight from the first 100 trials from next session
%relative_weight
%% 2. response time correlation with relative weight of stimulus
% multiple linear regression
%(bias/stickiness)
allRT = [];
allST = []; % sample time (center_out-center_in)
allDT = []; % decision time (side_in - center_out)
allITI = [];
allBias = [];
allStim = [];
allStick = [];
allGenotype = [];
allAnimal = [];

% cut off largest 5% of RT in every session, calculate the mean RT per
% session
meanRT_session = nan(numel(rt_all),1);
for rr = 1:numel(rt_all)
    tempRT = rt_all{rr};
    cutRT = prctile(tempRT,0.95);
    RTIncluded = tempRT(tempRT<cutRT);
    meanRT_session(rr) = nanmean(RTIncluded);
end



for iAnimal = 1:numel(rt_all)
    rt = rt_all{iAnimal}(:);  
    st = sample_time_all{iAnimal}(:);
    dt = decision_time_all{iAnimal}(:);
    iti = iti_all{iAnimal}(:);
    
    w = abs_weight{iAnimal};  % 3 x Ntrials
    nTrials = numel(rt);
    bias  = w(1,1:nTrials)';
    stim  = w(2,1:nTrials)';
    stick = w(3,1:nTrials)';
    
    % --- Keep only trials where all timings are >= 0 ---
    validIdx = (rt >= 0) & (st >= 0) & (dt >= 0) & (iti >= 0);
    
    % Append only valid trials
    allRT       = [allRT; rt(validIdx)];
    allST       = [allST; st(validIdx)];
    allDT       = [allDT; dt(validIdx)];
    allITI      = [allITI; iti(validIdx)];
    allBias     = [allBias; bias(validIdx)];
    allStim     = [allStim; stim(validIdx)];
    allStick    = [allStick; stick(validIdx)];
    allGenotype = [allGenotype; repmat(genotypes(iAnimal), sum(validIdx), 1)];
    allAnimal   = [allAnimal; repmat(iAnimal, sum(validIdx), 1)];
end
% Build table
[pathstr, name, ext] = fileparts(psymodelName);
savedatapath = pathstr;

% save the average RT and corresponding genotype
RTData.meanRT = meanRT_session;
RTData.geno = genotypes;
save(fullfile(savedatapath, ['meanRT', label, '.mat']), 'RTData');
 

% cutoff highest 5% of data
cutoff = prctile(allRT, 95);
normalRTMask = allRT<cutoff;
tbl = table(allRT(normalRTMask), allBias(normalRTMask), allStim(normalRTMask), allStick(normalRTMask), categorical(allGenotype(normalRTMask)), categorical(allAnimal(normalRTMask)), ...
    'VariableNames', {'ResponseTime','Bias','Stimulus','Stick','Genotype','Animal'});
tbl.Genotype = reordercats(tbl.Genotype, {'WT',mutGene});

linear_label = ['response time_', label];
mixed_linear_psyweight(tbl, linear_label, savefigpath, savedatapath, mutGene)

% cutoff lowest 2.5% and highest 2.5%
cutlow = prctile(allST,2.5);
cuthigh = prctile(allST, 97.5);
normalSTMask = allST>cutlow & allST<cuthigh;
tbl = table(allST(normalSTMask), allBias(normalSTMask), allStim(normalSTMask), allStick(normalSTMask), categorical(allGenotype(normalSTMask)), categorical(allAnimal(normalSTMask)), ...
    'VariableNames', {'ResponseTime','Bias','Stimulus','Stick','Genotype','Animal'});
tbl.Genotype = reordercats(tbl.Genotype, {'WT', mutGene});
linear_label = ['sample time_', label];
mixed_linear_psyweight(tbl, linear_label, savefigpath, savedatapath,mutGene)


cutoff = prctile(allDT, 95);
normalDTMask = allDT<cutoff;
tbl = table(allDT(normalDTMask), allBias(normalDTMask), allStim(normalDTMask), allStick(normalDTMask), categorical(allGenotype(normalDTMask)), categorical(allAnimal(normalDTMask)), ...
    'VariableNames', {'ResponseTime','Bias','Stimulus','Stick','Genotype','Animal'});
tbl.Genotype = reordercats(tbl.Genotype, {'WT', mutGene});
linear_label = ['decision time_', label];
mixed_linear_psyweight(tbl, linear_label, savefigpath, savedatapath,mutGene)

% remove too large ITIs
cutoff = prctile(allITI, 95);
normalITIMask = allITI<cutoff;
tbl = table(allITI(normalITIMask), allBias(normalITIMask), allStim(normalITIMask), allStick(normalITIMask), categorical(allGenotype(normalITIMask)), categorical(allAnimal(normalITIMask)), ...
    'VariableNames', {'ResponseTime','Bias','Stimulus','Stick','Genotype','Animal'});
tbl.Genotype = reordercats(tbl.Genotype, {'WT', mutGene});
linear_label = ['ITI_', label];
mixed_linear_psyweight(tbl, linear_label, savefigpath,savedatapath, mutGene)


%% fitting GLM on individual animals
% 
% unique animals
% 
% nAnimals = numel(subjects);
% 
% coefNames = {'intercept','bias','stim','stick'};
% dependentVars = {'rt', 'st', 'dt', 'iti'};
% 
% for dd =1:length(dependentVars)
%     dependent = dependentVars{dd};
%     nCoef = numel(coefNames);
% 
%     allGLMResults = struct;
%     allGLMResults.animals = subjects;
%     allGLMResults.coefNames = coefNames;
%     allGLMResults.coefficients = nan(nAnimals, nCoef);
%     allGLMResults.pvalues = nan(nAnimals, nCoef);
%     allGLMResults.genotype = cell(nAnimals,1);
% 
%     for iAnimal = 1:nAnimals
%         thisAnimal = subjects(iAnimal);
%         switch dependent
%             case 'rt'
%                 y = rt_all{iAnimal};
%             case 'st'
%                 y = sample_time_all{iAnimal};
%             case 'dt'
%                 y = decision_time_all{iAnimal};
%             case 'iti'
%                 y = iti_all{iAnimal};
%         end
% 
%         w = abs_weight{iAnimal};
%         nTrials = numel(y);
% 
%         bias = w(1,1:nTrials)';
%         stim = w(2,1:nTrials)';
%         stick = w(3,1:nTrials)';
% 
%         % remove possible negative values
%         % filter tiny values
%         bias(bias<1e-3) = 1e-3;
%         stim(stim<1e-3) = 1e-3;
%         stick(stick<1e-3) = 1e-3;
%         cleanMask = y>0 & ~isnan(y);
%         bias = bias(cleanMask,:);
%         stim = stim(cleanMask,:);
%         stick = stick(cleanMask,:);
%         y = y(cleanMask>0);
% 
%         % fit GLM (normal distribution, identity link since RT is continuous)
%         tbl = table(y,bias, stim, stick);
% 
%         % Fit GLM (gamma distribution for RTs)
%         mdl = fitglm(tbl, 'y ~ bias + stim + stick', ...
%             'Distribution','inverse gaussian','Link','log');
% 
%         % Save coefficients & p-values
%         allGLMResults.coefficients(iAnimal,:) = mdl.Coefficients.Estimate(1:end)'; % skip intercept
%         allGLMResults.pvalues(iAnimal,:) = mdl.Coefficients.pValue(1:end)';
%         allGLMResults.genotype(iAnimal) = genotypes(iAnimal);
% 
%         % Store the whole model too (optional, can take memory)
%         %allGLMResults.models{iAnimal} = mdl;
%     end
% 
%     % save the result
%     [pathstr, name, ext] = fileparts(psymodelName);
%     savedatapath = pathstr;
%     save(fullfile(savedatapath, [dependent,'_GLM_',label, '.mat']), 'allGLMResults');
% end