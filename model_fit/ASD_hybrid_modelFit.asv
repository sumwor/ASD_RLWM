function ASD_hybrid_modelFit(dataIndex, files,savedatapath, savefigpath)

% fit the models for different sessions
odors = {'AB', 'CD', 'DC', 'AB-CD', 'AB-DC'};
protocols = {'AB', 'AB-CD','AB-DC', 'AB-CD', 'AB-DC'};
trialsList = {'nAB', 'nCD', 'nDC', 'nAB', 'nAB'};
nSessions = [3, 3, 5, 3, 5];

% test if out first
nFiles = length(files);

for ii = 1:6
    prepDataPath = files{ii};
    
    % for AB2, AB3, CD1, CD3, add a forgetting parameter between sessions 
    % for Q value

    fit_result = fit_hybrid_bias_models(prepDataPath, savedatapath,dataIndex);



    %prepData = load(prepDataPath);
    %% make individual session plots
    filename = regexp(prepDataPath, '(?<=data4model).*', 'match');
    label = filename{1}(1:end-4);
    modelNum = 1;
    plot_hybrid_fit(fit_result, savefigpath, label, modelNum)

    % model estimation based on fitted parameters
    % P_correct
    % psychometric curve
    estimate_Latent_fit(prepDataPath, fit_result, dataIndex, label, modelNum, savefigpath)
    
    psymodelName = fullfile(savedatapath, ['PsyFit',label, '.csv']);

    %rename_file(psymodelName,label,dataIndex)
    compare_hybrid_psymodels(fit_result, psymodelName,label, savefigpath)
    plot_psy_fit(psymodelName, dataIndex, savefigpath, label)
    estimate_Latent_fit_psy(psymodelName, dataIndex, label, savefigpath)
    
    % scatter plot of last weight from last session and new weight from new
    % session

    % put mixed linear result in different sessions together
end


for ii = 1:6
    prepDataPath = files{ii};
        filename = regexp(prepDataPath, '(?<=data4model).*', 'match');
    label = filename{1}(1:end-4);
 psymodelName = fullfile(savedatapath, ['PsyFit',label, '.csv']);
estimate_Latent_fit_psy(psymodelName, dataIndex, label, savefigpath)
end


% 3 coefficient + baseline for 6 sessions in WT/HET genotypes 
vars = {'response time', 'sample time', 'decision time', 'ITI'};
structVar = {'rt', 'st', 'dt', 'iti'};

% also make a plot of average response time 

% load AIC and average response time as well
for ii = 1:6
    prepDataPath = files{ii};
    filename = regexp(prepDataPath, '(?<=data4model).*', 'match');
    label = filename{1}(1:end-4);
    psymodelName = fullfile(savedatapath, ['PsyFit',label, '.csv']);
    fit_result = fit_hybrid_bias_models(prepDataPath, savedatapath,dataIndex);
    if ii == 1
        psy_fit = readtable(psymodelName);
        subjects = psy_fit.Animal;
        nSubject = length(subjects);
        genotypes = fit_result.genotypes;
        weight_session = nan(nSubject, 6, 2, 3);
        meanRT_session = nan(nSubject, 6);
        meanRT_geno = cell(nSubject, 6);
        GLMr = struct;
        for vv = 1:length(vars)
            currVar = vars{vv};
            GLMr.(structVar{vv}) = struct;
            GLMr.(structVar{vv}).coeff_WT = nan(length(vars),6);
            GLMr.(structVar{vv}).SE_WT = nan(length(vars), 6);
            GLMr.(structVar{vv}).coeff_Mut = nan(length(vars),6);
            GLMr.(structVar{vv}).SE_Mut = nan(length(vars), 6);
            GLMr.(structVar{vv}).p_WT = nan(length(vars), 6);
            GLMr.(structVar{vv}).p_Mut = nan(length(vars), 6);
            GLMr.(structVar{vv}).p_geno = nan(length(vars), 6);
        end


    end

    % load weight
    for ss = 1:nSubject
        analysis = dataIndex.BehPath(strcmp(dataIndex.Animal,num2str(subjects(ss))));
        modelpath = fullfile(analysis{1},'latent',['psy_fit_',label,'.json']);
        if exist(modelpath)
            txt = fileread(modelpath);
            psy_latent = jsondecode(txt);
            if size(psy_latent.wMode,2) > 200
                weight_session(ss, ii, 1,:) = mean(psy_latent.wMode(:,1:60),2);
                weight_session(ss,ii,2,:) = mean(psy_latent.wMode(:,end-59:end),2);
            end
        end

    end

    % load GLM results

    for vv = 1:length(vars)
        currVar = vars{vv};

        [pathstr, name, ext] = fileparts(psymodelName);
        savedatapath = pathstr;
        GLM_results= fullfile(savedatapath, ['glm_coeff_summary_',currVar,'_', label]);
        load(GLM_results);
        % load baseline
        GLMr.(structVar{vv}).coeff_WT(:,ii) = coef_summary.WT_est;
        GLMr.(structVar{vv}).SE_WT(:,ii) = coef_summary.WT_SE;
        GLMr.(structVar{vv}).coeff_Mut(:,ii) = coef_summary.HET_est;
        GLMr.(structVar{vv}).SE_mut(:,ii) = coef_summary.HET_SE;
        GLMr.(structVar{vv}).p_WT(:,ii) = coef_summary.WT_p;
        GLMr.(structVar{vv}).p_Mut(:,ii) = coef_summary.HET_p;
        GLMr.(structVar{vv}).p_geno(:,ii) = coef_summary.diff_p;


    end
    
    % load average resposne time
    meanRTPath= fullfile(savedatapath, ['meanRT', label,'.mat']);
    load(meanRTPath);
    meanRT_session(1:length(RTData.meanRT),ii) = RTData.meanRT;
    [meanRT_geno{1:length(RTData.meanRT),ii}] = RTData.geno{:};
end

All_geno = unique(genotypes);

if sum(contains(All_geno, 'KO'))>0
    mutGene = 'KO';
elseif sum(contains(All_geno,'HET'))>0
    mutGene = 'HET';
elseif sum(contains(All_geno,'HEM'))>0
    mutGene = 'HEM';
end

%% plot average response time
figure;

%% ======== Subplot 1: AB sessions (1–3) ========
subplot(2,1,1)

% --- WT group ---
wt_idx = strcmp(meanRT_geno(:,1), 'WT');
wt_data = meanRT_session(wt_idx, 1:3);
wt_mean = nanmean(wt_data, 1);
wt_ste  = nanstd(wt_data, [], 1) ./ sqrt(sum(~isnan(wt_data), 1));

% --- Mutant group ---
mut_idx = strcmp(meanRT_geno(:,1), mutGene);
mut_data = meanRT_session(mut_idx, 1:3);
mut_mean = nanmean(mut_data, 1);
mut_ste  = nanstd(mut_data, [], 1) ./ sqrt(sum(~isnan(mut_data), 1));

% --- Plot means + error bars ---
hold on;
x = 1:3;
errorbar(x, wt_mean, wt_ste, '-o', 'LineWidth', 1.5, 'DisplayName', 'WT');
errorbar(x, mut_mean, mut_ste, '-o', 'LineWidth', 1.5, 'DisplayName', mutGene);

xlabel('Session');
ylabel('Mean RT');
lgd = legend('show');
set(lgd, 'box' , 'off')
xlim([0.5 3.5])
xticks(1:3)
box off;
title('Response time across AB 1–3 sessions (mean ± SEM)');

% --- Statistical test: Mann–Whitney (ranksum) for each session ---
pvals = nan(1,3);
for i = 1:3
    x1 = wt_data(:,i);
    x2 = mut_data(:,i);
    % remove NaNs
    x1 = x1(~isnan(x1)); x2 = x2(~isnan(x2));
    if ~isempty(x1) && ~isempty(x2)
        pvals(i) = ranksum(x1, x2);
    end
end

% --- Multiple comparison correction (FDR or Bonferroni) ---
pvals_corr = mafdr(pvals, 'BHFDR', true);  % Benjamini–Hochberg FDR

% --- Annotate p-values on the figure ---
ymax = max([wt_mean+wt_ste, mut_mean+mut_ste], [], 'all');
for i = 1:3
    text(i, ymax*0.95, sprintf('p=%.3f', pvals_corr(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 20);
    if pvals_corr(i) < 0.05
        text(i, ymax*0.65, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    end
end


%% ======== Subplot 2: CD sessions (4–6) ========
subplot(2,1,2)

% --- WT group ---
wt_idx = strcmp(meanRT_geno(:,4), 'WT');
wt_data = meanRT_session(wt_idx, 4:6);
wt_mean = nanmean(wt_data, 1);
wt_ste  = nanstd(wt_data, [], 1) ./ sqrt(sum(~isnan(wt_data), 1));

% --- Mutant group ---
mut_idx = strcmp(meanRT_geno(:,4), mutGene);
mut_data = meanRT_session(mut_idx, 4:6);
mut_mean = nanmean(mut_data, 1);
mut_ste  = nanstd(mut_data, [], 1) ./ sqrt(sum(~isnan(mut_data), 1));

% --- Plot means + error bars ---
hold on;
x = 1:3;
errorbar(x, wt_mean, wt_ste, '-o', 'LineWidth', 1.5, 'DisplayName', 'WT');
errorbar(x, mut_mean, mut_ste, '-o', 'LineWidth', 1.5, 'DisplayName', mutGene);

xlabel('Session');
ylabel('Mean RT');
lgd = legend('show');
set(lgd, 'box' , 'off')
xlim([0.5 3.5])
xticks(1:3)
box off;
title('Response time across CD 1–3 sessions (mean ± SEM)');

% --- Statistical test: Mann–Whitney (ranksum) for each session ---
pvals = nan(1,3);
for i = 1:3
    x1 = wt_data(:,i);
    x2 = mut_data(:,i);
    x1 = x1(~isnan(x1)); x2 = x2(~isnan(x2));
    if ~isempty(x1) && ~isempty(x2)
        pvals(i) = ranksum(x1, x2);
    end
end

% --- Multiple comparison correction ---
pvals_corr = mafdr(pvals, 'BHFDR', true);

% --- Annotate p-values ---
ymax = max([wt_mean+wt_ste, mut_mean+mut_ste], [], 'all');
for i = 1:3
    text(i, ymax*0.85, sprintf('p=%.3f', pvals_corr(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 20);
    if pvals_corr(i) < 0.05
        text(i, ymax*0.65, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    end
end
    print(gcf,'-dpng',fullfile(savefigpath,'Response time'));    %png format
    saveas(gcf, fullfile(savefigpath,'Response time'), 'fig');
    saveas(gcf, fullfile(savefigpath,'Response time'),'svg');

%% scatter plot of last sessions end weight and next session's start weight
WTMask = strcmp(genotypes,'WT');
mutMask = strcmp(genotypes,mutGene);

% get percentage weight (weight next session trial 1-100 / weight last session
% trial last 1-100
% session 1-2
% session 2-3
% weight_session: subject x session (1-6) x start(1)/end(2) x weight(3)

WT_weight1_rel_end = squeeze(abs(weight_session(WTMask, 1, 2, :))./sum(abs(weight_session(WTMask,1,2,:)),4));
WT_weight2_rel_end = squeeze(abs(weight_session(WTMask,2,2,:))./sum(abs(weight_session(WTMask,2,2,:)),4));
WT_weight3_rel_end = squeeze(abs(weight_session(WTMask,3,2,:))./sum(abs(weight_session(WTMask,3,2,:)),4));
WT_weight4_rel_end = squeeze(abs(weight_session(WTMask, 4, 2, :))./sum(abs(weight_session(WTMask,4,2,:)),4));
WT_weight5_rel_end = squeeze(abs(weight_session(WTMask,5,2,:))./sum(abs(weight_session(WTMask,5,2,:)),4));
WT_weight6_rel_end = squeeze(abs(weight_session(WTMask,6,2,:))./sum(abs(weight_session(WTMask,6,2,:)),4));

WT_weight1_rel_start = squeeze(abs(weight_session(WTMask, 1,1, :))./sum(abs(weight_session(WTMask,1,1,:)),4));
WT_weight2_rel_start = squeeze(abs(weight_session(WTMask,2,1,:))./sum(abs(weight_session(WTMask,2,1,:)),4));
WT_weight3_rel_start = squeeze(abs(weight_session(WTMask,3,1,:))./sum(abs(weight_session(WTMask,3,1,:)),4));
WT_weight4_rel_start = squeeze(abs(weight_session(WTMask, 4,1, :))./sum(abs(weight_session(WTMask,4,1,:)),4));
WT_weight5_rel_start = squeeze(abs(weight_session(WTMask,5,1,:))./sum(abs(weight_session(WTMask,5,1,:)),4));
WT_weight6_rel_start = squeeze(abs(weight_session(WTMask,6,1,:))./sum(abs(weight_session(WTMask,6,1,:)),4));

Mut_weight1_rel_end = squeeze(abs(weight_session(mutMask, 1, 2, :))./sum(abs(weight_session(mutMask,1,2,:)),4));
Mut_weight2_rel_end = squeeze(abs(weight_session(mutMask,2,2,:))./sum(abs(weight_session(mutMask,2,2,:)),4));
Mut_weight3_rel_end = squeeze(abs(weight_session(mutMask,3,2,:))./sum(abs(weight_session(mutMask,3,2,:)),4));
Mut_weight4_rel_end = squeeze(abs(weight_session(mutMask, 4, 2, :))./sum(abs(weight_session(mutMask,4,2,:)),4));
Mut_weight5_rel_end = squeeze(abs(weight_session(mutMask,5,2,:))./sum(abs(weight_session(mutMask,5,2,:)),4));
Mut_weight6_rel_end = squeeze(abs(weight_session(mutMask,6,2,:))./sum(abs(weight_session(mutMask,6,2,:)),4));

Mut_weight1_rel_start = squeeze(abs(weight_session(mutMask, 1,1, :))./sum(abs(weight_session(mutMask,1,1,:)),4));
Mut_weight2_rel_start = squeeze(abs(weight_session(mutMask,2,1,:))./sum(abs(weight_session(mutMask,2,1,:)),4));
Mut_weight3_rel_start = squeeze(abs(weight_session(mutMask,3,1,:))./sum(abs(weight_session(mutMask,3,1,:)),4));
Mut_weight4_rel_start = squeeze(abs(weight_session(mutMask, 4,1, :))./sum(abs(weight_session(mutMask,4,1,:)),4));
Mut_weight5_rel_start = squeeze(abs(weight_session(mutMask,5,1,:))./sum(abs(weight_session(mutMask,5,1,:)),4));Mut_weight6_rel_start = squeeze(abs(weight_session(mutMask,6,1,:))./sum(abs(weight_session(mutMask,6,1,:)),4));

%% absolute weight
WT_weight1_end = squeeze(abs(weight_session(WTMask, 1, 2, :)));
WT_weight2_end = squeeze(abs(weight_session(WTMask,2,2,:)));
WT_weight3_end = squeeze(abs(weight_session(WTMask,3,2,:)));
WT_weight4_end = squeeze(abs(weight_session(WTMask, 4, 2, :)));
WT_weight5_end = squeeze(abs(weight_session(WTMask,5,2,:)));
WT_weight6_end = squeeze(abs(weight_session(WTMask,6,2,:)));

WT_weight1_start = squeeze(abs(weight_session(WTMask, 1,1, :)));
WT_weight2_start = squeeze(abs(weight_session(WTMask,2,1,:)));
WT_weight3_start = squeeze(abs(weight_session(WTMask,3,1,:)));
WT_weight4_start = squeeze(abs(weight_session(WTMask, 4,1, :)));
WT_weight5_start = squeeze(abs(weight_session(WTMask,5,1,:)));
WT_weight6_start = squeeze(abs(weight_session(WTMask,6,1,:)));

Mut_weight1_end = squeeze(abs(weight_session(mutMask, 1, 2, :)));
Mut_weight2_end = squeeze(abs(weight_session(mutMask,2,2,:)));
Mut_weight3_end = squeeze(abs(weight_session(mutMask,3,2,:)));
Mut_weight4_end = squeeze(abs(weight_session(mutMask, 4, 2, :)));
Mut_weight5_end = squeeze(abs(weight_session(mutMask,5,2,:)));
Mut_weight6_end = squeeze(abs(weight_session(mutMask,6,2,:)));

Mut_weight1_start = squeeze(abs(weight_session(mutMask, 1,1, :)));
Mut_weight2_start = squeeze(abs(weight_session(mutMask,2,1,:)));
Mut_weight3_start = squeeze(abs(weight_session(mutMask,3,1,:)));
Mut_weight4_start = squeeze(abs(weight_session(mutMask, 4,1, :)));
Mut_weight5_start = squeeze(abs(weight_session(mutMask,5,1,:)));
Mut_weight6_start = squeeze(abs(weight_session(mutMask,6,1,:)));


colWT = [0.3 0.6 0.9];
colMut = [0.9 0.4 0.4];

% plot WT T
predictors = {'Bias', 'Stimulus', 'Stickiness'};
figure;
set(gcf, 'Position', [100 100 900 300]); % wider figure for 3 subplots

for p = 1:3  % loop over predictors
    subplot(1,3,p); hold on;
    
    % concatenate sessions for each group (subjects x sessions)
    WT_all = [WT_weight1_rel_end(:,p), WT_weight2_rel_end(:,p), WT_weight3_rel_end(:,p)];
    Mut_all = [Mut_weight1_rel_end(:,p), Mut_weight2_rel_end(:,p), Mut_weight3_rel_end(:,p)];
    
    % plot individual subjects (WT)
    plot(1:3, WT_all', '--', 'Color', [colWT 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    % plot individual subjects (Mut)
    plot(1:3, Mut_all', '--', 'Color', [colMut 0.5], 'LineWidth', 0.8, 'HandleVisibility', 'off');
    
    % calculate mean and SEM
    WT_mean = mean(WT_all,1,'omitnan');
    WT_sem  = std(WT_all,0,1,'omitnan') ./ sqrt(sum(~isnan(WT_all(:,1))));
    Mut_mean = mean(Mut_all,1,'omitnan');
    Mut_sem  = std(Mut_all,0,1,'omitnan') ./ sqrt(sum(~isnan(Mut_all(:,1))));
    
    % plot mean ± SEM
    h1=errorbar(1:3, WT_mean, WT_sem, 'o-', 'Color', colWT, ...
        'MarkerFaceColor', colWT, 'LineWidth', 2, 'CapSize', 8);
    h2=errorbar(1:3, Mut_mean, Mut_sem, 'o-', 'Color', colMut, ...
        'MarkerFaceColor', colMut, 'LineWidth', 2, 'CapSize', 8);
    
    % format
    xlim([0.8 3.2]);
    xticks(1:3);
    xticklabels({'1','2','3'});
    if p==2
        xlabel('Session')
    end
    ylabel('Relative weight');
    title([predictors{p}]);
    box off;
    
    if p == 3
        legend([h1 h2], {'WT', mutGene}, 'Location','best', 'Box','off');
    end
end

sgtitle('Average relative weight across AB sessions');

%% savings across days
figure;
set(gcf, 'Position', [100 100 900 300]);

% Store p-values for FDR correction
pvals = nan(3,2);

for p = 1:3  % predictors
    subplot(2,3,p); hold on;
    
    % Extract WT and Mut data
    % WT_end = {WT_weight1_rel_end(:,p), WT_weight4_rel_end(:,p)};
    % WT_start = {WT_weight2_rel_start(:,p), WT_weight5_rel_start(:,p)};
    % Mut_end = {Mut_weight1_rel_end(:,p), Mut_weight4_rel_end(:,p)};
    % Mut_start = {Mut_weight2_rel_start(:,p), Mut_weight5_rel_start(:,p)};
    
    WT_end = {WT_weight1_end(:,p),WT_weight2_end(:,p), WT_weight4_end(:,p),WT_weight5_end(:,p)};
    WT_start = {WT_weight2_start(:,p),WT_weight3_start(:,p), WT_weight5_start(:,p),WT_weight6_start(:,p)};
    Mut_end = {Mut_weight1_end(:,p),Mut_weight2_end(:,p), Mut_weight4_end(:,p),Mut_weight5_end(:,p)};
    Mut_start = {Mut_weight2_start(:,p),Mut_weight2_start(:,p), Mut_weight5_start(:,p),Mut_weight6_start(:,p)};


    % Combine sequentially for plotting
    WT_seq = [WT_end{1}, WT_start{1}, WT_end{2}, WT_start{2}];
    Mut_seq = [Mut_end{1}, Mut_start{1}, Mut_end{2}, Mut_start{2}];
    
    % Plot individual subjects (WT)
    for i = 1:size(WT_seq,1)
        plot([1 2], WT_seq(i,[1 2]), '--', 'Color', [colWT 0.4], 'LineWidth', 0.8, 'HandleVisibility','off');
        plot([3 4], WT_seq(i,[3 4]), '--', 'Color', [colWT 0.4], 'LineWidth', 0.8, 'HandleVisibility','off');
    end
    % Plot individual subjects (Mut)
    for i = 1:size(Mut_seq,1)
        plot([1 2], Mut_seq(i,[1 2]), '--', 'Color', [colMut 0.4], 'LineWidth', 0.8, 'HandleVisibility','off');
        plot([3 4], Mut_seq(i,[3 4]), '--', 'Color', [colMut 0.4], 'LineWidth', 0.8, 'HandleVisibility','off');
    end
    
    % --- Compute means and SEM ---
    WT_mean = mean(WT_seq,1,'omitnan');
    WT_sem  = std(WT_seq,0,1,'omitnan') ./ sqrt(sum(~isnan(WT_seq(:,1))));
    Mut_mean = mean(Mut_seq,1,'omitnan');
    Mut_sem  = std(Mut_seq,0,1,'omitnan') ./ sqrt(sum(~isnan(Mut_seq(:,1))));
    
    % --- Plot means and error bars ---
    h1 = errorbar([1 2], WT_mean([1 2]), WT_sem([1 2]), 'o-', ...
        'Color', colWT, 'MarkerFaceColor', colWT, 'LineWidth', 2, 'CapSize', 8);
    errorbar([3 4], WT_mean([3 4]), WT_sem([3 4]), 'o-', ...
        'Color', colWT, 'MarkerFaceColor', colWT, 'LineWidth', 2, 'CapSize', 8, 'HandleVisibility','off');
    
    h2 = errorbar([1 2], Mut_mean([1 2]), Mut_sem([1 2]), 'o-', ...
        'Color', colMut, 'MarkerFaceColor', colMut, 'LineWidth', 2, 'CapSize', 8);
    errorbar([3 4], Mut_mean([3 4]), Mut_sem([3 4]), 'o-', ...
        'Color', colMut, 'MarkerFaceColor', colMut, 'LineWidth', 2, 'CapSize', 8, 'HandleVisibility','off');
    
    % --- Format ---
    xlim([0.7 4.3]);
    xticks(1:4);
    xticklabels({'1 end','2 start','2 end','3 start'});
    ylabel('Relative weight');
    if p<4
        title(predictors{p});
    end
    box off;

    % Build wide-format data for each subject

    % segment 2: 2end - 3start
    pvals(p,2) = ranksum(WT_slope23, Mut_slope23);

    if p == 1
        legend([h1 h2], {'WT mean', 'Mut mean'}, 'Location', 'best', 'Box', 'off');
    end
end

% --- FDR correction ---
[~, ~, ~, p_fdr] = fdr_bh(pvals_all, 0.05, 'pdep', 'yes');

% --- Plot corrected p-values ---
for p = 1:3
    subplot(1,3,p); hold on;
    idx = p_idx{p};
    y_max = max(ylim);
    y_text = y_max * 0.95;
    
    % WT comparisons
    text(1.5, y_text, sprintf('WT: p_{1-2}=%.3f\nWT: p_{2-3}=%.3f', p_fdr(idx(1)), p_fdr(idx(2))), ...
        'Color', colWT, 'HorizontalAlignment', 'center', 'FontSize', 8);
    
    % Mut comparisons
    text(2.5, y_text*0.9, sprintf('Mut: p_{1-2}=%.3f\nMut: p_{2-3}=%.3f', p_fdr(idx(3)), p_fdr(idx(4))), ...
        'Color', colMut, 'HorizontalAlignment', 'center', 'FontSize', 8);
end

sgtitle('Session transitions (End → Start) with Wilcoxon tests and FDR correction');

% test for differnce
for p =1:3
WT_tbl = table( (1:size(WT_weight1_rel_end,1))', ...
    WT_weight1_rel_end(:,p), WT_weight2_rel_start(:,p), ...
    repmat({'WT'}, size(WT_weight1_rel_end,1),1), ...
    'VariableNames', {'Subject','w1_end','w2_start','Genotype'});

Mut_tbl = table( (1:size(Mut_weight1_rel_end,1))'+100, ...
    Mut_weight1_rel_end(:,p), Mut_weight2_rel_start(:,p), ...
    repmat({'Mut'}, size(Mut_weight1_rel_end,1),1), ...
    'VariableNames', {'Subject','w1_end','w2_start','Genotype'});

T = [WT_tbl; Mut_tbl];

% Define within-subject design (2 levels)
within = table({'1_end';'2_start'}, 'VariableNames', {'Session'});

% Fit repeated-measures ANOVA
rm = fitrm(T, 'w1_end,w2_start ~ Genotype', 'WithinDesign', within);

% Run repeated-measures ANOVA
ranovatbl = ranova(rm, 'WithinModel', 'Session');
disp(ranovatbl)
end

%% absolute weight
figure; 
set(gcf,'Position',[100 100 1200 600]);

% Groupings: first 3 predictors, second 3 predictors
for p = 1:3
    % Determine subplot
    
    subplot(3,1,p); hold on;

    
    % Extract weights
    WT_end = {eval(['WT_weight1_end']),eval(['WT_weight2_end']),eval(['WT_weight4_end']),eval(['WT_weight5_end'])};
    WT_start ={eval(['WT_weight2_start']),eval(['WT_weight3_start']),eval(['WT_weight4_start']),eval(['WT_weight5_start'])};
    Mut_end = {eval(['Mut_weight1_end']),eval(['Mut_weight2_end']),eval(['Mut_weight4_end']),eval(['Mut_weight5_end'])};
    Mut_start = {eval(['Mut_weight2_start']),eval(['Mut_weight3_start']),eval(['Mut_weight4_start']),eval(['Mut_weight5_start'])};
    pred = p;
    WT_seq = [WT_end{1}(:,pred), WT_start{1}(:,pred),WT_end{2}(:,pred), WT_start{2}(:,pred), WT_end{3}(:,pred), WT_start{3}(:,pred),WT_end{4}(:,pred), WT_start{4}(:,pred),];
    Mut_seq = [Mut_end{1}(:,pred), Mut_start{1}(:,pred),Mut_end{2}(:,pred), Mut_start{2}(:,pred),Mut_end{3}(:,pred), Mut_start{3}(:,pred),Mut_end{4}(:,pred), Mut_start{4}(:,pred)];
    % Compute slope per animal
    % WT_slope = WT_start_next - WT_end;
    % Mut_slope = Mut_start_next - Mut_end;
    % 
    % Plot individual subjects (WT)
    for i = 1:size(WT_seq,1)
        plot([1 2], WT_seq(i,[1 2]), '--', 'Color', [colWT 0.4], 'LineWidth', 0.8, 'HandleVisibility','off');
        plot([3 4], WT_seq(i,[3 4]), '--', 'Color', [colWT 0.4], 'LineWidth', 0.8, 'HandleVisibility','off');
        plot([5 6], WT_seq(i,[5 6]), '--', 'Color', [colWT 0.4], 'LineWidth', 0.8, 'HandleVisibility','off');
        plot([7 8], WT_seq(i,[7 8]), '--', 'Color', [colWT 0.4], 'LineWidth', 0.8, 'HandleVisibility','off');
    end
    % Plot individual subjects (Mut)
    for i = 1:size(Mut_seq,1)
        plot([1 2], Mut_seq(i,[1 2]), '--', 'Color', [colMut 0.4], 'LineWidth', 0.8, 'HandleVisibility','off');
        plot([3 4], Mut_seq(i,[3 4]), '--', 'Color', [colMut 0.4], 'LineWidth', 0.8, 'HandleVisibility','off');
        plot([5 6], Mut_seq(i,[5 6]), '--', 'Color', [colMut 0.4], 'LineWidth', 0.8, 'HandleVisibility','off');
        plot([7 8], Mut_seq(i,[7 8]), '--', 'Color', [colMut 0.4], 'LineWidth', 0.8, 'HandleVisibility','off');

    end
    
    
    % --- Compute means and SEM ---
    WT_mean = mean(WT_seq,1,'omitnan');
    WT_sem  = std(WT_seq,0,1,'omitnan') ./ sqrt(sum(~isnan(WT_seq(:,1))));
    Mut_mean = mean(Mut_seq,1,'omitnan');
    Mut_sem  = std(Mut_seq,0,1,'omitnan') ./ sqrt(sum(~isnan(Mut_seq(:,1))));
    
    % --- Plot means and error bars ---
    h1 = errorbar([1 2], WT_mean([1 2]), WT_sem([1 2]), 'o-', ...
        'Color', colWT, 'MarkerFaceColor', colWT, 'LineWidth', 2, 'CapSize', 8);
    errorbar([3 4], WT_mean([3 4]), WT_sem([3 4]), 'o-', ...
        'Color', colWT, 'MarkerFaceColor', colWT, 'LineWidth', 2, 'CapSize', 8, 'HandleVisibility','off');
    errorbar([5 6], WT_mean([5 6]), WT_sem([5 6]), 'o-', ...
        'Color', colWT, 'MarkerFaceColor', colWT, 'LineWidth', 2, 'CapSize', 8, 'HandleVisibility','off');
    errorbar([7 8], WT_mean([7 8]), WT_sem([7 8]), 'o-', ...
        'Color', colWT, 'MarkerFaceColor', colWT, 'LineWidth', 2, 'CapSize', 8, 'HandleVisibility','off');

    h2 = errorbar([1 2], Mut_mean([1 2]), Mut_sem([1 2]), 'o-', ...
        'Color', colMut, 'MarkerFaceColor', colMut, 'LineWidth', 2, 'CapSize', 8);
    errorbar([3 4], Mut_mean([3 4]), Mut_sem([3 4]), 'o-', ...
        'Color', colMut, 'MarkerFaceColor', colMut, 'LineWidth', 2, 'CapSize', 8, 'HandleVisibility','off');
    errorbar([5 6], Mut_mean([5 6]), Mut_sem([5 6]), 'o-', ...
        'Color', colMut, 'MarkerFaceColor', colMut, 'LineWidth', 2, 'CapSize', 8, 'HandleVisibility','off');
    errorbar([7 8], Mut_mean([7 8]), Mut_sem([7 8]), 'o-', ...
        'Color', colMut, 'MarkerFaceColor', colMut, 'LineWidth', 2, 'CapSize', 8, 'HandleVisibility','off');
    
    
    if p == 1 
        legend([h1 h2],{'WT',mutGene},'Location','best','Box','off');
    end
    
    % ---- Wilcoxon test within genotype ----
    % p_wt = signrank(WT_end, WT_start_next);
    % p_mut = signrank(Mut_end, Mut_start_next);
    % 
    % % Optional FDR correction can be applied across all 6 predictors later
    % 
    % % Display p-values
    % yl = ylim;
    % text(1.5, yl(2)-0.05*(yl(2)-yl(1)), sprintf('WT p=%.3f',p_wt),'Color',colWT,'HorizontalAlignment','center');
    % text(1.5, yl(2)-0.1*(yl(2)-yl(1)), sprintf('Mut p=%.3f',p_mut),'Color',colMut,'HorizontalAlignment','center');
    % 
    % Labels
    xticks([1 2 3 4 5 6 7 8]); xticklabels({'AB1', 'AB2', 'AB2', 'AB3', 'CD1', 'CD2', 'CD2', 'CD3' });
    ylabel(predictors{p}); box off; 
    %title(predictors{p});
    end
end

sgtitle('Session transitions (End → Start) with individual and mean ± SEM');

%% plot glm result as a function of sessions

% Assume you have loaded GLMr.(structVar{vv}) as described
sessions = 1:3;  % first 3 sessions only
predictors = {'Baseline','Bias','Stimulus','Stick'};

for vv = 1:length(vars)
figure;
sgtitle(['GLM coefficient of weights (AB) on ', vars{vv}])
for pp = 1:4
    subplot(2,2,pp); hold on;

    % Extract WT & Mut coefficients and SE
    coeff_WT  = GLMr.(structVar{vv}).coeff_WT(pp,sessions);
    se_WT     = GLMr.(structVar{vv}).SE_WT(pp,sessions);
    coeff_Mut = GLMr.(structVar{vv}).coeff_Mut(pp,sessions);
    se_Mut    = GLMr.(structVar{vv}).SE_mut(pp,sessions);

    % Plot with error bars
    errorbar(sessions, coeff_WT, se_WT, '-o', 'LineWidth', 3);
    errorbar(sessions, coeff_Mut, se_Mut, '-o', 'LineWidth', 3);
    xlim([0.5 3.5])
    % FDR-corrected p-values
    pvals0_corr = [GLMr.rt.p_WT(pp,sessions); GLMr.rt.p_Mut(pp,sessions)]';
    %pvals0_corr = mafdr(pvals0(:));   % FDR correction
    %pvals0_corr = reshape(pvals0_corr,[],2);

    pvalsDiff_corr = GLMr.rt.p_geno(pp,sessions);

    % Add significance markers
    for ss = 1:numel(sessions)
        y = max([coeff_WT(ss)+se_WT(ss), coeff_Mut(ss)+se_Mut(ss)]) * 1.1;

        % WT vs 0
        if pvals0_corr(ss,1) < 0.05
            text(sessions(ss)-0.1, y, '*','Color','b','FontSize',30);
        end
        % Mut vs 0
        if pvals0_corr(ss,2) < 0.05
            text(sessions(ss)+0.1, y, '*','Color','r','FontSize',30);
        end
        % WT vs Mut
        if pvalsDiff_corr(ss) < 0.05
            plot(sessions(ss), y*1.05, 'k*','MarkerSize',8);
        end
    end

    title(predictors{pp});
    if pp>2
        xlabel('Session');
    end
    if pp == 1 | pp == 3
    ylabel('Coefficient');
    end
    if pp == 4
        lgd = legend({'WT',mutGene},'Location','best');
        set(lgd,'Box','off','Color','none');
    end
end


% save the figure
    print(gcf,'-dpng',fullfile(savefigpath,'latent', ['GLM coefficient of weights (AB) on ', vars{vv}]));    %png format
    saveas(gcf, fullfile(savefigpath,'latent', ['GLM coefficient of weights (AB) on ', vars{vv}]), 'fig');
    saveas(gcf, fullfile(savefigpath, 'latent',['GLM coefficient of weights (AB) on ', vars{vv}]),'svg');
end

% Assume you have loaded GLMr.(structVar{vv}) as described
sessions = 4:6;  % CD sessions
predictors = {'Baseline','Bias','Stimulus','Stick'};

for vv = 1:length(vars)
figure;
sgtitle(['GLM coefficient of weights (CD) on ', vars{vv}])
for pp = 1:4
    subplot(2,2,pp); hold on;

    % Extract WT & Mut coefficients and SE
    coeff_WT  = GLMr.(structVar{vv}).coeff_WT(pp,sessions);
    se_WT     = GLMr.(structVar{vv}).SE_WT(pp,sessions);
    coeff_Mut = GLMr.(structVar{vv}).coeff_Mut(pp,sessions);
    se_Mut    = GLMr.(structVar{vv}).SE_mut(pp,sessions);

    % Plot with error bars
    errorbar(sessions, coeff_WT, se_WT, '-o', 'LineWidth', 3);
    errorbar(sessions, coeff_Mut, se_Mut, '-o', 'LineWidth', 3);
    xlim([3.5 6.5])
    % FDR-corrected p-values
    pvals0_corr = [GLMr.rt.p_WT(pp,sessions); GLMr.rt.p_Mut(pp,sessions)]';
    %pvals0_corr = mafdr(pvals0(:));   % FDR correction
    %pvals0_corr = reshape(pvals0_corr,[],2);

    pvalsDiff_corr = GLMr.rt.p_geno(pp,sessions);

    % Add significance markers
    for ss = 1:numel(sessions)
        y = max([coeff_WT(ss)+se_WT(ss), coeff_Mut(ss)+se_Mut(ss)]) * 1.1;

        % WT vs 0
        if pvals0_corr(ss,1) < 0.05
            text(sessions(ss)-0.1, y, '*','Color','b','FontSize',30);
        end
        % Mut vs 0
        if pvals0_corr(ss,2) < 0.05
            text(sessions(ss)+0.1, y, '*','Color','r','FontSize',30);
        end
        % WT vs Mut
        if pvalsDiff_corr(ss) < 0.05
            plot(sessions(ss), y*1.05, 'k*','MarkerSize',8);
        end
    end

    title(predictors{pp});
    if pp>2
        xlabel('Session');
    end
    if pp == 1 | pp == 3
    ylabel('Coefficient');
    end
    if pp == 4
        lgd = legend({'WT',mutGene},'Location','best');
        set(lgd,'Box','off','Color','none');
    end
end


% save the figure
    print(gcf,'-dpng',fullfile(savefigpath,'latent', ['GLM coefficient of weights (CD) on ', vars{vv}]));    %png format
    saveas(gcf, fullfile(savefigpath,'latent', ['GLM coefficient of weights (CD) on ', vars{vv}]), 'fig');
    saveas(gcf, fullfile(savefigpath, 'latent',['GLM coefficient of weights (CD) on ', vars{vv}]),'svg');
end

end
