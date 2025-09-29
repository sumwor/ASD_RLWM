function mixed_linear_psyweight(tbl, label, savefigpath, savedatapath, mutGene)

lme_inter = fitglm(tbl, ...
    'ResponseTime ~ Bias*Genotype + Stimulus*Genotype + Stick*Genotype', ...
    'Distribution','inverse gaussian');

% Coefficients
beta = lme_inter.Coefficients.Estimate;
names = lme_inter.CoefficientNames;
SE   = lme_inter.Coefficients.SE;
pVal = lme_inter.Coefficients.pValue;

% --- Critical t-value for 95% CI ---
df = lme_inter.DFE;
tcrit = tinv(0.975, df);

% Predictors of interest
predictors = {'Bias','Stimulus','Stick'};
nPred = numel(predictors);

% Initialize arrays
WT_est = zeros(1,nPred);
WT_CI  = zeros(2,nPred);
HET_est = zeros(1,nPred);
HET_CI  = zeros(2,nPred);

for k = 1:nPred
    % --- WT slope ---
    idx_main = find(strcmp(names, predictors{k}));
    est_WT = beta(idx_main);
    se_WT  = SE(idx_main);
    WT_est(k) = est_WT;
    WT_CI(:,k) = est_WT + [-1 1]'*tcrit*se_WT;

    % --- HET slope = WT + interaction ---
    idx_int = find(strcmp(names, [predictors{k} ':Genotype_' mutGene]));
    if ~isempty(idx_int)
        est_HET = est_WT + beta(idx_int);
        var_HET = se_WT^2 + SE(idx_int)^2;
        se_HET  = sqrt(var_HET);
    else
        est_HET = est_WT;
        se_HET  = se_WT;
    end
    HET_est(k) = est_HET;
    HET_CI(:,k) = est_HET + [-1 1]'*tcrit*se_HET;
end


% --- Pseudo R² ---
yhat = predict(lme_inter, tbl);

% Observed values
y = tbl.ResponseTime;

% Null model (intercept only)
glm_null = fitglm(tbl, 'ResponseTime ~ 1', 'Distribution','inverse gaussian');
yhat_null = predict(glm_null, tbl);

% Deviance-based pseudo-R²
pseudoR2 = 1 - sum((y - yhat).^2) / sum((y - yhat_null).^2);

%% --- Plot predictor slopes ---
catPred = categorical(predictors, predictors, 'Ordinal',true);
slopes = [WT_est(:), HET_est(:)];
figure; hold on;
b = bar(catPred, slopes, 'grouped');
b(1).FaceColor = [0.3 0.6 0.9];  % WT
b(2).FaceColor = [0.9 0.4 0.4];  % HET

% Error bars
for i = 1:nPred
    errorbar(b(1).XEndPoints(i), WT_est(i), WT_est(i)-WT_CI(1,i), WT_CI(2,i)-WT_est(i), ...
        'k','LineStyle','none','LineWidth',1.2,'CapSize',8);
    errorbar(b(2).XEndPoints(i), HET_est(i), HET_est(i)-HET_CI(1,i), HET_CI(2,i)-HET_est(i), ...
        'k','LineStyle','none','LineWidth',1.2,'CapSize',8);
end

ylabel('Slope (Effect on RT)');
legendHandle = legend({'WT', mutGene}, 'Location','best'); 
set(legendHandle, 'Box', 'off', 'Color', 'none');
title(['Predictor Effects on ', label, ' by Genotype'], 'Interpreter', 'none');

% --- Add genotype × weight p-values ---
p_weight = [];
for i = 1:length(predictors)
    idx_int = find(strcmp(names, [predictors{i} ':Genotype_' mutGene]));
    if ~isempty(idx_int)
        pVal_int = pVal(idx_int);
        p_weight = [p_weight, pVal_int];
        xMid = mean([b(1).XEndPoints(i), b(2).XEndPoints(i)]);
        yTop = mean([WT_est(i); HET_est(i)]) + 0.02;
        text(xMid, yTop, sprintf('p=%.3g', pVal_int), ...
            'HorizontalAlignment','center','FontSize',12,'Color','k');
    end
end

% Save figure
print(gcf,'-dpng',fullfile(savefigpath,'latent', ['glm_weight_', label]));
saveas(gcf, fullfile(savefigpath,'latent', ['glm_weight_', label]), 'fig');
saveas(gcf, fullfile(savefigpath,'latent', ['glm_weight_', label]), 'svg');

% save the data
coef_summary = struct();
coef_summary.predictors = predictors;
coef_summary.WT_est  = WT_est;
coef_summary.WT_CI  = WT_CI; % back out SE from CI
coef_summary.HET_est = HET_est;
coef_summary.HET_CI  = HET_CI;


%% --- Plot baseline/intercept ---
idx_intercept = find(strcmp(names, '(Intercept)'));
idx_geno      = find(strcmp(names, ['Genotype_' mutGene]));

est_WT_base = beta(idx_intercept);
se_WT_base  = SE(idx_intercept);
WT_CI       = est_WT_base + [-1 1]'*tcrit*se_WT_base;

if ~isempty(idx_geno)
    est_HET_base = est_WT_base + beta(idx_geno);
    se_HET_base  = sqrt(se_WT_base^2 + SE(idx_geno)^2);  % approx
    HET_CI       = est_HET_base + [-1 1]'*tcrit*se_HET_base;
else
    est_HET_base = est_WT_base;
    se_HET_base  = se_WT_base;
    HET_CI       = WT_CI;
end

figure; hold on;
x = [1 2];
b1 = bar(x(1), est_WT_base, 0.5, 'FaceColor', [0.3 0.6 0.9]);
b2 = bar(x(2), est_HET_base, 0.5, 'FaceColor', [0.9 0.4 0.4]);

errorbar(x(1), est_WT_base, est_WT_base - WT_CI(1), WT_CI(2) - est_WT_base, ...
    'k', 'LineStyle', 'none', 'LineWidth',1.2, 'CapSize',8);
errorbar(x(2), est_HET_base, est_HET_base - HET_CI(1), HET_CI(2) - est_HET_base, ...
    'k', 'LineStyle', 'none', 'LineWidth',1.2, 'CapSize',8);

set(gca, 'XTick', x, 'XTickLabel', {'WT', mutGene});
lg = legend([b1 b2], {'WT','HET'}, 'Box','off','Color','none');
ylabel('Baseline RT (Intercept)');
title(['Baseline ', label, ' by Genotype'], 'Interpreter','none');

% p-value
if ~isempty(idx_geno)
    pVal_geno = pVal(idx_geno);
    yl = ylim;
    yText = max(WT_CI(2), HET_CI(2)) + 0.05*(yl(2)-yl(1));
    text(mean(x), yText, sprintf('p=%.3g', pVal_geno), 'HorizontalAlignment','center', 'FontSize',12);
end

box off;
hold off;



print(gcf,'-dpng',fullfile(savefigpath,'latent', ['glm_baseline_', label]));
saveas(gcf, fullfile(savefigpath,'latent', ['glm_baseline_', label]), 'fig');
saveas(gcf, fullfile(savefigpath,'latent', ['glm_baseline_', label]), 'svg');

% save the coefficient


%% --- Extract coefficients and p-values directly ---
names = lme_inter.CoefficientNames;
beta  = lme_inter.Coefficients.Estimate;
SE    = lme_inter.Coefficients.SE;
pVal  = lme_inter.Coefficients.pValue;
df    = lme_inter.DFE;

predictors = {'Bias','Stimulus','Stick'};
nPred = numel(predictors);

% Initialize output structure
coef_summary = struct();
coef_summary.predictors = [{'Baseline'},predictors];

% WT slopes
WT_est = zeros(1,nPred);
WT_SE  = zeros(1,nPred);
WT_p   = zeros(1,nPred);

% HET slopes
HET_est = zeros(1,nPred);
HET_SE  = zeros(1,nPred);
HET_p   = zeros(1,nPred);

for k = 1:nPred
    % --- WT slope ---
    idx_main = find(strcmp(names, predictors{k}));
    WT_est(k) = beta(idx_main);
    WT_SE(k)  = SE(idx_main);
    WT_p(k)   = pVal(idx_main);

    % --- HET slope (WT + interaction) ---
    idx_int = find(strcmp(names, [predictors{k} ':Genotype_' mutGene]));
    if ~isempty(idx_int)
        HET_est(k) = beta(idx_main) + beta(idx_int);
        HET_SE(k)  = sqrt(SE(idx_main)^2 + SE(idx_int)^2); % approximate
        t_HET      = HET_est(k) / HET_SE(k);
        HET_p(k)   = 2*(1 - tcdf(abs(t_HET), df)); % two-tailed
    else
        HET_est(k) = WT_est(k);
        HET_SE(k)  = WT_SE(k);
        HET_p(k)   = WT_p(k);
    end
end


%% --- Baseline / intercept ---
idx_intercept = find(strcmp(names,'(Intercept)'));
WT_base_est = beta(idx_intercept);
WT_base_SE  = SE(idx_intercept);
WT_base_p   = pVal(idx_intercept);

idx_geno = find(strcmp(names,['Genotype_' mutGene]));
if ~isempty(idx_geno)
    HET_base_est = beta(idx_intercept) + beta(idx_geno);
    HET_base_SE  = sqrt(SE(idx_intercept)^2 + SE(idx_geno)^2);
    t_HET_base   = HET_base_est / HET_base_SE;
    HET_base_p   = 2*(1 - tcdf(abs(t_HET_base), df));
else
    HET_base_est = WT_base_est;
    HET_base_SE  = WT_base_SE;
    HET_base_p   = WT_base_p;
end

coef_summary.WT_est = [WT_base_est,WT_est];
coef_summary.WT_SE  = [WT_base_SE, WT_SE];
coef_summary.WT_p   = [WT_base_p,WT_p];

coef_summary.HET_est = [HET_base_est,HET_est];
coef_summary.HET_SE  = [HET_base_SE,HET_SE];
coef_summary.HET_p   = [HET_base_p,HET_p];


diff_p = zeros(1,nPred);

for k = 1:nPred
    idx_int = find(strcmp(names, [predictors{k} ':Genotype_' mutGene]));
    if ~isempty(idx_int)
        % Interaction p-value = difference between HET and WT
        diff_p(k) = pVal(idx_int);
    else
        diff_p(k) = NaN; % no interaction term -> no difference
    end
end

% Baseline difference (intercept)
if ~isempty(idx_geno)
    baseline_diff_p = pVal(idx_geno);
else
    baseline_diff_p = NaN;
end

% Add to coef_summary
coef_summary.diff_p = [baseline_diff_p,diff_p];


%% --- Save structure ---
save(fullfile(savedatapath, ['glm_coeff_summary_', label, '.mat']), 'coef_summary');

% Save as .mat
save(fullfile(savedatapath,['glm_coeff_summary_',label,'.mat']), 'coef_summary');
 