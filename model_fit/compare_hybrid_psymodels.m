function compare_hybrid_psymodels(hybrid_fit, psymodelName,label, savefigpath)

% compare the AIC result of hybrid model and psy model

% average normalized AIC

% number of animals that are better fitted by hybrid and psymodel
psy_fit = readtable(psymodelName);

% figure 1: plot AIC and negLL

figure;

sgtitle(label)

AICs = zeros(length(psy_fit.Animal),2);
AICs(:,1) = hybrid_fit.All_fits(:,3,1);
AICs(:,2) = psy_fit.AIC;
    %ttAICs = [AICs, fit_psy.AIC];
    mAICs = AICs - repmat(mean(AICs, 2), 1, size(AICs, 2));

subplot(1, 2, 1)
hold on

% Bar plot (means)
barH = bar(mean(mAICs), 'FaceColor',[0.7 0.7 0.9]);
errorbar(1:2, mean(mAICs), std(mAICs)./sqrt(size(mAICs,1)), ...
    'k','LineStyle','none','LineWidth',1.2);

% Overlay individual data points (with jitter for visibility)
scatter(ones(size(mAICs,1),1)*1, mAICs(:,1), 40, 'k','filled', ...
    'jitter','on','jitterAmount',0.1);
scatter(ones(size(mAICs,1),1)*2, mAICs(:,2), 40, 'k','filled', ...
    'jitter','on','jitterAmount',0.1);

% Axis settings
xticks(1:2);
xticklabels({'Hybrid','Psy'});
set(gca, 'TickLabelInterpreter', 'none')
ylabel('\Delta AIC')
title('AIC')

% Sign-rank test
[p,~] = signrank(mAICs(:,1), mAICs(:,2));

% Display p-value at (1.5, maxY-2)
maxY = max(mAICs(:));
text(1.5, maxY-2, sprintf('p = %.3f', p), ...
    'HorizontalAlignment','center','FontSize',20);

%% Subplot 2: Genotype preference counts
subplot(1, 2, 2)

% Find which model has smaller AIC for each subject
[~, bestModel] = min(mAICs, [], 2); % 1 = Hybrid, 2 = Psy

% Extract genotypes
genotypes = hybrid_fit.genotypes;  % assume string/cell array of 'WT' or 'HET'

% mutGene
All_geno = unique(genotypes);

if sum(contains(All_geno, 'KO'))>0
    mutGene = 'KO';
elseif sum(contains(All_geno,'HET'))>0
    mutGene = 'HET';
elseif sum(contains(All_geno, 'HEM'))>0
    mutGene = 'HEM';
end

% Count preferences by genotype
nWT_Hybrid = sum(strcmp(genotypes,'WT') & bestModel==1);
nWT_Psy    = sum(strcmp(genotypes,'WT') & bestModel==2);
nHET_Hybrid = sum(strcmp(genotypes,mutGene) & bestModel==1);
nHET_Psy    = sum(strcmp(genotypes,mutGene) & bestModel==2);

counts = [nWT_Hybrid, nWT_Psy;
          nHET_Hybrid, nHET_Psy];

% Plot grouped bar
bar(counts,'stacked')
xticklabels({'WT',mutGene})
ylabel('Number of animals')
lgd = legend({'Hybrid better','Psy better'},'Location','eastoutside');
set(lgd,'Box','off','Color','none')
title('Model preference by genotype')
box off

het_counts = [nHET_Hybrid, nHET_Psy];
nHET_total = sum(het_counts);

% Expected counts under chance = 50/50 split
expected = [nHET_total/2, nHET_total/2];

% Chi-square goodness-of-fit test
[h_chi2, p_chi2, stats_chi2] = chi2gof(1:2, 'Freq', het_counts, 'Expected', expected);

text(1.5, 10, sprintf('p = %.3f', p_chi2), ...
    'HorizontalAlignment','center','FontSize',20);

if ~exist(fullfile(savefigpath,'model'))
    mkdir(fullfile(savefigpath,'model'))
end
savefigfile= fullfile(savefigpath,'model', ['modelComparison_',label]);
print(gcf,'-dpng',savefigfile);    %png format
saveas(gcf, savefigfile, 'fig');
saveas(gcf, savefigfile,'svg');

