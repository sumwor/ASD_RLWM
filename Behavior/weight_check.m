function weight_check(weightIndex, label, strain, savefigpath)

% plot pre- and post-weight

if strcmp(label, 'AB')
    nPlot = 3;
elseif strcmp(label, 'all')
    nPlot = 20;
end

if strcmp(strain, 'Nlgn3')
    geno_included = {'HEM','WT'};
    het_geno= 'HEM';
elseif strcmp(strain, 'Cntnap2_KO')
    geno_included = { 'KO','WT'};
    het_geno = 'KO';
elseif strcmp(strain, 'TSC2') | strcmp(strain, 'ChD8')
    geno_included = { 'HET','WT'};
    het_geno = 'HET';
end

animalList = unique(weightIndex.Animal);
nAnimals = length(animalList);

pre_weight = nan(nAnimals, nPlot);
post_weight = nan(nAnimals, nPlot);
genotype = cell(nAnimals,1);

for aa = 1:nAnimals
    currAnimal = animalList(aa);
    subIndex = weightIndex(weightIndex.Animal==currAnimal,:);
    if strcmp(label, 'AB')
        sesMask = strcmp(subIndex.Protocol,'AB') & (subIndex.ProtocolDay<=3);
        pre_weight(aa,:) = subIndex.PreWeight(sesMask);
        post_weight(aa,:) = subIndex.PostWeight(sesMask);
        genotype{aa} = subIndex.Genotype{1};
    end
end

% plot the pre and post weight by genotype

color_list = { 'red', 'blue','magenta'};

figure;

ax1=subplot(1,2,1)


for gg = 1:length(geno_included)
    mean_perf = nanmean(pre_weight(strcmp(genotype, geno_included{gg}),:),1);

    ste_perf = nanstd(pre_weight(strcmp(genotype, geno_included{gg}),:),1)/sum(strcmp(genotype,geno_included{gg}));


    %plot(mean_perf)
    x_plot = [1:nPlot];

    errorbar(x_plot, mean_perf, ste_perf, 'LineWidth', 2,'Color', color_list{gg});
    hold on;
    tempData= pre_weight(strcmp(genotype, geno_included{gg}),:);
    for tt =1:size(tempData,1)
        plot(x_plot, tempData(tt,:), 'LineWidth', 0.5, 'Color', color_list{gg},'LineStyle',':','HandleVisibility', 'off')
    end
end
xticks([1,2,3]);
xlim([0 4])
%legend(geno_included)
%legend('box', 'off')
ylabel('Body weight pre (g)');
set(gca,'box','off')
title('Weight pre-sessioin')
% ANOVA to test significance


Hetdata = pre_weight(strcmp(genotype, het_geno),:)';
nHet = sum(strcmp(genotype, het_geno));

WTdata = pre_weight(strcmp(genotype, 'WT'),:)';
nWT = sum(strcmp(genotype, 'WT'));

dat = [WTdata(:); Hetdata(:)];

% Create grouping variables
group = [repmat({'WT'}, numel(WTdata), 1); repmat({'HET'}, numel(Hetdata), 1)]; % Group labels
trials= [repmat((1:3)', size(WTdata, 2), 1); repmat((1:3)', size(Hetdata, 2), 1)];   % Quantile labels

% Perform 2-way ANOVA
[p, tbl, stats] = anovan(dat, {group, trials}, ...
    'model', 'interaction', ...
    'varnames', {'Group', 'Trial'});
% display p-value for group in the figure

text(0.2, max(dat), ['P(geno):',num2str(p(1), '%.3g')], 'FontSize', 20 )
text(0.2,max(dat)-1, ['P(growth):',num2str(p(2), '%.3g')], 'FontSize', 20 )

%text(1, 40, ['WT ', num2str(nWT)], 'FontSize', 20 );
%text(1, 38, ['HET ', num2str(nHet)], 'FontSize', 20 );
% save the
ax2=subplot(1,2,2)


for gg = 1:length(geno_included)
    mean_perf = nanmean(post_weight(strcmp(genotype, geno_included{gg}),:),1);

    ste_perf = nanstd(post_weight(strcmp(genotype, geno_included{gg}),:),1)/sum(strcmp(genotype,geno_included{gg}));


    %plot(mean_perf)
    x_plot = [1:nPlot];

    errorbar(x_plot, mean_perf, ste_perf, 'LineWidth', 2,'Color', color_list{gg});
    hold on;
    tempData= post_weight(strcmp(genotype, geno_included{gg}),:);
    for tt =1:size(tempData,1)
        plot(x_plot, tempData(tt,:), 'LineWidth', 0.5, 'Color', color_list{gg},'LineStyle',':','HandleVisibility', 'off')
    end
end
xticks([1,2,3]);
xlim([0 4])
legend(geno_included,'location','southeast')
legend('box', 'off')
ylabel('Body weight pre (g)');
set(gca,'box','off')
title('Weight post-sessioin')
% ANOVA to test significance


Hetdata = post_weight(strcmp(genotype, het_geno),:)';
nHet = sum(strcmp(genotype, het_geno));

WTdata = post_weight(strcmp(genotype, 'WT'),:)';
nWT = sum(strcmp(genotype, 'WT'));

dat = [WTdata(:); Hetdata(:)];

% Create grouping variables
group = [repmat({'WT'}, numel(WTdata), 1); repmat({'HET'}, numel(Hetdata), 1)]; % Group labels
trials= [repmat((1:3)', size(WTdata, 2), 1); repmat((1:3)', size(Hetdata, 2), 1)];   % Quantile labels

% Perform 2-way ANOVA
[p, tbl, stats] = anovan(dat, {group, trials}, ...
    'model', 'interaction', ...
    'varnames', {'Group', 'Trial'});
% display p-value for group in the figure

text(0.2, max(dat), ['P(geno):',num2str(p(1), '%.3g')], 'FontSize', 20 )
text(0.2,max(dat)-1, ['P(growth):',num2str(p(2), '%.3g')], 'FontSize', 20 )
linkaxes([ax1, ax2], 'y');
print(gcf,'-dpng',fullfile(savefigpath, ['Weight ', label]));    %png format
saveas(gcf, fullfile(savefigpath, ['Weight ', label]), 'fig');
saveas(gcf, fullfile(savefigpath, ['Weight ', label]),'svg');
