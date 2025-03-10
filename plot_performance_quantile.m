function plot_performance_quantile(data, genotype, savefigpath, tlabel)

if contains(tlabel,'AB') || contains(tlabel, 'CD')
    nPlot = 3;
elseif contains(tlabel, 'DC')
    nPlot = 6;
end

figure;
sgtitle(['Performance in the ', tlabel, ' sessions in quantile'])

for pp = 1:nPlot
    subplot(2,3,pp)

    genotype_list = unique(genotype);

    color_list = {'red','blue', 'magenta'};
    for gg = 1:length(genotype_list)
        mean_perf = nanmean(data(strcmp(genotype, genotype_list{gg}),:,pp),1);
        if sum(strcmp(genotype, genotype_list{gg}))>1
            ste_perf = nanstd(data(strcmp(genotype, genotype_list{gg}),:,pp),1)/sum(strcmp(genotype, genotype_list{gg}));
        else
            ste_perf = [0, 0, 0, 0];
        end
        %plot(mean_perf)
        hold on;
        errorbar([1:4], mean_perf, ste_perf, 'LineWidth', 2, 'Color',color_list{gg});

        % plot each session in thinner lines
        tempData= data(strcmp(genotype, genotype_list{gg}),:,pp);
        for tt =1:size(tempData,1)
            plot([1:4], tempData(tt,:), 'LineWidth', 0.5, 'Color', color_list{gg},'LineStyle',':','HandleVisibility', 'off')
        end
    end
    ylim([0, 1]);
    xticks([1,2,3,4]);
    %legend(genotype_list)
    title(['Session ', num2str(pp)])
    ylabel('Performance');

    % two-way anova
    if ismember('HEM', genotype_list)
        Hetdata = data(strcmp(genotype, 'HEM'),:,pp);
        nHet = sum(strcmp(genotype, 'HEM'));
    else
        Hetdata = data(strcmp(genotype, 'HET'),:,pp);
        nHet = sum(strcmp(genotype, 'HET'));
    end
    WTdata = data(strcmp(genotype, 'WT'),:,pp);
    nWT = sum(strcmp(genotype, 'WT'));

    dat = [WTdata(:); Hetdata(:)];

    % Create grouping variables
    group = [repmat({'WT'}, numel(WTdata), 1); repmat({'HET'}, numel(Hetdata), 1)]; % Group labels
    quantile = [repmat((1:4)', size(WTdata, 1), 1); repmat((1:4)', size(Hetdata, 1), 1)];   % Quantile labels

    % Perform 2-way ANOVA
    [p, tbl, stats] = anovan(dat, {group, quantile}, ...
        'model', 'interaction', ...
        'varnames', {'Group', 'Quantile'});
    % display p-value for group in the figure

    text(3, 1, ['P(group):',num2str(p(1))], 'FontSize', 10 )
    
    
    if pp==nPlot
        legend(genotype_list, 'Location','southeast')
        legend('Box','off')
        text(3, 0.1, ['WT ', num2str(nWT)], 'FontSize', 10 );
        text(3, 0.2, ['HET ', num2str(nHet)], 'FontSize', 10 );
    end
end
%% savefigs


print(gcf,'-dpng',fullfile(savefigpath,['Performance in the  ', tlabel, ' sessions in quantile']));    %png format
savefig(fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in quantile.fig']));
% saveas(gcf, fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in quantile']), 'fig');
saveas(gcf, fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in quantile']),'svg');


% close all including anova windows
delete(findall(0, 'Type', 'figure'))
