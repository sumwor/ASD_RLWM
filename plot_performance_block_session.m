function plot_performance_block_session(data, genotype, savefigpath, tlabel)

if contains(tlabel,'AB') || contains(tlabel, 'CD')
    nPlot = 3;
elseif contains(tlabel, 'DC')
    nPlot = 6;
end

x_plot = 1:20;

figure;
sgtitle(['Performance in the ', tlabel, ' sessions in block-session'])

for pp = 1:nPlot
    subplot(2,3,pp)

    genotype_list = unique(genotype);
    color_list = { 'red', 'blue','magenta'};

    for gg = 1:length(genotype_list)
        % plot with blocks that have more than 3 animals
        not_nan_counts = sum(~isnan(data(strcmp(genotype, genotype_list{gg}),:,pp)));

% Find indices of columns with more than 3 NaNs
        columns_with_nans = find(not_nan_counts  >= 3);

        mean_perf = nanmean(data(strcmp(genotype, genotype_list{gg}),columns_with_nans,pp),1);
        if sum(strcmp(genotype, genotype_list{gg}))>1
            ste_perf = nanstd(data(strcmp(genotype, genotype_list{gg}),columns_with_nans,pp),1)/sum(strcmp(genotype, genotype_list{gg}));
        else
            ste_perf = nan(length(columns_with_nans),1);
        end
        %plot(mean_perf)
        hold on;
        errorbar(columns_with_nans, mean_perf, ste_perf, 'LineWidth', 2,'Color', color_list{gg})

        % plot each session in thinner lines
        tempData= data(strcmp(genotype, genotype_list{gg}),:,pp);
        for tt =1:size(tempData,1)
            plot(x_plot, tempData(tt,:), 'LineWidth', 0.5, 'Color', color_list{gg},'LineStyle',':','HandleVisibility', 'off')
        end

    end
    ylim([0, 1]);
    %xticks([1,2,3,4]);
    %legend(genotype_list)
    title(['Session ', num2str(pp)])
    ylabel('Performance');
    % two-way anova
    if ismember('HEM', genotype_list) % nlgn3
        Hetdata = data(strcmp(genotype, 'HEM'),columns_with_nans,pp);
        nHet = sum(strcmp(genotype, 'HEM'));
    elseif length(genotype_list) == 3 % if it's cntnap, compare KO with WT
        Hetdata = data(strcmp(genotype, 'KO'),columns_with_nans,pp);
        nHet = sum(strcmp(genotype, 'KO'));
    else
        Hetdata = data(strcmp(genotype, 'HET'),columns_with_nans,pp);
        nHet = sum(strcmp(genotype, 'HET'));
    end
    WTdata = data(strcmp(genotype, 'WT'),columns_with_nans,pp);
    nWT = sum(strcmp(genotype, 'WT'));
    
    % check if there is empty data (not enough sessions/animals)
    if ~isempty(WTdata) && ~isempty(Hetdata)
        dat = [WTdata(:); Hetdata(:)];
    
        % Create grouping variables
        group = [repmat({'WT'}, numel(WTdata), 1); repmat({'HET'}, numel(Hetdata), 1)]; % Group labels
        block = [repmat((columns_with_nans)', size(WTdata, 1), 1); repmat((columns_with_nans)', size(Hetdata, 1), 1)];   % Quantile labels
    
        % Perform 2-way ANOVA
        [p, tbl, stats] = anovan(dat, {group, block}, ...
            'model', 'interaction', ...
            'varnames', {'Group', 'Block'});
        % display p-value for group in the figure
    
        text(3, 1, ['P(group):',num2str(p(1))], 'FontSize', 10 )
    end

    if pp==nPlot
        legend(genotype_list, 'Location', 'southeast')
        legend('Box','off')
        text(3, 0.1, ['WT ', num2str(nWT)], 'FontSize', 10 );
        text(3, 0.2, ['HET ', num2str(nHet)], 'FontSize', 10 );
    end
end

%% savefigs
print(gcf,'-dpng',fullfile(savefigpath,['Performance in the  ', tlabel, ' sessions in block-session']));    %png format
saveas(gcf, fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in block-session']), 'fig');
saveas(gcf, fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in block-session']),'svg');

% close all including anova windows
delete(findall(0, 'Type', 'figure'))

%close;
