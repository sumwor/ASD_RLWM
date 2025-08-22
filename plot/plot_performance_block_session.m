function plot_performance_block_session(data, strain, genotype, savefigpath, tlabel)

if strcmp(tlabel,'AB') | strcmp(tlabel, 'CD') | strcmp(tlabel, 'AB-CD')
    nPlot = 3;
elseif strcmp(tlabel, 'DC') | strcmp(tlabel, 'AB-DC')
    nPlot = 6;
end
if contains(tlabel, 'reversed') | contains(tlabel, 'notreversed')
    if contains(tlabel, 'AB') | contains(tlabel, 'CD')
        nPlot = 3;
    else
        nPlot= 6;
    end
end

% control for 1000 trials maximum
x_plot = 1:10;
data = data(:,x_plot,:);

figure;
sgtitle(['Performance in the ', tlabel, ' sessions in block-session'])

for pp = 1:nPlot
    subplot(2,3,pp)

    genotype_list = unique(genotype);
   if strcmp(strain, 'Cntnap2_KO')  % only look at WT and KO
        genotype_list = {'KO', 'WT'};
   elseif strcmp(strain, 'Shank3B')
       genotype_list={'HET', 'WT'};
    end
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
        nHet = sum(~all(isnan(Hetdata), 2));
    elseif ismember('KO', genotype_list) & length(genotype_list)==3 % if it's cntnap, compare KO with WT
        Hetdata = data(strcmp(genotype, 'HET'),columns_with_nans,pp);
        nHet = sum(~all(isnan(Hetdata), 2));
        KOdata = data(strcmp(genotype, 'KO'), columns_with_nans, pp);
        nKO = sum(~all(isnan(KOdata), 2));
    elseif ismember('KO', genotype_list) & length(genotype_list)==2
        Hetdata = data(strcmp(genotype, 'KO'),columns_with_nans,pp);
        nHet = sum(~all(isnan(Hetdata), 2));
    else
        Hetdata = data(strcmp(genotype, 'HET'),columns_with_nans,pp);
        nHet = sum(~all(isnan(Hetdata), 2));
    end
    WTdata = data(strcmp(genotype, 'WT'),columns_with_nans,pp);
    nWT = sum(~all(isnan(WTdata), 2));
    
    % check if there is empty data (not enough sessions/animals)
    if ~isempty(WTdata) && ~isempty(Hetdata)
        if length(genotype_list) == 2
            WTdata = WTdata'; Hetdata = Hetdata';
            dat = [WTdata(:); Hetdata(:)];
    
        % Create grouping variables
            group = [repmat({'WT'}, numel(WTdata), 1); repmat({'HET'}, numel(Hetdata), 1)]; % Group labels
            block = [repmat((columns_with_nans)', size(WTdata, 2), 1); repmat((columns_with_nans)', size(Hetdata, 2), 1)];   % Quantile labels
        elseif length(genotype_list) == 3
            dat = [WTdata(:); Hetdata(:);KOdata(:)];
            group = [repmat('WT')];
        end
        % Perform 2-way ANOVA

        % removing NaNs
        validIdx = ~isnan(dat);
        dat = dat(validIdx);
        group = group(validIdx);
        block = block(validIdx);

        [p, tbl, stats] = anovan(dat, {group, block}, ...
            'model', 'linear', ...
            'varnames', {'Group', 'Block'});
        % display p-value for group in the figure
        
%         tblData = table(dat, categorical(group), categorical(block), ...
%         'VariableNames', {'Y', 'Group', 'Block'});
%         lme = fitlme(tblData, 'Y ~ Group*Block + (1|Block)');
%         disp(anova(lme))

       text(1, 0.25, ['P(geno): ', num2str(p(1), '%.3g')], 'FontSize', 15)
        text(1, 0.15, ['P(learn): ', num2str(p(2), '%.3g')], 'FontSize', 15)
    end
    if pp==nPlot
        legend(genotype_list, 'Location', 'east','FontSize', 14)
        legend('Box','off')
        text(1, 1, ['WT ', num2str(nWT)], 'FontSize', 10 );
        if ismember('KO', genotype_list) & length(genotype_list)==2
            text(5, 1, ['KO ', num2str(nHet)], 'FontSize', 12 );
        else
        text(5, 1, ['HET ', num2str(nHet)], 'FontSize', 12 );
        end
    end
end



%% savefigs
print(gcf,'-dpng',fullfile(savefigpath,['Performance in the  ', tlabel, ' sessions in block-session']));    %png format
saveas(gcf, fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in block-session']), 'fig');
saveas(gcf, fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in block-session']),'svg');

% close all including anova windows
delete(findall(0, 'Type', 'figure'))

%% compute AUC to compare the performance difference
% normalized over the trial blocks, and chance level (0.5)

perf_AUC = squeeze(nanmean((data-0.5), 2));
figure
    for gg = 1:length(genotype_list)
        % plot with blocks that have more than 3 animals

        mean_perf = nanmean(perf_AUC(strcmp(genotype, genotype_list{gg}),:),1);
       
        ste_perf = nanstd(perf_AUC(strcmp(genotype, genotype_list{gg}),:),1)/sum(strcmp(genotype, genotype_list{gg}));

        %plot(mean_perf)
        hold on;
        errorbar(1:nPlot, mean_perf, ste_perf, 'LineWidth', 2,'Color', color_list{gg})

        % plot each session in thinner lines
        tempData= perf_AUC(strcmp(genotype, genotype_list{gg}),:);
        for tt =1:size(tempData,1)
            plot(1:nPlot, tempData(tt,:), 'LineWidth', 0.5, 'Color', color_list{gg},'LineStyle',':','HandleVisibility', 'off')
        end

    end
    ylim([-0.5, 0.5]);
    %xticks([1,2,3,4]);
    %legend(genotype_list)
    ylabel('Performance (AUC)');
    xticks(0:nPlot);
xticklabels(0:nPlot);
xlim([0.75,nPlot+0.25]);
xlabel('Sessions')
title(['AUC for ', tlabel]);
        legend(genotype_list, 'Location', 'southeast')
        legend('Box','off')

% ANOVA
    if ismember('HEM', genotype_list) % nlgn3
        Hetdata = perf_AUC(strcmp(genotype, 'HEM'),:)';
             nHet = sum(~all(isnan(Hetdata), 2));
    elseif ismember('KO', genotype_list) && length(genotype_list)==3 % if it's cntnap, compare KO with WT
        Hetdata = perf_AUC(strcmp(genotype, 'HET'),:)';
             nHet = sum(~all(isnan(Hetdata), 2));
        KOdata = perf_AUC(strcmp(genotype, 'KO'), :)';
             nHet = sum(~all(isnan(KOdata), 2));
    elseif ismember('KO', genotype_list) & length(genotype_list)==2
        Hetdata = perf_AUC(strcmp(genotype, 'KO'),:)';
             nHet = sum(~all(isnan(Hetdata), 2));
    else
        Hetdata = perf_AUC(strcmp(genotype, 'HET'),:)';
             nHet = sum(~all(isnan(Hetdata), 2));
    end
    WTdata = perf_AUC(strcmp(genotype, 'WT'),:)';
    nWT = sum(strcmp(genotype, 'WT'));
    
    % check if there is empty data (not enough sessions/animals)
    if ~isempty(WTdata) && ~isempty(Hetdata)
        if length(genotype_list) == 2
            dat = [WTdata(:); Hetdata(:)];
    
        % Create grouping variables
            group = [repmat({'WT'}, numel(WTdata), 1); repmat({'HET'}, numel(Hetdata), 1)]; % Group labels
            session = [repmat((1:nPlot)', size(WTdata, 2), 1); repmat((1:nPlot)', size(Hetdata, 2), 1)];   % Quantile labels
        elseif length(genotype_list) == 3
            dat = [WTdata(:); Hetdata(:);KOdata(:)];
            session = [repmat((1:nPlot)', size(WTdata, 2), 1); 
                repmat((1:nPlot)', size(Hetdata, 2), 1);
                 repmat((1:nPlot)', size(KOdata, 2), 1)];
        end
        % Perform 2-way ANOVA
        [p, tbl, stats] = anovan(dat, {group, session}, ...
            'model', 'interaction', ...
            'varnames', {'Group', 'Session'});
        % display p-value for group in the figure
    
        text(1.5, 0.5, ['P(geno): ', num2str(p(1), '%.3g')], 'FontSize', 20 )
        text(1.5, 0.45, ['P(reversal): ', num2str(p(2), '%.3g')], 'FontSize', 20 )
    end

print(gcf,'-dpng',fullfile(savefigpath,['Performance AUC in the  ', tlabel, ' sessions in block-session']));    %png format
saveas(gcf, fullfile(savefigpath, ['Performance in the AUC', tlabel, ' sessions in block-session']), 'fig');
saveas(gcf, fullfile(savefigpath, ['Performance in the AUC', tlabel, ' sessions in block-session']),'svg');

% close all including anova windows
delete(findall(0, 'Type', 'figure'))



