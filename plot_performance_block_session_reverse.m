function plot_performance_block_session_reverse(data, strain, genotype, savefigpath, tlabel, reversed, not_reversed)

%% compare performance for animals that reversed and not reversed
% in WT and HET group

if contains(tlabel,'AB') || contains(tlabel, 'CD')
    nPlot = 3;
elseif contains(tlabel, 'DC')
    nPlot = 6;
end

% control for 1000 trials maximum
x_plot = 1:10;
data = data(:,x_plot,:);
genotype_list = unique(genotype);
if strcmp(strain, 'Cntnap2')  % only look at WT and KO
    genotype_list = {'KO', 'WT'};
elseif strcmp(strain, 'Shank3B')
    genotype_list={'KO', 'WT', 'HET'};
end
color_list = { 'red', 'blue','magenta'};

%% compute AUC to compare the performance difference
% normalized over the trial blocks, and chance level (0.5)

perf_AUC = squeeze(nanmean((data-0.5), 2));


%% a mixed effect ANOVA to test for difference
WT_reversed = perf_AUC(intersect(find(strcmp(genotype, 'WT')), reversed),:)';
WT_not_reversed = perf_AUC(intersect(find(strcmp(genotype, 'WT')), not_reversed),:)';
if strcmp(strain, 'Nlgn3')
    HetGene = 'HEM';
elseif strcmp(strain, 'TSC2')
    HetGene = 'HET';
elseif strcmp(strain, 'ChD8')
    HetGene = 'HET';
elseif strcmp(strain, 'Cntnap2')
    HetGene = 'HET';
end
Het_reversed = perf_AUC(intersect(find(strcmp(genotype, HetGene)), reversed),:)';
Het_not_reversed = perf_AUC(intersect(find(strcmp(genotype, HetGene)), not_reversed),:)';

n1 = size(WT_reversed, 2);
n2 = size(WT_not_reversed, 2);
n3 = size(Het_reversed, 2);
n4 = size(Het_not_reversed, 2);

% Convert matrices into column vectors (long format)
performance = [WT_reversed(:); WT_not_reversed(:); Het_reversed(:); Het_not_reversed(:)];

% Create group labels
genotype_ANOVA = [repmat({'WT'}, n1*nPlot, 1); repmat({'WT'}, n2*nPlot, 1); ...
    repmat({'Het'}, n3*nPlot, 1); repmat({'Het'}, n4*nPlot, 1)];

reversal = [repmat({'Reversed'}, n1*nPlot, 1); repmat({'Not Reversed'}, n2*nPlot, 1); ...
    repmat({'Reversed'}, n3*nPlot, 1); repmat({'Not Reversed'}, n4*nPlot, 1)];

session = repmat([1:nPlot]', sum([n1, n2, n3, n4]), 1); % Session (within-subject factor)

% Assign unique animal IDs
animal_id = [repelem(1:n1,nPlot)'; repelem(n1 + (1:n2),nPlot)'; ...
    repelem(n1+n2 + (1:n3),nPlot)'; repelem(n1+n2+n3 + (1:n4),nPlot)'];

% Create a table
tbl = table(performance, genotype_ANOVA, reversal, session, animal_id, ...
    'VariableNames', {'Performance', 'Genotype', 'Reversal', 'Session', 'Animal'});

% Convert categorical variables
tbl.Genotype = categorical(tbl.Genotype);
tbl.Reversal = categorical(tbl.Reversal);
tbl.Session = categorical(tbl.Session);
tbl.Animal = categorical(tbl.Animal);

%% examine reveral in WT and reversal effect in HET
% and also genotype effect in reversed/not reversed

WT_tbl = tbl(tbl.Genotype == 'WT', :);
[p_revEffect_WT, tbl_anova, stats] = anovan(WT_tbl.Performance, {WT_tbl.Reversal, WT_tbl.Session}, ...
    'model', 'interaction', 'varnames', {'Reversal', 'Session'}, 'random', 2);

HET_tbl = tbl(tbl.Genotype == 'Het', :);
[p_revEffect_Het, tbl_anova, stats] = anovan(HET_tbl.Performance, {HET_tbl.Reversal, HET_tbl.Session}, ...
    'model', 'interaction', 'varnames', {'Reversal', 'Session'}, 'random', 2);

rev_tbl = tbl(tbl.Reversal == 'Reversed', :);
[p_genoEffect_reversed, tbl_anova, stats] = anovan(rev_tbl.Performance, {rev_tbl.Genotype, rev_tbl.Session}, ...
    'model', 'interaction', 'varnames', {'Genotype', 'Session'}, 'random', 2);

notrev_tbl = tbl(tbl.Reversal == 'Not Reversed', :);
[p_genoEffect_notreversed, tbl_anova, stats] = anovan(notrev_tbl.Performance, {notrev_tbl.Genotype, notrev_tbl.Session}, ...
    'model', 'interaction', 'varnames', {'Genotype', 'Session'}, 'random', 2);

%% missing NaNs - to be solved?

% lme = fitlme(tbl, 'Performance ~ Genotype * Reversal * Session + (1|Animal)', 'FitMethod', 'ML');
% % Display results
% disp(anova(lme));

%% make the plot
figure
sgtitle(['AUC for ', tlabel]);
for pp=1:2
    subplot(1,2,pp)
    for gg = 1:length(genotype_list)
        % plot with blocks that have more than 3 animals
        if pp==1
            animalMask = reversed;
        else
            animalMask = not_reversed;
        end
        matching_indices = intersect(find(strcmp(genotype, genotype_list{gg})), animalMask);
        mean_perf = nanmean(perf_AUC(matching_indices,:),1);

        ste_perf = nanstd(perf_AUC(matching_indices,:),1)/length(matching_indices);

        %plot(mean_perf)
        hold on;
        errorbar(1:nPlot, mean_perf, ste_perf, 'LineWidth', 2,'Color', color_list{gg})

        % plot each session in thinner lines
        tempData= perf_AUC(matching_indices,:);
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
    if pp==1
        title('reversed');
    else
        title('not reversed')
    end


    if pp==1
        text(2, -0.1, ['P(geno):',num2str(p_genoEffect_reversed(1))], 'FontSize', 10 )
    else
            legend(genotype_list, 'Location', 'southeast')
    legend('Box','off')
         text(2, -0.1, ['P(geno):',num2str(p_genoEffect_notreversed(1))], 'FontSize', 10 )
          text(2, -0.12, ['P(reversal), WT:',num2str(p_revEffect_WT(1))], 'FontSize', 10 )
            text(2, -0.14, ['P(reversal), Het:',num2str(p_revEffect_Het(1))], 'FontSize', 10 )
    end
end

print(gcf,'-dpng',fullfile(savefigpath,['Performance AUC in the  ', tlabel, ' sessions in block-session']));    %png format
saveas(gcf, fullfile(savefigpath, ['Performance in the AUC', tlabel, ' sessions in block-session']), 'fig');
saveas(gcf, fullfile(savefigpath, ['Performance in the AUC', tlabel, ' sessions in block-session']),'svg');


% close all including anova windows
delete(findall(0, 'Type', 'figure'))



