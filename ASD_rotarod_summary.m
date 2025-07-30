function ASD_rotarod_summary(rot_data, saverotpath)

%% simple plot of rotarod performance

Genotype = unique(rot_data.Genotype);
Animals = unique(rot_data.ASDID);
Trials = unique(rot_data.Trial);

perf = nan(length(Animals), length(Trials));
geno_list = cell(length(Animals),1);

for aa = 1:length(Animals)
    geno_list(aa)= unique(rot_data.Genotype(strcmp(rot_data.ASDID,Animals{aa})));
    tempResult = rot_data.TimeSpentOnRod(strcmp(rot_data.ASDID,Animals{aa}));
    tempFall = rot_data.fallByTurning(strcmp(rot_data.ASDID,Animals{aa}));
    tempResult(boolean(tempFall)) = NaN;
    perf(aa,:) = tempResult;
end

% ignore the performance when fall by turning

% convert time on rod to terminal speed
%perf(:,1:6) = 5 + ((40-5)/300) * perf(:,1:6);
%perf(:,7:12) = 10 + ((80-10)/300) * perf(:, 7:12);

perf= 5 + ((80-5)/300) * perf;

%geno_included = {'WT', 'HET'};
geno_included = {'WT', 'HEM'};
figure;

mean_perf_WT = nanmean(perf(strcmp(geno_list, 'WT'),:),1);

ste_perf_WT = nanstd(perf(strcmp(geno_list, 'WT'),:),1)/sum(strcmp(geno_list, 'WT'));

mean_perf_HET = nanmean(perf(strcmp(geno_list, 'HEM'),:),1);

ste_perf_HET = nanstd(perf(strcmp(geno_list, 'HEM'),:),1)/sum(strcmp(geno_list, 'WT'));


    %plot(mean_perf)
x_plot = [1:12];

errorbar(x_plot, mean_perf_WT, ste_perf_WT, 'LineWidth', 2);
hold on;
errorbar(x_plot, mean_perf_HET, ste_perf_HET, 'LineWidth', 2);

%xticks([1,2,3,4]);
legend('WT', 'HEM')
legend('box', 'off')
ylabel('Performance');
set(gca,'box','off')

% ANOVA to test significance
    if ismember('HEM', Genotype)
        Hetdata = perf(strcmp(geno_list, 'HEM'),:)';
        nHet = sum(strcmp(geno_list, 'HEM'));
    else
        Hetdata = perf(strcmp(geno_list, 'HET'),:)';
        nHet = sum(strcmp(geno_list, 'HET'));
    end
    WTdata = perf(strcmp(geno_list, 'WT'),:)';
    nWT = sum(strcmp(geno_list, 'WT'));

    dat = [WTdata(:); Hetdata(:)];

    % Create grouping variables
    group = [repmat({'WT'}, numel(WTdata), 1); repmat({'HET'}, numel(Hetdata), 1)]; % Group labels
    trials= [repmat((1:12)', size(WTdata, 2), 1); repmat((1:12)', size(Hetdata, 2), 1)];   % Quantile labels

    % Perform 2-way ANOVA
    [p, tbl, stats] = anovan(dat, {group, trials}, ...
        'model', 'interaction', ...
        'varnames', {'Group', 'Trial'});
    % display p-value for group in the figure

    text(10, 80, ['P(group):',num2str(p(1))], 'FontSize', 10 )
         

    text(10, 30, ['WT ', num2str(nWT)], 'FontSize', 10 );
    text(10, 28, ['HET ', num2str(nHet)], 'FontSize', 10 );
% save the 

print(gcf,'-dpng',fullfile(saverotpath, 'Rotarod performance'));    %png format
saveas(gcf, fullfile(saverotpath, 'Rotarod performance'), 'fig');
saveas(gcf, fullfile(saverotpath, 'Rotarod performance'),'svg');

