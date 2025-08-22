function ASD_rotarod_summary(rot_data, saverotpath, strain)

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
if strcmp(strain, 'Nlgn3')
    geno_included = {'WT', 'HEM'};
elseif strcmp(strain, 'Cntnap2_KO')
    geno_included = {'WT', 'KO'};
elseif strcmp(strain, 'TSC2') | strcmp(strain, 'ChD8')
    geno_included = {'WT', 'HET'};
end
figure;

mean_perf_WT = nanmean(perf(strcmp(geno_list, 'WT'),:),1);

ste_perf_WT = nanstd(perf(strcmp(geno_list, 'WT'),:),1)/sum(strcmp(geno_list, 'WT'));

het_geno = geno_included{2};
mean_perf_HET = nanmean(perf(strcmp(geno_list, het_geno),:),1);

ste_perf_HET = nanstd(perf(strcmp(geno_list, het_geno),:),1)/sum(strcmp(geno_list, 'WT'));


    %plot(mean_perf)
x_plot = [1:12];

errorbar(x_plot, mean_perf_WT, ste_perf_WT, 'LineWidth', 2);
hold on;
errorbar(x_plot, mean_perf_HET, ste_perf_HET, 'LineWidth', 2);

%xticks([1,2,3,4]);
legend('WT', het_geno)
legend('box', 'off')
ylabel('Performance');
set(gca,'box','off')

% ANOVA to test significance

        Hetdata = perf(strcmp(geno_list, het_geno),:)';
        nHet = sum(strcmp(geno_list, het_geno));

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

    text(1, 55, ['P(geno):',num2str(p(1), '%.3g')], 'FontSize', 20 )
         text(1, 53, ['P(learn):',num2str(p(2), '%.3g')], 'FontSize', 20 )

    text(1, 40, ['WT ', num2str(nWT)], 'FontSize', 20 );
    text(1, 38, ['HET ', num2str(nHet)], 'FontSize', 20 );
% save the 

print(gcf,'-dpng',fullfile(saverotpath, 'Rotarod performance'));    %png format
saveas(gcf, fullfile(saverotpath, 'Rotarod performance'), 'fig');
saveas(gcf, fullfile(saverotpath, 'Rotarod performance'),'svg');

