function plot_performance_block(data, genotype, savefigpath, tlabel)

figure;
sgtitle(['Performance in the first 2 ', tlabel, ' sessions in 100 blocks'])
subplot(2,1,1)

genotype_list = unique(genotype);

for gg = 1:length(genotype_list)
    mean_perf = nanmean(data(strcmp(genotype, genotype_list{gg}),:,1));
    ste_perf = nanstd(data(strcmp(genotype, genotype_list{gg}),:,1))/sum(strcmp(genotype, genotype_list{gg}));
    %plot(mean_perf)
    hold on;
    errorbar([1:4], mean_perf, ste_perf, 'LineWidth', 2)
end
ylim([0, 1]);
xticks([1,2,3,4]);
%legend(genotype_list)
title('Session 1')
ylabel('Performance');

subplot(2,1,2)

for gg = 1:length(genotype_list)
    mean_perf = nanmean(data(strcmp(genotype, genotype_list{gg}),:,2));
    ste_perf = nanstd(data(strcmp(genotype, genotype_list{gg}),:,2))/sum(strcmp(genotype, genotype_list{gg}));
    %plot(mean_perf)
    hold on;
    errorbar([1:4], mean_perf, ste_perf, 'LineWidth', 2)
end
ylim([0, 1])
legend(genotype_list)
legend('Box','off')
title('Session 2')
xticks([1,2,3,4]);
xlabel('Quantile');
ylabel('Performance');

%% savefigs
print(gcf,'-dpng',fullfile(savefigpath,['Performance in the first 2 ', tlabel, ' sessions in quantile']));    %png format
saveas(gcf, fullfile(savefigpath, ['Performance in the first 2 ', tlabel, ' sessions in quantile']), 'fig');
saveas(gcf, fullfile(savefigpath, ['Performance in the first 2 ', tlabel, ' sessions in quantile']),'svg');

close;
