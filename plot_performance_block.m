function plot_performance_block(data, blockLength,genotype, savefigpath, tlabel)


figure;
title(['Performance in the ', tlabel, ' sessions in block'])


genotype_list = unique(genotype);
x_plot = [1:blockLength:50*blockLength];
for gg = 1:length(genotype_list)
    mean_perf = nanmean(data(strcmp(genotype, genotype_list{gg}),:),1);
    if sum(strcmp(genotype, genotype_list{gg}))>1
        ste_perf = nanstd(data(strcmp(genotype, genotype_list{gg}),:),1)/sum(strcmp(genotype, genotype_list{gg}));
    else
        ste_perf = zeros(1, length(x_plot));
    end
    %plot(mean_perf)
    hold on;
    errorbar(x_plot, mean_perf, ste_perf, 'LineWidth', 2)
end
ylim([0, 1]);
%xticks([1,2,3,4]);
legend(genotype_list)
legend('box', 'off')
ylabel('Performance');


%% savefigs
print(gcf,'-dpng',fullfile(savefigpath,['Performance in the  ', tlabel, ' sessions in block']));    %png format
saveas(gcf, fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in block']), 'fig');
saveas(gcf, fullfile(savefigpath, ['Performance in the ', tlabel, ' sessions in block']),'svg');

close;
