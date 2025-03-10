function ASD_odor_WT_summary(strain_list, root_dir)

%% go over every strains and compare the WT performance

savesumpath = fullfile(root_dir, 'Summary_allstrains');
if ~exist(savesumpath)
    mkdir(savesumpath);
end

data_block_session = cell(0);
data_quantile = cell(0);
for ii =1:length(strain_list)
    savedatafolder = fullfile(root_dir,strain_list{ii},'Summary','Results');
    % load the results
    data_block_session{ii} = load(fullfile(savedatafolder,'perf_block_session.mat' ));
    data_quantile{ii} = load(fullfile(savedatafolder,'perf_block_session.mat' ));
end


% look for wild type
% number of trials performed
% performance
% intertrial interval
% response time

%% performance in blocks per session

geno = 'WT';
phase = {'AB', 'CD', 'DC'};
for pp =1:length(phase)
    curr_phase = phase{pp};
    trials.(curr_phase) = cell(length(strain_list),1);
    if strcmp(curr_phase, 'DC')
            trials.(curr_phase){ss} = nan(nAnimals,20,6);
        else
            trials.(curr_phase){ss} = nan(nAnimals,20,3);
        end
end

for ss =1:length(strain_list)
    nAnimals = sum(strcmp(data_block_session{ss}.genotype,geno));
    for pp =1:length(phase)
        curr_phase=phase{pp};
         struct_inneed = ['perf_',curr_phase,'_block_session'];
        WT_perf = data_block_session{ss}.(struct_inneed)(strcmp(data_block_session{ss}.genotype,geno),:,:);
    % load the data
        for aa =1:nAnimals  
            trials.(curr_phase){ss}(aa,:,:) = WT_perf(aa,:,:);
        end
    end
end

for pphase = 1:length(phase)
    curr_phase = phase{pphase};
    figure;
    sgtitle(['Performance comparison in ', curr_phase, ' between strains'])

    if curr_phase=='DC'
        nplot = 6;
    else
        nplot = 3;
    end
    for pplot = 1:nplot
        subplot(2,3,pplot)

        x_plot = 1:20;
        for gg = 1:length(strain_list)
            mean_perf = nanmean(trials.(curr_phase){gg}(:,:,pplot),1);
            ste_perf = nanstd(trials.(curr_phase){gg}(:,:,pplot),1)/size(trials.(curr_phase){gg},1);

            %plot(mean_perf)
            hold on;
            errorbar(x_plot, mean_perf, ste_perf, 'LineWidth', 2)
        end
        ylim([0 1])
        title(['Session ',num2str(pplot)])

    end
    lgd = legend(strain_list,'Location', 'bestoutside');
    lgd.Position=[0.9 0.25 0.1 0.1];
            legend('box', 'off')
            ylabel('Performance');

     % save figures
     print(gcf,'-dpng',fullfile(savesumpath,['Performance in the  ', curr_phase, ' sessions in block']));    %png format
    saveas(gcf, fullfile(savesumpath, ['Performance in the ', curr_phase, ' sessions in block']), 'fig');
    saveas(gcf, fullfile(savesumpath, ['Performance in the ', curr_phase, ' sessions in block']),'svg');

end
%xticks([1,2,3,4]);


%% savefigs
print(gcf,'-dpng',fullfile(savefigpath,['Number of trials performed']));    %png format
saveas(gcf, fullfile(savefigpath, ['Number of trials performed']), 'fig');
%savefig(fullfile(savefigpath, 'Number of trials performed.fig'));
saveas(gcf, fullfile(savefigpath, ['Number of trials performed']),'svg');

close;

