function plot_error_summary(input, triallabel, dataIndex, edges,tlabel, savefigpath)
% data: error rate
% data1: error_reinforce

setup_figprop;

% get trial Mask
trialTypes = {'left', 'right'};
Genotypes = unique(dataIndex.Genotype);
if strcmp(tlabel, 'AB')
    nSessions = 3;
    ABMask1 = strcmp(dataIndex.Protocol, 'AB') & (cellfun(@(x) x == 1, dataIndex.ProtocolDay));
    ABMask2 = strcmp(dataIndex.Protocol, 'AB') & (cellfun(@(x) x == 2, dataIndex.ProtocolDay));
    ABMask3 = strcmp(dataIndex.Protocol, 'AB') & (cellfun(@(x) x == 3, dataIndex.ProtocolDay));
    data1 = {}; data2 = {}; data3 = {};
    for fd = 1:length(trialTypes)
        data1.(trialTypes{fd}) = input.AB(:,fd,ABMask1);
        data2.(trialTypes{fd}) = input.AB(:,fd,ABMask2);
        data3.(trialTypes{fd}) = input.AB(:,fd, ABMask3);
    end
    genotype1 = dataIndex.Genotype(ABMask1);
    genotype2 = dataIndex.Genotype(ABMask2);
    genotype3 = dataIndex.Genotype(ABMask3);

elseif strcmp(tlabel, 'AB-CD-AB') % AB trials in AB-CD session
    nSessions = 3;
    ABMask1 = strcmp(dataIndex.Protocol, 'AB-CD') & (cellfun(@(x) x == 1, dataIndex.ProtocolDay));
    ABMask2 = strcmp(dataIndex.Protocol, 'AB-CD') & (cellfun(@(x) x == 2, dataIndex.ProtocolDay));
    ABMask3 = strcmp(dataIndex.Protocol, 'AB-CD') & (cellfun(@(x) x == 3, dataIndex.ProtocolDay));
    data1 = {}; data2 = {}; data3 = {};
    for fd = 1:length(trialTypes)
        data1.(trialTypes{fd}) = input.AB(:,fd,ABMask1);
        data2.(trialTypes{fd}) = input.AB(:,fd,ABMask2);
        data3.(trialTypes{fd}) = input.AB(:,fd,ABMask3);
    end
    genotype1 = dataIndex.Genotype(ABMask1);
    genotype2 = dataIndex.Genotype(ABMask2);
    genotype3 = dataIndex.Genotype(ABMask3);

elseif strcmp(tlabel, 'AB-CD') % CD trials in AB-CD session
    nSessions = 3;
    ABMask1 = strcmp(dataIndex.Protocol, 'AB-CD') & (cellfun(@(x) x == 1, dataIndex.ProtocolDay));
    ABMask2 = strcmp(dataIndex.Protocol, 'AB-CD') & (cellfun(@(x) x == 2, dataIndex.ProtocolDay));
    ABMask3 = strcmp(dataIndex.Protocol, 'AB-CD') & (cellfun(@(x) x == 3, dataIndex.ProtocolDay));
    data1 = {}; data2 = {}; data3 = {};
    for fd = 1:length(trialTypes)
        data1.(trialTypes{fd}) = input.CD(:,fd, ABMask1);
        data2.(trialTypes{fd}) = input.CD(:,fd, ABMask2);
        data3.(trialTypes{fd}) = input.CD(:,fd, ABMask3);
    end
    genotype1 = dataIndex.Genotype(ABMask1);
    genotype2 = dataIndex.Genotype(ABMask2);
    genotype3 = dataIndex.Genotype(ABMask3);

elseif strcmp(tlabel, 'AB-DC-AB') % AB trials in AB-DC session
    nSessions = 5;
    ABMask1 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 1, dataIndex.ProtocolDay));
    ABMask2 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 2, dataIndex.ProtocolDay));
    ABMask3 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 3, dataIndex.ProtocolDay));
    ABMask4 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 4, dataIndex.ProtocolDay));
    ABMask5 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 5, dataIndex.ProtocolDay));

    data1 = {}; data2 = {}; data3 = {};data4 = {}; data5 = {};
    for fd = 1:length(trialTypes)
        data1.(trialTypes{fd}) = input.AB(:,fd,ABMask1);
        data2.(trialTypes{fd}) = input.AB(:,fd,ABMask2);
        data3.(trialTypes{fd}) = input.AB(:,fd,ABMask3);
        data4.(trialTypes{fd}) = input.AB(:,fd,ABMask4);
        data5.(trialTypes{fd}) = input.AB(:,fd,ABMask5);

    end
    genotype1 = dataIndex.Genotype(ABMask1);
    genotype2 = dataIndex.Genotype(ABMask2);
    genotype3 = dataIndex.Genotype(ABMask3);
    genotype4 = dataIndex.Genotype(ABMask4);
    genotype5 = dataIndex.Genotype(ABMask5);

elseif strcmp(tlabel, 'AB-DC') % AB trials in AB-DC session
    nSessions = 5;
    ABMask1 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 1, dataIndex.ProtocolDay));
    ABMask2 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 2, dataIndex.ProtocolDay));
    ABMask3 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 3, dataIndex.ProtocolDay));
    ABMask4 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 4, dataIndex.ProtocolDay));
    ABMask5 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 5, dataIndex.ProtocolDay));

    data1 = {}; data2 = {}; data3 = {};data4 = {}; data5 = {};
    for fd = 1:length(trialTypes)
        data1.(trialTypes{fd}) = input.DC(:,fd,ABMask1);
        data2.(trialTypes{fd}) = input.DC(:,fd,ABMask2);
        data3.(trialTypes{fd}) = input.DC(:,fd,ABMask3);
        data4.(trialTypes{fd}) = input.DC(:,fd,ABMask4);
        data5.(trialTypes{fd}) = input.DC(:,fd,ABMask5);

    end
    genotype1 = dataIndex.Genotype(ABMask1);
    genotype2 = dataIndex.Genotype(ABMask2);
    genotype3 = dataIndex.Genotype(ABMask3);
    genotype4 = dataIndex.Genotype(ABMask4);
    genotype5 = dataIndex.Genotype(ABMask5);


end

if length(Genotypes) == 2
    colors = {'red', 'black'};
elseif length(Genotypes) == 3
    colors = {'blue', 'red', 'black'};
end

%%
for nSes = 1:nSessions
    switch nSes
        case 1
            data = data1;
            genotype = genotype1;
        case 2
            data = data2;
            genotype = genotype2;
        case 3
            data = data3;
            genotype = genotype3;
        case 4
            data = data4;
            genotype = genotype4;
        case 5
            data = data5;
            genotype = genotype5;
        otherwise
            error('Invalid nSes value');
    end
    figure;
    figname = "Error rate with consecutive same " + triallabel + ' '+ tlabel + " session " + nSes + " (s)";
    sgtitle(figname)
    subplot(1,2,1);
    for geno = 1:length(Genotypes)
        hold on;
        y = nanmean(data.left(:,strcmp(genotype,Genotypes{geno})),2);
        sem = nanstd(data.left(:,strcmp(genotype,Genotypes{geno})),0,2)/sqrt(sum(strcmp(genotype,Genotypes{geno})));
        plot(edges(1:end-1),y,'Color',colors{geno});
        errorshade(edges(1:end-1), y+sem, y-sem, colors{geno},0.3);
    end
    xlim([0 6])
    ylim([0, 1.0])
    title('Left choice');
    set(gca,'box','off')
    subplot(1,2,2);
    line_handles = gobjects(length(Genotypes),1);
    for geno = 1:length(Genotypes)
        hold on;
        y = nanmean(data.right(:,strcmp(genotype,Genotypes{geno})),2);
        sem = nanstd(data.right(:,strcmp(genotype,Genotypes{geno})),0,2)/sum(strcmp(genotype,Genotypes{geno}));
        line_handles(geno) = plot(edges(1:end-1),y,'Color',colors{geno});
        errorshade(edges(1:end-1), y+sem, y-sem, colors{geno},0.3);
    end
    xlim([0 6])
    ylim([0, 1.0])
    title('Right choice');
    set(gca,'box','off')

    legend(line_handles, Genotypes, 'box', 'off', 'Color', 'none');

    print(gcf,'-dpng',fullfile(savefigpath, figname));    %png format
    %saveas(gcf, fullfile(BehPath, ['response-time_', tlabel]), 'fig');
    %saveas(gcf, fullfile(BehPath, ['response-time_', tlabel]),'svg');
end


end


