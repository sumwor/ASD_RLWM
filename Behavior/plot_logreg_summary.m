function plot_logreg_summary(input, logLabel, dataIndex, tLabel, savefigpath)

setup_figprop;
Genotypes = unique(dataIndex.Genotype);
% get trial Mask
Genotypes = unique(dataIndex.Genotype);
if strcmp(tLabel, 'AB')
    nSessions = 3;
    ABMask1 = strcmp(dataIndex.Protocol, 'AB') & (cellfun(@(x) x == 1, dataIndex.ProtocolDay));
    ABMask2 = strcmp(dataIndex.Protocol, 'AB') & (cellfun(@(x) x == 2, dataIndex.ProtocolDay));
    ABMask3 = strcmp(dataIndex.Protocol, 'AB') & (cellfun(@(x) x == 3, dataIndex.ProtocolDay));
    data1= input(ABMask1); data1 = cellfun(@(s) s.AB, data1, 'UniformOutput', false);
    data2 = input(ABMask2);data2 = cellfun(@(s) s.AB, data2, 'UniformOutput', false);
    data3 = input(ABMask3);data3 = cellfun(@(s) s.AB, data3, 'UniformOutput', false);
    genotype1 = dataIndex.Genotype(ABMask1);
    genotype2 = dataIndex.Genotype(ABMask2);
    genotype3 = dataIndex.Genotype(ABMask3);

elseif strcmp(tLabel, 'AB-CD-AB') % AB trials in AB-CD session
    nSessions = 3;
    ABMask1 = strcmp(dataIndex.Protocol, 'AB-CD') & (cellfun(@(x) x == 1, dataIndex.ProtocolDay));
    ABMask2 = strcmp(dataIndex.Protocol, 'AB-CD') & (cellfun(@(x) x == 2, dataIndex.ProtocolDay));
    ABMask3 = strcmp(dataIndex.Protocol, 'AB-CD') & (cellfun(@(x) x == 3, dataIndex.ProtocolDay));
    data1= input(ABMask1); data1 = cellfun(@(s) s.AB, data1, 'UniformOutput', false);
    data2 = input(ABMask2);data2 = cellfun(@(s) s.AB, data2, 'UniformOutput', false);
    data3 = input(ABMask3);data3 = cellfun(@(s) s.AB, data3, 'UniformOutput', false);
    genotype1 = dataIndex.Genotype(ABMask1);
    genotype2 = dataIndex.Genotype(ABMask2);
    genotype3 = dataIndex.Genotype(ABMask3);

elseif strcmp(tLabel, 'AB-CD') % CD trials in AB-CD session
    nSessions = 3;
    ABMask1 = strcmp(dataIndex.Protocol, 'AB-CD') & (cellfun(@(x) x == 1, dataIndex.ProtocolDay));
    ABMask2 = strcmp(dataIndex.Protocol, 'AB-CD') & (cellfun(@(x) x == 2, dataIndex.ProtocolDay));
    ABMask3 = strcmp(dataIndex.Protocol, 'AB-CD') & (cellfun(@(x) x == 3, dataIndex.ProtocolDay));
    data1= input(ABMask1); data1 = cellfun(@(s) s.CD, data1, 'UniformOutput', false);
    data2 = input(ABMask2);data2 = cellfun(@(s) s.CD, data2, 'UniformOutput', false);
    data3 = input(ABMask3);data3 = cellfun(@(s) s.CD, data3, 'UniformOutput', false);
    genotype1 = dataIndex.Genotype(ABMask1);
    genotype2 = dataIndex.Genotype(ABMask2);
    genotype3 = dataIndex.Genotype(ABMask3);

elseif strcmp(tLabel, 'AB-DC-AB') % AB trials in AB-DC session
    nSessions = 5;
    ABMask1 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 1, dataIndex.ProtocolDay));
    ABMask2 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 2, dataIndex.ProtocolDay));
    ABMask3 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 3, dataIndex.ProtocolDay));
    ABMask4 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 4, dataIndex.ProtocolDay));
    ABMask5 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 5, dataIndex.ProtocolDay));

    data1= input(ABMask1); data1 = cellfun(@(s) s.AB, data1, 'UniformOutput', false);
    data2 = input(ABMask2);data2 = cellfun(@(s) s.AB, data2, 'UniformOutput', false);
    data3 = input(ABMask3);data3 = cellfun(@(s) s.AB, data3, 'UniformOutput', false);
    data4= input(ABMask1); data4 = cellfun(@(s) s.AB, data4, 'UniformOutput', false);
    data5 = input(ABMask2);data5 = cellfun(@(s) s.AB, data5, 'UniformOutput', false);

    genotype1 = dataIndex.Genotype(ABMask1);
    genotype2 = dataIndex.Genotype(ABMask2);
    genotype3 = dataIndex.Genotype(ABMask3);
    genotype4 = dataIndex.Genotype(ABMask4);
    genotype5 = dataIndex.Genotype(ABMask5);

elseif strcmp(tLabel, 'AB-DC') % AB trials in AB-DC session
    nSessions = 5;
    ABMask1 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 1, dataIndex.ProtocolDay));
    ABMask2 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 2, dataIndex.ProtocolDay));
    ABMask3 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 3, dataIndex.ProtocolDay));
    ABMask4 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 4, dataIndex.ProtocolDay));
    ABMask5 = strcmp(dataIndex.Protocol, 'AB-DC') & (cellfun(@(x) x == 5, dataIndex.ProtocolDay));

    data1= input(ABMask1); data1 = cellfun(@(s) s.DC, data1, 'UniformOutput', false);
    data2 = input(ABMask2);data2 = cellfun(@(s) s.DC, data2, 'UniformOutput', false);
    data3 = input(ABMask3);data3 = cellfun(@(s) s.DC, data3, 'UniformOutput', false);
    data4= input(ABMask1); data4 = cellfun(@(s) s.DC, data4, 'UniformOutput', false);
    data5 = input(ABMask2);data5 = cellfun(@(s) s.DC, data5, 'UniformOutput', false);

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

%% plot logistic regression result
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
    figname = [tLabel,num2str(nSes),' ', logLabel];
    plot_logreg_genotype(data,genotype, savefigpath, figname, colors)

end

