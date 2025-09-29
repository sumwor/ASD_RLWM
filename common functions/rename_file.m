function rename_file(psymodelName, label,dataIndex)

% edit file name


psy_fit = readtable(psymodelName);
if strcmp(label,'AB-AB1')
    wrongLabel = 'AB1';
    protocol = 'AB';
    protocolDay = 1;
elseif strcmp(label, 'AB-AB2')
    wrongLabel = 'AB2';
    protocol = 'AB';
    protocolDay = 2;
elseif strcmp(label, 'AB-AB3')
    wrongLabel = 'AB3';
    protocol = 'AB';
    protocolDay = 3;
elseif strcmp(label, 'AB-CD-CD1')
    protocol = 'AB-CD';
    protocolDay = 1;
    elseif strcmp(label, 'AB-CD-CD2')
    protocol = 'AB-CD';
    protocolDay = 2;
    elseif strcmp(label, 'AB-CD-CD3')
    protocol = 'AB-CD';
    protocolDay = 3;
end

subjects = psy_fit.Animal;
nSubject = length(subjects);

genotypes = cell(nSubject,1);


for k = 1:length(subjects)
    % load csv file for response time and intertrial interval
    tempgeno = dataIndex.Genotype(strcmp(dataIndex.Animal,num2str(subjects(k))));
    genotypes{k} = tempgeno{1};
    csvfilemask = strcmp(dataIndex.Animal, num2str(subjects(k))) & strcmp(dataIndex.Protocol,protocol) & cell2mat(dataIndex.ProtocolDay) == protocolDay;
    csvfile = dataIndex.BehCSV{csvfilemask};
    

    analysis = dataIndex.BehPath(strcmp(dataIndex.Animal,num2str(subjects(k))));

    
    wrongmodelpath = fullfile(analysis{1},'latent',['psy_fit_',wrongLabel,'.json']);

    modelpath = fullfile(analysis{1},'latent',['psy_fit_',label,'.json']);
    movefile(wrongmodelpath, modelpath)

end
