function master_ASD_RLWM

% masterfile to process the behavior log files of ASD_RLWM behavior

root_dir = 'Z:\HongliWang\Juvi_ASD Deterministic';

strain_list = {'ChD8'; 'Nlgn3'; 'Cntnap2'; 'Shank3B'; 'TSC2'};

%% add for loop later

animalList = readtable(fullfile(root_dir, 'Cntnap2','Data','AnimalList.csv'));

logfilepath = fullfile(root_dir, 'Cntnap2','Data');
analysispath = fullfile(root_dir, 'Cntnap2','Analysis');
dataIndex = makeDataIndex_ASD(logfilepath, analysispath);

nFiles = size(dataIndex,1);
% parse every .mat file, generate .csv files
for ii = 1:nFiles
    resultdf = extract_behavior_df(fullfile(dataIndex.LogFilePath{ii}, dataIndex.LogFileName{ii}), ...
        fullfile(dataIndex.BehPath{ii}, dataIndex.LogFileName{ii}));
    outfname = fullfile(dataIndex.BehPath{ii}, sprintf('%s_%s_behaviorDF.csv', dataIndex.Animal{ii}, dataIndex.Session{ii}));
    if ii==1
        dataIndex.BehCSV = cell(nFiles,1);
    end
    dataIndex.BehCSV{ii} = sprintf('%s_%s_behaviorDF.csv', dataIndex.Animal{ii}, dataIndex.Session{ii});
    if ~isfile(outfname)
        writetable(resultdf, outfname);
    end

    %% process every session, generate single session analysis
    % 1. behavior session summary
    % 2. response time
    % 3. intertial interval
    ASD_session(resultdf,dataIndex.Animal{ii}, dataIndex.Session{ii},dataIndex.BehPath{ii} );

end