function master_ASD_RLWM

%% to do:
% check the timestamp of trial start to make sure there is no carrying over
% from last session under same protocol


%% masterfile to process the behavior log files of ASD_RLWM behavior

root_dir = 'Z:\HongliWang\Juvi_ASD Deterministic';

strain_list = {'ChD8'; 'Nlgn3'; 'Cntnap2'; 'TSC2';'Shank3B'};

%% add for loop later
strainNum =2;
animalList = readtable(fullfile(root_dir, strain_list{strainNum},'Data','AnimalList.csv'));

%dataIndex = makeDataIndex_ASD(logfilepath);

logfilepath = fullfile(root_dir, strain_list{strainNum},'Data');
analysispath = fullfile(root_dir, strain_list{strainNum},'Analysis');
dataIndex = makeDataIndex_ASD(logfilepath, analysispath);

nFiles = size(dataIndex,1);
% parse every .mat file, generate .csv files

ErrorList = table([],[],'VariableNames',{'Session','ErrorMessage'});

odors = cell(nFiles,1);
for ii = 1:nFiles

    if ii==1
        dataIndex.BehCSV = cell(nFiles,1);
        dataIndex.OdorPresented = cell(nFiles,1);
    end

    outfname = fullfile(dataIndex.BehPath{ii}, sprintf('%s_%s_behaviorDF.csv', dataIndex.Animal{ii}, dataIndex.Session{ii}));
    dataIndex.BehCSV{ii} = sprintf('%s_%s_behaviorDF.csv', dataIndex.Animal{ii}, dataIndex.Session{ii});

    try
        resultdf = extract_behavior_df(fullfile(dataIndex.LogFilePath{ii}, dataIndex.LogFileName{ii}));
        writetable(resultdf, outfname);
        %   end

        % parse through the result to check what odor presented in the behavior
        % file
        odor_presented = unique(resultdf.schedule);
        dataIndex.OdorPresented{ii} = odor_presented;

    catch ME
        newEntry = {[dataIndex.Animal{ii},'_', dataIndex.Session{ii}], ME.message}
        ErrorList = [ErrorList; newEntry];
    end

end

% save the error table
writetable(ErrorList, fullfile(root_dir,strain_list{strainNum}));

%% process every session, generate single session analysis, move this
% to a independent for loop
% 1. behavior session summary
% 2. number of AB/CD/DC trials experienced per session
% 2. response time
% 3. intertial interval
buggedFileList = [];
for ii = 1:nFiles
    resultdf = readtable(fullfile(dataIndex.BehPath{ii},dataIndex.BehCSV{ii}));
    stim = unique(resultdf.schedule);
    if length(stim)==4
        if ismember(3, stim)
            numPorts = length(unique(resultdf.port_side(resultdf.schedule==3)));
        elseif ismember(5, stim)
            numPorts = length(unique(resultdf.port_side(resultdf.schedule==5)));
        end
    elseif length(stim) == 6
        numPorts = max(length(unique(resultdf.port_side(resultdf.schedule==3))),length(unique(resultdf.port_side(resultdf.schedule==5))));
    else
        numPorts = 0;
    end

    if numPorts >= 2
        buggedFileList = [buggedFileList, ii];
    end
end

%ASD_session(resultdf,dataIndex.Animal{ii}, dataIndex.Session{ii},dataIndex.BehPath{ii} );


%% check how trial number build up over time

%% go over every session to plot performance in blocks for
% 1. first two sessions of AB
% 2. first two sessions of CD
% 3. DC reverse

savefigpath = fullfile(root_dir,strain_list{strainNum},'Summary','BehPlot');
if ~exist(savefigpath)
    mkdir(savefigpath)
end
savedatapath = fullfile(root_dir,strain_list{strainNum},'Summary','Results');
if ~exist(savedatapath)
    mkdir(savedatapath)
end
ASD_odor_summary(dataIndex, savefigpath, savedatapath);

%% compare WT performance across strain

ASD_odor_WT_summary(strain_list, root_dir);


%% rotarod performance
strainNum =2;
rotarod_path = fullfile(root_dir, strain_list{strainNum}, 'Data','rotarod.csv');
saverotpath = fullfile(root_dir, strain_list{strainNum},'Summary', 'Rotarod');
if ~exist(saverotpath)
    mkdir(saverotpath);
end
rot_data = readtable(rotarod_path);
ASD_rotarod_summary(rot_data, saverotpath);

%% correlate performance between odor and rotarod
ASD_odor_rotarod(strain_list,strainNum, root_dir);