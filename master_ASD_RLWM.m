function master_ASD_RLWM

% masterfile to process the behavior log files of ASD_RLWM behavior

root_dir = 'Z:\HongliWang\Juvi_ASD Deterministic';

strain_list = {'ChD8'; 'Nlgn3'; 'Cntnap2'; 'Shank3B'; 'TSC2'};

%% add for loop later

animalList = readtable(fullfile(root_dir, 'Cntnap2','Data','AnimalList.csv'));

logfilepath = fullfile(root_dir, 'Nlgn3','Data');
analysispath = fullfile(root_dir, 'Nlgn3','Analysis');
dataIndex = makeDataIndex_ASD(logfilepath, analysispath);

nFiles = size(dataIndex,1);
% parse every .mat file, generate .csv files

odors = cell(nFiles,1);
for ii = 1:nFiles

    if ii==1
        dataIndex.BehCSV = cell(nFiles,1);
        dataIndex.OdorPresented = cell(nFiles,1);
    end
        
    outfname = fullfile(dataIndex.BehPath{ii}, sprintf('%s_%s_behaviorDF.csv', dataIndex.Animal{ii}, dataIndex.Session{ii}));
    dataIndex.BehCSV{ii} = sprintf('%s_%s_behaviorDF.csv', dataIndex.Animal{ii}, dataIndex.Session{ii});
    
    resultdf = extract_behavior_df(fullfile(dataIndex.LogFilePath{ii}, dataIndex.LogFileName{ii}), ...
            fullfile(dataIndex.BehPath{ii}, dataIndex.LogFileName{ii}));

%    if ~isfile(outfname)

        writetable(resultdf, outfname);
 %   end
    
    % parse through the result to check what odor presented in the behavior
    % file
    odor_presented = unique(resultdf.schedule);
    dataIndex.OdorPresented{ii} = odor_presented;
    % calculate raw performance with AB, CD, and DC, add to dataIndex table
%     switch len(odor_presented)
%         case 2
%             % only AB odor
%         case 4
%             % AB,CD or AB,DC
%         case 6
%             % AB, CD, and DC
%     end


    %% process every session, generate single session analysis, move this
    % to a independent for loop
    % 1. behavior session summary
    % 2. response time
    % 3. intertial interval
    %ASD_session(resultdf,dataIndex.Animal{ii}, dataIndex.Session{ii},dataIndex.BehPath{ii} );

end

%% go over every session to plot performance in blocks for
% 1. first two sessions of AB
% 2. first two sessions of CD
% 3. DC reverse

savefigpath = fullfile(root_dir,'Nlgn3','Summary','BehPlot');
if ~exist(savefigpath)
    mkdir(savefigpath)
end
ASD_odor_summary(dataIndex, savefigpath);

