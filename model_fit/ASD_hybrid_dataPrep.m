function fileNames = ASD_hybrid_dataPrep(dataIndex,savedatafolder)

%% prepare the data for model fitting
% Chase, Lee, 2025

% data structure:
% counted_trial: trial number per session
% schedule: odor presented
% reward1: reward
% animal_ids: .mat file name  Not used
% session: session number for odor pairs (AB1, AB2, etc)   not used
% genotype: 0WT; 1Het
% Animal_ID: ASDID
% action: 0 left; 1 right

odors = {'AB', 'CD', 'DC', 'AB-CD', 'AB-DC'};
protocols = {'AB', 'AB-CD','AB-DC', 'AB-CD', 'AB-DC'};
trialsList = {'nAB', 'nCD', 'nDC', 'nAB', 'nAB'};
nSessions = [3, 3, 5, 3, 5];
nTrials = nan(length(protocols),1);
fileNames = cell(0);
fileCount = 1;
% count the trials first
for pp = 1:length(protocols)
    for ses = 1:nSessions(pp)
        savefilename = fullfile(savedatafolder,['data4model',protocols{pp},'-',odors{pp},num2str(ses),'.mat']);
        fileNames{fileCount} = savefilename;
        if ~exist(savefilename)


            subdataIndex = dataIndex((strcmp(dataIndex.Protocol,protocols{pp}) & cellfun(@(x) x == ses, dataIndex.ProtocolDay)),:);
            totalTrials = sum(subdataIndex.(trialsList{pp}));

            % initiate data table
            X = table(...
                nan(totalTrials,1),...
                nan(totalTrials,1),...
                nan(totalTrials,1),...
                nan(totalTrials,1),...
                cell(totalTrials,1),...
                nan(totalTrials,1)...
                );

            X.Properties.VariableNames = {...
                'counted_trials',...
                'schedule',...
                'reward1',...
                'action',...
                'genotype',...
                'Animal_ID'...
                };

            % load the results
            nAnimals = size(subdataIndex,1);
            trialCount = 1;
            for aa = 1:nAnimals
                resultdf = readtable(fullfile(subdataIndex.BehPath{aa},subdataIndex.BehCSV{aa}));
                X.counted_trials(trialCount:trialCount+subdataIndex.(trialsList{pp})(aa)-1) = 1:subdataIndex.(trialsList{pp})(aa);
                if strcmp(odors{pp}, 'AB') | strcmp(odors{pp}, 'AB-CD') | strcmp(odors{pp}, 'AB-DC')
                    startTrial = 1; endTrial = subdataIndex.nAB(aa);
                elseif strcmp(odors{pp}, 'CD')
                    startTrial = subdataIndex.nAB(aa)+1; endTrial = subdataIndex.nCD(aa)+subdataIndex.nAB(aa);
                elseif strcmp(odors{pp}, 'DC')
                    startTrial = subdataIndex.nAB(aa)+1; endTrial = subdataIndex.nDC(aa)+subdataIndex.nAB(aa);
                end
                tempReward = resultdf.reward; % convert NaN to 0
                tempReward(isnan(resultdf.reward)) = 0;
                tempAction = resultdf.actions; % convert 0,1 to 1,2
                tempAction(resultdf.actions==0) = 1;
                tempAction(resultdf.actions==1) = 2;
                X.schedule(trialCount:trialCount+subdataIndex.(trialsList{pp})(aa)-1) = resultdf.schedule(startTrial:endTrial);
                X.reward1(trialCount:trialCount+subdataIndex.(trialsList{pp})(aa)-1) = tempReward(startTrial:endTrial);
                X.action(trialCount:trialCount+subdataIndex.(trialsList{pp})(aa)-1) = tempAction(startTrial: endTrial);
                X.genotype(trialCount:trialCount+subdataIndex.(trialsList{pp})(aa)-1) = repmat(subdataIndex.Genotype(aa), endTrial-startTrial+1,1);
                X.Animal_ID(trialCount:trialCount+subdataIndex.(trialsList{pp})(aa)-1) = repmat(str2num(subdataIndex.Animal{aa}), endTrial-startTrial+1,1);
                trialCount = trialCount + subdataIndex.(trialsList{pp})(aa);
                %X.schedule()


            end

            save(savefilename, 'X');

        end
        fileCount = fileCount + 1;
    end
end

end


