function perf = perf_in_block(results, protocol, blockLength)

% AB: calculate 3 sessions
% CD: calculate 3 sessions
% DC: calculate 6 sessions

if strcmp(protocol, 'AB') || strcmp(protocol,'CD')
    nSessions = 3;
elseif strcmp(protocol, 'DC')
    nSessions = 6;
end

perf = nan(1,50);

for nn = 1:nSessions
    if size(results{nn}) > 1
    if nn==1
        combined_result = results{nn};
    else
        combined_result = [combined_result; results{nn}];
    end
    end
end
nTrials = size(combined_result,1);

for qq = 1:50
    startTrial = blockLength*(qq-1) + 1;
    endTrial = blockLength*qq;
    if endTrial < nTrials && startTrial < nTrials
        perf(1,qq) = sum(~isnan(combined_result.reward(startTrial:endTrial)))/blockLength;
    elseif endTrial > nTrials && startTrial < nTrials
        perf(1,qq) = sum(~isnan(combined_result.reward(startTrial:end)))/(nTrials-startTrial+1);
    end
end
