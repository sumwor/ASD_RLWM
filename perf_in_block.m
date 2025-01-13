function perf = perf_in_block(result1, result2, blockLength)


perf = nan(1,50);
combined_result = [result1; result2];
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
