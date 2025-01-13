function perf = perf_in_quantile(result1, result2)


perf = zeros(1,4, 2);
nTrials1 = size(result1,1); nTrials2 = size(result2, 1);
for qq = 1:4
    startTrial1 = floor(nTrials1/4)*(qq-1) + 1;
    endTrial1 = floor(nTrials1/4)*qq;
    startTrial2 = floor(nTrials2/4)*(qq-1) + 1;
    endTrial2 = floor(nTrials2/4)*qq;
    perf(1,qq, 1) = sum(~isnan(result1.reward(startTrial1:endTrial1)))/(endTrial1-startTrial1+1);
    perf(1,qq, 2) = sum(~isnan(result2.reward(startTrial2:endTrial2)))/(endTrial2-startTrial2+1);

end
