
function negloglike=test_RPE(xpar,dat)


beta=xpar(2);
alpha=xpar(1);
nt=size(dat,1);
negloglike=0;

state_Id = unique(dat(:,3));
nStates = length(state_Id);
originalState = dat(:,3);
if nStates ==4 & ismember(5, state_Id)
    convertedStates = dat(:,3);
    convertedStates(dat(:,3)==5) = 3;
    convertedStates(dat(:,3)==6) = 4;

    dat(:,3) = convertedStates;
end

% initialize action values here
Q = 0.5*ones(nStates, 2); % two values (left/right for each state)  
Q_values = 0.5*ones(nStates,2, nt);

for k=1:nt
    currState = dat(k,3);
    currState_ori = originalState(k);
    if ismember(currState_ori, [1,3,6])
        correctAction = 0;
    elseif ismember(currState_ori, [2,4,5])
        correctAction = 1;
    end

    pCorrect=exp(beta*Q(currState, correctAction+1))/(exp(beta*Q(currState, correctAction+1))+exp(beta*Q(currState, correctAction+1)));
    pIncorrect=1-pCorrect;
        
    if pCorrect==0 
        pCorrect=realmin;   % Smallest positive normalized floating point number, because otherwise log(zero) is -Inf
    end        
    if pIncorrect==0 
        pIncorrect=realmin;
    end            
  
    %compare with actual choice to calculate log-likelihood
    if ~isnan(dat(k,1))
        if dat(k,2)>1
            logp=log(pCorrect);
        else
            logp=log(pIncorrect);
        end  
    else
        logp = 0;
    end
    negloglike=negloglike-logp;  % calculate log likelihood
    
    % update value for the performed action     %chose right
    if ~isnan(dat(k,1))
        Q(currState, dat(k,1)+1) = Q(currState, dat(k,1)+1)+alpha*(dat(k,2)-Q(currState, dat(k,1)+1));   
    end

    %% update Q values
    for qq = 1:nStates
        Q_values(:, :, k) = Q;
    end

end

%% plot code
figure;