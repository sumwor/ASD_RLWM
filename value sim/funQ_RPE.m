function negloglike=funQ_RPE(xpar,dat, v_init)
% % funQ_RPE %
%PURPOSE:   Function for maximum likelihood estimation, called by
%           fit_fun().
%
%INPUT ARGUMENTS
%   xpar:       alpha, beta
%   dat:        data
%               dat(:,1) = choice vector
%               dat(:,2) = reward vector
%
%OUTPUT ARGUMENTS
%   negloglike:      the negative log-likelihood to be minimized

%%
beta=xpar(2);
alpha=xpar(1);
nt=size(dat,1);
negloglike=0;

state_Id = unique(dat(:,3));
nStates = length(state_Id);

% initialize action values here
Q = v_init; % two values (left/right for each state)  

for k=1:nt
    deltaQA = Q(1,1)-Q(1,2);
    deltaQB = Q(2,1) - Q(2,2);

    currState = dat(k,3);
    if ismember(currState, [1,3,6])
        correctAction = 0;
    elseif ismember(currState, [2,4,5])
        correctAction = 1;
    end
    
    if ismember(currState, [1,3,5])
        pLeft=1/(1+exp(-beta*deltaQA));
    elseif ismember(currState, [2,4,6])
        pLeft=1/(1+exp(-beta*deltaQB));
    end
    pRight=1-pLeft;
        
    if pLeft==0 
        pLeft=realmin;   % Smallest positive normalized floating point number, because otherwise log(zero) is -Inf
    end        
    if pRight==0 
        pRight=realmin;
    end            
  
    %compare with actual choice to calculate log-likelihood
    if dat(k,1)==1
        logp=log(pRight);
    elseif dat(k,1)==0
        logp=log(pLeft);
    else
        logp = 0;
    end
    negloglike=negloglike-logp;  % calculate log likelihood
    
    % update value for the performed action     %chose right
    if dat(k,3) == 1  % odor A
        if dat(k,1)==0      %chose left
            Q(1,1) =Q(1,1)+alpha*(dat(k,2)-Q(1,1));
        elseif dat(k,1)==1 %chose right
            Q(1,2) =Q(1,2)+alpha*(dat(k,2)-Q(1,2));
        end
    elseif dat(k,3) == 2 % odor B
        if dat(k,1)==0      %chose left
            Q(2,1) =Q(2,1)+alpha*(dat(k,2)-Q(2,1));
        elseif dat(k,1)==1 %chose right
            Q(2,2) =Q(2,2)+alpha*(dat(k,2)-Q(2,2));
        end
    end


end

end
