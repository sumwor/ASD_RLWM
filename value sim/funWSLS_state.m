function negloglike=funWSLS_state(xpar,dat)
% % funWSLS %
%PURPOSE:   Function for maximum likelihood estimation, called by
%           fit_fun().
%
%INPUT ARGUMENTS
%   xpar:       p (probability of following win-stay, lose-switch
%   dat:        data
%               dat(:,1) = choice vector
%               dat(:,2) = reward vector
%               dat(:,3) = state vector
%
%OUTPUT ARGUMENTS
%   negloglike:      the negative log-likelihood to be minimized

%% for sessions with more than 1 pair of odors
% skip the readiness check
% skip the AB-CD-DC sessions as well
%%
nStates = length(xpar);
state_Id = unique(dat(:,3));
% 
% CD sessions
% convert state 3,4 to 1,2
% convert state 5,6 to 2,1
if nStates ==4 & ismember(5, state_Id)
    convertedStates = dat(:,3);
    convertedStates(dat(:,3)==5) = 3;
    convertedStates(dat(:,3)==6) = 4;

    dat(:,3) = convertedStates;
end

p1=xpar(1);  % prob to use WSLS in state 1
p2 = xpar(2);  % prob to use WSLS in state 2
nt=size(dat,1);
negloglike=0;

for k=1:nt
    currState = dat(k,3);
    %probability of choosing right
    if k==1  %first trial
        pCorrect=0.5;
    else % find the outcome of the previous nearest trial with the same odor
        prevTrial = find(dat(1:k-1,3) == currState);
        if ~isempty(prevTrial)
            %display(prevTrial)
            if dat(prevTrial(end),2)>0   %last trial = win + right
                pCorrect=xpar(currState);
            elseif dat(prevTrial(end),2)<0  %last trial = win + left
                pCorrect=1-xpar(currState);
            end
        else
            pCorrect = 0.5;
        end
    end
        
    if pCorrect==0
        pCorrect=realmin;
    elseif pCorrect==1
        pCorrect = 1-realmin;% Smallest positive normalized floating point number, because otherwise log(zero) is -Inf
    end        

    %compare with actual choice to calculate log-likelihood
    if dat(k,2)>0
        logp=log(pCorrect);
    else
        logp=log(1-pCorrect);
    end   
    negloglike=negloglike-logp;  % calculate log likelihood
    
end
