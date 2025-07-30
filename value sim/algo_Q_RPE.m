function stats=algo_Q_RPE(stats,params,init_v)
% % algo_Q_RPE %
%PURPOSE:   Simulate player based on q-learning with reward predictione
%           errors (see Sul et at., Neuron 2010)
%AUTHORS:   AC Kwan 180423
%
%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%       params(1) = a - learning rate
%       params(2) = b - inverse temperature
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step

a = params(1);
b = params(2);
stats.qA = nan(length(stats.c),2);
stats.qB = nan(length(stats.c),2);
stats.pCorrect = nan(length(stats.c),1);
stats.rpe = nan(length(stats.c),1);

for trial = 1:length(stats.c)
    stats.currTrial = trial;
    if stats.currTrial == 1  %if this is the first trial
        stats.qA(stats.currTrial,1:2) = init_v(1,:);
        stats.qB(stats.currTrial,1:2) = init_v(2,:);
        stats.rpe(stats.currTrial) = NaN;
        
        deltaQA = stats.qA(stats.currTrial,1) - stats.qA(stats.currTrial,2);
        deltaQB = stats.qB(stats.currTrial,1) - stats.qB(stats.currTrial,2);


    else
        %% update action values
        updateState = stats.s(stats.currTrial-1);
        if ismember(updateState, [1,3,5])
            stats.qB(stats.currTrial,:) = stats.qB(stats.currTrial-1,:);
            if stats.c(stats.currTrial-1)==0   % if chose left on last trial
                stats.rpe(stats.currTrial-1)=stats.r(stats.currTrial-1)-stats.qA(stats.currTrial-1,1);
                stats.qA(stats.currTrial,1)=stats.qA(stats.currTrial-1,1)+a*stats.rpe(stats.currTrial-1);
                stats.qA(stats.currTrial,2)=stats.qA(stats.currTrial-1,2);
            elseif stats.c(stats.currTrial-1)==1   % else, chose right
                stats.rpe(stats.currTrial-1)=stats.r(stats.currTrial-1)-stats.qA(stats.currTrial-1,2);
                stats.qA(stats.currTrial,1)=stats.qA(stats.currTrial-1,1);
                stats.qA(stats.currTrial,2)=stats.qA(stats.currTrial-1,2)+a*stats.rpe(stats.currTrial-1);
            else  %miss trials
                stats.rpe(stats.currTrial-1) = 0;
                stats.qA(stats.currTrial,1)=stats.qA(stats.currTrial-1,1);
                stats.qA(stats.currTrial,2)=stats.qA(stats.currTrial-1,2);
            end
        elseif ismember(updateState, [2,4,6])
            stats.qA(stats.currTrial,:) = stats.qA(stats.currTrial-1,:);
            if stats.c(stats.currTrial-1)==0   % if chose left on last trial
                stats.rpe(stats.currTrial-1)=stats.r(stats.currTrial-1)-stats.qB(stats.currTrial-1,1);
                stats.qB(stats.currTrial,1)=stats.qB(stats.currTrial-1,1)+a*stats.rpe(stats.currTrial-1);
                stats.qB(stats.currTrial,2)=stats.qB(stats.currTrial-1,2);
            elseif stats.c(stats.currTrial-1)==1   % else, chose right
                stats.rpe(stats.currTrial-1)=stats.r(stats.currTrial-1)-stats.qB(stats.currTrial-1,2);
                stats.qB(stats.currTrial,1)=stats.qB(stats.currTrial-1,1);
                stats.qB(stats.currTrial,2)=stats.qB(stats.currTrial-1,2)+a*stats.rpe(stats.currTrial-1);
            else  %miss trials
                stats.rpe(stats.currTrial-1) = 0;
                stats.qB(stats.currTrial,1)=stats.qB(stats.currTrial-1,1);
                stats.qB(stats.currTrial,2)=stats.qB(stats.currTrial-1,2);
            end
        end

        %% softmax rule for action selection
        deltaQA = stats.qA(stats.currTrial,1) - stats.qA(stats.currTrial,2);
        deltaQB = stats.qB(stats.currTrial,1) - stats.qB(stats.currTrial,2);


    end
    currState = stats.s(stats.currTrial);
            if ismember(currState, [1,3])
            stats.pCorrect(stats.currTrial)=1/(1+exp(-b*deltaQA));
        elseif currState == 5
            stats.pCorrect(stats.currTrial)=1-1/(1+exp(-b*deltaQA));
        elseif ismember(currState, [2,4])
            stats.pCorrect(stats.currTrial)=1-1/(1+exp(-b*deltaQB));
        elseif currState==6
            stats.pCorrect(stats.currTrial)=1/(1+exp(-b*deltaQB));
        end

end

end
