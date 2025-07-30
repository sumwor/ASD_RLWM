function stats=algo_FQ_RPE_CK(stats,params,init_v)
% % algo_FQ_RPE_CK %
%PURPOSE:   Simulate player based on q-learning with reward prediction
%           errors, with differential learning rates with choice tendency
%AUTHORS:   H Atilgan & AC Kwan 04302020
%
%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%       params(1) = a - learning rate for choice yielding reward
%       params(2) = beta - inverse temperature
%       params(3) = a CK- learning rate for choice tendency
%       params(4) = beta CK - inverse temperature for choice tendency
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step
%%
a = params(1);
b= params(2);
alpha_c = params(3);
beta_c = params(4);

stats.qA = nan(length(stats.c),2);
stats.qB = nan(length(stats.c),2);
stats.ck = nan(length(stats.c),2);
stats.pCorrect = nan(length(stats.c),1);
stats.rpe = nan(length(stats.c),1);

for trial = 1:length(stats.c)
    stats.currTrial = trial;

    if stats.currTrial == 1  %if this is the first trial
        stats.qA(stats.currTrial,1:2) = init_v.Q(1,:);
        stats.qB(stats.currTrial,1:2) = init_v.Q(2,:);
        stats.ck(stats.currTrial,:) = init_v.C;
        stats.rpe(stats.currTrial) = NaN;

        deltaQA = stats.qA(stats.currTrial,1) - stats.qA(stats.currTrial,2);
        deltaQB = stats.qB(stats.currTrial,1) - stats.qB(stats.currTrial,2);
        deltaCK = stats.ck(stats.currTrial,1)-stats.ck(stats.currTrial,2);


    else
        %% update choice kernel
        if stats.c(stats.currTrial-1)==1   % else, chose right
            stats.ck(stats.currTrial,2) = stats.ck(stats.currTrial-1,2)+ alpha_c * (1-stats.ck(stats.currTrial-1,2));
            stats.ck(stats.currTrial,1) = (1-alpha_c)*stats.ck(stats.currTrial-1,1);
        elseif stats.c(stats.currTrial-1)==0
            stats.ck(stats.currTrial,2) = (1-alpha_c)*stats.ck(stats.currTrial-1,2);
            stats.ck(stats.currTrial,1) = stats.ck(stats.currTrial-1,1) + alpha_c * (1-stats.ck(stats.currTrial-1,1)); % left
        else
            stats.ck(stats.currTrial,:) = stats.ck(stats.currTrial-1,:);

        end

        %% update action values
        updateState = stats.s(stats.currTrial-1);
        if ismember(updateState, [1,3,5])
            stats.qB(stats.currTrial,:) =(1-a)* stats.qB(stats.currTrial-1,:);
            if stats.c(stats.currTrial-1)==0   % if chose left on last trial
                stats.rpe(stats.currTrial-1)=stats.r(stats.currTrial-1)-stats.qA(stats.currTrial-1,1);
                stats.qA(stats.currTrial,1)=stats.qA(stats.currTrial-1,1)+a*stats.rpe(stats.currTrial-1);
                stats.qA(stats.currTrial,2)=(1-a)*stats.qA(stats.currTrial-1,2);
            elseif stats.c(stats.currTrial-1)==1   % else, chose right
                stats.rpe(stats.currTrial-1)=stats.r(stats.currTrial-1)-stats.qA(stats.currTrial-1,2);
                stats.qA(stats.currTrial,1)=(1-a)*stats.qA(stats.currTrial-1,1);
                stats.qA(stats.currTrial,2)=stats.qA(stats.currTrial-1,2)+a*stats.rpe(stats.currTrial-1);
            else  %miss trials
                stats.rpe(stats.currTrial-1) = 0;
                stats.qA(stats.currTrial,1)=(1-a)*stats.qA(stats.currTrial-1,1);
                stats.qA(stats.currTrial,2)=(1-a)*stats.qA(stats.currTrial-1,2);
            end
        elseif ismember(updateState,[2,4,6])
            stats.qA(stats.currTrial,:) = (1-a)*stats.qA(stats.currTrial-1,:);
            if stats.c(stats.currTrial-1)==0   % if chose left on last trial
                stats.rpe(stats.currTrial-1)=stats.r(stats.currTrial-1)-stats.qB(stats.currTrial-1,1);
                stats.qB(stats.currTrial,1)=stats.qB(stats.currTrial-1,1)+a*stats.rpe(stats.currTrial-1);
                stats.qB(stats.currTrial,2)=(1-a)*stats.qB(stats.currTrial-1,2);
            elseif stats.c(stats.currTrial-1)==1   % else, chose right
                stats.rpe(stats.currTrial-1)=stats.r(stats.currTrial-1)-stats.qB(stats.currTrial-1,2);
                stats.qB(stats.currTrial,1)=(1-a)*stats.qB(stats.currTrial-1,1);
                stats.qB(stats.currTrial,2)=stats.qB(stats.currTrial-1,2)+a*stats.rpe(stats.currTrial-1);
            else  %miss trials
                stats.rpe(stats.currTrial-1) = 0;
                stats.qB(stats.currTrial,1)=(1-a)*stats.qB(stats.currTrial-1,1);
                stats.qB(stats.currTrial,2)=(1-a)*stats.qB(stats.currTrial-1,2);
            end
        end

        %% softmax with choice kernel
        %% softmax rule for action selection
        deltaQA = stats.qA(stats.currTrial,1) - stats.qA(stats.currTrial,2);
        deltaQB = stats.qB(stats.currTrial,1) - stats.qB(stats.currTrial,2);
        deltaCK = stats.ck(stats.currTrial,1) - stats.ck(stats.currTrial,2);

    end

    currState =  stats.s(stats.currTrial);
    if ismember(currState, [1,3])
        stats.pCorrect(stats.currTrial)  = (1/(1+exp(-(b*deltaQA+beta_c*deltaCK))));
    elseif currState == 5
        stats.pCorrect(stats.currTrial)  = 1-(1/(1+exp(-(b*deltaQA+beta_c*deltaCK))));
    elseif ismember(currState, [2,4])
        stats.pCorrect(stats.currTrial) = 1-(1/(1+exp(-(b*deltaQB+beta_c*deltaCK))));
    elseif currState==6
        stats.pCorrect(stats.currTrial) = 1/(1+exp(-(b*deltaQB+beta_c*deltaCK)));
    end
end
end