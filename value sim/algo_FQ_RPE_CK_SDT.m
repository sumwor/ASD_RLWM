function stats=algo_FQ_RPE_CK_SDT(stats,params,init_v)
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
sigma = params(5);

stats.qA = nan(length(stats.c),2);
stats.qB = nan(length(stats.c),2);
stats.ck = nan(length(stats.c),2);
stats.pCorrect = nan(length(stats.c),1);
stats.rpe = nan(length(stats.c),1);
stats.perceiveds = nan(length(stats.c),1);

for trial = 1:length(stats.c)
    stats.currTrial = trial;

    if stats.currTrial == 1  %if this is the first trial
        stats.qA(stats.currTrial,1:2) = init_v.Q(1,:);
        stats.qB(stats.currTrial,1:2) = init_v.Q(2,:);
        stats.ck(stats.currTrial,:) = init_v.C;
        stats.rpe(stats.currTrial) = NaN;

        % determine the perceived state
        currTrueState = stats.s(stats.currTrial);
        if ismember(currTrueState, [1,3,5])
            O = 0 + sigma*randn();
        elseif ismember(currTrueState, [2,4,6])
            O = 1 + sigma*randn();
        end

        PA_over_PB = normpdf(O, 0, sigma)/normpdf(O, 1, sigma);
        p_B = 1/(1+PA_over_PB);
        p_A = PA_over_PB/(1+PA_over_PB);

        % determine the perceived state
        if rand() <= p_A %
            currPerceivedState = 1;
        else
            currPerceivedState = 2;
        end

        stats.perceiveds(stats.currTrial) = currPerceivedState;

        deltaQA = stats.qA(stats.currTrial,1) - stats.qA(stats.currTrial,2);
        deltaQB = stats.qB(stats.currTrial,1) - stats.qB(stats.currTrial,2);
         deltaCK = stats.ck(stats.currTrial,1) - stats.ck(stats.currTrial,2);

       
        if ismember(currTrueState, [1,3,6])
            stats.pCorrect(stats.currTrial)  =  p_A*(1/(1+exp(-(b*deltaQA+beta_c*deltaCK))))+ p_B*(1/(1+exp(-(b*deltaQB+beta_c*deltaCK))));
        elseif ismember(currTrueState, [2,4,5]) 
            stats.pCorrect(stats.currTrial)  = p_A*(1-1/(1+exp(-(b*deltaQA+beta_c*deltaCK))))+ p_B*(1-1/(1+exp(-(b*deltaQB+beta_c*deltaCK))));
        end


    else
        %% update choice kernels
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
        
        updateState = stats.perceiveds(stats.currTrial-1);
        if updateState == 1
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
        elseif updateState == 2
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
        currTrueState = stats.s(stats.currTrial);
        if ismember(currTrueState, [1,3,5])
            O = 0 + sigma*randn();
        elseif ismember(currTrueState, [2,4,6])
            O = 1 + sigma*randn();
        end


        PA_over_PB = normpdf(O, 0, sigma)/normpdf(O, 1, sigma);
        if PA_over_PB == inf
            p_A = 1-realmin;
        else
            p_A = PA_over_PB/(1+PA_over_PB);
        end
        p_B = 1/(1+PA_over_PB);


        % determine the perceived state
        if rand() <= p_A %
            currPerceivedState = 1;
        else
            currPerceivedState = 2;
        end

        stats.perceiveds(stats.currTrial) = currPerceivedState;

        deltaQA = stats.qA(stats.currTrial,1) - stats.qA(stats.currTrial,2);
        deltaQB = stats.qB(stats.currTrial,1) - stats.qB(stats.currTrial,2);
        deltaCK = stats.ck(stats.currTrial,1) - stats.ck(stats.currTrial,2);

       
        if ismember(currTrueState, [1,3,6])
            stats.pCorrect(stats.currTrial)  =  p_A*(1/(1+exp(-(b*deltaQA+beta_c*deltaCK))))+ p_B*(1/(1+exp(-(b*deltaQB+beta_c*deltaCK))));
        elseif ismember(currTrueState, [2,4,5]) 
            stats.pCorrect(stats.currTrial)  = p_A*(1-1/(1+exp(-(b*deltaQA+beta_c*deltaCK))))+ p_B*(1-1/(1+exp(-(b*deltaQB+beta_c*deltaCK))));
        end

    end
end
end