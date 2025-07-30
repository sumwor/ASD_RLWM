function varargout=funFQ_RPE_CK_SDT(xpar,dat,init_v)
% % funFQ_RPE_CK %
%PURPOSE:   Function for maximum likelihood estimation, called by fit_fun().
%
%INPUT ARGUMENTS
%   xpar:       alpha, beta, alpha CK, beta CK, 
%   stats:      info about animal/agent's performance
%
%OUTPUT ARGUMENTS
%   negloglike:      the negative log-likelihood to be minimized
% updated: 04/30/2020
%%
if nargin == 0 % Used for optimizer testing
    varargout{1} = [0,0,0];               % LB
    varargout{2} = [1,10,2];             % UB
    varargout{3} = [0,3000];                % Zlim
    varargout{4} = [-30 75];                % view
    varargout{5} = [0.3 5 0.2 3];           % xmin
    varargout{6} = 'FQ_RPE_SDT function';    % name
    varargout{7} = [0.3 5 0.2 3];           % default x0
    return;
end

alpha=xpar(1);
beta=xpar(2);
alpha_c = xpar(3);
beta_c  = xpar(4);
sigma_sensory = xpar(5);

nt=size(dat,1);
ChoiceProb = nan(nt,1);
%CK =[0 0];  % right, left choice kernel, will be updated in each trial.

Q = init_v.Q; % 2x2 matrix for Q(A,L), Q(A,R), Q(B,L), Q(B,R)
CK = init_v.C; % 2x2 matrix for CK(L), CK(R)

for k=1:nt
    % generate noise
    stimulus = dat(k, 3);
    if ismember(stimulus, [1,3,5])
        O = sigma_sensory*randn();
    elseif ismember(stimulus, [2,4,6])
        O = 1+ sigma_sensory*randn();
    end
    
    
    % calculate posterior distribution P(A|O) and P(B|O)
    PA_over_PB = normpdf(O, 0, sigma_sensory)/normpdf(O, 1, sigma_sensory);
    p_B = 1/(1+PA_over_PB);
    p_A = PA_over_PB/(1+PA_over_PB);

    % given P_A and p_B, generating p_left with softmax funtion
    %% softmax with choice kernel
    Q_A = Q(1,:);
    Q_B = Q(2,:);
    deltaQA = Q(1,1)-Q(1,2);
    deltaQB = Q(2,1)-Q(2,2);


    deltaCK = CK(1)-CK(2);

    %V = (beta*Q) + (beta_c*CK);
    
    pleft  = p_A*(1/(1+exp(-(beta*deltaQA+beta_c*deltaCK)))) + p_B*(1/(1+exp(-(beta*deltaQB+beta_c*deltaCK))));
    pright = 1-pleft;
    
    if pright==0
        pright=realmin;   % Smallest positive normalized floating point number, because otherwise log(zero) is -Inf
    end
    if pleft==0
        pleft=realmin;
    end
    
    % calculate choice probability for actual choice
    if dat(k,1)==1
        ChoiceProb(k) = pright;
    elseif dat(k,1)==0
        ChoiceProb(k) = pleft;
    else
        ChoiceProb(k) = nan;
    end
    
    %% choice kernel update
    CK = (1-alpha_c) * CK;
    if dat(k,1)==1    % if chosen right
        CK(1) = CK(1)+ (alpha_c *1);
    elseif dat(k,1)==0
        CK(2) = CK(2)+ (alpha_c *1); % left
    else % for miss - no action, no update
        CK = CK;
    end
    
    %% update value for the performed action
    if dat(k,3) == 1  % odor A
        if dat(k,1)==0      %chose left
            Q(1,1) =Q(1,1)+alpha*(dat(k,2)-Q(1,1));
            Q(1,2)  =(1-alpha)*Q(1,2);
            Q(2,1) = (1-alpha)*Q(2,1);
            Q(2,2) = (1-alpha)*Q(2,2);
        elseif dat(k,1)==1 %chose right
            Q(1,2) =Q(1,2)+alpha*(dat(k,2)-Q(1,2));
            Q(1,1)  =(1-alpha)*Q(1,1);
            Q(2,1) = (1-alpha)*Q(2,1);
            Q(2,2) = (1-alpha)*Q(2,2);
        end
    elseif dat(k,3) == 2 % odor B
        if dat(k,1)==0      %chose left
            Q(2,1) =Q(2,1)+alpha*(dat(k,2)-Q(2,1));
            Q(1,2)  =(1-alpha)*Q(1,2);
            Q(1,1) = (1-alpha)*Q(1,1);
            Q(2,2) = (1-alpha)*Q(2,2);
        elseif dat(k,1)==1 %chose right
            Q(2,2) =Q(2,2)+alpha*(dat(k,2)-Q(2,2));
            Q(1,1)  =(1-alpha)*Q(1,1);
            Q(2,1) = (1-alpha)*Q(2,1);
            Q(1,2) = (1-alpha)*Q(1,2);
    end
end
negloglike = - sum(log(ChoiceProb(~isnan(ChoiceProb))));
varargout{1} = negloglike;
end