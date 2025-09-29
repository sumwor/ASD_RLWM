function neglogpost = actor_critic_llh(theta, data)
% params = [alpha_actor, alpha_critic]
alpha_actor_stimuli  = theta(1); % constrain to (0,1) via logistic
alpha_actor_stick = theta(2);
alpha_actor_bias = theta(3);
alpha_critic = theta(4);
%bias = theta(5);
if length(theta) == 5
    forget = [theta(5),1;
        1, theta(5)];
else
    forget = 1;
end
% fixed policy inverse temperature
nTrials = length(data.s);

stimuli = data.s;
choices = data.c;
rewards = data.r;

Q = (data.Q-0.5).*forget+0.5;

% Initialize
W = zeros(3,1);
prev_choice = 0;

loglik = 0;

for t = 1:nTrials
    stim = (stimuli(t)-1.5)*2;        % -1 or 1
    state = stimuli(t);      % 1=odorA, 2=odorB
    action = choices(t);         % 1=Left, 2=Right
    reward = rewards(t);         % 0/1

    % Policy
    x = [stim; prev_choice; 1];
    logit = W' * x;
    pR = 1 / (1 + exp(logit));
    pL = 1 - pR;

    if action == 2
        p = pR;
    else
        p = pL;
    end
    loglik = loglik + log(p + eps);  % accumulate log-likelihood

    % Critic update
    if reward == 1
        Q(state, action) = Q(state, action) + ...
            alpha_critic * (reward - Q(state, action));
    end

    % Actor update
    delta = reward - Q(state, action);
    if action == 2
        grad = (1 - pR) * x;
    else
        grad = -pR * x;
    end
    alpha_vec = [alpha_actor_stimuli; alpha_actor_stick; alpha_actor_bias];
    W = W + alpha_vec .* (delta * grad);
   

    % Stickiness update
    if reward == 1
        if action == 1
            prev_choice = 1;
        else
            prev_choice = -1;
        end
    else
        prev_choice = 0;
    end
end

neglogpost = -loglik;
end