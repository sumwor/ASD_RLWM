function llh = a0bsnoise_hybrid_llh(theta,data)

% use stickness with noise sensory input

beta = theta(1);
alpha = [1e-6 theta(2)];%.2;
stick = theta(3);
epsilon = 1e-6;%theta(7);
lapse = theta(4);%theta(2);%.05;
ret = theta(5);%theta(3);%.1;
bias = theta(6);
sigma = theta(7);

nSamples = 50;
% state 2: RL
% state 1: random
% lapse: probability of lapsing from 2 to 1
% ret: probability of returning from random to RL
T = [1-ret lapse; ret 1-lapse];

stimuli = data(:,1);
choices = data(:,2);
rewards = data(:,3);

Q = [.5 .5; .5 .5];
lt = log(.5);

llh_total = 0; % accumulate likelihood over samples

for sampleIdx = 1:nSamples
    p = [lapse 1-lapse];
    s = stimuli(1);
    b(2) = epsilon/2 + (1-epsilon)/(1+exp(beta*(Q(s,1)-Q(s,2))));
    b(1) = 1-b(2);
    b1 = [.5 .5];
    firstS = stimuli(1);
    lt_sample = lt;
    llh_total = llh_total+lt_sample;
for k = 2:length(choices)
    s = stimuli(k-1);
    
    %% perceptual noise

    if s==1
        O = sigma*randn();
    elseif s==2
        O = 1+ sigma*randn();
    end
        % calculate posterior distribution P(A|O) and P(B|O)
logPA = log(normpdf(O,0,sigma));
logPB = log(normpdf(O,1,sigma));
log_ratio = logPA - logPB;   % log( PA / PB )
PA_over_PB = exp(log_ratio);

p_B = 1 ./ (1 + PA_over_PB);
p_A = 1 ./ (1 + exp(-log_ratio));

%     PA_over_PB = normpdf(O, 0, sigma)/normpdf(O, 1, sigma);
%     p_B = 1/(1+PA_over_PB);
%     p_A = PA_over_PB/(1+PA_over_PB);

    choice = choices(k-1);
    r = rewards(k-1);
    % probablity of engaged/disengaged state
    p = (b1(choice)*p(1)*T(:,1) + b(choice)*p(2)*T(:,2)) / exp(lt_sample);
    
    Q(s,choice) = Q(s,choice) + alpha(r+1)*(r-Q(s,choice));

    if choices(k-1) == 2 % if previously choose right
        I = 1;
    else
        I = -1;
    end
    
    % probability to choose right in perceived state A or perceived state B
    tempb_A = epsilon/2 + (1-epsilon)/(1+exp(beta*(Q(1,1)-Q(1,2))-stick*I));
    tempb_B = epsilon/2 + (1-epsilon)/(1+exp(beta*(Q(2,1)-Q(2,2))-stick*I));
    b(2) = p_A*tempb_A + p_B*tempb_B;
    b(1) = 1-b(2);
    mQ = mean(Q);
    b1(firstS) = bias;
    b1(3-firstS) = 1-b1(firstS);
    
    lt_sample = log(b1(choices(k))*p(1) + b(choices(k))*p(2));
    llh_total = llh_total + lt_sample;
end
end
llh = -llh_total / nSamples;

%llh = -llh;

end