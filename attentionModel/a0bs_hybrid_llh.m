function llh = a0bs_hybrid_llh(theta,data)

% use choice kernel instead of stickness

beta = theta(1);
alpha = [1e-6 theta(2)];%.2;
stick = theta(3);
epsilon = 1e-6;%theta(7);
lapse = theta(4);%theta(2);%.05;
ret = theta(5);%theta(3);%.1;
bias = theta(6);

% state 2: RL
% state 1: random
% lapse: probability of lapsing from 2 to 1
% ret: probability of returning from random to RL

T = [1-ret lapse;ret 1-lapse];

stimuli = data(:,1);
choices = data(:,2);
rewards = data(:,3);

Q = [.5 .5;.5 .5];
lt = log(.5);
llh = lt;
p = [lapse 1-lapse];
s = stimuli(1);
b(2) = epsilon/2 + (1-epsilon)/(1+exp(beta*(Q(s,1)-Q(s,2))));
b(1) = 1-b(2);
b1 = [1-bias bias];

firstS = stimuli(1);

for k = 2:length(choices)
    s = stimuli(k-1);
    choice = choices(k-1);
    r = rewards(k-1);
    p = (b1(choice)*p(1)*T(:,1) + b(choice)*p(2)*T(:,2)) / exp(lt);
    
    Q(s,choice) = Q(s,choice) + alpha(r+1)*(r-Q(s,choice));
    %CK(choice) = CK(choice) + alpha_CK*(1-CK(choice));
%     side=zeros(4,1);
%     if stimuli(k-1) == stimuli(k) && rewards(k-1) == 0
%         side(1) = 2 * (1.5 - choices(k-1)); % 1 for A1 and -1 for A2
%     end
%     if stimuli(k-1) ~= stimuli(k) && rewards(k-1) == 0
%         side(2) = 2 * (1.5 - choices(k-1));
%     end
%     if stimuli(k-1) == stimuli(k) && rewards(k-1) > 0
%         side(3) = 2 * (1.5 - choices(k-1));
%     end
%     if stimuli(k-1) ~= stimuli(k) && rewards(k-1) > 0
%         side(4) = 2 * (1.5 - choices(k-1));
%     end
    if choices(k-1) == 2 % if previously choose right
        I = 1;
    else
        I = -1;
    end
    b(2) = epsilon/2 + (1-epsilon)/(1+exp(beta*(Q(stimuli(k),1)-Q(stimuli(k),2))-stick*I));
    b(1) = 1-b(2);
    mQ = mean(Q);
    b1(2) = bias;
    b1(1) = 1-b1(2);
       
    lt = log(b1(choices(k))*p(1) + b(choices(k))*p(2));
    llh = llh + lt;
end

llh = -llh;

end