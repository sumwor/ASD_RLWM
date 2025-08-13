function entropy_rate = get_stimulus_entropy(seq)

% entropy rate via Markov chain

if ismember(3, seq) & ~ismember(5, seq)
    seq = seq-2;
    seq = seq(seq>0);
elseif ismember(5,seq)
    seq = seq-4;
    seq = seq(seq>0);
end
% Count transitions
transitions = zeros(2,2);
for i = 1:length(seq)-1
    transitions(seq(i), seq(i+1)) = transitions(seq(i), seq(i+1)) + 1;
end

% Transition probabilities (conditional probabilities)
trans_probs = transitions ./ sum(transitions, 2);

% Joint probabilities of states (estimate stationary distribution)
state_counts = histcounts(seq, 0.5:1:2.5, 'Normalization','probability');

% Calculate entropy rate
entropy_rate = 0;
for i = 1:2
    for j = 1:2
        if trans_probs(i,j) > 0
            entropy_rate = entropy_rate - state_counts(i) * trans_probs(i,j) * log2(trans_probs(i,j));
        end
    end
end

%disp(['Entropy rate (bits per symbol): ', num2str(entropy_rate)]);