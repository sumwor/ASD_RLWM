function error_rate = plot_errors(resultdf, protocol, edges, tlabel, savefigpath);

% plot the error rate as a function of previous consecutive choice

% separate AB-CD-DC part

prev_trials = 0:5;
err_choice = struct;
err_reinforce=struct;
err_choice.L = struct;
err_choice.R = struct;
err_reinforce.L = struct;
err_reinforce.R = struct;

schedules = {'AB', 'CD', 'DC'};

for pp = 1:length(schedules)

    err_choice.L.(schedules{pp}) = nan(length(prev_trials),2);
    err_choice.R.(schedules{pp}) = nan(length(prev_trials),2);
    err_reinforce.L.(schedules{pp}) = nan(length(prev_trials),2);
    err_reinforce.R.(schedules{pp}) = nan(length(prev_trials),2);
end

%% determine the current protocol
if strcmp(protocol,'AB')
    ABEnd = size(resultdf,1);
    CDEnd = NaN;
    DCEnd = NaN;
elseif strcmp(protocol,'AB-CD')
    ABEnd = find(list == 3 | list == 4, 1, 'first')-1;
    CDEnd = size(resultdf,1);
    DCEnd = NaN;
elseif strcmp(protocol, 'AB-CD-DC')
    ABEnd = find(list == 3 | list == 4, 1, 'first')-1;
    CDEnd = find(list == 5 | list == 6, 1, 'first')-1;
    DCEnd = size(resultdf,1);
elseif strcmp(protocol, 'AB-DC')
    ABEnd = find(list == 5 | list == 6, 1, 'first')-1;
    CDEnd = NaN;
    DCEnd = size(resultdf,1);
end

for ii = 1:length(prev_trials)
    % count trials and correct choice for previous choice

    for pp = 1:length(schedules)
        nCorrect= NaN; nTrials= NaN;
        if strcmp(schedules{pp}, 'AB')
            [nCorrect, nTrials] = get_err_rate(resultdf(1:ABEnd-1,:),ii);

        elseif strcmp(schedules{pp}, 'CD') && ~isnan(CDEnd)
            [nCorrect, nTrials] = get_err_rate(resultdf(ABEnd+1:CDEnd-1,:),ii);
        elseif strcmp(schedules{pp}, 'DC') && ~isnan(DCEnd)
            if isnan(CDEnd)
                [nCorrect, nTrials] = get_err_rate(resultdf(ABEnd+1:DCEnd-1,:),ii);
            else
                [nCorrect, nTrials] = get_err_rate(resultdf(CDEnd+1:DCEnd-1,:),ii);
            end
        end

        if isstruct(nCorrect)
            err_choice.L.(schedules{pp})(ii,1) = nTrials.L-nCorrect.L; err_choice.L.(schedules{pp})(ii,2) = nTrials.L;
            err_choice.R.(schedules{pp})(ii,1) = nTrials.R-nCorrect.R; err_choice.R.(schedules{pp})(ii,2) = nTrials.R;

            err_reinforce.L.(schedules{pp})(ii,1) = nTrials.L_r-nCorrect.L_r; err_reinforce.L.(schedules{pp})(ii,2) = nTrials.L_r;
            err_reinforce.R.(schedules{pp})(ii,1) = nTrials.R_r-nCorrect.R_r; err_reinforce.R.(schedules{pp})(ii,2) = nTrials.R_r;
        end
    end
end

figure;
for pp = 1:length(schedules)

    subplot(3,2,2*pp-1)
    plot(prev_trials, err_choice.L.(schedules{pp})(:,1)./err_choice.L.(schedules{pp})(:,2))
    hold on; plot(prev_trials, err_choice.R.(schedules{pp})(:,1)./err_choice.R.(schedules{pp})(:,2))
    set(gca, 'box', 'off');
    if pp==1
        title('Choice');
         ylabel('Error rate AB');
    elseif pp==2
        ylabel('Error rate CD');
    elseif pp==3
        ylabel('Error rate DC')
    end
   
    if pp==3
        xlabel('Number of consecutive previous choices')
    end
    ylim([0 1])

    subplot(3,2,2*pp)
    plot(prev_trials, err_choice.R.(schedules{pp})(:,1)./err_choice.R.(schedules{pp})(:,2))
    hold on; plot(prev_trials, err_reinforce.R.(schedules{pp})(:,1)./err_reinforce.R.(schedules{pp})(:,2))
    set(gca, 'box', 'off');
    if pp==3
    lgd = legend('Left', 'Right');
    set(lgd, 'Box', 'off');            % Remove the border box
    set(lgd, 'Color', 'none');
    end
    if pp==1
    title('Rewarded choice')
    end
    ylim([0 1])
end

    function [nCorrect, nTrials] = get_err_rate(result, nPrev)

    nCorrect = struct; nTrials = struct;
    nTrials.L = 0; nCorrect.L = 0;
    nTrials.R = 0; nCorrect.R = 0;
    nTrials.L_r = 0; nCorrect.L_r = 0;
    nTrials.R_r = 0; nCorrect.R_r = 0;

    for tt = (1+nPrev):size(result,1)-1
        prev_choice = result.actions(tt-nPrev:tt-1);
        prev_reward = result.reward(tt-nPrev:tt-1);
        if all(prev_choice == 0)
            nTrials.L = nTrials.L + 1;
            if result.reward(tt+1) > 0 && result.actions(tt+1) == 1
                nCorrect.L = nCorrect.L + 1;
            end
            if all(prev_reward>0)
                nTrials.L_r = nTrials.L_r + 1;
                if result.reward(tt+1) > 0 && result.actions(tt+1) == 1
                    nCorrect.L_r = nCorrect.L_r + 1;
                end
            end
        elseif all(prev_choice == 1)
            nTrials.R = nTrials.R + 1;
            if result.reward(tt+1) > 0 && result.actions(tt+1) == 0
                nCorrect.R = nCorrect.R + 1;
            end
            if all(prev_reward>0)
                nTrials.R_r = nTrials.R_r + 1;
                if result.reward(tt+1) > 0 && result.actions(tt+1) == 0
                    nCorrect.R_r = nCorrect.R_r + 1;
                end
            end
        end
    end

    end

end
