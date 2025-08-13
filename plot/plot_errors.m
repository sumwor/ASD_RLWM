function [error_rate, error_rate_reinforce] = plot_errors(resultdf, protocol, edges, tlabel, savefigpath);

% plot the error rate as a function of previous consecutive choice

% separate AB-CD-DC part

prev_trials = 0:5;
% left and right error rate

err_choice = struct;
err_reinforce=struct;
err_choice.L = struct;
err_choice.R = struct;
err_reinforce.L = struct;
err_reinforce.R = struct;

schedules = {'AB', 'CD', 'DC'};
error_rate = struct;
error_rate_reinforce = struct;
for pp = 1:length(schedules)
    error_rate.(schedules{pp}) = nan(length(prev_trials),2);
    error_rate_reinforce.(schedules{pp}) = nan(length(prev_trials),2);
end

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
    ABEnd = find(resultdf.schedule == 3 | resultdf.schedule == 4, 1, 'first')-1;
    CDEnd = size(resultdf,1);
    DCEnd = NaN;
elseif strcmp(protocol, 'AB-CD-DC')
    ABEnd = find(resultdf.schedule == 3 | resultdf.schedule == 4, 1, 'first')-1;
    CDEnd = find(resultdf.schedule == 5 | resultdf.schedule == 6, 1, 'first')-1;
    DCEnd = size(resultdf,1);
elseif strcmp(protocol, 'AB-DC')
    ABEnd = find(resultdf.schedule == 5 | resultdf.schedule == 6, 1, 'first')-1;
    CDEnd = NaN;
    DCEnd = size(resultdf,1);
end

for ii = 1:length(prev_trials)
    % count trials and correct choice for previous choice

    for pp = 1:length(schedules)
        nError= NaN; nTrials= NaN;
        if strcmp(schedules{pp}, 'AB')
            [nError, nTrials] = get_err_rate(resultdf(1:ABEnd-1,:),ii-1);

        elseif strcmp(schedules{pp}, 'CD') && ~isnan(CDEnd)
            [nError, nTrials] = get_err_rate(resultdf(ABEnd+1:CDEnd-1,:),ii-1);
        elseif strcmp(schedules{pp}, 'DC') && ~isnan(DCEnd)
            if isnan(CDEnd)
                [nError, nTrials] = get_err_rate(resultdf(ABEnd+1:DCEnd-1,:),ii-1);
            else
                [nError, nTrials] = get_err_rate(resultdf(CDEnd+1:DCEnd-1,:),ii-1);
            end
        end

        if isstruct(nError)
            err_choice.L.(schedules{pp})(ii,1) = nError.L; err_choice.L.(schedules{pp})(ii,2) = nTrials.L;
            err_choice.R.(schedules{pp})(ii,1) = nError.R; err_choice.R.(schedules{pp})(ii,2) = nTrials.R;

            err_reinforce.L.(schedules{pp})(ii,1) = nError.L_r; err_reinforce.L.(schedules{pp})(ii,2) = nTrials.L_r;
            err_reinforce.R.(schedules{pp})(ii,1) = nError.R_r; err_reinforce.R.(schedules{pp})(ii,2) = nTrials.R_r;
            error_rate.(schedules{pp})(ii,1) =  nError.L/nTrials.L;
            error_rate.(schedules{pp})(ii,2) = nError.R/nTrials.R;
            error_rate_reinforce.(schedules{pp})(ii,1) = nError.L_r/nTrials.L;
            error_rate_reinforce.(schedules{pp})(ii,2) = nError.R_r/nTrials.R;

        end

    end
end

figure;
for pp = 1:length(schedules)

    subplot(3,2,2*pp-1)
    plot(prev_trials, err_choice.L.(schedules{pp})(:,1)./err_choice.L.(schedules{pp})(:,2))
    hold on; plot(prev_trials, err_reinforce.L.(schedules{pp})(:,1)./err_reinforce.L.(schedules{pp})(:,2))
    set(gca, 'box', 'off');
    if pp==1
        title('Left');
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
        lgd = legend('Choice', 'Rewarded choice');
        set(lgd, 'Box', 'off');            % Remove the border box
        set(lgd, 'Color', 'none');
    end
    if pp==1
        title('Right')
    end
    ylim([0 1])
end

% save the plot
print(gcf,'-dpng',fullfile(savefigpath, ['error_rate_', tlabel]));    %png format
%saveas(gcf, fullfile(BehPath, ['error_rate_', tlabel]), 'fig');
%saveas(gcf, fullfile(BehPath, ['error_rate_', tlabel]),'svg');

close()


    function [nError, nTrials] = get_err_rate(result, nPrev)

        nError = struct; nTrials = struct;
        nTrials.L = 0; nError.L = 0;
        nTrials.R = 0; nError.R = 0;
        nTrials.L_r = 0; nError.L_r = 0;
        nTrials.R_r = 0; nError.R_r = 0;

        if nPrev == 0
            nTrials.L = sum(result.actions==0);
            nTrials.R = sum(result.actions==1);
            nError.L = sum(result.actions==0 & isnan(result.reward));
            nError.R = sum(result.actions==1 & isnan(result.reward));
            nError.L_r = NaN; nError.R_r = NaN;
        else
            for tt = (1+nPrev):size(result,1)-nPrev
                prev_choice = result.actions(tt-nPrev:tt-1);
                prev_reward = result.reward(tt-nPrev:tt-1);
                if all(prev_choice == 0)
                    nTrials.L = nTrials.L + 1;
                    if isnan(result.reward(tt+1)) && result.actions(tt+1) == 0
                        nError.L = nError.L + 1;
                    end
                    if all(prev_reward>0)
                        nTrials.L_r = nTrials.L_r + 1;
                        if isnan(result.reward(tt+1)) && result.actions(tt+1) == 0
                            nError.L_r = nError.L_r + 1;
                        end
                    end
                elseif all(prev_choice == 1)
                    nTrials.R = nTrials.R + 1;
                    if isnan(result.reward(tt+1)) && result.actions(tt+1) == 1
                        nError.R = nError.R + 1;
                    end
                    if all(prev_reward>0)
                        nTrials.R_r = nTrials.R_r + 1;
                        if isnan(result.reward(tt+1)) && result.actions(tt+1) == 1
                            nError.R_r = nError.R_r + 1;
                        end
                    end
                end
            end
        end

    end

end
