function mixed_logreg_RCUC_indOdor(dataIndex,tlabel, nSession, step_back, savefigpath, savedatapath)

% mixed logstic regression
% tlabel: odor pairs to checkk
% nSession: which session to look (1,2,3 for AB/AB-CD; 12345 for AB-DC)
% separate odors

setup_figprop;
odor_list={'A', 'B', 'C', 'D', 'C', 'D'};
% Preallocate
allData = [];
switch tlabel
    case 'AB'
        subDataIndex = dataIndex(strcmp(dataIndex.Protocol, 'AB') & (cellfun(@(x) x == nSession, dataIndex.ProtocolDay)),:);
        odor_pair = [1,2];
    case 'AB-CD-AB'
        subDataIndex = dataIndex(strcmp(dataIndex.Protocol,'AB-CD') & (cellfun(@(x) x == nSession, dataIndex.ProtocolDay)),:);
        odor_pair = [1,2];
    case 'AB-CD'
        subDataIndex = dataIndex(strcmp(dataIndex.Protocol,'AB-CD') & (cellfun(@(x) x == nSession, dataIndex.ProtocolDay)),:);
        odor_pair = [3,4];
    case 'AB-DC-AB'
        subDataIndex = dataIndex(strcmp(dataIndex.Protocol,'AB-DC') & (cellfun(@(x) x == nSession, dataIndex.ProtocolDay)),:);
        odor_pair = [1,2];
    case 'AB-DC'
        subDataIndex = dataIndex(strcmp(dataIndex.Protocol,'AB-DC') & (cellfun(@(x) x == nSession, dataIndex.ProtocolDay)),:);
        odor_pair = [5,6];

        
end

maxLag = step_back;
nSessions = size(subDataIndex,1);

for oo = 1:length(odor_pair)

    for s = 1:nSessions
        csvfilepath = fullfile(subDataIndex.BehPath{s}, subDataIndex.BehCSV{s});
        resultdf = readtable(csvfilepath);
        protocol = subDataIndex.Protocol{s};
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

        switch tlabel
            case {'AB', 'AB-CD-AB', 'AB-DC-AB'}
                actions = resultdf.actions(1:ABEnd); % column vector
                rewards = resultdf.reward(1:ABEnd);
                schedule = resultdf.schedule(1:ABEnd);
            case 'AB-CD'
                actions = resultdf.actions(ABEnd+1:CDEnd); % column vector
                rewards = resultdf.reward(ABEnd+1:CDEnd);
                schedule = resultdf.schedule(ABEnd+1:CDEnd);
            case {'AB-DC'}
                actions = resultdf.actions(ABEnd+1:DCEnd); % column vector
                rewards = resultdf.reward(ABEnd+1:DCEnd);
                schedule = resultdf.schedule(ABEnd+1:DCEnd);
        end


        currOdor = odor_pair(oo);
        temp_r = rewards;
        rewards(isnan(temp_r))=0;
        rewards(~isnan(temp_r)) = 1;

        temp_c = actions;
        actions(actions==0) = -1;

        rewards = rewards(schedule==currOdor);
        actions = actions(schedule==currOdor);
        geno = subDataIndex.Genotype{s};
        nTrials = length(actions);

        % Skip sessions with too few trials
        if nTrials <= maxLag
            continue;
        end

        % Build predictors



        %% create regressor vectors: rewarded/right = 1; rewarded/left = -1

        YR=zeros(length(actions),1);
        YR=rewards.*((actions==1) & (rewards>0)) + (-1*rewards).*((actions==-1) & (rewards>0));  %allows for reward size r!=1
        % YR=1*((c==1) & (r==1)) + (-1).*((c==-1) & (r==1));

        NR=zeros(length(actions),1);   % unrewarded/right = 1; unrewarded/left = -1
        NR=1*((actions==1) & (rewards==0)) + (-1)*((actions==-1) & (rewards==0));

        rmat=zeros(length(NR)-step_back,2*step_back);
        for i=1+step_back:length(NR)
            for j=1:step_back
                rmat(i-step_back,j)=YR(i-j);
                rmat(i-step_back,j+step_back)=NR(i-j);
            end
        end


        %allData = [allData; row]; % append row
        currChoice = actions(1+maxLag:end);

        % Store in struct array for later table conversion
        for t = 1:length(currChoice)
            if ~isnan(currChoice(t))
                row = struct();
                row.Choice = (currChoice(t)+1)/2;
                row.Session = s;
                row.Genotype = geno;
                for lag = 1:maxLag
                    row.(sprintf('Rwd_m%d',lag)) = rmat(t, lag);
                    row.(sprintf('Unrwd_m%d',lag)) = rmat(t, lag+maxLag);
                end
                allData = [allData; row];
            end
        end
    end
        % Convert to table
        tbl = struct2table(allData);

        % Make Session a categorical variable
        tbl.Session = categorical(tbl.Session);
        tbl.Genotype = categorical(tbl.Genotype);
        otherGenos = setdiff(categories(tbl.Genotype), {'WT'}, 'stable')';  % force row vector
        tbl.Genotype = setcats(tbl.Genotype, [{'WT'}, otherGenos]);

        % Define formula
        predictorTerms = [];
        for lag = 1:maxLag
            predictorTerms = [predictorTerms, ...
                sprintf(' + Rwd_m%d*Genotype + Unrwd_m%d*Genotype', lag, lag)];
        end

        formula = ['Choice ~ 1 + Genotype' predictorTerms ' + (1|Session)'];
        % Fit mixed-effects logistic regression
        glme = fitglme(tbl, formula, ...
            'Distribution', 'Binomial', 'Link', 'Logit');

        % View results
        disp(glme);

        coefTable = glme.Coefficients;
        coefNames = coefTable.Name;
        coefVals  = coefTable.Estimate;
        coefSE    = coefTable.SE;        % Standard error
        coefPvals = coefTable.pValue;

        genoLevels = categories(tbl.Genotype);
        nGeno = numel(genoLevels);

        lags = 0:maxLag; % include bias at lag=0

        % Preallocate
        betas_Rwd   = zeros(nGeno, maxLag+1);
        betas_Unrwd = zeros(nGeno, maxLag+1);
        SE_Rwd      = zeros(nGeno, maxLag+1);
        SE_Unrwd    = zeros(nGeno, maxLag+1);
        pvals_Rwd   = ones(nGeno, maxLag+1);
        pvals_Unrwd = ones(nGeno, maxLag+1);

        % Extract intercept (bias) and genotype effect SEs
        interceptIdx = strcmp(coefNames, '(Intercept)');
        intercept = coefVals(interceptIdx);
        interceptSE = coefSE(interceptIdx);
        interceptP = coefPvals(interceptIdx);

        for g = 1:nGeno
            thisGeno = genoLevels{g};
            if g == 1
                betas_Rwd(g,1) = intercept;
                betas_Unrwd(g,1) = intercept;
                SE_Rwd(g,1) = interceptSE;
                SE_Unrwd(g,1) = interceptSE;
                pvals_Rwd(g,1) = interceptP;
                pvals_Unrwd(g,1) = interceptP;
            else
                genoTerm = sprintf('Genotype_%s', thisGeno);
                idx = strcmp(coefNames, genoTerm);
                if any(idx)
                    genoEffect = coefVals(idx);
                    genoSE = coefSE(idx);
                    genoP = coefPvals(idx);
                    betas_Rwd(g,1) = intercept + genoEffect;
                    betas_Unrwd(g,1) = intercept + genoEffect;
                    % Approximate SE by sum in quadrature (assuming independence)
                    SE_Rwd(g,1) = sqrt(interceptSE^2 + genoSE^2);
                    SE_Unrwd(g,1) = sqrt(interceptSE^2 + genoSE^2);
                    pvals_Rwd(g,1) = max(interceptP, genoP);
                    pvals_Unrwd(g,1) = max(interceptP, genoP);
                else
                    betas_Rwd(g,1) = intercept;
                    betas_Unrwd(g,1) = intercept;
                    SE_Rwd(g,1) = interceptSE;
                    SE_Unrwd(g,1) = interceptSE;
                    pvals_Rwd(g,1) = interceptP;
                    pvals_Unrwd(g,1) = interceptP;
                end
            end
        end

        % Extract betas, SEs, and pvals for lags 1:maxLag
        for g = 1:nGeno
            thisGeno = genoLevels{g};
            for lag = 1:maxLag
                baseRwdIdx = strcmp(coefNames, sprintf('Rwd_m%d', lag));
                baseUnrwdIdx = strcmp(coefNames, sprintf('Unrwd_m%d', lag));

                baseRwd = coefVals(baseRwdIdx);
                baseUnrwd = coefVals(baseUnrwdIdx);
                baseRwdSE = coefSE(baseRwdIdx);
                baseUnrwdSE = coefSE(baseUnrwdIdx);
                baseRwdP = coefPvals(baseRwdIdx);
                baseUnrwdP = coefPvals(baseUnrwdIdx);

                if g == 1
                    betas_Rwd(g, lag+1) = baseRwd;
                    betas_Unrwd(g, lag+1) = baseUnrwd;
                    SE_Rwd(g, lag+1) = baseRwdSE;
                    SE_Unrwd(g, lag+1) = baseUnrwdSE;
                    pvals_Rwd(g, lag+1) = baseRwdP;
                    pvals_Unrwd(g, lag+1) = baseUnrwdP;
                else
                    interRwdIdx = strcmp(coefNames, sprintf('Genotype_%s:Rwd_m%d', thisGeno, lag));
                    interUnrwdIdx = strcmp(coefNames, sprintf('Genotype_%s:Unrwd_m%d', thisGeno, lag));

                    interRwd = 0; interUnrwd = 0;
                    interRwdSE = 0; interUnrwdSE = 0;
                    interRwdP = 1; interUnrwdP = 1;

                    if any(interRwdIdx)
                        interRwd = coefVals(interRwdIdx);
                        interRwdSE = coefSE(interRwdIdx);
                        interRwdP = coefPvals(interRwdIdx);
                    end
                    if any(interUnrwdIdx)
                        interUnrwd = coefVals(interUnrwdIdx);
                        interUnrwdSE = coefSE(interUnrwdIdx);
                        interUnrwdP = coefPvals(interUnrwdIdx);
                    end

                    betas_Rwd(g, lag+1) = baseRwd + interRwd;
                    betas_Unrwd(g, lag+1) = baseUnrwd + interUnrwd;

                    % SE sum in quadrature
                    SE_Rwd(g, lag+1) = sqrt(baseRwdSE^2 + interRwdSE^2);
                    SE_Unrwd(g, lag+1) = sqrt(baseUnrwdSE^2 + interUnrwdSE^2);

                    pvals_Rwd(g, lag+1) = min(baseRwdP, interRwdP);
                    pvals_Unrwd(g, lag+1) = min(baseUnrwdP, interUnrwdP);
                end
            end
        end

        colors = lines(nGeno);
        sigThresh = 0.05;
        CIfactor = 1.96; % 95% CI multiplier

        figure('Position', [100, 100, 1200, 400]);
        sgtitle(['mixed RCUC ', tlabel,num2str(nSession),' odor',odor_list{currOdor}]);
        % Rewarded history + bias plot
        h1=subplot(1,2,1);
        hold on;
        for g = 1:nGeno
            % Plot line
            plot(lags, betas_Rwd(g,:), '-o', 'Color', colors(g,:), 'MarkerFaceColor', colors(g,:), 'MarkerSize', 8);
            plot([-1 maxLag+1], [0,0], 'k', 'LineWidth', 0.5)
            % Plot error bars (95% CI)
            %errorbar(lags, betas_Rwd(g,:), SE_Rwd(g,:)*CIfactor,  ...
            %'Color', colors(g,:), 'CapSize', 0);    % Plot significance stars
            for i = 1:length(lags)
                x = lags(i);
                y = betas_Rwd(g,i);
                err = SE_Rwd(g,i)*CIfactor;
                plot([x x], [y-err, y+err], 'Color', colors(g,:));
            end

            sigIdx = pvals_Rwd(g,:) < sigThresh;
            yStar = betas_Rwd(g,sigIdx) + 0.1*range(betas_Rwd(g,:));
            text(lags(sigIdx)-0.05, yStar, '*', 'FontSize', 30, 'Color', colors(g,:), ...
                'HorizontalAlignment', 'center');
        end
        xlabel('Lag ');
        ylabel('Log-odds beta');
        title('Rewarded Choice');
        xlim([-0.1 maxLag+1])
        %ylim([-1 1])
        %legend(genoLevels, 'Location', 'best');
        %grid on;
        set(gca, 'XDir', 'reverse');
        hold off;

        % Unrewarded history + bias plot
        h2=subplot(1,2,2);
        hold on;
        for g = 1:nGeno
            line_handles(g)=plot(lags, betas_Unrwd(g,:), '-o', 'Color', colors(g,:), 'MarkerFaceColor', colors(g,:), 'MarkerSize', 8);
            plot([-1 maxLag+1], [0,0], 'k', 'LineWidth', 0.5)
            for i = 1:length(lags)
                x = lags(i);
                y = betas_Unrwd(g,i);
                err = SE_Unrwd(g,i)*CIfactor;
                plot([x x], [y-err, y+err], 'Color', colors(g,:));
            end
            %errorbar(lags, betas_Unrwd(g,:), SE_Unrwd(g,:)*CIfactor, 'LineStyle', 'none', 'Color', colors(g,:));
            sigIdx = pvals_Unrwd(g,:) < sigThresh;
            yStar = betas_Unrwd(g,sigIdx) + 0.1*range(betas_Unrwd(g,:));
            text(lags(sigIdx)-0.05, yStar, '*', 'FontSize', 30, 'Color', colors(g,:), ...
                'HorizontalAlignment', 'center');
        end
        xlabel('Lag ');
        %ylabel('Log-odds beta');
        xlim([-0.1 maxLag+1])
        %ylim([-1 1])
        title('Unrewarded Choice');
        legend(line_handles,genoLevels, 'Location', 'best','box','off');
        set(gca, 'XDir', 'reverse');
        hold off;
        linkaxes([h1, h2], 'y');

        print(gcf,'-dpng',fullfile(savefigpath, ['mixed_RCUC_', tlabel,num2str(nSession),'_odor',odor_list{currOdor}]));    %png format
        saveas(gcf, fullfile(savefigpath, ['mixed_RCUC_', tlabel,num2str(nSession),'_odor',odor_list{currOdor}]), 'fig');
        saveas(gcf, fullfile(savefigpath, ['mixed_RCUC_', tlabel,num2str(nSession),'_odor',odor_list{currOdor}]),'svg');

        % save the logstic regression result
        savelogname = fullfile(savedatapath, ['mixed_RCUC_', tlabel,num2str(nSession),' odor',odor_list{currOdor},'.mat']);
        save(savelogname, 'glme');
    end
end