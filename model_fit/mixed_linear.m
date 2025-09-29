function mixed_linear_psy_weight(tbl, label, savefigpath)

lme_inter = fitlme(tbl, ...
    'ResponseTime ~ Bias*Genotype + Stimulus*Genotype + Stick*Genotype + (1|Animal)');
disp(lme_inter)

predictors = {'Bias','Stimulus','Stick'};
catPred = categorical(predictors, predictors, 'Ordinal',true);  % force order
slopes = [WT_est(:), HET_est(:)];  % rows = predictors, cols = WT/HET
figure; hold on;
b = bar(catPred, slopes, 'grouped');
b(1).FaceColor = [0.3 0.6 0.9];  % WT
b(2).FaceColor = [0.9 0.4 0.4];  % HET

% --- Error bars ---
for i = 1:length(predictors)
    % WT
    xWT = b(1).XEndPoints(i);
    yWT = WT_est(i);
    errorbar(xWT, yWT, yWT-WT_CI(1,i), WT_CI(2,i)-yWT, 'k','LineStyle','none','LineWidth',1.2,'CapSize',8);
    
    % HET
    xHET = b(2).XEndPoints(i);
    yHET = HET_est(i);
    errorbar(xHET, yHET, yHET-HET_CI(1,i), HET_CI(2,i)-yHET, 'k','LineStyle','none','LineWidth',1.2,'CapSize',8);
end

ylabel('Slope (Effect on RT)');
legendHandle = legend({'WT','HET'}, 'Location','best'); 
set(legendHandle, 'Box', 'off', 'Color', 'none');  % make legend box off and transparent
title('Predictor Effects on Response Time by Genotype');

% --- Add genotype × weight p-values between bars ---
for i = 1:length(predictors)
    % Find the interaction term
    idx_int = find(strcmp(lme_inter.CoefficientNames, [predictors{i} ':Genotype_HET']));
    if ~isempty(idx_int)
        pVal = lme_inter.Coefficients.pValue(idx_int);
        
        % Midpoint between the two bars
        xMid = mean([b(1).XEndPoints(i), b(2).XEndPoints(i)]);
        yTop = max(WT_est(i), HET_est(i)) + 0.001;  % position above the taller bar
        
        % Draw a line connecting the two bars
        %plot([b(1).XEndPoints(i), b(2).XEndPoints(i)], [yTop yTop], 'k-', 'LineWidth',1.2);
        
        % Add p-value text
        text(xMid, yTop + 0.02, sprintf('p=%.3g', pVal), ...
            'HorizontalAlignment','center','FontSize',20,'Color','k');
    end
end