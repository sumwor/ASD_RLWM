function plot_hybrid_fit(fit_result, savefigpath, label, modelNum)

% plot beta, alpha+, s1, lapse, rec, bias for different genotype
%params = fit_result.All_params{modelNum};         
params = fit_result.All_params{1};
genos  = fit_result.genotypes;        
paramNames = {'\beta','log(\alpha^+)', 's1','lapse','recover','normalized bias'};

% Apply log transform to alpha+ (col 2)
params(:,2) = log10(params(:,2));

% normalize bias to absolute bias from 50%
params(:,6) = abs(params(:,6)-0.5);
% Unique genotypes
genoList = unique(genos);

figure('Position',[100 100 1600 1200]);
for p = 1:6
    subplot(2,3,p); hold on;
    
    % Plot bars first (group mean)
    for g = 1:numel(genoList)
        this_idx = strcmp(genos, genoList{g});
        y = params(this_idx,p);
        
        % bar
        bar(g, mean(y), 'FaceAlpha',0.5); 
        
        % errorbar (SEM)
        errorbar(g, mean(y), std(y)/sqrt(numel(y)), ...
                 'k','LineStyle','none','LineWidth',1.2);
        
        % scatter individual dots (with jitter)
        xJitter = (rand(size(y))-0.5)*0.2;  % small jitter
        scatter(g + xJitter, y, 50, 'filled', ...
                'MarkerFaceAlpha',0.7, 'MarkerEdgeColor','k');
    end
    
    % Formatting
    set(gca,'XTick',1:numel(genoList),'XTickLabel',genoList);
    ylabel(paramNames{p});
    title(paramNames{p});
    
    % Mann–Whitney U test if exactly 2 genotypes
    if numel(genoList)==2
        y1 = params(strcmp(genos,genoList{1}),p);
        y2 = params(strcmp(genos,genoList{2}),p);
        pval = ranksum(y1,y2);

        % Find y-axis limits and adjust
        yl = ylim;
        yText = yl(2) - 0.1*(yl(2)-yl(1));  % place 10% below top
        text(1.5, yText, ...
             ['p = ' num2str(pval,'%.3g')], ...
             'HorizontalAlignment','center', 'FontSize',15);

        % Optionally expand ylim so bars+dots+p-value fit nicely
        ylim([yl(1), yl(2)+(yl(2)-yl(1))*0.2]);
    end
end
sgtitle(['Fitted parameters by genotype ', label]);

print(gcf,'-dpng',fullfile(savefigpath, ['Fitted parameters by genotype ', label]));    %png format
saveas(gcf, fullfile(savefigpath, ['Fitted parameters by genotype ', label]), 'fig');
saveas(gcf, fullfile(savefigpath, ['Fitted parameters by genotype ', label]),'svg');