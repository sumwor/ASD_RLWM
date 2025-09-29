function plot_psy_fit(psymodelName, dataIndex, savefigpath, label)

psy_fit = readtable(psymodelName);

if strcmp(label,'AB-AB1')
    protocol = 'AB';
    protocolDay = 1;
end

subjects = psy_fit.Animal;
nSubject = length(subjects);
genotypes = cell(nSubject,1);

%% load the data
for k = 1:length(subjects)
    tempgeno = dataIndex.Genotype(strcmp(dataIndex.Animal,num2str(subjects(k))));
    genotypes{k} = tempgeno{1};
    
    analysis = dataIndex.BehPath(strcmp(dataIndex.Animal,num2str(subjects(k))));
    modelpath = fullfile(analysis{1},'latent',['psy_fit_',label,'.json']);

    txt = fileread(modelpath);
    psy_latent = jsondecode(txt);
    
    if k==1
        nParams = length(psy_latent.opt_hyper.alpha) + length(psy_latent.opt_hyper.alpha);
        paramNames = {'bias', 'Stimulus', 'stick'};
        params = zeros(nSubject,nParams);
    end
    params(k,1:length(paramNames)) = log10(psy_latent.opt_hyper.alpha);
    params(k,length(paramNames)+1:end) = log10(psy_latent.opt_hyper.sigma);
end


% plot beta, alpha+, s1, lapse, rec, bias for different genotype
%params = fit_result.All_params{modelNum};         

genos  = genotypes;        


% Unique genotypes
genoList = unique(genos);

figure('Position',[100 100 1600 1200]);
for p = 1:6
    % plot learning rate
    subplot(2,3,p); hold on;
    
    % Plot bars first (group mean)
    if p <=3
        title(paramNames{p});
    end
    if p == 1
        ylabel('Log_{10}\alpha')
    end
    if p == 4
        ylabel('Log_{10}\sigma')
    end
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
    %ylabel(paramNames{p});
    %title(paramNames{p});
    
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
sgtitle(['Fitted psy parameters by genotype ', label]);

print(gcf,'-dpng',fullfile(savefigpath, ['Fitted psy parameters by genotype ', label]));    %png format
saveas(gcf, fullfile(savefigpath, ['Fitted psy parameters by genotype ', label]), 'fig');
saveas(gcf, fullfile(savefigpath, ['Fitted psy parameters by genotype ', label]),'svg');