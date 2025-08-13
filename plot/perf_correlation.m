function perf_correlation(perf_AB, perf_CD, perf_DC, genotype, savefigpath)

%% plot performacne correlation

figure;

% 
DC_ave_perform = nanmean(perf_DC(:,5:6),2);
notnanMask_DC = ~isnan(DC_ave_perform);
notnanMask_AB = ~isnan(perf_AB);
notnanMask_CD = ~isnan(perf_CD);

coeff_WT = nan(2,3);
coeff_HET = nan(2,3);
coeff_KO = nan(2,3);
p_WT = nan(2,3);
p_HET = nan(2,3);
p_KO = nan(2,3);

WT_mask = strcmp(genotype, 'WT');
if ismember('HEM',genotype)
    HET_mask = strcmp(genotype,'HEM');
elseif ismember('KO', genotype) & length(unique(genotype))==3
    KO_mask = strcmp(genotype,'KO');
    HET_mask = strcmp(genotype,'HET');
elseif ismember('KO', genotype) & length(unique(genotype))==2
    HET_mask = strcmp(genotype, 'KO');
else
    HET_mask = strcmp(genotype,'HET');
    KO_mask = nan;
end

for ii = 1:3  % AB sessions
    
    subplot(2,4,ii)
    title(['AB session', num2str(ii)])
    hold on;
    scatter(perf_AB(WT_mask,ii), DC_ave_perform(WT_mask), 100, 'filled', 'blue');
    scatter(perf_AB(HET_mask,ii), DC_ave_perform(HET_mask), 100, 'filled','red');
    xlim([-0.5,0.5])
    ylim([-0.5, 0.5])
    % calculate correlation

    [tempcorr,tempP] = corrcoef(perf_AB(WT_mask&notnanMask_DC&notnanMask_AB(:,ii), ii), ...
    DC_ave_perform(WT_mask&notnanMask_DC&notnanMask_AB(:,ii)));
    coeff_WT(1,ii) = tempcorr(1,2);
    p_WT(1,ii) = tempP(1,2);

    [tempcorr, tempP] = corrcoef(perf_AB(HET_mask&notnanMask_DC&notnanMask_AB(:,ii), ii), ...
    DC_ave_perform(HET_mask&notnanMask_DC&notnanMask_AB(:,ii)));
    coeff_HET(1,ii) = tempcorr(1,2);
    p_HET(1,ii) = tempP(1,2);

    if ii == 1
        ylabel('DC performance (AUC)')
    end
end

% plot correlation coefficient and p-value
subplot(2,4,4)
plot(1:3, coeff_WT(1,:), 'blue')
hold on; plot(1:3, coeff_HET(1,:), 'red')
xlim([0.8, 3.2])
ylim([-1 1])
set(gca, 'box','off')
title('Pearson correlation')
for pp = 1:3
    if p_WT(1,pp) < 0.05
        text(pp, coeff_WT(1,pp),num2str(p_WT(1,pp)), 'FontSize',15);
    end
    if p_HET(1,pp) < 0.05
        text(pp, coeff_HET(1,pp),num2str(p_HET(1,pp)), 'FontSize',15);
    end
end

for jj = 1:3  % AB sessions
    
    subplot(2,4,jj+4)
    title(['CD session', num2str(jj)])
    hold on;
    scatter(perf_CD(WT_mask,jj), DC_ave_perform(WT_mask), 100, 'filled', 'blue');
    scatter(perf_CD(HET_mask,jj), DC_ave_perform(HET_mask), 100, 'filled','red');
    xlim([-0.5,0.5])
    ylim([-0.5, 0.5])
    if ii == 1
        ylabel('DC performance (AUC)')
    end

        [tempcorr,tempP] = corrcoef(perf_CD(WT_mask&notnanMask_DC&notnanMask_CD(:,jj), jj), ...
        DC_ave_perform(WT_mask&notnanMask_DC&notnanMask_CD(:,jj)));
    coeff_WT(2,jj) = tempcorr(1,2);
    p_WT(2,jj) = tempP(1,2);
         [tempcorr,tempP] = corrcoef(perf_CD(HET_mask&notnanMask_DC&notnanMask_CD(:,jj), jj), ...
        DC_ave_perform(HET_mask&notnanMask_DC&notnanMask_CD(:,jj)));
    coeff_HET(2,jj) = tempcorr(1,2);
    p_HET(2,jj) = tempP(1,2);

end

% plot correlation coefficient and p-value
subplot(2,4,8)
plot(1:3, coeff_WT(2,:), 'blue')
hold on; plot(1:3, coeff_HET(2,:), 'red')
xlim([0.8, 3.2])
ylim([-1 1])
set(gca, 'box','off')
title('Pearson correlation')
for pp = 1:3
    if p_WT(2,pp) < 0.05
        text(pp, coeff_WT(2,pp),num2str(p_WT(2,pp)), 'FontSize',15);
    end
    if p_HET(2,pp) < 0.05
        text(pp, coeff_HET(2,pp),num2str(p_HET(2,pp)), 'FontSize',15);
    end
end

print(gcf,'-dpng',fullfile(savefigpath,['Performance correlation (AUC)']));    %png format
saveas(gcf, fullfile(savefigpath, ['Performance correlation (AUC)']), 'fig');
saveas(gcf, fullfile(savefigpath, ['Performance correlation (AUC)']),'svg');


