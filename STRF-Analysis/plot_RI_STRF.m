
%% Plot results

coordLabels_x = num2str([-7:7]');
coordLabels_y = flipud(coordLabels_x);

titleStr{1} = sprintf('%% Proportion of Clusters with significant RI',monkeyName{m},numel(SigRISTRF),p_crit);
titleStr{2} = 'Median p value of RI (log scale)';
titleStr{3} = 'Best Frequency';
titleStr{4} = 'Peak Delay';
titleStr{5} = sprintf('T_{10}','Interpreter', 'none');
titleStr{6} = sprintf('T_{50}','Interpreter', 'none');
titleStr{7} = 'Duration';
titleStr{8} = 'PLI';
titleStr{9} = 'DSI';

% Map
h_map(m) = figure;

for ldx_plt = 1:8
	subplot_tight(2,4,ldx_plt,[.1 .06])
	h_n = heatmap(flipud(htbl(:,:,ldx_plt,m)));
	h_n.CellLabelFormat = '%.1f';
	set(gca,'XDisplayLabels',coordLabels_x)
	set(gca,'YDisplayLabels',coordLabels_y)
	colormap jet;
	if strcmp(monkeyName{m},'Cassius')
		xlabel('Medial <----> Lateral')
	else
		xlabel('Lateral <----> Medial')
	end
	ylabel('Posterior <----> Anterior')
	title(titleStr{ldx_plt})
end

sgtitle(sprintf('%s (N=%i, p_{thr}=%.0d)',monkeyName{m},numel(SigRISTRF),p_crit),'FontSize',24)
currPos = get(h_map(m),'Position'); set(h_map(m),'Position',[currPos(1),currPos(2),1600,800]);


% Distributions
idx_a1 = find(cellfun(@(x) strcmp(x,'A1'), {SigRISTRF.area}));
idx_r = find(cellfun(@(x) strcmp(x,'R'), {SigRISTRF.area}));
idx_bm = find(cellfun(@(x) strcmp(x,'BeltM'), {SigRISTRF.area}));
idx_bl = find(cellfun(@(x) strcmp(x,'BeltL'), {SigRISTRF.area}));
idx_core = sort([idx_a1 idx_r]);
idx_belt = sort([idx_bm idx_bl]);

[h_p, p_p, ci_p, stat_p] = ttest2([SigRISTRF(idx_core).p],[SigRISTRF(idx_belt).p]);
if p_p <.05
    sig_p = '*';
else
    sig_p = '';
end
[h_pd, p_pd, ci_pd, stat_pd] = ttest2([SigRISTRF(idx_core).PD],[SigRISTRF(idx_belt).PD]);
if p_pd <.05
    sig_pd = '*';
else
    sig_pd = '';
end
[h_t10, p_t10, ci_t10, stat_t10] = ttest2([SigRISTRF(idx_core).t10],[SigRISTRF(idx_belt).t10]);
if p_t10 <.05
    sig_t10 = '*';
else
    sig_t10 = '';
end
[h_t50, p_t50, ci_t50, stat_t50] = ttest2([SigRISTRF(idx_core).PD],[SigRISTRF(idx_belt).PD]);
if p_t50 <.05
    sig_t50 = '*';
else
    sig_t50 = '';
end
[h_dur, p_dur, ci_dur, stat_dur] = ttest2([SigRISTRF(idx_core).Dur],[SigRISTRF(idx_belt).Dur]);
if p_dur <.05
    sig_dur = '*';
else
    sig_dur = '';
end
[h_pli, p_pli, ci_pli, stat_pli] = ttest2([SigRISTRF(idx_core).PLI],[SigRISTRF(idx_belt).PLI]);
if p_pli <.05
    sig_pli = '*';
else
    sig_pli = '';
end
[h_dsi, p_dsi, ci_dsi, stat_dsi] = ttest2([SigRISTRF(idx_core).DSI],[SigRISTRF(idx_belt).DSI]);
if p_dsi <.05
    sig_dsi = '*';
else
    sig_dsi = '';
end

Edges_p = [-35:1:max(max(log([SigRISTRF(idx_core).p])),max(log([SigRISTRF(idx_belt).p])))];
Edges_BF = [0:1000:max(max([SigRISTRF(idx_core).BF]),max([SigRISTRF(idx_belt).BF]))];
Edges_PD = [0:5:max(max([SigRISTRF(idx_core).PD]),max([SigRISTRF(idx_belt).PD]))];
Edges_t10 = [0:5:max(max([SigRISTRF(idx_core).t10]),max([SigRISTRF(idx_belt).t10]))];
Edges_t50 = [0:5:max(max([SigRISTRF(idx_core).t50]),max([SigRISTRF(idx_belt).t50]))];
Edges_Dur = [0:5:max(max([SigRISTRF(idx_core).Dur]),max([SigRISTRF(idx_belt).Dur]))];
Edges_PLI = [0:.02:max(max([SigRISTRF(idx_core).PLI]),max([SigRISTRF(idx_belt).PLI]))];
Edges_DSI = [min(min([SigRISTRF(idx_core).DSI]),min([SigRISTRF(idx_belt).DSI])):.05:max(max([SigRISTRF(idx_core).DSI]),max([SigRISTRF(idx_belt).DSI]))];

legendLabel = {'Core','Belt'};

h_dist(m) = figure;

subplot_tight(2,4,1,[.1 .06])
h1 = histogram(log10([SigRISTRF(idx_core).p]),Edges_p,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
hold on
h2 = histogram(log10([SigRISTRF(idx_belt).p]),Edges_p,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
xline(nanmedian(log([SigRISTRF(idx_core).p])),'b-.',sprintf('Median(Core) = %.1f',nanmedian(log([SigRISTRF(idx_core).p]))),'LineWidth',2,'FontSize',14)
xline(nanmedian(log([SigRISTRF(idx_belt).p])),'b-.',sprintf('Median(Belt) = %.1f',nanmedian(log([SigRISTRF(idx_belt).p]))),'LineWidth',2,'FontSize',14)
legend([h1 h2],legendLabel)
ylabel('Proportion of Clusters')
title([titleStr{2},sprintf('\np%s = %.2d (two-sample t-test)',sig_p,p_p)])

subplot_tight(2,4,2,[.1 .06])
h1 = histogram([SigRISTRF(idx_core).BF],Edges_BF,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
hold on
h2 = histogram([SigRISTRF(idx_belt).BF],Edges_BF,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
legend([h1 h2],legendLabel)
xlabel('[Hz]')
ylabel('Proportion of Clusters')
title(titleStr{3})

subplot_tight(2,4,3,[.1 .06])
h1 = histogram([SigRISTRF(idx_core).PD],Edges_PD,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
hold on
h2 = histogram([SigRISTRF(idx_belt).PD],Edges_PD,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
xline(nanmedian([SigRISTRF(idx_core).PD]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([SigRISTRF(idx_core).PD])),'LineWidth',2,'FontSize',14)
xline(nanmedian([SigRISTRF(idx_belt).PD]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([SigRISTRF(idx_belt).PD])),'LineWidth',2,'FontSize',14)
legend([h1 h2],legendLabel)
xlabel('[ms]')
ylabel('Proportion of Clusters')
title([titleStr{4},sprintf('\np%s = %.2d (two-sample t-test)',sig_pd,p_pd)])

subplot_tight(2,4,4,[.1 .06])
h1 = histogram([SigRISTRF(idx_core).t10],Edges_t10,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
hold on
h2 = histogram([SigRISTRF(idx_belt).t10],Edges_t10,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
xline(nanmedian([SigRISTRF(idx_core).t10]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([SigRISTRF(idx_core).t10])),'LineWidth',2,'FontSize',14)
xline(nanmedian([SigRISTRF(idx_belt).t10]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([SigRISTRF(idx_belt).t10])),'LineWidth',2,'FontSize',14)
legend([h1 h2],legendLabel)
xlabel('[ms]')
ylabel('Proportion of Clusters')
title([titleStr{5},sprintf('\np%s = %.2d (two-sample t-test)',sig_t10,p_t10)])

subplot_tight(2,4,5,[.1 .06])
h1 = histogram([SigRISTRF(idx_core).t50],Edges_t50,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
hold on
h2 = histogram([SigRISTRF(idx_belt).t50],Edges_t50,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
xline(nanmedian([SigRISTRF(idx_core).t50]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([SigRISTRF(idx_core).t50])),'LineWidth',2,'FontSize',14)
xline(nanmedian([SigRISTRF(idx_belt).t50]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([SigRISTRF(idx_belt).t50])),'LineWidth',2,'FontSize',14)
legend([h1 h2],legendLabel)
xlabel('[ms]')
ylabel('Proportion of Clusters')
title([titleStr{6},sprintf('\np%s = %.2d (two-sample t-test)',sig_t50,p_t50)])

subplot_tight(2,4,6,[.1 .06])
h1 = histogram([SigRISTRF(idx_core).Dur],Edges_Dur,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
hold on
h2 = histogram([SigRISTRF(idx_belt).Dur],Edges_Dur,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
xline(nanmedian([SigRISTRF(idx_core).Dur]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([SigRISTRF(idx_core).Dur])),'LineWidth',2,'FontSize',14)
xline(nanmedian([SigRISTRF(idx_belt).Dur]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([SigRISTRF(idx_belt).Dur])),'LineWidth',2,'FontSize',14)
legend([h1 h2],legendLabel)
xlabel('[ms]')
ylabel('Proportion of Clusters')
title([titleStr{7},sprintf('\np%s = %.2d (two-sample t-test)',sig_dur,p_dur)])

subplot_tight(2,4,7,[.1 .06])
h1 = histogram([SigRISTRF(idx_core).PLI],Edges_PLI,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
hold on
h2 = histogram([SigRISTRF(idx_belt).PLI],Edges_PLI,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
xline(nanmedian([SigRISTRF(idx_core).PLI]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([SigRISTRF(idx_core).PLI])),'LineWidth',2,'FontSize',14)
xline(nanmedian([SigRISTRF(idx_belt).PLI]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([SigRISTRF(idx_belt).PLI])),'LineWidth',2,'FontSize',14)
legend([h1 h2],legendLabel)
ylabel('Proportion of Clusters')
title([titleStr{8},sprintf('\np%s = %.2d (two-sample t-test)',sig_pli,p_pli)])

subplot_tight(2,4,8,[.1 .06])
h1 = histogram([SigRISTRF(idx_core).DSI],Edges_DSI,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
hold on
h2 = histogram([SigRISTRF(idx_belt).DSI],Edges_DSI,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
xline(nanmedian([SigRISTRF(idx_core).DSI]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([SigRISTRF(idx_core).DSI])),'LineWidth',2,'FontSize',14)
xline(nanmedian([SigRISTRF(idx_belt).DSI]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([SigRISTRF(idx_belt).DSI])),'LineWidth',2,'FontSize',14)
legend([h1 h2],legendLabel)
ylabel('Proportion of Clusters')
title([titleStr{9},sprintf('\np%s = %.2d (two-sample t-test)',sig_dsi,p_dsi)])

sgtitle(sprintf('%s (N=%i, p_{thr}=%.0d)',monkeyName{m},numel(SigRISTRF),p_crit),'FontSize',24)
currPos = get(h_dist(m),'Position'); set(h_dist(m),'Position',[currPos(1),currPos(2),1600,800]);

%% Create reference area map
%htbl_area = ones(15,15)*-1;
%if strcmp(monkeyName{m},'Miyagi')
%	coord_core = [coord_ac(2:20,1); coord_ac(2:25,2)];	
%	coord_belt = [coord_ac(2:31,3); coord_ac(2:16,4)];
%else
%	coord_core = [coord_ac(2:23,1); coord_ac(2:19,2)];	
%	coord_belt = [coord_ac(2:36,3); coord_ac(2:24,4)];
%end

%for ldx_core = 1:numel(coord_core)
%	htbl_area(coord_core{ldx_core}(2)+8,coord_core{ldx_core}(1)+8) = 1;
%end
%for ldx_belt = 1:numel(coord_belt)
%	htbl_area(coord_belt{ldx_belt}(2)+8,coord_belt{ldx_belt}(1)+8) = 0;
%end

%cmap = [0 0 0; 0 1 0; 1 0 0];
%figure;
%heatmap(flipud(htbl_area),'CellLabelColor','none')
%grid off
%colormap(cmap);
%colorbar('off')
%set(gca,'XDisplayLabels',coordLabels_x)
%set(gca,'YDisplayLabels',coordLabels_y)
%if strcmp(monkeyName{m},'Cassius')
%	xlabel('Medial <----> Lateral')
%else
%	xlabel('Lateral <----> Medial')
%end
%ylabel('Posterior <----> Anterior')
%sgtitle(sprintf('Area Map for %s',monkeyName{m}),'FontSize',14)
%currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),412,372]);
