%% Plot correlations
h_corr{m} = figure;
% Bandwidth(Hz) vs Best Frequency
subplot_tight(2,3,1,[.1 .06])
scatter([SigRISTRF(idx_core).BFHz],[SigRISTRF(idx_core).BWHz],'filled','r')
hold on
scatter([SigRISTRF(idx_belt).BFHz],[SigRISTRF(idx_belt).BWHz],'filled','g')
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('Best Frequency (Hz)')
ylabel('Bandwidth (Hz)')
[corr1_core, p1_core] = corr([SigRISTRF(idx_core).BFHz]',[SigRISTRF(idx_core).BWHz]','rows','complete');
[corr1_belt, p1_belt] = corr([SigRISTRF(idx_belt).BFHz]',[SigRISTRF(idx_belt).BWHz]','rows','complete');
if p1_core < .05 && p1_belt < .05
	p_corrcoeff1 = compare_correlation_coefficients(corr1_core,corr1_belt,length(idx_core),length(idx_belt));
	if p_corrcoeff1 < .05
		sig_r_label1 = '*';
	else
		sig_r_label1 = '';
	end
else
	sig_r_label1 = '';
end
legend({sprintf('Core (N=%i,r=%.2f)',numel([SigRISTRF(idx_core).BWHz]),corr1_core),sprintf('Belt (N=%i,r=%.2f)%s',numel([SigRISTRF(idx_belt).BWHz]),corr1_belt,sig_r_label1)})

% Bandwidth(octave) vs Best Frequency
subplot_tight(2,3,2,[.1 .06])
scatter([SigRISTRF(idx_core).BFHz],[SigRISTRF(idx_core).BW],'filled','r')
hold on
scatter([SigRISTRF(idx_belt).BFHz],[SigRISTRF(idx_belt).BW],'filled','g')
set(gca,'XScale','log')
%set(gca,'YScale','log')
xlabel('Best Frequency (Hz)')
ylabel('Bandwidth (octave)')
[corr2_core, p2_core] = corr([SigRISTRF(idx_core).BFHz]',[SigRISTRF(idx_core).BW]','rows','complete');
[corr2_belt, p2_belt] = corr([SigRISTRF(idx_belt).BFHz]',[SigRISTRF(idx_belt).BW]','rows','complete');
if p2_core < .05 && p2_belt < .05
	p_corrcoeff2 = compare_correlation_coefficients(corr2_core,corr2_belt,length(idx_core),length(idx_belt));
	if p_corrcoeff2 < .05
		sig_r_label2 = '*';
	else
		sig_r_label2 = '';
	end
else
	sig_r_label2 = '';
end
legend({sprintf('Core (N=%i,r=%.2f)',numel([SigRISTRF(idx_core).BW]),corr2_core),sprintf('Belt (N=%i,r=%.2f)%s',numel([SigRISTRF(idx_belt).BW]),corr2_belt,sig_r_label2)})

% Best Spectral Modulation Frequency vs Bandwidth
subplot_tight(2,3,3,[.1 .06])
BW_core = [];
for j = 1:numel(idx_core)
	if isempty(SigRISTRF(idx_core(j)).BW)
		BW_core(end+1,1) = nan;
	else
		BW_core(end+1,1) = SigRISTRF(idx_core(j)).BW;
	end
end
bSMF_core = [];
for j = 1:numel(idx_core)
	if isempty(SigRISTRF(idx_core(j)).bSMF)
		bSMF_core(end+1,1) = nan;
	else
		bSMF_core(end+1,1) = SigRISTRF(idx_core(j)).bSMF(1);
	end
end
BW_belt = [];
for j = 1:numel(idx_belt)
	if isempty(SigRISTRF(idx_belt(j)).BW)
		BW_belt(end+1,1) = nan;
	else
		BW_belt(end+1,1) = SigRISTRF(idx_belt(j)).BW;
	end
end
bSMF_belt = [];
for j = 1:numel(idx_belt)
	if isempty(SigRISTRF(idx_belt(j)).bSMF)
		bSMF_belt(end+1,1) = nan;
	else
		bSMF_belt(end+1,1) = SigRISTRF(idx_belt(j)).bSMF(1);
	end
end
%scatter([SigRISTRF(idx_core).BW],[SigRISTRF(idx_core).bSMF],'filled','r')
scatter(BW_core,bSMF_core,'filled','r')
hold
%scatter([SigRISTRF(idx_belt).BW],[SigRISTRF(idx_belt).bSMF],'filled','g')
scatter(BW_belt,bSMF_belt,'filled','g')
%set(gca,'XScale','log')
xlabel('Bandwidth (octave)')
ylabel('Best Spectral Mod. Freq.')
[corr3_core, p3_core] = corr(BW_core,bSMF_core,'rows','complete');
[corr3_belt, p3_belt] = corr(BW_belt,bSMF_belt,'rows','complete');
if p3_core < .05 && p3_belt < .05
	p_corrcoeff3 = compare_correlation_coefficients(corr3_core,corr3_belt,length(idx_core),length(idx_belt));
	if p_corrcoeff3 < .05
		sig_r_label3 = '*';
	else
		sig_r_label3 = '';
	end
else
	sig_r_label3 = '';
end
legend({sprintf('Core (N=%i,r=%.2f)',sum(~isnan(bSMF_core)),corr3_core),sprintf('Belt (N=%i,r=%.2f)%s',sum(~isnan(bSMF_belt)),corr3_belt,sig_r_label3)})

% Delay(t1_10 or simply Delay) vs Best Frequency
subplot_tight(2,3,4,[.1 .06])
scatter([SigRISTRF(idx_core).BFHz],[SigRISTRF(idx_core).t10],'filled','r')
hold on
scatter([SigRISTRF(idx_belt).BFHz],[SigRISTRF(idx_belt).t10],'filled','g')
set(gca,'XScale','log')
xlabel('Best Frequency (Hz)')
ylabel('Delay (ms)')
[corr4_core, p4_core] = corr([SigRISTRF(idx_core).BFHz]',[SigRISTRF(idx_core).Delay]','rows','complete');
[corr4_belt, p4_belt] = corr([SigRISTRF(idx_belt).BFHz]',[SigRISTRF(idx_belt).Delay]','rows','complete');
if p4_core < .05 && p4_belt < .05
	p_corrcoeff4 = compare_correlation_coefficients(corr4_core,corr4_belt,length(idx_core),length(idx_belt));
	if p_corrcoeff4 < .05
		sig_r_label4 = '*';
	else
		sig_r_label4 = '';
	end
else
	sig_r_label4 = '';
end
legend({sprintf('Core (N=%i,r=%.2f)',numel([SigRISTRF(idx_core).Delay]),corr4_core),sprintf('Belt (N=%i,r=%.2f)%s',numel([SigRISTRF(idx_belt).Delay]),corr4_belt,sig_r_label4)})

% Temporal Cutoff Frequency vs Integration Time
subplot_tight(2,3,5,[.1 .06])
scatter([SigRISTRF(idx_core).Dur],[SigRISTRF(idx_core).FmUpperCutoff],'filled','r')
hold on
scatter([SigRISTRF(idx_belt).Dur],[SigRISTRF(idx_belt).FmUpperCutoff],'filled','g')
xlabel('Integration Time (ms)')
ylabel('Temporal Cutoff Freq. (Hz)')
[corr5_core, p5_core] = corr([SigRISTRF(idx_core).Dur]',[SigRISTRF(idx_core).FmUpperCutoff]','rows','complete');
[corr5_belt, p5_belt] = corr([SigRISTRF(idx_belt).Dur]',[SigRISTRF(idx_belt).FmUpperCutoff]','rows','complete');
if p5_core < .05 && p5_belt < .05
	p_corrcoeff5 = compare_correlation_coefficients(corr5_core,corr5_belt,length(idx_core),length(idx_belt));
	if p_corrcoeff5 < .05
		sig_r_label5 = '*';
	else
		sig_r_label5 = '';
	end
else
	sig_r_label5 = '';
end
legend({sprintf('Core (N=%i,r=%.2f)',numel([SigRISTRF(idx_core).Dur]),corr5_core),sprintf('Belt (N=%i,r=%.2f)%s',numel([SigRISTRF(idx_belt).Dur]),corr5_belt,sig_r_label5)})

% Delay vs Integration Time
subplot_tight(2,3,6,[.1 .06])
scatter([SigRISTRF(idx_core).Dur],[SigRISTRF(idx_core).Delay],'filled','r')
hold on
scatter([SigRISTRF(idx_belt).Dur],[SigRISTRF(idx_belt).Delay],'filled','g')
xlabel('Integration Time (ms)')
ylabel('Delay (ms)')
[corr6_core, p6_core] = corr([SigRISTRF(idx_core).Dur]',[SigRISTRF(idx_core).Delay]','rows','complete');
[corr6_belt, p6_belt] = corr([SigRISTRF(idx_belt).Dur]',[SigRISTRF(idx_belt).Delay]','rows','complete');
if p6_core < .05 && p6_belt < .05
	p_corrcoeff6 = compare_correlation_coefficients(corr6_core,corr6_belt,length(idx_core),length(idx_belt));
	if p_corrcoeff6 < .05
		sig_r_label6 = '*';
	else
		sig_r_label6 = '';
	end
else
	sig_r_label6 = '';
end
legend({sprintf('Core (N=%i,r=%.2f)',numel([SigRISTRF(idx_core).Dur]),corr6_core),sprintf('Belt (N=%i,r=%.2f)%s',numel([SigRISTRF(idx_belt).Dur]),corr6_belt,sig_r_label6)})

%sgtitle(sprintf('N=%i, p_{thr}=%.0d',numel(SigRISTRF),p_crit),'FontSize',24)
currPos = get(h_corr{m},'Position'); set(gcf,'Position',[currPos(1),currPos(2),900,600]);	
saveas(h_corr{m},fullfile(dir_RISTRF_Fig,sprintf('RISTRF_Corr_Comb_%s.png',tag)))
%save(fullfile(dir_RISTRF,sprintf('RISTRF_Comb_%s.mat',tag)))

return


for m = 1:2

	%% Plot results
	
	coordLabels_x = num2str([-7:7]');
	coordLabels_y = flipud(coordLabels_x);
	
	titleStr{1} = sprintf('%% Proportion of Clusters with significant RI',monkeyName{m},numel(SigRISTRF),p_crit);
	titleStr{2} = 'Median p value of RI (log scale)';
	titleStr{3} = 'Best Frequency';
	titleStr{4} = 'Best Frequency (Hz)';
	titleStr{5} = 'Delay';
	titleStr{6} = 'Peak Delay';
	titleStr{7} = sprintf('T_{10}','Interpreter', 'none');
	titleStr{8} = sprintf('T_{50}','Interpreter', 'none');
	titleStr{9} = 'Duration';
	titleStr{10} = 'PLI';
	titleStr{11} = 'Bandwidth (octave)';
	titleStr{12} = 'Bandwidth (Hz)';
	titleStr{13} = 'Spectral Mod. Freq.';
	titleStr{14} = 'Temporal Cutoff Freq.';

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

	h_map2(m) = figure;

	for ldx_plt = 9:14
		subplot_tight(2,3,ldx_plt-8,[.1 .06])
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
	currPos = get(h_map2(m),'Position'); set(h_map2(m),'Position',[currPos(1),currPos(2),1200,800]);

	return

	% Distributions
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
	h11 = histogram(log10([SigRISTRF(idx_core).p]),Edges_p,'Normalization','count','FaceColor','r','FaceAlpha',.5);
	hold on
	h12 = histogram(log10([SigRISTRF(idx_belt).p]),Edges_p,'Normalization','count','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian(log([SigRISTRF(idx_core).p])),'b-.',sprintf('Median(Core) = %.1f',nanmedian(log([SigRISTRF(idx_core).p]))),'LineWidth',2,'FontSize',12)
	xline(nanmedian(log([SigRISTRF(idx_belt).p])),'b-.',sprintf('Median(Belt) = %.1f',nanmedian(log([SigRISTRF(idx_belt).p]))),'LineWidth',2,'FontSize',12)
	legend([h11 h12],legendLabel)
	ylabel('Number of Clusters')
	title([titleStr{2},sprintf('\np%s = %.2d (two-sample t-test)',sig_p,p_p)])
	
	subplot_tight(2,4,2,[.1 .06])
	h21 = histogram([SigRISTRF(idx_core).BF],Edges_BF,'Normalization','count','FaceColor','r','FaceAlpha',.5);
	hold on
	h22 = histogram([SigRISTRF(idx_belt).BF],Edges_BF,'Normalization','count','FaceColor','g','FaceAlpha',.5);
	legend([h21 h22],legendLabel)
	xlabel('[Hz]')
	ylabel('Number of Clusters')
	title(titleStr{3})
	
	subplot_tight(2,4,3,[.1 .06])
	h31 = histogram([SigRISTRF(idx_core).PD],Edges_PD,'Normalization','count','FaceColor','r','FaceAlpha',.5);
	hold on
	h32 = histogram([SigRISTRF(idx_belt).PD],Edges_PD,'Normalization','count','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian([SigRISTRF(idx_core).PD]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([SigRISTRF(idx_core).PD])),'LineWidth',2,'FontSize',12)
	xline(nanmedian([SigRISTRF(idx_belt).PD]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([SigRISTRF(idx_belt).PD])),'LineWidth',2,'FontSize',12)
	legend([h31 h32],legendLabel)
	xlabel('[ms]')
	ylabel('Number of Clusters')
	title([titleStr{4},sprintf('\np%s = %.2d (two-sample t-test)',sig_pd,p_pd)])
	
	subplot_tight(2,4,4,[.1 .06])
	h41 = histogram([SigRISTRF(idx_core).t10],Edges_t10,'Normalization','count','FaceColor','r','FaceAlpha',.5);
	hold on
	h42 = histogram([SigRISTRF(idx_belt).t10],Edges_t10,'Normalization','count','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian([SigRISTRF(idx_core).t10]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([SigRISTRF(idx_core).t10])),'LineWidth',2,'FontSize',12)
	xline(nanmedian([SigRISTRF(idx_belt).t10]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([SigRISTRF(idx_belt).t10])),'LineWidth',2,'FontSize',12)
	legend([h41 h42],legendLabel)
	xlabel('[ms]')
	ylabel('Number of Clusters')
	title([titleStr{5},sprintf('\np%s = %.2d (two-sample t-test)',sig_t10,p_t10)])
	
	subplot_tight(2,4,5,[.1 .06])
	h51 = histogram([SigRISTRF(idx_core).t50],Edges_t50,'Normalization','count','FaceColor','r','FaceAlpha',.5);
	hold on
	h52 = histogram([SigRISTRF(idx_belt).t50],Edges_t50,'Normalization','count','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian([SigRISTRF(idx_core).t50]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([SigRISTRF(idx_core).t50])),'LineWidth',2,'FontSize',12)
	xline(nanmedian([SigRISTRF(idx_belt).t50]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([SigRISTRF(idx_belt).t50])),'LineWidth',2,'FontSize',12)
	legend([h51 h52],legendLabel)
	xlabel('[ms]')
	ylabel('Number of Clusters')
	title([titleStr{6},sprintf('\np%s = %.2d (two-sample t-test)',sig_t50,p_t50)])
	
	subplot_tight(2,4,6,[.1 .06])
	h61 = histogram([SigRISTRF(idx_core).Dur],Edges_Dur,'Normalization','count','FaceColor','r','FaceAlpha',.5);
	hold on
	h62 = histogram([SigRISTRF(idx_belt).Dur],Edges_Dur,'Normalization','count','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian([SigRISTRF(idx_core).Dur]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([SigRISTRF(idx_core).Dur])),'LineWidth',2,'FontSize',12)
	xline(nanmedian([SigRISTRF(idx_belt).Dur]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([SigRISTRF(idx_belt).Dur])),'LineWidth',2,'FontSize',12)
	legend([h61 h62],legendLabel)
	xlabel('[ms]')
	ylabel('Number of Clusters')
	title([titleStr{7},sprintf('\np%s = %.2d (two-sample t-test)',sig_dur,p_dur)])
	
	subplot_tight(2,4,7,[.1 .06])
	h71 = histogram([SigRISTRF(idx_core).PLI],Edges_PLI,'Normalization','count','FaceColor','r','FaceAlpha',.5);
	hold on
	h72 = histogram([SigRISTRF(idx_belt).PLI],Edges_PLI,'Normalization','count','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian([SigRISTRF(idx_core).PLI]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([SigRISTRF(idx_core).PLI])),'LineWidth',2,'FontSize',12)
	xline(nanmedian([SigRISTRF(idx_belt).PLI]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([SigRISTRF(idx_belt).PLI])),'LineWidth',2,'FontSize',12)
	legend([h71 h72],legendLabel)
	ylabel('Number of Clusters')
	title([titleStr{8},sprintf('\np%s = %.2d (two-sample t-test)',sig_pli,p_pli)])
	
	subplot_tight(2,4,8,[.1 .06])
	h81 = histogram([SigRISTRF(idx_core).DSI],Edges_DSI,'Normalization','count','FaceColor','r','FaceAlpha',.5);
	hold on
	h82 = histogram([SigRISTRF(idx_belt).DSI],Edges_DSI,'Normalization','count','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian([SigRISTRF(idx_core).DSI]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([SigRISTRF(idx_core).DSI])),'LineWidth',2,'FontSize',12)
	xline(nanmedian([SigRISTRF(idx_belt).DSI]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([SigRISTRF(idx_belt).DSI])),'LineWidth',2,'FontSize',12)
	legend([h81 h82],legendLabel)
	ylabel('Number of Clusters')
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
	%sgtitle(sprintf('Area Map for %s',monkeyName{m}),'FontSize',24)
	%currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),412,372]);
	save(fullfile(dir_RISTRF,sprintf('RISTRF_%s_%s.mat',tag,monkeyName{m})))
end

