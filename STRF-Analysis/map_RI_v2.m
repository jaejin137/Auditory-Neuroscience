monkeyName = {'Miyagi','Cassius'};
%monkeyName = {'Cassius','Miyagi'};
p_crit = .01;


htbl = nan(15,15,9,2);

for m = 1:2	
	dir_Casmiya = fullfile('~','Casmiya');
	dir_STRF = fullfile('~','STRF');
	dir_RI = fullfile(dir_STRF,'RI');
	dir_Param = fullfile(dir_STRF,'STRFParams');
	
	
	% Load coordinates data
	load(fullfile(dir_STRF,sprintf('coord_ac_%s.mat',monkeyName{m})))
	% Load laminar profile data
	load(fullfile(dir_Casmiya,sprintf('lamProf_%s.mat',monkeyName{m})))
	
	list_RI = dir(fullfile(dir_RI,sprintf('RI_%s*_8-100.mat',monkeyName{m})));
	%list_Param = dir(fullfile(dir_Param,sprintf('STRFParams_%s_*.mat',monkeyName{m})));
	
	
	UberRIStat = struct([]);
	j = 1;
	% Construct RI stats
	for idx_RI = 1:length(list_RI)
		% session date
		nameTok = split(list_RI(idx_RI).name,'_');
		dateTok = split(nameTok{2},'-');
		RIStat(idx_RI).date = dateTok{2};
		% coordinate
		RIStat(idx_RI).coord = lamProf(find([lamProf.date]==str2num(RIStat(idx_RI).date))).coord;
		% area name
		try
			[row, col] = ind2sub(size(coord_ac),find(cellfun(@(x) isequal(x,RIStat(idx_RI).coord), coord_ac)));
			RIStat(idx_RI).area = coord_ac{1,col};
		end

		% Load RI
		% number of clusters and median p value thresholded by p_crit
		load(fullfile(dir_RI,list_RI(idx_RI).name),'STRFSig');
		idx_valid = find(cellfun(@(x) ~isempty(x), {STRFSig.clustNum}));
		idx_clust = find([STRFSig(idx_valid).p] <= p_crit);
		RIStat(idx_RI).numClust = numel(idx_clust);
		RIStat(idx_RI).p_median = nanmedian([STRFSig(find([STRFSig.p] <= p_crit)).p]');

		% Load STRF parameters
		try
			load(fullfile(dir_Param,sprintf('STRFParams_%s_%s_%s_%s_%s.mat',monkeyName{m},dateTok{2},nameTok{3:5})),'UberSTRF')
		catch
			fprintf('\n!!!Warning: UberSTRF not available for %s_%s_%s-%s-%s!!!\n',monkeyName{m},dateTok{2},nameTok{3:5})
		end
			

		for i = 1:numel(idx_valid)
			if STRFSig(idx_valid(i)).p <= p_crit
				UberRIStat(j).date = RIStat(idx_RI).date;
				UberRIStat(j).coord = RIStat(idx_RI).coord;
				UberRIStat(j).area = RIStat(idx_RI).area;
				UberRIStat(j).clust = STRFSig(idx_valid(i)).clustNum;
				UberRIStat(j).p = STRFSig(idx_valid(i)).p;
				try
					idx_param = find(cellfun(@(x) isequal(x,STRFSig(idx_valid(i)).clustNum),{UberSTRF.ClustNum}));
					if ~isempty(idx_param)
						UberRIStat(j).BF = UberSTRF(idx_param).RFParam.BFHz;
						UberRIStat(j).PD = UberSTRF(idx_param).RFParam.PeakDelay;
						UberRIStat(j).t10 = UberSTRF(idx_param).RFParam.t1_10;
						UberRIStat(j).t50 = UberSTRF(idx_param).RFParam.t1_50;
						UberRIStat(j).Dur = UberSTRF(idx_param).RFParam.Duration;
						UberRIStat(j).PLI = UberSTRF(idx_param).RFParam.PLI;
						UberRIStat(j).DSI = UberSTRF(idx_param).RFParam.DSI;
					end
				end
				j = j+1;
			end
		end
	end
	
	
	% Unique coordinates
	list_coord = {};
	for i = 1:numel(UberRIStat)
		if i == 1
			list_coord{1,1} = UberRIStat(i).coord;
		else
			if ~any(cellfun(@(x) isequal(x,UberRIStat(i).coord), list_coord))
				list_coord{end+1,1} = UberRIStat(i).coord;
			end
		end
	end
	
	
	for k = 1:numel(list_coord)
		idx_sameCoord{k} = find(cellfun(@(x) isequal(x,list_coord{k}), {UberRIStat.coord}));
		RIMap(k).coord = list_coord{k};
		RIMap(k).numClust = numel(idx_sameCoord{k});
		RIMap(k).p_median = nanmedian([UberRIStat(idx_sameCoord{k}).p]);
		RIMap(k).BF = nanmean([UberRIStat(idx_sameCoord{k}).BF]);
		RIMap(k).PD = nanmean([UberRIStat(idx_sameCoord{k}).PD]);
		RIMap(k).t10 = nanmean([UberRIStat(idx_sameCoord{k}).t10]);
		RIMap(k).t50 = nanmean([UberRIStat(idx_sameCoord{k}).t50]);
		RIMap(k).Dur = nanmean([UberRIStat(idx_sameCoord{k}).Dur]);
		RIMap(k).PLI = nanmean([UberRIStat(idx_sameCoord{k}).PLI]);
		RIMap(k).DSI = nanmean([UberRIStat(idx_sameCoord{k}).DSI]);

		% Construct tables for heatmap
		htbl(list_coord{k}(2)+8,list_coord{k}(1)+8,1,m) = RIMap(k).numClust/numel(UberRIStat)*100;
		%htbl(list_coord{k}(2)+8,list_coord{k}(1)+8,1) = RIMap(k).numClust;
		htbl(list_coord{k}(2)+8,list_coord{k}(1)+8,2,m) = log10(RIMap(k).p_median);
		%htbl(list_coord{k}(2)+8,list_coord{k}(1)+8,2) = RIMap(k).p_median;

		try
			htbl(list_coord{k}(2)+8,list_coord{k}(1)+8,3,m) = RIMap(k).BF;
			htbl(list_coord{k}(2)+8,list_coord{k}(1)+8,4,m) = RIMap(k).PD;
			htbl(list_coord{k}(2)+8,list_coord{k}(1)+8,5,m) = RIMap(k).t10;
			htbl(list_coord{k}(2)+8,list_coord{k}(1)+8,6,m) = RIMap(k).t50;
			htbl(list_coord{k}(2)+8,list_coord{k}(1)+8,7,m) = RIMap(k).Dur;
			htbl(list_coord{k}(2)+8,list_coord{k}(1)+8,8,m) = RIMap(k).PLI;
			htbl(list_coord{k}(2)+8,list_coord{k}(1)+8,9,m) = RIMap(k).DSI;
		end
	end
	
	coordLabels_x = num2str([-7:7]');
	coordLabels_y = flipud(coordLabels_x);
	
	titleStr{1} = sprintf('%% Proportion of Clusters with significant RI',monkeyName{m},numel(UberRIStat),p_crit);
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

	for idx_plt = 1:8
		subplot_tight(2,4,idx_plt,[.08 .04])
		h_n = heatmap(flipud(htbl(:,:,idx_plt,m)));
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
		title(titleStr{idx_plt})
	end

	sgtitle(sprintf('%s (N=%i, p_{thr}=%.0d)',monkeyName{m},numel(UberRIStat),p_crit),'FontSize',32)
	currPos = get(h_map(m),'Position'); set(h_map(m),'Position',[currPos(1),currPos(2),1600,800]);


	% Distributions
	idx_a1 = find(cellfun(@(x) strcmp(x,'A1'), {UberRIStat.area}));
	idx_r = find(cellfun(@(x) strcmp(x,'R'), {UberRIStat.area}));
	idx_bm = find(cellfun(@(x) strcmp(x,'BeltM'), {UberRIStat.area}));
	idx_bl = find(cellfun(@(x) strcmp(x,'BeltL'), {UberRIStat.area}));
	idx_core = sort([idx_a1 idx_r]);
	idx_belt = sort([idx_bm idx_bl]);

	[h_p, p_p, ci_p, stat_p] = ttest2([UberRIStat(idx_core).p],[UberRIStat(idx_belt).p]);
    if p_p <.05
        sig_p = '*';
    else
        sig_p = '';
    end
	[h_pd, p_pd, ci_pd, stat_pd] = ttest2([UberRIStat(idx_core).PD],[UberRIStat(idx_belt).PD]);
    if p_pd <.05
        sig_pd = '*';
    else
        sig_pd = '';
    end
	[h_t10, p_t10, ci_t10, stat_t10] = ttest2([UberRIStat(idx_core).t10],[UberRIStat(idx_belt).t10]);
    if p_t10 <.05
        sig_t10 = '*';
    else
        sig_t10 = '';
    end
	[h_t50, p_t50, ci_t50, stat_t50] = ttest2([UberRIStat(idx_core).PD],[UberRIStat(idx_belt).PD]);
    if p_t50 <.05
        sig_t50 = '*';
    else
        sig_t50 = '';
    end
	[h_dur, p_dur, ci_dur, stat_dur] = ttest2([UberRIStat(idx_core).Dur],[UberRIStat(idx_belt).Dur]);
    if p_dur <.05
        sig_dur = '*';
    else
        sig_dur = '';
    end
	[h_pli, p_pli, ci_pli, stat_pli] = ttest2([UberRIStat(idx_core).PLI],[UberRIStat(idx_belt).PLI]);
    if p_pli <.05
        sig_pli = '*';
    else
        sig_pli = '';
    end
	[h_dsi, p_dsi, ci_dsi, stat_dsi] = ttest2([UberRIStat(idx_core).DSI],[UberRIStat(idx_belt).DSI]);
    if p_dsi <.05
        sig_dsi = '*';
    else
        sig_dsi = '';
    end

	Edges_p = [-35:1:max(max(log([UberRIStat(idx_core).p])),max(log([UberRIStat(idx_belt).p])))];
	Edges_BF = [0:1000:max(max([UberRIStat(idx_core).BF]),max([UberRIStat(idx_belt).BF]))];
	Edges_PD = [0:5:max(max([UberRIStat(idx_core).PD]),max([UberRIStat(idx_belt).PD]))];
	Edges_t10 = [0:5:max(max([UberRIStat(idx_core).t10]),max([UberRIStat(idx_belt).t10]))];
	Edges_t50 = [0:5:max(max([UberRIStat(idx_core).t50]),max([UberRIStat(idx_belt).t50]))];
	Edges_Dur = [0:5:max(max([UberRIStat(idx_core).Dur]),max([UberRIStat(idx_belt).Dur]))];
	Edges_PLI = [0:.02:max(max([UberRIStat(idx_core).PLI]),max([UberRIStat(idx_belt).PLI]))];
	Edges_DSI = [min(min([UberRIStat(idx_core).DSI]),min([UberRIStat(idx_belt).DSI])):.05:max(max([UberRIStat(idx_core).DSI]),max([UberRIStat(idx_belt).DSI]))];




	legendLabel = {'Core','Belt'};

	h_dist(m) = figure;

	subplot_tight(2,4,1,[.08 .04])
	h1 = histogram(log10([UberRIStat(idx_core).p]),Edges_p,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
	hold on
	h2 = histogram(log10([UberRIStat(idx_belt).p]),Edges_p,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian(log([UberRIStat(idx_core).p])),'b-.',sprintf('Median(Core) = %.1f',nanmedian(log([UberRIStat(idx_core).p]))),'LineWidth',4,'FontSize',16)
	xline(nanmedian(log([UberRIStat(idx_belt).p])),'b-.',sprintf('Median(Belt) = %.1f',nanmedian(log([UberRIStat(idx_belt).p]))),'LineWidth',4,'FontSize',16)
	legend([h1 h2],legendLabel)
	ylabel('Proportion of Clusters')
	title([titleStr{2},sprintf('\np%s = %.2d (two-sample t-test)',sig_p,p_p)])
	
	subplot_tight(2,4,2,[.08 .04])
	h1 = histogram([UberRIStat(idx_core).BF],Edges_BF,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
	hold on
	h2 = histogram([UberRIStat(idx_belt).BF],Edges_BF,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
	legend([h1 h2],legendLabel)
	xlabel('[Hz]')
	ylabel('Proportion of Clusters')
	title(titleStr{3})
	
	subplot_tight(2,4,3,[.08 .04])
	h1 = histogram([UberRIStat(idx_core).PD],Edges_PD,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
	hold on
	h2 = histogram([UberRIStat(idx_belt).PD],Edges_PD,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian([UberRIStat(idx_core).PD]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([UberRIStat(idx_core).PD])),'LineWidth',4,'FontSize',16)
	xline(nanmedian([UberRIStat(idx_belt).PD]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([UberRIStat(idx_belt).PD])),'LineWidth',4,'FontSize',16)
	legend([h1 h2],legendLabel)
	xlabel('[ms]')
	ylabel('Proportion of Clusters')
	title([titleStr{4},sprintf('\np%s = %.2d (two-sample t-test)',sig_pd,p_pd)])
	
	subplot_tight(2,4,4,[.08 .04])
	h1 = histogram([UberRIStat(idx_core).t10],Edges_t10,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
	hold on
	h2 = histogram([UberRIStat(idx_belt).t10],Edges_t10,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian([UberRIStat(idx_core).t10]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([UberRIStat(idx_core).t10])),'LineWidth',4,'FontSize',16)
	xline(nanmedian([UberRIStat(idx_belt).t10]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([UberRIStat(idx_belt).t10])),'LineWidth',4,'FontSize',16)
	legend([h1 h2],legendLabel)
	xlabel('[ms]')
	ylabel('Proportion of Clusters')
	title([titleStr{5},sprintf('\np%s = %.2d (two-sample t-test)',sig_t10,p_t10)])
	
	subplot_tight(2,4,5,[.08 .04])
	h1 = histogram([UberRIStat(idx_core).t50],Edges_t50,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
	hold on
	h2 = histogram([UberRIStat(idx_belt).t50],Edges_t50,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian([UberRIStat(idx_core).t50]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([UberRIStat(idx_core).t50])),'LineWidth',4,'FontSize',16)
	xline(nanmedian([UberRIStat(idx_belt).t50]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([UberRIStat(idx_belt).t50])),'LineWidth',4,'FontSize',16)
	legend([h1 h2],legendLabel)
	xlabel('[ms]')
	ylabel('Proportion of Clusters')
	title([titleStr{6},sprintf('\np%s = %.2d (two-sample t-test)',sig_t50,p_t50)])
	
	subplot_tight(2,4,6,[.08 .04])
	h1 = histogram([UberRIStat(idx_core).Dur],Edges_Dur,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
	hold on
	h2 = histogram([UberRIStat(idx_belt).Dur],Edges_Dur,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian([UberRIStat(idx_core).Dur]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([UberRIStat(idx_core).Dur])),'LineWidth',4,'FontSize',16)
	xline(nanmedian([UberRIStat(idx_belt).Dur]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([UberRIStat(idx_belt).Dur])),'LineWidth',4,'FontSize',16)
	legend([h1 h2],legendLabel)
	xlabel('[ms]')
	ylabel('Proportion of Clusters')
	title([titleStr{7},sprintf('\np%s = %.2d (two-sample t-test)',sig_dur,p_dur)])
	
	subplot_tight(2,4,7,[.08 .04])
	h1 = histogram([UberRIStat(idx_core).PLI],Edges_PLI,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
	hold on
	h2 = histogram([UberRIStat(idx_belt).PLI],Edges_PLI,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian([UberRIStat(idx_core).PLI]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([UberRIStat(idx_core).PLI])),'LineWidth',4,'FontSize',16)
	xline(nanmedian([UberRIStat(idx_belt).PLI]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([UberRIStat(idx_belt).PLI])),'LineWidth',4,'FontSize',16)
	legend([h1 h2],legendLabel)
	ylabel('Proportion of Clusters')
	title([titleStr{8},sprintf('\np%s = %.2d (two-sample t-test)',sig_pli,p_pli)])
	
	subplot_tight(2,4,8,[.08 .04])
	h1 = histogram([UberRIStat(idx_core).DSI],Edges_DSI,'Normalization','probability','FaceColor','r','FaceAlpha',.5);
	hold on
	h2 = histogram([UberRIStat(idx_belt).DSI],Edges_DSI,'Normalization','probability','FaceColor','g','FaceAlpha',.5);
	xline(nanmedian([UberRIStat(idx_core).DSI]),'b-.',sprintf('Median(Core) = %.1f',nanmedian([UberRIStat(idx_core).DSI])),'LineWidth',4,'FontSize',16)
	xline(nanmedian([UberRIStat(idx_belt).DSI]),'b-.',sprintf('Median(Belt) = %.1f',nanmedian([UberRIStat(idx_belt).DSI])),'LineWidth',4,'FontSize',16)
	legend([h1 h2],legendLabel)
	ylabel('Proportion of Clusters')
	title([titleStr{9},sprintf('\np%s = %.2d (two-sample t-test)',sig_dsi,p_dsi)])
	
	sgtitle(sprintf('%s (N=%i, p_{thr}=%.0d)',monkeyName{m},numel(UberRIStat),p_crit),'FontSize',32)
	currPos = get(h_dist(m),'Position'); set(h_dist(m),'Position',[currPos(1),currPos(2),1600,800]);

	% Create reference area map
	htbl_area = ones(15,15)*-1;
	if strcmp(monkeyName{m},'Miyagi')
		coord_core = [coord_ac(2:18,1); coord_ac(2:25,2)];	
		coord_belt = [coord_ac(2:31,3); coord_ac(2:16,4)];
	else
		coord_core = [coord_ac(2:23,1); coord_ac(2:19,2)];	
		coord_belt = [coord_ac(2:35,3); coord_ac(2:24,4)];
	end

	for i = 1:numel(coord_core)
		htbl_area(coord_core{i}(2)+8,coord_core{i}(1)+8) = 1;
	end
	for i = 1:numel(coord_belt)
		htbl_area(coord_belt{i}(2)+8,coord_belt{i}(1)+8) = 0;
	end

	cmap = [0 0 0; 0 1 0; 1 0 0];
	figure;
	heatmap(flipud(htbl_area),'CellLabelColor','none')
	grid off
	colormap(cmap);
	colorbar('off')
	set(gca,'XDisplayLabels',coordLabels_x)
	set(gca,'YDisplayLabels',coordLabels_y)
	if strcmp(monkeyName{m},'Cassius')
		xlabel('Medial <----> Lateral')
	else
		xlabel('Lateral <----> Medial')
	end
	ylabel('Posterior <----> Anterior')
	sgtitle(sprintf('Area Map for %s',monkeyName{m}),'FontSize',20)
	currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),412,372]);
end

