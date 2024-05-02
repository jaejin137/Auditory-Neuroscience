monkeyName = {'Miyagi','Cassius'};
%monkeyName = {'Cassius','Miyagi'};
p_crit = .01;	% Default(=.01)
%p_crit = 1.e-6;

dir_RISTRF = fullfile('~','STRF','RISTRF');
if ~exist(dir_RISTRF,'dir')
	mkdir(dir_RISTRF)
end

tag = datestr(now,'yymmddHHMM');

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
	
	
	SigRISTRF = struct([]);
	UberRISTRF = struct([]);
	gldx_uber = 1;
	gldx_sig = 1;
	% Construct RI stats
	for ldx_RI = 1:length(list_RI)
		% session date
		nameTok = split(list_RI(ldx_RI).name,'_');
		dateTok = split(nameTok{2},'-');
		RIStat(ldx_RI).date = dateTok{2};
		% coordinate
		RIStat(ldx_RI).coord = lamProf(find([lamProf.date]==str2num(RIStat(ldx_RI).date))).coord;
		% area name
		try
			[row, col] = ind2sub(size(coord_ac),find(cellfun(@(x) isequal(x,RIStat(ldx_RI).coord), coord_ac)));
			RIStat(ldx_RI).area = coord_ac{1,col};
		end

		% Load RI
		% number of clusters and median p value thresholded by p_crit
		load(fullfile(dir_RI,list_RI(ldx_RI).name),'STRFSig');
		idx_valid = find(cellfun(@(x) ~isempty(x), {STRFSig.clustNum}))';
		idx_clust_sig = find([STRFSig(idx_valid).p] <= p_crit)';
		RIStat(ldx_RI).numClust = numel(idx_clust_sig);
		RIStat(ldx_RI).p_median = nanmedian([STRFSig(find([STRFSig.p] <= p_crit)).p]');

		% Load STRF parameters
		try
			load(fullfile(dir_Param,sprintf('STRFParams_%s_%s_%s_%s_%s.mat',monkeyName{m},dateTok{2},nameTok{3:5})),'UberSTRF')
		catch
			fprintf('\n!!!Warning: UberSTRF not available for %s_%s_%s-%s-%s!!!\n',monkeyName{m},dateTok{2},nameTok{3:5})
		end
			
		% Process stats for all clusters regardless of RI
		for ldx_valid = 1:numel(idx_valid)
			UberRISTRF(gldx_uber).date = RIStat(ldx_RI).date;
			UberRISTRF(gldx_uber).coord = RIStat(ldx_RI).coord;
			UberRISTRF(gldx_uber).area = RIStat(ldx_RI).area;
			UberRISTRF(gldx_uber).clust = STRFSig(idx_valid(ldx_valid)).clustNum;
			UberRISTRF(gldx_uber).p = STRFSig(idx_valid(ldx_valid)).p;
			try
				idx_param = find(cellfun(@(x) isequal(x,STRFSig(idx_valid(ldx_valid)).clustNum),{UberSTRF.ClustNum}));
				if ~isempty(idx_param)
					UberRISTRF(gldx_uber).BF = UberSTRF(idx_param).RFParam.BFHz;
					UberRISTRF(gldx_uber).PD = UberSTRF(idx_param).RFParam.PeakDelay;
					UberRISTRF(gldx_uber).t10 = UberSTRF(idx_param).RFParam.t1_10;
					UberRISTRF(gldx_uber).t50 = UberSTRF(idx_param).RFParam.t1_50;
					UberRISTRF(gldx_uber).Dur = UberSTRF(idx_param).RFParam.Duration;
					UberRISTRF(gldx_uber).PLI = UberSTRF(idx_param).RFParam.PLI;
					UberRISTRF(gldx_uber).DSI = UberSTRF(idx_param).RFParam.DSI;
				end
			end
			gldx_uber = gldx_uber+1;
		end

		% Process stats for clusters with significant RI
		for ldx_valid = 1:numel(idx_valid)
			if STRFSig(idx_valid(ldx_valid)).p <= p_crit
				SigRISTRF(gldx_sig).date = RIStat(ldx_RI).date;
				SigRISTRF(gldx_sig).coord = RIStat(ldx_RI).coord;
				SigRISTRF(gldx_sig).area = RIStat(ldx_RI).area;
				SigRISTRF(gldx_sig).clust = STRFSig(idx_valid(ldx_valid)).clustNum;
				SigRISTRF(gldx_sig).p = STRFSig(idx_valid(ldx_valid)).p;
				try
					idx_param = find(cellfun(@(x) isequal(x,STRFSig(idx_valid(ldx_valid)).clustNum),{UberSTRF.ClustNum}));
					if ~isempty(idx_param)
						SigRISTRF(gldx_sig).BF = UberSTRF(idx_param).RFParam.BFHz;
						SigRISTRF(gldx_sig).PD = UberSTRF(idx_param).RFParam.PeakDelay;
						SigRISTRF(gldx_sig).t10 = UberSTRF(idx_param).RFParam.t1_10;
						SigRISTRF(gldx_sig).t50 = UberSTRF(idx_param).RFParam.t1_50;
						SigRISTRF(gldx_sig).Dur = UberSTRF(idx_param).RFParam.Duration;
						SigRISTRF(gldx_sig).PLI = UberSTRF(idx_param).RFParam.PLI;
						SigRISTRF(gldx_sig).DSI = UberSTRF(idx_param).RFParam.DSI;
					end
				end
				gldx_sig = gldx_sig+1;
			end
		end
	end
	
	
	% Unique coordinates for all clusters
	coord_uniq_uber = {};
	for ldx_uber = 1:numel(UberRISTRF)
		if ldx_uber == 1
			coord_uniq_uber{1,1} = UberRISTRF(ldx_uber).coord;
		else
			if ~any(cellfun(@(x) isequal(x,UberRISTRF(ldx_uber).coord), coord_uniq_uber))
				coord_uniq_uber{end+1,1} = UberRISTRF(ldx_uber).coord;
			end
		end
	end

	% Unique coordinates for clusters with significant RI
	coord_uniq_sig = {};
	for ldx_sig = 1:numel(SigRISTRF)
		if ldx_sig == 1
			coord_uniq_sig{1,1} = SigRISTRF(ldx_sig).coord;
		else
			if ~any(cellfun(@(x) isequal(x,SigRISTRF(ldx_sig).coord), coord_uniq_sig))
				coord_uniq_sig{end+1,1} = SigRISTRF(ldx_sig).coord;
			end
		end
	end
	
	% Map RI and STRF
	for k = 1:numel(coord_uniq_sig)
		coord_match_uber{k} = find(cellfun(@(x) isequal(x,coord_uniq_uber{k}), {UberRISTRF.coord}));
		coord_match_sig{k} = find(cellfun(@(x) isequal(x,coord_uniq_sig{k}), {SigRISTRF.coord}));
		[row,col] = ind2sub(size(coord_ac),find(cellfun(@(x) isequal(x,coord_uniq_uber{k}), coord_ac)));
		RISTRFMap(k).coord = coord_uniq_uber{k};
		if any(strcmp(coord_ac{1,col},{'A1','R'}))
			RISTRFMap(k).area = "core";
		elseif any(strcmp(coord_ac{1,col},{'BeltL','BeltM'}))
			RISTRFMap(k).area = "belt";
		else
			RISTRFMap(k).area = "unknown";
		end
		RISTRFMap(k).numAllClust = numel(coord_match_uber{k});
		RISTRFMap(k).numSigClust = numel(coord_match_sig{k});
		RISTRFMap(k).p_median = nanmedian([UberRISTRF(coord_match_uber{k}).p]);
		RISTRFMap(k).BF = nanmean([UberRISTRF(coord_match_sig{k}).BF]);
		RISTRFMap(k).PD = nanmean([UberRISTRF(coord_match_sig{k}).PD]);
		RISTRFMap(k).t10 = nanmean([UberRISTRF(coord_match_sig{k}).t10]);
		RISTRFMap(k).t50 = nanmean([UberRISTRF(coord_match_sig{k}).t50]);
		RISTRFMap(k).Dur = nanmean([UberRISTRF(coord_match_sig{k}).Dur]);
		RISTRFMap(k).PLI = nanmean([UberRISTRF(coord_match_sig{k}).PLI]);
		RISTRFMap(k).DSI = nanmean([UberRISTRF(coord_match_sig{k}).DSI]);

		% Construct tables for heatmap
		%htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,1,m) = RISTRFMap(k).numClust/numel(UberRISTRF)*100;
		htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,1,m) = RISTRFMap(k).numSigClust/RISTRFMap(k).numAllClust*100;
		htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,2,m) = log10(RISTRFMap(k).p_median);

		try
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,3,m) = RISTRFMap(k).BF;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,4,m) = RISTRFMap(k).PD;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,5,m) = RISTRFMap(k).t10;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,6,m) = RISTRFMap(k).t50;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,7,m) = RISTRFMap(k).Dur;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,8,m) = RISTRFMap(k).PLI;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,9,m) = RISTRFMap(k).DSI;
		end
	end

	%% Basic statistics
	idx_loc_core = find([RISTRFMap.area]=="core");
	idx_loc_belt = find([RISTRFMap.area]=="belt");
	if m==1
		ratioSigClust_Core_Miyagi = nanmean([RISTRFMap(idx_loc_core).numSigClust]./[RISTRFMap(idx_loc_core).numAllClust]);
		ratioSigClust_Belt_Miyagi = nanmean([RISTRFMap(idx_loc_belt).numSigClust]./[RISTRFMap(idx_loc_belt).numAllClust]);
		pMedian_Core_Miyagi = nanmean([RISTRFMap(idx_loc_core).p_median]);
		pMedian_Belt_Miyagi = nanmean([RISTRFMap(idx_loc_belt).p_median]);
		t10_Core_Miyagi = nanmean([RISTRFMap(idx_loc_core).t10]);
		t10_Belt_Miyagi = nanmean([RISTRFMap(idx_loc_belt).t10]);
	else
		ratioSigClust_Core_Cassius = nanmean([RISTRFMap(idx_loc_core).numSigClust]./[RISTRFMap(idx_loc_core).numAllClust]);
		ratioSigClust_Belt_Cassius = nanmean([RISTRFMap(idx_loc_belt).numSigClust]./[RISTRFMap(idx_loc_belt).numAllClust]);
		pMedian_Core_Cassius = nanmean([RISTRFMap(idx_loc_core).p_median]);
		pMedian_Belt_Cassius = nanmean([RISTRFMap(idx_loc_belt).p_median]);
		t10_Core_Cassius = nanmean([RISTRFMap(idx_loc_core).t10]);
		t10_Belt_Cassius = nanmean([RISTRFMap(idx_loc_belt).t10]);
	end


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

