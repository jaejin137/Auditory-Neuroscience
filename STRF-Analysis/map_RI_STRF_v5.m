%monkeyName = {'Miyagi','Cassius'};
monkeyName = {'Cassius','Miyagi'};
p_crit = .01;	% Default(=.01)
%p_crit = 1.e-6;

dir_RISTRF = fullfile('~','STRF','RISTRF');
if ~exist(dir_RISTRF,'dir')
	mkdir(dir_RISTRF)
end

dir_Casmiya = fullfile('~','Casmiya');
dir_STRF = fullfile('~','STRF');
dir_RI = fullfile(dir_STRF,'RI');
dir_Param = fullfile(dir_STRF,'STRFParams');
	
tag = datestr(now,'yymmddHHMM');

htbl = zeros(15,15,14,2);
htbl_supra = zeros(15,15,14,2);
htbl_infra = zeros(15,15,14,2);

for m = 1:2
	
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

		% Load channel information file
		fname_clustInfo = fullfile('~','kiloSorted_DMR',sprintf('Mr%s-%s/%s_%s_%s/KS2_7_AC/ClusterInfo/cluster_info_new.tsv',monkeyName{m},dateTok{2},nameTok{3:5}));
		clustInfo = importfile_clustinfo(fname_clustInfo);
			
		% Process stats for all clusters regardless of RI
		for ldx_valid = 1:numel(idx_valid)
			UberRISTRF(gldx_uber).date = RIStat(ldx_RI).date;
			UberRISTRF(gldx_uber).coord = RIStat(ldx_RI).coord;
			UberRISTRF(gldx_uber).area = RIStat(ldx_RI).area;
			UberRISTRF(gldx_uber).clust = STRFSig(idx_valid(ldx_valid)).clustNum;
			% Assign channel number to each cluster.
			UberRISTRF(gldx_uber).chan = clustInfo.ch(find(clustInfo.id==STRFSig(idx_valid(ldx_valid)).clustNum));
			% Classify each cluster either as supra or infra layer.
			if UberRISTRF(gldx_uber).chan < 12
				UberRISTRF(gldx_uber).layer = 'supra';
			elseif UberRISTRF(gldx_uber).chan >= 12
				UberRISTRF(gldx_uber).layer = 'infra';
			else
				UberRISTRF(gldx_uber).layer = 'unknown';
			end
			UberRISTRF(gldx_uber).p = STRFSig(idx_valid(ldx_valid)).p;
			try
				idx_param = find(cellfun(@(x) isequal(x,STRFSig(idx_valid(ldx_valid)).clustNum),{UberSTRF.ClustNum}));
				if ~isempty(idx_param)
					UberRISTRF(gldx_uber).BF = UberSTRF(idx_param).RFParam.BF;
					UberRISTRF(gldx_uber).BFHz = UberSTRF(idx_param).RFParam.BFHz;
					UberRISTRF(gldx_uber).Delay = UberSTRF(idx_param).RFParam.Delay;
					UberRISTRF(gldx_uber).PD = UberSTRF(idx_param).RFParam.PeakDelay;
					UberRISTRF(gldx_uber).t10 = UberSTRF(idx_param).RFParam.t1_10;
					UberRISTRF(gldx_uber).t50 = UberSTRF(idx_param).RFParam.t1_50;
					UberRISTRF(gldx_uber).Dur = UberSTRF(idx_param).RFParam.Duration;
					UberRISTRF(gldx_uber).PLI = UberSTRF(idx_param).RFParam.PLI;
					UberRISTRF(gldx_uber).DSI = UberSTRF(idx_param).RFParam.DSI;
					UberRISTRF(gldx_uber).BW = UberSTRF(idx_param).RFParam.BW;
					UberRISTRF(gldx_uber).BWHz = UberSTRF(idx_param).RFParam.BWHz;
					UberRISTRF(gldx_uber).bSMF = UberSTRF(idx_param).RFParam.bSMF;
					UberRISTRF(gldx_uber).FmUpperCutoff = UberSTRF(idx_param).RFParam.FmUpperCutoff;
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
				% Assign channel number to each cluster.
				SigRISTRF(gldx_sig).chan = clustInfo.ch(find(clustInfo.id==STRFSig(idx_valid(ldx_valid)).clustNum));
				% Classify each cluster either as supra or infra layer.
				if SigRISTRF(gldx_sig).chan < 12
					SigRISTRF(gldx_sig).layer = 'supra';
				elseif SigRISTRF(gldx_sig).chan >= 12
					SigRISTRF(gldx_sig).layer = 'infra';
				else
					SigRISTRF(gldx_sig).layer = 'unknown';
				end
				SigRISTRF(gldx_sig).p = STRFSig(idx_valid(ldx_valid)).p;
				try
					idx_param = find(cellfun(@(x) isequal(x,STRFSig(idx_valid(ldx_valid)).clustNum),{UberSTRF.ClustNum}));
					if ~isempty(idx_param)
						SigRISTRF(gldx_sig).BF = UberSTRF(idx_param).RFParam.BF;
						SigRISTRF(gldx_sig).BFHz = UberSTRF(idx_param).RFParam.BFHz;
						SigRISTRF(gldx_sig).Delay = UberSTRF(idx_param).RFParam.Delay;
						SigRISTRF(gldx_sig).PD = UberSTRF(idx_param).RFParam.PeakDelay;
						SigRISTRF(gldx_sig).t10 = UberSTRF(idx_param).RFParam.t1_10;
						SigRISTRF(gldx_sig).t50 = UberSTRF(idx_param).RFParam.t1_50;
						SigRISTRF(gldx_sig).Dur = UberSTRF(idx_param).RFParam.Duration;
						SigRISTRF(gldx_sig).PLI = UberSTRF(idx_param).RFParam.PLI;
						SigRISTRF(gldx_sig).DSI = UberSTRF(idx_param).RFParam.DSI;
						SigRISTRF(gldx_sig).BW = UberSTRF(idx_param).RFParam.BW;
						SigRISTRF(gldx_sig).BWHz = UberSTRF(idx_param).RFParam.BWHz;
						SigRISTRF(gldx_sig).bSMF = UberSTRF(idx_param).RFParam.bSMF;
						SigRISTRF(gldx_sig).FmUpperCutoff = UberSTRF(idx_param).RFParam.FmUpperCutoff;
					end
				end
				gldx_sig = gldx_sig+1;
			end
		end
	end
	
	% Indices of each area
	idx_a1 = find(cellfun(@(x) strcmp(x,'A1'), {SigRISTRF.area}));
	idx_r = find(cellfun(@(x) strcmp(x,'R'), {SigRISTRF.area}));
	idx_bm = find(cellfun(@(x) strcmp(x,'BeltM'), {SigRISTRF.area}));
	idx_bl = find(cellfun(@(x) strcmp(x,'BeltL'), {SigRISTRF.area}));
	idx_core = sort([idx_a1 idx_r])';
	idx_belt = sort([idx_bm idx_bl])';
	idx_supra = find(cellfun(@(x) strcmp(x,'supra'), {SigRISTRF.layer}))';
	idx_infra = find(cellfun(@(x) strcmp(x,'infra'), {SigRISTRF.layer}))';
	idx_core_supra = intersect(idx_core,idx_supra);
	idx_core_infra = intersect(idx_core,idx_infra);
	idx_belt_supra = intersect(idx_belt,idx_supra);
	idx_belt_infra = intersect(idx_belt,idx_infra);

	
	% Unique coordinates for all clusters
	coord_uniq_uber = {};
	for ldx_uber = 1:numel(UberRISTRF)
		if ldx_uber == 1 && ~isempty(UberRISTRF(ldx_uber).coord) 
			coord_uniq_uber{1,1} = UberRISTRF(ldx_uber).coord;
		else
			if ~isempty(UberRISTRF(ldx_uber).coord) && ~any(cellfun(@(x) isequal(x,UberRISTRF(ldx_uber).coord), coord_uniq_uber))
				coord_uniq_uber{end+1,1} = UberRISTRF(ldx_uber).coord;
			end
		end
	end

	% Unique coordinates for clusters with significant RI
	coord_uniq_sig = {};
	for ldx_sig = 1:numel(SigRISTRF)
		if ldx_sig == 1 && ~isempty(SigRISTRF(ldx_sig).coord)
			coord_uniq_sig{1,1} = SigRISTRF(ldx_sig).coord;
		else
			if ~isempty(SigRISTRF(ldx_sig).coord) && ~any(cellfun(@(x) isequal(x,SigRISTRF(ldx_sig).coord), coord_uniq_sig))
				coord_uniq_sig{end+1,1} = SigRISTRF(ldx_sig).coord;
			end
		end
	end
	
	% Map RI and STRF
	for k = 1:numel(coord_uniq_sig)
		coord_match_uber{k,1} = find(cellfun(@(x) isequal(x,coord_uniq_uber{k}), {UberRISTRF.coord}))';
		coord_match_sig{k,1} = find(cellfun(@(x) isequal(x,coord_uniq_sig{k}), {SigRISTRF.coord}))';
		[row,col] = ind2sub(size(coord_ac),find(cellfun(@(x) isequal(x,coord_uniq_uber{k}), coord_ac)));
		RISTRFMap(k).coord = coord_uniq_uber{k};
		RISTRFMap_Supra(k).coord = coord_uniq_uber{k};
		RISTRFMap_Infra(k).coord = coord_uniq_uber{k};
		if any(strcmp(coord_ac{1,col},{'A1','R'}))
			RISTRFMap(k).area = "core";
			RISTRFMap_Supra(k).area = "core";
			RISTRFMap_Infra(k).area = "core";
		elseif any(strcmp(coord_ac{1,col},{'BeltL','BeltM'}))
			RISTRFMap(k).area = "belt";
			RISTRFMap_Supra(k).area = "belt";
			RISTRFMap_Infra(k).area = "belt";
		else
			RISTRFMap(k).area = "unknown";
			RISTRFMap_Supra(k).area = "unknown";
			RISTRFMap_Infra(k).area = "unknown";
		end
		RISTRFMap(k).numAllClust = numel(coord_match_uber{k});
		RISTRFMap(k).numSigClust = numel(coord_match_sig{k});
		RISTRFMap(k).p_median = nanmedian([UberRISTRF(coord_match_sig{k}).p]);
		RISTRFMap(k).BF = nanmean([UberRISTRF(coord_match_sig{k}).BF]);
		RISTRFMap(k).BFHz = nanmean([UberRISTRF(coord_match_sig{k}).BFHz]);
		RISTRFMap(k).Delay = nanmean([UberRISTRF(coord_match_sig{k}).Delay]);
		RISTRFMap(k).PD = nanmean([UberRISTRF(coord_match_sig{k}).PD]);
		RISTRFMap(k).t10 = nanmean([UberRISTRF(coord_match_sig{k}).t10]);
		RISTRFMap(k).t50 = nanmean([UberRISTRF(coord_match_sig{k}).t50]);
		RISTRFMap(k).Dur = nanmean([UberRISTRF(coord_match_sig{k}).Dur]);
		RISTRFMap(k).PLI = nanmean([UberRISTRF(coord_match_sig{k}).PLI]);
		RISTRFMap(k).DSI = nanmean([UberRISTRF(coord_match_sig{k}).DSI]);
		RISTRFMap(k).BW = nanmean([UberRISTRF(coord_match_sig{k}).BW]);
		RISTRFMap(k).BWHz = nanmean([UberRISTRF(coord_match_sig{k}).BWHz]);
		RISTRFMap(k).bSMF = nanmean([UberRISTRF(coord_match_sig{k}).bSMF]);
		RISTRFMap(k).FmUpperCutoff = nanmean([UberRISTRF(coord_match_sig{k}).FmUpperCutoff]);

		RISTRFMap_Supra(k).numAllClust = numel(coord_match_uber{k});
		RISTRFMap_Supra(k).numSigClust = numel(intersect(coord_match_sig{k},idx_supra));
		RISTRFMap_Supra(k).p_median = nanmedian([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).p]);
		RISTRFMap_Supra(k).BF = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).BF]);
		RISTRFMap_Supra(k).BFHz = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).BFHz]);
		RISTRFMap_Supra(k).Delay = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).Delay]);
		RISTRFMap_Supra(k).PD = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).PD]);
		RISTRFMap_Supra(k).t10 = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).t10]);
		RISTRFMap_Supra(k).t50 = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).t50]);
		RISTRFMap_Supra(k).Dur = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).Dur]);
		RISTRFMap_Supra(k).PLI = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).PLI]);
		RISTRFMap_Supra(k).DSI = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).DSI]);
		RISTRFMap_Supra(k).BW = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).BW]);
		RISTRFMap_Supra(k).BWHz = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).BWHz]);
		RISTRFMap_Supra(k).bSMF = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).bSMF]);
		RISTRFMap_Supra(k).FmUpperCutoff = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_supra)).FmUpperCutoff]);

		RISTRFMap_Infra(k).numAllClust = numel(coord_match_uber{k});
		RISTRFMap_Infra(k).numSigClust = numel(intersect(coord_match_sig{k},idx_infra));
		RISTRFMap_Infra(k).p_median = nanmedian([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).p]);
		RISTRFMap_Infra(k).BF = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).BF]);
		RISTRFMap_Infra(k).BFHz = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).BFHz]);
		RISTRFMap_Infra(k).Delay = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).Delay]);
		RISTRFMap_Infra(k).PD = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).PD]);
		RISTRFMap_Infra(k).t10 = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).t10]);
		RISTRFMap_Infra(k).t50 = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).t50]);
		RISTRFMap_Infra(k).Dur = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).Dur]);
		RISTRFMap_Infra(k).PLI = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).PLI]);
		RISTRFMap_Infra(k).DSI = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).DSI]);
		RISTRFMap_Infra(k).BW = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).BW]);
		RISTRFMap_Infra(k).BWHz = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).BWHz]);
		RISTRFMap_Infra(k).bSMF = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).bSMF]);
		RISTRFMap_Infra(k).FmUpperCutoff = nanmean([UberRISTRF(intersect(coord_match_sig{k},idx_infra)).FmUpperCutoff]);

		%% Construct tables for heatmap
		% All layers
		htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,1,m) = RISTRFMap(k).numSigClust/RISTRFMap(k).numAllClust*100;
		htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,2,m) = log10(RISTRFMap(k).p_median);
		try
			htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,3,m) = RISTRFMap(k).BF;
			htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,4,m) = RISTRFMap(k).BFHz;
			htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,5,m) = RISTRFMap(k).Delay;
			htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,6,m) = RISTRFMap(k).PD;
			htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,7,m) = RISTRFMap(k).t10;
			htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,8,m) = RISTRFMap(k).t50;
			htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,9,m) = RISTRFMap(k).Dur;
			htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,10,m) = RISTRFMap(k).PLI;
			htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,11,m) = RISTRFMap(k).BW;
			htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,12,m) = RISTRFMap(k).BWHz;
			htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,13,m) = RISTRFMap(k).bSMF;
			htbl(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,14,m) = RISTRFMap(k).FmUpperCutoff;
		end

		% Supra-layers
		htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,1,m) = RISTRFMap_Supra(k).numSigClust/RISTRFMap_Supra(k).numAllClust*100;
		htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,2,m) = log10(RISTRFMap_Supra(k).p_median);
		try
			htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,3,m) = RISTRFMap_Supra(k).BF;
			htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,4,m) = RISTRFMap_Supra(k).BFHz;
			htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,5,m) = RISTRFMap_Supra(k).Delay;
			htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,6,m) = RISTRFMap_Supra(k).PD;
			htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,7,m) = RISTRFMap_Supra(k).t10;
			htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,8,m) = RISTRFMap_Supra(k).t50;
			htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,9,m) = RISTRFMap_Supra(k).Dur;
			htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,10,m) = RISTRFMap_Supra(k).PLI;
			htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,11,m) = RISTRFMap_Supra(k).BW;
			htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,12,m) = RISTRFMap_Supra(k).BWHz;
			htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,13,m) = RISTRFMap_Supra(k).bSMF;
			htbl_supra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,14,m) = RISTRFMap_Supra(k).FmUpperCutoff;
		end

		% Infra-layers
		htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,1,m) = RISTRFMap_Infra(k).numSigClust/RISTRFMap_Infra(k).numAllClust*100;
		htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,2,m) = log10(RISTRFMap_Infra(k).p_median);
		try
			htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,3,m) = RISTRFMap_Infra(k).BF;
			htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,4,m) = RISTRFMap_Infra(k).BFHz;
			htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,5,m) = RISTRFMap_Infra(k).Delay;
			htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,6,m) = RISTRFMap_Infra(k).PD;
			htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,7,m) = RISTRFMap_Infra(k).t10;
			htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,8,m) = RISTRFMap_Infra(k).t50;
			htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,9,m) = RISTRFMap_Infra(k).Dur;
			htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,10,m) = RISTRFMap_Infra(k).PLI;
			htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,11,m) = RISTRFMap_Infra(k).BW;
			htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,12,m) = RISTRFMap_Infra(k).BWHz;
			htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,13,m) = RISTRFMap_Infra(k).bSMF;
			htbl_infra(-coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,14,m) = RISTRFMap_Infra(k).FmUpperCutoff;
		end

	end

	%% Basic statistics
	idx_loc_core = find([RISTRFMap.area]=="core");
	idx_loc_belt = find([RISTRFMap.area]=="belt");
	if m==1
		ratioSigClust_Core_Cassius = nanmedian([RISTRFMap(idx_loc_core).numSigClust]./[RISTRFMap(idx_loc_core).numAllClust]);
		ratioSigClust_Belt_Cassius = nanmedian([RISTRFMap(idx_loc_belt).numSigClust]./[RISTRFMap(idx_loc_belt).numAllClust]);
		pMedian_Core_Cassius = nanmedian([RISTRFMap(idx_loc_core).p_median]);
		pMedian_Belt_Cassius = nanmedian([RISTRFMap(idx_loc_belt).p_median]);
		t10_Core_Cassius = nanmedian([RISTRFMap(idx_loc_core).t10]);
		t10_Belt_Cassius = nanmedian([RISTRFMap(idx_loc_belt).t10]);
	else
		ratioSigClust_Core_Miyagi = nanmedian([RISTRFMap(idx_loc_core).numSigClust]./[RISTRFMap(idx_loc_core).numAllClust]);
		ratioSigClust_Belt_Miyagi = nanmedian([RISTRFMap(idx_loc_belt).numSigClust]./[RISTRFMap(idx_loc_belt).numAllClust]);
		pMedian_Core_Miyagi = nanmedian([RISTRFMap(idx_loc_core).p_median]);
		pMedian_Belt_Miyagi = nanmedian([RISTRFMap(idx_loc_belt).p_median]);
		t10_Core_Miyagi = nanmedian([RISTRFMap(idx_loc_core).t10]);
		t10_Belt_Miyagi = nanmedian([RISTRFMap(idx_loc_belt).t10]);
	end

	%% Plot correlations
	h_corr{m} = figure;
	% Bandwidth(Hz) vs Best Frequency
	subplot_tight(2,3,1,[.1 .06])
	%scatter([SigRISTRF(idx_core).BFHz],[SigRISTRF(idx_core).BWHz],'filled','r')
	%hold on
	%scatter([SigRISTRF(idx_core).BFHz],[SigRISTRF(idx_core).BWHz],'filled','b')
	%corr1_core = corr([SigRISTRF(idx_core).BFHz]',[SigRISTRF(idx_core).BWHz]','rows','complete');
	%corr1_belt = corr([SigRISTRF(idx_belt).BFHz]',[SigRISTRF(idx_belt).BWHz]','rows','complete');
	%p_corrcoeff1 = compare_correlation_coefficients(corr1_core,corr1_belt,length(idx_core),length(idx_belt));
	%if p_corrcoeff1 < .05
	%	sig_r_label1 = '*';
	%else
	%	sig_r_label1 = '';
	%end
	%legend({sprintf('Core (N=%i,r=%.2f)',numel([SigRISTRF(idx_core).BWHz]),corr1_core),sprintf('Belt (N=%i,r=%.2f)%s',numel([SigRISTRF(idx_belt).BWHz]),corr1_belt,sig_r_label1)},'Location','best')
	scatter([SigRISTRF(idx_core_supra).BFHz],[SigRISTRF(idx_core_supra).BWHz],'filled','r')
	hold on
	scatter([SigRISTRF(idx_core_infra).BFHz],[SigRISTRF(idx_core_infra).BWHz],'filled','b')
	scatter([SigRISTRF(idx_belt_supra).BFHz],[SigRISTRF(idx_belt_supra).BWHz],'filled','g')
	scatter([SigRISTRF(idx_belt_infra).BFHz],[SigRISTRF(idx_belt_infra).BWHz],'filled','c')
	NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
	SW = [min(xlim) min(ylim)]+[diff(xlim) diff(ylim)]*0.05;
	NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
	SE = [max(xlim) min(ylim)]+[-diff(xlim) diff(ylim)]*0.05;
	set(gca,'XScale','log')
	set(gca,'YScale','log')
	xlabel('Best Frequency (Hz)')
	ylabel('Bandwidth (Hz)')
	corr1_core_supra = corr([SigRISTRF(idx_core_supra).BFHz]',[SigRISTRF(idx_core_supra).BWHz]','rows','complete');
	corr1_belt_supra = corr([SigRISTRF(idx_belt_supra).BFHz]',[SigRISTRF(idx_belt_supra).BWHz]','rows','complete');
	corr1_core_infra = corr([SigRISTRF(idx_core_infra).BFHz]',[SigRISTRF(idx_core_infra).BWHz]','rows','complete');
	corr1_belt_infra = corr([SigRISTRF(idx_belt_infra).BFHz]',[SigRISTRF(idx_belt_infra).BWHz]','rows','complete');
	p_corrcoeff1_supra = compare_correlation_coefficients(corr1_core_supra,corr1_belt_supra,length(idx_core_supra),length(idx_belt_supra));
	p_corrcoeff1_infra = compare_correlation_coefficients(corr1_core_infra,corr1_belt_infra,length(idx_core_infra),length(idx_belt_infra));
	p_corrcoeff1_core = compare_correlation_coefficients(corr1_core_supra,corr1_core_infra,length(idx_core_supra),length(idx_core_infra));
	p_corrcoeff1_belt = compare_correlation_coefficients(corr1_belt_supra,corr1_belt_infra,length(idx_belt_supra),length(idx_belt_infra));
	if p_corrcoeff1_supra < .05
		sig_r_label1_supra = 'S';
	else
		sig_r_label1_supra = 'NS';
	end
	if p_corrcoeff1_infra < .05
		sig_r_label1_infra = 'S';
	else
		sig_r_label1_infra = 'NS';
	end
	if p_corrcoeff1_core < .05
		sig_r_label1_core = 'S';
	else
		sig_r_label1_core = 'NS';
	end
	if p_corrcoeff1_belt < .05
		sig_r_label1_belt = 'S';
	else
		sig_r_label1_belt = 'NS';
	end
	legend({sprintf('Core\\_sup (N=%i,r=%.2f)',numel([SigRISTRF(idx_core_supra).BWHz]),corr1_core_supra), sprintf('Belt\\_sup (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt_supra).BWHz]),corr1_belt_supra), sprintf('Core\\_inf (N=%i,r=%.2f)',numel([SigRISTRF(idx_core_infra).BWHz]),corr1_core_infra), sprintf('Belt\\_inf (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt_infra).BWHz]),corr1_belt_infra)})
	%text(0,0,{sprintf('Supra-layer\\_Core-vs-Belt:%s',sig_r_label1_supra), sprintf('Infra-layer\\_Core-vs-Belt:%s',sig_r_label1_infra), sprintf('Core\\_Supra-vs-Infra:%s',sig_r_label1_core), sprintf('Belt\\_Supra-vs-Infra:%s',sig_r_label1_belt)})
	title(sprintf('Supra-layer\\_Core-vs-Belt:%s, Infra-layer\\_Core-vs-Belt:%s\nCore\\_Supra-vs-Infra:%s, Belt\\_Supra-vs-Infra:%s', sig_r_label1_supra, sig_r_label1_infra, sig_r_label1_core, sig_r_label1_belt))

	% Bandwidth(octave) vs Best Frequency
	subplot_tight(2,3,2,[.1 .06])
	%scatter([SigRISTRF(idx_core).BFHz],[SigRISTRF(idx_core).BW],'filled','r')
	%hold on
	%scatter([SigRISTRF(idx_belt).BFHz],[SigRISTRF(idx_belt).BW],'filled','b')
	%corr2_core = corr([SigRISTRF(idx_core).BFHz]',[SigRISTRF(idx_core).BW]','rows','complete');
	%corr2_belt = corr([SigRISTRF(idx_belt).BFHz]',[SigRISTRF(idx_belt).BW]','rows','complete');
	%legend({sprintf('Core (N=%i,r=%.2f)',numel([SigRISTRF(idx_core).BW]),corr2_core),sprintf('Belt (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt).BW]),corr2_belt)})
	scatter([SigRISTRF(idx_core_supra).BFHz],[SigRISTRF(idx_core_supra).BW],'filled','r')
	hold on
	scatter([SigRISTRF(idx_core_infra).BFHz],[SigRISTRF(idx_core_infra).BW],'filled','b')
	scatter([SigRISTRF(idx_belt_supra).BFHz],[SigRISTRF(idx_belt_supra).BW],'filled','g')
	scatter([SigRISTRF(idx_belt_infra).BFHz],[SigRISTRF(idx_belt_infra).BW],'filled','c')
	set(gca,'XScale','log')
	xlabel('Best Frequency (Hz)')
	ylabel('Bandwidth (octave)')
	corr2_core_supra = corr([SigRISTRF(idx_core_supra).BFHz]',[SigRISTRF(idx_core_supra).BW]','rows','complete');
	corr2_belt_supra = corr([SigRISTRF(idx_belt_supra).BFHz]',[SigRISTRF(idx_belt_supra).BW]','rows','complete');
	corr2_core_infra = corr([SigRISTRF(idx_core_infra).BFHz]',[SigRISTRF(idx_core_infra).BW]','rows','complete');
	corr2_belt_infra = corr([SigRISTRF(idx_belt_infra).BFHz]',[SigRISTRF(idx_belt_infra).BW]','rows','complete');
	p_corrcoeff2_supra = compare_correlation_coefficients(corr2_core_supra,corr2_belt_supra,length(idx_core_supra),length(idx_belt_supra));
	p_corrcoeff2_infra = compare_correlation_coefficients(corr2_core_infra,corr2_belt_infra,length(idx_core_infra),length(idx_belt_infra));
	p_corrcoeff2_core = compare_correlation_coefficients(corr2_core_supra,corr2_core_infra,length(idx_core_supra),length(idx_core_infra));
	p_corrcoeff2_belt = compare_correlation_coefficients(corr2_belt_supra,corr2_belt_infra,length(idx_belt_supra),length(idx_belt_infra));
	if p_corrcoeff2_supra < .05
		sig_r_label2_supra = 'S';
	else
		sig_r_label2_supra = 'NS';
	end
	if p_corrcoeff2_infra < .05
		sig_r_label2_infra = 'S';
	else
		sig_r_label2_infra = 'NS';
	end
	if p_corrcoeff2_core < .05
		sig_r_label2_core = 'S';
	else
		sig_r_label2_core = 'NS';
	end
	if p_corrcoeff2_belt < .05
		sig_r_label2_belt = 'S';
	else
		sig_r_label2_belt = 'NS';
	end
	legend({sprintf('Core\\_sup (N=%i,r=%.2f)',numel([SigRISTRF(idx_core_supra).BW]),corr2_core_supra), sprintf('Belt\\_sup (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt_supra).BW]),corr2_belt_supra), sprintf('Core\\_inf (N=%i,r=%.2f)',numel([SigRISTRF(idx_core_infra).BW]),corr2_core_infra), sprintf('Belt\\_inf (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt_infra).BW]),corr2_belt_infra)})
	%text(0,0,{sprintf('Supra-layer\\_Core-vs-Belt:%s',sig_r_label2_supra), sprintf('Infra-layer\\_Core-vs-Belt:%s',sig_r_label2_infra), sprintf('Core\\_Supra-vs-Infra:%s',sig_r_label2_core), sprintf('Belt\\_Supra-vs-Infra:%s',sig_r_label2_belt)})
	title(sprintf('Supra-layer\\_Core-vs-Belt:%s, Infra-layer\\_Core-vs-Belt:%s\nCore\\_Supra-vs-Infra:%s, Belt\\_Supra-vs-Infra:%s', sig_r_label2_supra, sig_r_label2_infra, sig_r_label2_core, sig_r_label2_belt))

	% Best Spectral Modulation Frequency vs Bandwidth
	subplot_tight(2,3,3,[.1 .06])
	BW_core_supra = [];
	for j = 1:numel(idx_core_supra)
		if isempty(SigRISTRF(idx_core_supra(j)).BW)
			BW_core_supra(end+1,1) = nan;
		else
			BW_core_supra(end+1,1) = SigRISTRF(idx_core_supra(j)).BW;
		end
	end
	bSMF_core_supra = [];
	for j = 1:numel(idx_core_supra)
		if isempty(SigRISTRF(idx_core_supra(j)).bSMF)
			bSMF_core_supra(end+1,1) = nan;
		else
			bSMF_core_supra(end+1,1) = SigRISTRF(idx_core_supra(j)).bSMF(1);
		end
	end
	BW_core_infra = [];
	for j = 1:numel(idx_core_infra)
		if isempty(SigRISTRF(idx_core_infra(j)).BW)
			BW_core_infra(end+1,1) = nan;
		else
			BW_core_infra(end+1,1) = SigRISTRF(idx_core_infra(j)).BW;
		end
	end
	bSMF_core_infra = [];
	for j = 1:numel(idx_core_infra)
		if isempty(SigRISTRF(idx_core_infra(j)).bSMF)
			bSMF_core_infra(end+1,1) = nan;
		else
			bSMF_core_infra(end+1,1) = SigRISTRF(idx_core_infra(j)).bSMF(1);
		end
	end
	BW_belt_supra = [];
	for j = 1:numel(idx_belt_supra)
		if isempty(SigRISTRF(idx_belt_supra(j)).BW)
			BW_belt_supra(end+1,1) = nan;
		else
			BW_belt_supra(end+1,1) = SigRISTRF(idx_belt_supra(j)).BW;
		end
	end
	bSMF_belt_supra = [];
	for j = 1:numel(idx_belt_supra)
		if isempty(SigRISTRF(idx_belt_supra(j)).bSMF)
			bSMF_belt_supra(end+1,1) = nan;
		else
			bSMF_belt_supra(end+1,1) = SigRISTRF(idx_belt_supra(j)).bSMF(1);
		end
	end
	BW_belt_infra = [];
	for j = 1:numel(idx_belt_infra)
		if isempty(SigRISTRF(idx_belt_infra(j)).BW)
			BW_belt_infra(end+1,1) = nan;
		else
			BW_belt_infra(end+1,1) = SigRISTRF(idx_belt_infra(j)).BW;
		end
	end
	bSMF_belt_infra = [];
	for j = 1:numel(idx_belt_infra)
		if isempty(SigRISTRF(idx_belt_infra(j)).bSMF)
			bSMF_belt_infra(end+1,1) = nan;
		else
			bSMF_belt_infra(end+1,1) = SigRISTRF(idx_belt_infra(j)).bSMF(1);
		end
	end
	%scatter([SigRISTRF(idx_core).BW],[SigRISTRF(idx_core).bSMF],'filled','r')
	scatter(BW_core_supra,bSMF_core_supra,'filled','r')
	hold
	scatter(BW_core_infra,bSMF_core_infra,'filled','b')
	%scatter([SigRISTRF(idx_belt).BW],[SigRISTRF(idx_belt).bSMF],'filled','g')
	scatter(BW_belt_supra,bSMF_belt_supra,'filled','g')
	scatter(BW_belt_infra,bSMF_belt_infra,'filled','c')
	%set(gca,'XScale','log')
	xlabel('Bandwidth (octave)')
	ylabel('Best Spectral Mod. Freq.')
	%corr3_core = corr(BW_core,bSMF_core,'rows','complete');
	%corr3_belt = corr(BW_belt,bSMF_belt,'rows','complete');
	%legend({sprintf('Core (N=%i,r=%.2f)',sum(~isnan(bSMF_core)),corr3_core),sprintf('Belt (N=%i,r=%.2f)',sum(~isnan(bSMF_belt)),corr3_belt)})
	%legend({'core, supra','core, infra','belt, supra','belt, infra'},'Location','best')
	corr3_core_supra = corr(BW_core_supra,bSMF_core_supra,'rows','complete');
	corr3_belt_supra = corr(BW_belt_supra,bSMF_belt_supra,'rows','complete');
	corr3_core_infra = corr(BW_core_infra,bSMF_core_infra,'rows','complete');
	corr3_belt_infra = corr(BW_belt_infra,bSMF_belt_infra,'rows','complete');
	p_corrcoeff3_supra = compare_correlation_coefficients(corr3_core_supra,corr3_belt_supra,length(idx_core_supra),length(idx_belt_supra));
	p_corrcoeff3_infra = compare_correlation_coefficients(corr3_core_infra,corr3_belt_infra,length(idx_core_infra),length(idx_belt_infra));
	p_corrcoeff3_core = compare_correlation_coefficients(corr3_core_supra,corr3_core_infra,length(idx_core_supra),length(idx_core_infra));
	p_corrcoeff3_belt = compare_correlation_coefficients(corr3_belt_supra,corr3_belt_infra,length(idx_belt_supra),length(idx_belt_infra));
	if p_corrcoeff3_supra < .05
		sig_r_label3_supra = 'S';
	else
		sig_r_label3_supra = 'NS';
	end
	if p_corrcoeff3_infra < .05
		sig_r_label3_infra = 'S';
	else
		sig_r_label3_infra = 'NS';
	end
	if p_corrcoeff3_core < .05
		sig_r_label3_core = 'S';
	else
		sig_r_label3_core = 'NS';
	end
	if p_corrcoeff3_belt < .05
		sig_r_label3_belt = 'S';
	else
		sig_r_label3_belt = 'NS';
	end
	legend({sprintf('Core\\_sup (N=%i,r=%.2f)',numel(BW_core_supra),corr3_core_supra), sprintf('Belt\\_sup (N=%i,r=%.2f)',numel(BW_belt_supra),corr3_belt_supra), sprintf('Core\\_inf (N=%i,r=%.2f)',numel(BW_core_infra),corr3_core_infra), sprintf('Belt\\_inf (N=%i,r=%.2f)',numel(BW_belt_infra),corr3_belt_infra)})
	%text(0,0,{sprintf('Supra-layer\\_Core-vs-Belt:%s',sig_r_label3_supra), sprintf('Infra-layer\\_Core-vs-Belt:%s',sig_r_label3_infra), sprintf('Core\\_Supra-vs-Infra:%s',sig_r_label3_core), sprintf('Belt\\_Supra-vs-Infra:%s',sig_r_label3_belt)})
	title(sprintf('Supra-layer\\_Core-vs-Belt:%s, Infra-layer\\_Core-vs-Belt:%s\nCore\\_Supra-vs-Infra:%s, Belt\\_Supra-vs-Infra:%s', sig_r_label3_supra, sig_r_label3_infra, sig_r_label3_core, sig_r_label3_belt))

	% Delay(t1_10 or simply Delay) vs Best Frequency
	subplot_tight(2,3,4,[.1 .06])
	%scatter([SigRISTRF(idx_core).BFHz],[SigRISTRF(idx_core).Delay],'filled','r')
	%hold on
	%scatter([SigRISTRF(idx_belt).BFHz],[SigRISTRF(idx_belt).Delay],'filled','g')
	%corr4_core = corr([SigRISTRF(idx_core).BFHz]',[SigRISTRF(idx_core).Delay]','rows','complete');
	%corr4_belt = corr([SigRISTRF(idx_belt).BFHz]',[SigRISTRF(idx_belt).Delay]','rows','complete');
	%legend({sprintf('Core (N=%i,r=%.2f)',numel([SigRISTRF(idx_core).Delay]),corr4_core),sprintf('Belt (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt).Delay]),corr4_belt)})
	scatter([SigRISTRF(idx_core_supra).BFHz],[SigRISTRF(idx_core_supra).Delay],'filled','r')
	hold on
	scatter([SigRISTRF(idx_core_infra).BFHz],[SigRISTRF(idx_core_infra).Delay],'filled','b')
	scatter([SigRISTRF(idx_belt_supra).BFHz],[SigRISTRF(idx_belt_supra).Delay],'filled','g')
	scatter([SigRISTRF(idx_belt_infra).BFHz],[SigRISTRF(idx_belt_infra).Delay],'filled','c')
	set(gca,'XScale','log')
	set(gca,'XScale','log')
	xlabel('Best Frequency (Hz)')
	ylabel('Delay (ms)')
	corr4_core_supra = corr([SigRISTRF(idx_core_supra).BFHz]',[SigRISTRF(idx_core_supra).Delay]','rows','complete');
	corr4_belt_supra = corr([SigRISTRF(idx_belt_supra).BFHz]',[SigRISTRF(idx_belt_supra).Delay]','rows','complete');
	corr4_core_infra = corr([SigRISTRF(idx_core_infra).BFHz]',[SigRISTRF(idx_core_infra).Delay]','rows','complete');
	corr4_belt_infra = corr([SigRISTRF(idx_belt_infra).BFHz]',[SigRISTRF(idx_belt_infra).Delay]','rows','complete');
	p_corrcoeff4_supra = compare_correlation_coefficients(corr4_core_supra,corr4_belt_supra,length(idx_core_supra),length(idx_belt_supra));
	p_corrcoeff4_infra = compare_correlation_coefficients(corr4_core_infra,corr4_belt_infra,length(idx_core_infra),length(idx_belt_infra));
	p_corrcoeff4_core = compare_correlation_coefficients(corr4_core_supra,corr4_core_infra,length(idx_core_supra),length(idx_core_infra));
	p_corrcoeff4_belt = compare_correlation_coefficients(corr4_belt_supra,corr4_belt_infra,length(idx_belt_supra),length(idx_belt_infra));
	if p_corrcoeff4_supra < .05
		sig_r_label4_supra = 'S';
	else
		sig_r_label4_supra = 'NS';
	end
	if p_corrcoeff4_infra < .05
		sig_r_label4_infra = 'S';
	else
		sig_r_label4_infra = 'NS';
	end
	if p_corrcoeff4_core < .05
		sig_r_label4_core = 'S';
	else
		sig_r_label4_core = 'NS';
	end
	if p_corrcoeff4_belt < .05
		sig_r_label4_belt = 'S';
	else
		sig_r_label4_belt = 'NS';
	end
	legend({sprintf('Core\\_sup (N=%i,r=%.2f)',numel([SigRISTRF(idx_core_supra).Delay]),corr4_core_supra), sprintf('Belt\\_sup (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt_supra).Delay]),corr4_belt_supra), sprintf('Core\\_inf (N=%i,r=%.2f)',numel([SigRISTRF(idx_core_infra).Delay]),corr4_core_infra), sprintf('Belt\\_inf (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt_infra).Delay]),corr4_belt_infra)})
	%text(0,0,{sprintf('Supra-layer\\_Core-vs-Belt:%s',sig_r_label4_supra), sprintf('Infra-layer\\_Core-vs-Belt:%s',sig_r_label4_infra), sprintf('Core\\_Supra-vs-Infra:%s',sig_r_label4_core), sprintf('Belt\\_Supra-vs-Infra:%s',sig_r_label4_belt)})
	title(sprintf('Supra-layer\\_Core-vs-Belt:%s, Infra-layer\\_Core-vs-Belt:%s\nCore\\_Supra-vs-Infra:%s, Belt\\_Supra-vs-Infra:%s', sig_r_label4_supra, sig_r_label4_infra, sig_r_label4_core, sig_r_label4_belt))

	% Temporal Cutoff Frequency vs Integration Time
	subplot_tight(2,3,5,[.1 .06])
	%scatter([SigRISTRF(idx_core).Dur],[SigRISTRF(idx_core).FmUpperCutoff],'filled','r')
	%hold on
	%scatter([SigRISTRF(idx_belt).Dur],[SigRISTRF(idx_belt).FmUpperCutoff],'filled','g')
	%corr5_core = corr([SigRISTRF(idx_core).Dur]',[SigRISTRF(idx_core).FmUpperCutoff]','rows','complete');
	%corr5_belt = corr([SigRISTRF(idx_belt).Dur]',[SigRISTRF(idx_belt).FmUpperCutoff]','rows','complete');
	%legend({sprintf('Core (N=%i,r=%.2f)',numel([SigRISTRF(idx_core).Dur]),corr5_core),sprintf('Belt (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt).Dur]),corr5_belt)})
	scatter([SigRISTRF(idx_core_supra).Dur],[SigRISTRF(idx_core_supra).FmUpperCutoff],'filled','r')
	hold on
	scatter([SigRISTRF(idx_core_infra).Dur],[SigRISTRF(idx_core_infra).FmUpperCutoff],'filled','b')
	scatter([SigRISTRF(idx_belt_supra).Dur],[SigRISTRF(idx_belt_supra).FmUpperCutoff],'filled','g')
	scatter([SigRISTRF(idx_belt_infra).Dur],[SigRISTRF(idx_belt_infra).FmUpperCutoff],'filled','c')
	xlabel('Integration Time (ms)')
	ylabel('Temporal Cutoff Freq. (Hz)')
	corr5_core_supra = corr([SigRISTRF(idx_core_supra).Dur]',[SigRISTRF(idx_core_supra).FmUpperCutoff]','rows','complete');
	corr5_belt_supra = corr([SigRISTRF(idx_belt_supra).Dur]',[SigRISTRF(idx_belt_supra).FmUpperCutoff]','rows','complete');
	corr5_core_infra = corr([SigRISTRF(idx_core_infra).Dur]',[SigRISTRF(idx_core_infra).FmUpperCutoff]','rows','complete');
	corr5_belt_infra = corr([SigRISTRF(idx_belt_infra).Dur]',[SigRISTRF(idx_belt_infra).FmUpperCutoff]','rows','complete');
	p_corrcoeff5_supra = compare_correlation_coefficients(corr5_core_supra,corr5_belt_supra,length(idx_core_supra),length(idx_belt_supra));
	p_corrcoeff5_infra = compare_correlation_coefficients(corr5_core_infra,corr5_belt_infra,length(idx_core_infra),length(idx_belt_infra));
	p_corrcoeff5_core = compare_correlation_coefficients(corr5_core_supra,corr5_core_infra,length(idx_core_supra),length(idx_core_infra));
	p_corrcoeff5_belt = compare_correlation_coefficients(corr5_belt_supra,corr5_belt_infra,length(idx_belt_supra),length(idx_belt_infra));
	if p_corrcoeff5_supra < .05
		sig_r_label5_supra = 'S';
	else
		sig_r_label5_supra = 'NS';
	end
	if p_corrcoeff5_infra < .05
		sig_r_label5_infra = 'S';
	else
		sig_r_label5_infra = 'NS';
	end
	if p_corrcoeff5_core < .05
		sig_r_label5_core = 'S';
	else
		sig_r_label5_core = 'NS';
	end
	if p_corrcoeff5_belt < .05
		sig_r_label5_belt = 'S';
	else
		sig_r_label5_belt = 'NS';
	end
	legend({sprintf('Core\\_sup (N=%i,r=%.2f)',numel([SigRISTRF(idx_core_supra).FmUpperCutoff]),corr5_core_supra), sprintf('Belt\\_sup (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt_supra).FmUpperCutoff]),corr5_belt_supra), sprintf('Core\\_inf (N=%i,r=%.2f)',numel([SigRISTRF(idx_core_infra).FmUpperCutoff]),corr5_core_infra), sprintf('Belt\\_inf (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt_infra).FmUpperCutoff]),corr5_belt_infra)})
	%text(0,0,{sprintf('Supra-layer\\_Core-vs-Belt:%s',sig_r_label5_supra), sprintf('Infra-layer\\_Core-vs-Belt:%s',sig_r_label5_infra), sprintf('Core\\_Supra-vs-Infra:%s',sig_r_label5_core), sprintf('Belt\\_Supra-vs-Infra:%s',sig_r_label5_belt)})
	title(sprintf('Supra-layer\\_Core-vs-Belt:%s, Infra-layer\\_Core-vs-Belt:%s\nCore\\_Supra-vs-Infra:%s, Belt\\_Supra-vs-Infra:%s', sig_r_label5_supra, sig_r_label5_infra, sig_r_label5_core, sig_r_label5_belt))

	% Delay vs Integration Time
	subplot_tight(2,3,6,[.1 .06])
	%scatter([SigRISTRF(idx_core).Dur],[SigRISTRF(idx_core).Delay],'filled','r')
	%hold on
	%scatter([SigRISTRF(idx_belt).Dur],[SigRISTRF(idx_belt).Delay],'filled','g')
	%corr6_core = corr([SigRISTRF(idx_core).Dur]',[SigRISTRF(idx_core).Delay]','rows','complete');
	%corr6_belt = corr([SigRISTRF(idx_belt).Dur]',[SigRISTRF(idx_belt).Delay]','rows','complete');
	%legend({sprintf('Core (N=%i,r=%.2f)',numel([SigRISTRF(idx_core).Dur]),corr6_core),sprintf('Belt (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt).Dur]),corr6_belt)})
	scatter([SigRISTRF(idx_core_supra).Dur],[SigRISTRF(idx_core_supra).Delay],'filled','r')
	hold on
	scatter([SigRISTRF(idx_core_infra).Dur],[SigRISTRF(idx_core_infra).Delay],'filled','b')
	scatter([SigRISTRF(idx_belt_supra).Dur],[SigRISTRF(idx_belt_supra).Delay],'filled','g')
	scatter([SigRISTRF(idx_belt_infra).Dur],[SigRISTRF(idx_belt_infra).Delay],'filled','c')
	xlabel('Integration Time (ms)')
	ylabel('Delay (ms)')
	corr6_core_supra = corr([SigRISTRF(idx_core_supra).Dur]',[SigRISTRF(idx_core_supra).BWHz]','rows','complete');
	corr6_belt_supra = corr([SigRISTRF(idx_belt_supra).Dur]',[SigRISTRF(idx_belt_supra).BWHz]','rows','complete');
	corr6_core_infra = corr([SigRISTRF(idx_core_infra).Dur]',[SigRISTRF(idx_core_infra).BWHz]','rows','complete');
	corr6_belt_infra = corr([SigRISTRF(idx_belt_infra).Dur]',[SigRISTRF(idx_belt_infra).Delay]','rows','complete');
	p_corrcoeff6_supra = compare_correlation_coefficients(corr6_core_supra,corr6_belt_supra,length(idx_core_supra),length(idx_belt_supra));
	p_corrcoeff6_infra = compare_correlation_coefficients(corr6_core_infra,corr6_belt_infra,length(idx_core_infra),length(idx_belt_infra));
	p_corrcoeff6_core = compare_correlation_coefficients(corr6_core_supra,corr6_core_infra,length(idx_core_supra),length(idx_core_infra));
	p_corrcoeff6_belt = compare_correlation_coefficients(corr6_belt_supra,corr6_belt_infra,length(idx_belt_supra),length(idx_belt_infra));
	if p_corrcoeff6_supra < .05
		sig_r_label6_supra = 'S';
	else
		sig_r_label6_supra = 'NS';
	end
	if p_corrcoeff6_infra < .05
		sig_r_label6_infra = 'S';
	else
		sig_r_label6_infra = 'NS';
	end
	if p_corrcoeff6_core < .05
		sig_r_label6_core = 'S';
	else
		sig_r_label6_core = 'NS';
	end
	if p_corrcoeff6_belt < .05
		sig_r_label6_belt = 'S';
	else
		sig_r_label6_belt = 'NS';
	end
	legend({sprintf('Core\\_sup (N=%i,r=%.2f)',numel([SigRISTRF(idx_core_supra).Delay]),corr6_core_supra), sprintf('Belt\\_sup (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt_supra).Delay]),corr6_belt_supra), sprintf('Core\\_inf (N=%i,r=%.2f)',numel([SigRISTRF(idx_core_infra).Delay]),corr6_core_infra), sprintf('Belt\\_inf (N=%i,r=%.2f)',numel([SigRISTRF(idx_belt_infra).Delay]),corr6_belt_infra)})
	%text(0,0,{sprintf('Supra-layer\\_Core-vs-Belt:%s',sig_r_label6_supra), sprintf('Infra-layer\\_Core-vs-Belt:%s',sig_r_label6_infra), sprintf('Core\\_Supra-vs-Infra:%s',sig_r_label6_core), sprintf('Belt\\_Supra-vs-Infra:%s',sig_r_label6_belt)})
	title(sprintf('Supra-layer\\_Core-vs-Belt:%s, Infra-layer\\_Core-vs-Belt:%s\nCore\\_Supra-vs-Infra:%s, Belt\\_Supra-vs-Infra:%s', sig_r_label6_supra, sig_r_label6_infra, sig_r_label6_core, sig_r_label6_belt))

	sgtitle(sprintf('%s (N=%i, p_{thr}=%.0d)',monkeyName{m},numel(SigRISTRF),p_crit),'FontSize',24)
	currPos = get(h_corr{m},'Position'); set(gcf,'Position',[currPos(1),currPos(2),1200,800]);	
	saveas(h_corr{m},sprintf('~/STRF/RISTRF_Figures/RISTRF_Corr_v5_%s_%s.png',monkeyName{m},tag))
	%save(fullfile(dir_RISTRF,sprintf('RISTRF_v5_%s_%s.mat',tag,monkeyName{m})))
%end


% Plot maps

%for m = 1:2
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

	% All significant clusters
	h_map(m) = figure;
	for ldx_plt = 1:8
		subplot_tight(2,4,ldx_plt,[.1 .06])
		%h_n = heatmap(flipud(htbl(:,:,ldx_plt,m)));
		h_n = heatmap(htbl(:,:,ldx_plt,m));
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
	currPos = get(h_map(m),'Position'); set(h_map(m),'Position',[currPos(1),currPos(2),1740,800]);

	h_map2(m) = figure;
	for ldx_plt = 9:14
		subplot_tight(2,3,ldx_plt-8,[.1 .06])
		%h_n = heatmap(flipud(htbl(:,:,ldx_plt,m)));
		h_n = heatmap(htbl(:,:,ldx_plt,m));
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

	% Significant clusters in supra-layers
	h_map_supra(m) = figure;
	for ldx_plt = 1:8
		subplot_tight(2,4,ldx_plt,[.1 .06])
		%h_n = heatmap(flipud(htbl_supra(:,:,ldx_plt,m)));
		h_n = heatmap(htbl_supra(:,:,ldx_plt,m));
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
	sgtitle(sprintf('%s (N=%i, p_{thr}=%.0d), Supra-layer',monkeyName{m},numel(idx_supra),p_crit),'FontSize',24)
	currPos = get(h_map_supra(m),'Position'); set(h_map_supra(m),'Position',[currPos(1),currPos(2),1740,800]);

	h_map2_supra(m) = figure;
	for ldx_plt = 9:14
		subplot_tight(2,3,ldx_plt-8,[.1 .06])
		%h_n = heatmap(flipud(htbl_supra(:,:,ldx_plt,m)));
		h_n = heatmap(htbl_supra(:,:,ldx_plt,m));
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
	sgtitle(sprintf('%s (N=%i, p_{thr}=%.0d), Supra-layer',monkeyName{m},numel(idx_supra),p_crit),'FontSize',24)
	currPos = get(h_map2_supra(m),'Position'); set(h_map2_supra(m),'Position',[currPos(1),currPos(2),1200,800]);

	% Significant clusters in infra-layers
	h_map_infra(m) = figure;
	for ldx_plt = 1:8
		subplot_tight(2,4,ldx_plt,[.1 .06])
		%h_n = heatmap(flipud(htbl_infra(:,:,ldx_plt,m)));
		h_n = heatmap(htbl_infra(:,:,ldx_plt,m));
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
	sgtitle(sprintf('%s (N=%i, p_{thr}=%.0d), Infra-layer',monkeyName{m},numel(idx_infra),p_crit),'FontSize',24)
	currPos = get(h_map_infra(m),'Position'); set(h_map_infra(m),'Position',[currPos(1),currPos(2),1740,800]);

	h_map2_infra(m) = figure;
	for ldx_plt = 9:14
		subplot_tight(2,3,ldx_plt-8,[.1 .06])
		%h_n = heatmap(flipud(htbl_infra(:,:,ldx_plt,m)));
		h_n = heatmap(htbl_infra(:,:,ldx_plt,m));
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
	sgtitle(sprintf('%s (N=%i, p_{thr}=%.0d), Infra-layer',monkeyName{m},numel(idx_infra),p_crit),'FontSize',24)
	currPos = get(h_map2_infra(m),'Position'); set(h_map2_infra(m),'Position',[currPos(1),currPos(2),1200,800]);
	
	% Difference between supra-layer and infra-layer
	h_map_diff(m) = figure;
	for ldx_plt = 1:8
		ax_n1 = subplot_tight(2,4,ldx_plt,[.1 .06]);
		h_n1 = heatmap(htbl_supra(:,:,ldx_plt,m)-htbl_infra(:,:,ldx_plt,m),'MissingDataLabel','');
		cmin = min(min(htbl_supra(:,:,ldx_plt,m)-htbl_infra(:,:,ldx_plt,m)));
		cmax = max(max(htbl_supra(:,:,ldx_plt,m)-htbl_infra(:,:,ldx_plt,m)));
		caxis([cmin cmax])
		inc = (cmax-cmin)/256;
		cmap = [ linspace(1,0,round(cmax/inc)) zeros(1,abs(round(cmin/inc))); zeros(1,256); zeros(1,round(cmax/inc)) linspace(0,1,abs(round(cmin/inc)))]';
		cmap = flipud(cmap);
		colormap(gca, cmap)
		h_n1.CellLabelFormat = '%.1f';
		set(gca,'XDisplayLabels',coordLabels_x)
		set(gca,'YDisplayLabels',coordLabels_y)
		%colormap jet;
		if strcmp(monkeyName{m},'Cassius')
			xlabel('Medial <----> Lateral')
		else
			xlabel('Lateral <----> Medial')
		end
		ylabel('Posterior <----> Anterior')
		title(titleStr{ldx_plt})
	end
	sgtitle(sprintf('%s, Diff. btw Supra-Infra',monkeyName{m}),'FontSize',24)
	currPos = get(h_map_diff(m),'Position'); set(h_map_diff(m),'Position',[currPos(1),currPos(2),1740,800]);

	h_map2_diff(m) = figure;
	for ldx_plt = 9:14
		ax_n2 = subplot_tight(2,3,ldx_plt-8,[.1 .06]);
		h_n1 = heatmap(htbl_supra(:,:,ldx_plt,m)-htbl_infra(:,:,ldx_plt,m),'MissingDataLabel','');
		cmin = min(min(htbl_supra(:,:,ldx_plt,m)-htbl_infra(:,:,ldx_plt,m)));
		cmax = max(max(htbl_supra(:,:,ldx_plt,m)-htbl_infra(:,:,ldx_plt,m)));
		caxis([cmin cmax])
		inc = (cmax-cmin)/256;
		cmap = [ linspace(1,0,round(cmax/inc)) zeros(1,abs(round(cmin/inc))); zeros(1,256); zeros(1,round(cmax/inc)) linspace(0,1,abs(round(cmin/inc)))]';
		cmap = flipud(cmap);
		colormap(gca, cmap)
		h_n2.CellLabelFormat = '%.1f';
		set(gca,'XDisplayLabels',coordLabels_x)
		set(gca,'YDisplayLabels',coordLabels_y)
		%colormap jet;
		if strcmp(monkeyName{m},'Cassius')
			xlabel('Medial <----> Lateral')
		else
			xlabel('Lateral <----> Medial')
		end
		ylabel('Posterior <----> Anterior')
		title(titleStr{ldx_plt})
	end
	sgtitle(sprintf('%s, Diff. btw Supra-Infa',monkeyName{m}),'FontSize',24)
	currPos = get(h_map2_diff(m),'Position'); set(h_map2_diff(m),'Position',[currPos(1),currPos(2),1200,800]);

	saveas(h_map(m),sprintf('~/STRF/RISTRF_Figures/RISTRF_Map1_v5_%s_%s.png',monkeyName{m},tag))
	saveas(h_map2(m),sprintf('~/STRF/RISTRF_Figures/RISTRF_Map2_v5_%s_%s.png',monkeyName{m},tag))
	saveas(h_map_supra(m),sprintf('~/STRF/RISTRF_Figures/RISTRF_Map1_v5_%s_Supra_%s.png',monkeyName{m},tag))
	saveas(h_map2_supra(m),sprintf('~/STRF/RISTRF_Figures/RISTRF_Map2_v5_%s_Supra_%s.png',monkeyName{m},tag))
	saveas(h_map_infra(m),sprintf('~/STRF/RISTRF_Figures/RISTRF_Map1_v5_%s_Infra_%s.png',monkeyName{m},tag))
	saveas(h_map2_infra(m),sprintf('~/STRF/RISTRF_Figures/RISTRF_Map2_v5_%s_Infra_%s.png',monkeyName{m},tag))
	saveas(h_map_diff(m),sprintf('~/STRF/RISTRF_Figures/RISTRF_Map1_v5_%s_Diff_%s.png',monkeyName{m},tag))
	saveas(h_map2_diff(m),sprintf('~/STRF/RISTRF_Figures/RISTRF_Map2_v5_%s_Diff_%s.png',monkeyName{m},tag))


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
	Edges_BF = [0:1:max(max([SigRISTRF(idx_core).BF]),max([SigRISTRF(idx_belt).BF]))];
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
	currPos = get(h_dist(m),'Position'); set(h_dist(m),'Position',[currPos(1),currPos(2),1740,800]);

	%if strcmp(monkeyName{m},'Cassius')
	%	coord_core = [coord_ac(2:23,1); coord_ac(2:19,2)];	
	%	coord_belt = [coord_ac(2:36,3); coord_ac(2:24,4)];
	%else
	%	coord_core = [coord_ac(2:20,1); coord_ac(2:25,2)];	
	%	coord_belt = [coord_ac(2:31,3); coord_ac(2:16,4)];
	%end

	create_ref_map_v5

	save(fullfile(dir_RISTRF,sprintf('RISTRF_v5_%s_%s.mat',tag,monkeyName{m})))
end

