%monkeyName = {'Miyagi','Cassius'};
monkeyName = {'Cassius','Miyagi'};
p_crit = .01;	% Default(=.01)
%p_crit = 1.e-6;

datatag = input('Enter the data tag (including underscore if needed): ','s');

dir_STRF = fullfile('~','STRF');
dir_RI = fullfile(dir_STRF,'RI');
dir_Param = fullfile(dir_STRF,sprintf('STRFParams%s',datatag));
dir_RISTRF = fullfile(dir_STRF,'RISTRF');
dir_RISTRF_Fig = fullfile(dir_STRF,'RISTRF_Figures');
%dir_Casmiya = fullfile('~','Casmiya');
	
if ~exist(dir_RISTRF,'dir')
	mkdir(dir_RISTRF)
end

%savetag = datestr(now,'yymmddHHMM');
savetag = input('Enter a save tag: ','s');

htbl = nan(15,15,14,2);

UberRISTRF = struct([]);
SigRISTRF = struct([]);
RISTRFMap = struct([]);

gldx_uber = 1;
gldx_sig = 1;
gldx_map = 1;
split_uber = 1;
split_sig = 1;

coord_ac = {};

for m = 1:2
	
	% Load coordinates data
	load(fullfile(dir_STRF,sprintf('coord_ac_%s.mat',monkeyName{m})))
	% Load laminar profile data
	load(fullfile(dir_STRF,sprintf('lamProf_%s.mat',monkeyName{m})))
	
	list_RI = dir(fullfile(dir_RI,sprintf('RI_%s*_8-100.mat',monkeyName{m})));
	%list_Param = dir(fullfile(dir_Param,sprintf('STRFParams_%s_*.mat',monkeyName{m})));
	
	
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
			UberRISTRF(gldx_uber).name = monkeyName{m};
			UberRISTRF(gldx_uber).date = RIStat(ldx_RI).date;
			UberRISTRF(gldx_uber).coord = RIStat(ldx_RI).coord;
			UberRISTRF(gldx_uber).area = RIStat(ldx_RI).area;
			UberRISTRF(gldx_uber).clust = STRFSig(idx_valid(ldx_valid)).clustNum;
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
				SigRISTRF(gldx_sig).name = monkeyName{m};
				SigRISTRF(gldx_sig).date = RIStat(ldx_RI).date;
				SigRISTRF(gldx_sig).coord = RIStat(ldx_RI).coord;
				SigRISTRF(gldx_sig).area = RIStat(ldx_RI).area;
				SigRISTRF(gldx_sig).clust = STRFSig(idx_valid(ldx_valid)).clustNum;
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
	% Unique coordinates for all clusters
	coord_uniq_uber = {};
	for ldx_uber = split_uber:numel(UberRISTRF)
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
	for ldx_sig = split_sig:numel(SigRISTRF)
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
		coord_match_uber{k} = find(cellfun(@(x) isequal(x,coord_uniq_uber{k}), {UberRISTRF.coord}));
		coord_match_sig{k} = find(cellfun(@(x) isequal(x,coord_uniq_sig{k}), {SigRISTRF.coord}));
		%[row,col] = ind2sub(size(coord_ac),find(cellfun(@(x) isequal(x,coord_uniq_uber{k}), coord_ac)));
		RISTRFMap(gldx_map).name = monkeyName{m};
		RISTRFMap(gldx_map).coord = coord_uniq_uber{k};
		[row, col] = ind2sub(size(coord_ac),find(cellfun(@(x) isequal(x, RISTRFMap(gldx_map).coord), coord_ac)));
		if any(strcmp(coord_ac{1,col},{'A1','R'}))
			RISTRFMap(gldx_map).area = "core";
		elseif any(strcmp(coord_ac{1,col},{'BeltL','BeltM'}))
			RISTRFMap(gldx_map).area = "belt";
		else
			RISTRFMap(gldx_map).area = "unknown";
		end
		RISTRFMap(gldx_map).numAllClust = numel(coord_match_uber{k});
		RISTRFMap(gldx_map).numSigClust = numel(coord_match_sig{k});
		RISTRFMap(gldx_map).p_median = nanmedian([UberRISTRF(coord_match_uber{k}).p]);
		RISTRFMap(gldx_map).BF = nanmean([UberRISTRF(coord_match_sig{k}).BF]);
		RISTRFMap(gldx_map).BFHz = nanmean([UberRISTRF(coord_match_sig{k}).BFHz]);
		RISTRFMap(gldx_map).Delay = nanmean([UberRISTRF(coord_match_sig{k}).Delay]);
		RISTRFMap(gldx_map).PD = nanmean([UberRISTRF(coord_match_sig{k}).PD]);
		RISTRFMap(gldx_map).t10 = nanmean([UberRISTRF(coord_match_sig{k}).t10]);
		RISTRFMap(gldx_map).t50 = nanmean([UberRISTRF(coord_match_sig{k}).t50]);
		RISTRFMap(gldx_map).Dur = nanmean([UberRISTRF(coord_match_sig{k}).Dur]);
		RISTRFMap(gldx_map).PLI = nanmean([UberRISTRF(coord_match_sig{k}).PLI]);
		RISTRFMap(gldx_map).DSI = nanmean([UberRISTRF(coord_match_sig{k}).DSI]);
		RISTRFMap(gldx_map).BW = nanmean([UberRISTRF(coord_match_sig{k}).BW]);
		RISTRFMap(gldx_map).BWHz = nanmean([UberRISTRF(coord_match_sig{k}).BWHz]);
		RISTRFMap(gldx_map).bSMF = nanmean([UberRISTRF(coord_match_sig{k}).bSMF]);
		RISTRFMap(gldx_map).FmUpperCutoff = nanmean([UberRISTRF(coord_match_sig{k}).FmUpperCutoff]);
	
		% Construct tables for heatmap
		%htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,1,m) = RISTRFMap(gldx_map).numClust/numel(UberRISTRF)*100;
		htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,1,m) = RISTRFMap(gldx_map).numSigClust/RISTRFMap(k).numAllClust*100;
		htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,2,m) = log10(RISTRFMap(gldx_map).p_median);
	
		try
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,3,m) = RISTRFMap(gldx_map).BF;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,4,m) = RISTRFMap(gldx_map).BFHz;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,5,m) = RISTRFMap(gldx_map).Delay;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,6,m) = RISTRFMap(gldx_map).PD;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,7,m) = RISTRFMap(gldx_map).t10;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,8,m) = RISTRFMap(gldx_map).t50;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,9,m) = RISTRFMap(gldx_map).Dur;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,10,m) = RISTRFMap(gldx_map).PLI;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,11,m) = RISTRFMap(gldx_map).BW;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,12,m) = RISTRFMap(gldx_map).BWHz;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,13,m) = RISTRFMap(gldx_map).bSMF;
			htbl(coord_uniq_uber{k}(2)+8,coord_uniq_uber{k}(1)+8,14,m) = RISTRFMap(gldx_map).FmUpperCutoff;
		end
		gldx_map = gldx_map +1;
	end
	split_uber = gldx_uber-1;
	split_sig = gldx_sig-1;
end


%% Basic statistics
idx_loc_core = find([RISTRFMap.area]=="core");
idx_loc_belt = find([RISTRFMap.area]=="belt");
if m==1
	ratioSigClust_Core_Cassius = nanmean([RISTRFMap(idx_loc_core).numSigClust]./[RISTRFMap(idx_loc_core).numAllClust]);
	ratioSigClust_Belt_Cassius = nanmean([RISTRFMap(idx_loc_belt).numSigClust]./[RISTRFMap(idx_loc_belt).numAllClust]);
	pMedian_Core_Cassius = nanmean([RISTRFMap(idx_loc_core).p_median]);
	pMedian_Belt_Cassius = nanmean([RISTRFMap(idx_loc_belt).p_median]);
	t10_Core_Cassius = nanmean([RISTRFMap(idx_loc_core).t10]);
	t10_Belt_Cassius = nanmean([RISTRFMap(idx_loc_belt).t10]);
else
	ratioSigClust_Core_Miyagi = nanmean([RISTRFMap(idx_loc_core).numSigClust]./[RISTRFMap(idx_loc_core).numAllClust]);
	ratioSigClust_Belt_Miyagi = nanmean([RISTRFMap(idx_loc_belt).numSigClust]./[RISTRFMap(idx_loc_belt).numAllClust]);
	pMedian_Core_Miyagi = nanmean([RISTRFMap(idx_loc_core).p_median]);
	pMedian_Belt_Miyagi = nanmean([RISTRFMap(idx_loc_belt).p_median]);
	t10_Core_Miyagi = nanmean([RISTRFMap(idx_loc_core).t10]);
	t10_Belt_Miyagi = nanmean([RISTRFMap(idx_loc_belt).t10]);
end

% Indices of each area
idx_a1 = find(cellfun(@(x) strcmp(x,'A1'), {SigRISTRF.area}));
idx_r = find(cellfun(@(x) strcmp(x,'R'), {SigRISTRF.area}));
idx_bm = find(cellfun(@(x) strcmp(x,'BeltM'), {SigRISTRF.area}));
idx_bl = find(cellfun(@(x) strcmp(x,'BeltL'), {SigRISTRF.area}));
idx_core = sort([idx_a1 idx_r])';
idx_belt = sort([idx_bm idx_bl])';

%% Eliminate weirdly big BW
%for i = 1:length(SigRISTRF)
%	if SigRISTRF(i).BW > 4.
%		SigRISTRF(i).BW = [];
%	end
%end


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
	%[h_ks_1, p_ks_1, k_ks_1] = kstest_2s_2d([[SigRISTRF(idx_core).BFHz]' [SigRISTRF(idx_core).BWHz]'], [[SigRISTRF(idx_belt).BFHz]' [SigRISTRF(idx_belt).BWHz]']);
	if p_corrcoeff1 < .05
	%if h_ks_1
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
	%[h_ks_2, p_ks_2, k_ks_2] = kstest_2s_2d([[SigRISTRF(idx_core).BFHz]' [SigRISTRF(idx_core).BW]'], [[SigRISTRF(idx_belt).BFHz]' [SigRISTRF(idx_belt).BW]']);
	if p_corrcoeff2 < .05
	%if h_ks_2
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
	%[h_ks_3, p_ks_3, k_ks_3] = kstest_2s_2d([BW_core bSMF_core], [BW_belt bSMF_belt]);
	if p_corrcoeff3 < .05
	%if h_ks_3
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
	%[h_ks_4, p_ks_4, k_ks_4] = kstest_2s_2d([[SigRISTRF(idx_core).BFHz]' [SigRISTRF(idx_core).Delay]'], [[SigRISTRF(idx_belt).BFHz]' [SigRISTRF(idx_belt).Delay]']);
	if p_corrcoeff4 < .05
	%if h_ks_4
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
	%[h_ks_5, p_ks_5, k_ks_5] = kstest_2s_2d([[SigRISTRF(idx_core).Dur]' [SigRISTRF(idx_core).FmUpperCutoff]'], [[SigRISTRF(idx_belt).Dur]' [SigRISTRF(idx_belt).FmUpperCutoff]']);
	if p_corrcoeff5 < .05
	%if h_ks_5
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
	%[h_ks_6, p_ks_1, k_ks_1] = kstest_2s_2d([[SigRISTRF(idx_core).Dur]' [SigRISTRF(idx_core).Delay]'], [[SigRISTRF(idx_belt).Dur]' [SigRISTRF(idx_belt).Delay]']);
	if p_corrcoeff6 < .05
	%if h_ks_6
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
saveas(h_corr{m},fullfile(dir_RISTRF_Fig,sprintf('RISTRF_Corr_Comb_%s.png',savetag)))
save(fullfile(dir_RISTRF,sprintf('RISTRF_Comb_%s.mat',savetag)))

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
	save(fullfile(dir_RISTRF,sprintf('RISTRF_%s_%s.mat',savetag,monkeyName{m})))
end

