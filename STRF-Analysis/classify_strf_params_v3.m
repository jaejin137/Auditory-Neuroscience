monkeyName = input('Enter monkey name: ','s');
%monkeyName = 'Cassius';

clustType = input('Which unit type to process?[(s)ingle, (m)ua, (b)oth, or (c)ustom]: ','s');
%clustType = 's';

%ifPlot = input('Plot the result? (y/n) ','s');
ifPlot = 'y';

dir_result = fullfile('~','STRF','Result');
dir_kilosorted = fullfile('~','kiloSorted_DMR');

list_result = dir(fullfile(dir_result,sprintf('STRFParams_%s*.mat',monkeyName)));

%load('targetInfo_lamProf.mat')
load('areaColor.mat')
sessionColor = [[0,0,0];[1,1,0];[1,0,1];[0,1,1];[1,0,0];[0,1,0];[0,0,1];[.5,0,0];[0,.5,0];[0,0,.5]];


%if clustType == 's'
%	clustTypeName = 'Single Unit';
%elseif clustType == 'm'
%	clustTypeName = 'MUA';
%elseif clustType == 'b'
%	clustTypeName = 'Single Unit & MUA';
%else
%	display('!!!Error: Wrong cluster type given!!!')
%	return
%end

Area = {};
idx_core = [];
idx_belt = [];
for idx_data = 1:length(list_result)
	load(fullfile(dir_result,list_result(idx_data).name),'area');
	Area{idx_data,1} = area;
	if area == 'core'
		idx_core(end+1) = idx_data;
	else
		idx_belt(end+1) = idx_data;
	end
end

chanArr = {};


%parfor idx_data = 1:length(list_result)
for idx_data = 1:length(list_result)
	%try

		result = load(fullfile(dir_result,list_result(idx_data).name));
		area = result.area;
		UberSTRF = result.UberSTRF;

		basename = split(list_result(idx_data).name,'.');
		basename = split(list_result(idx_data).name,'.');
		dataTok = split(basename{1},'_');
		sessionDate = dataTok{3};
		driveID = sprintf('%s_%s_%s',dataTok{4:6});
		
		fprintf('\n\nData for %s loaded. Area: %s\n',sessionDate,area)

		path_clustInfo = fullfile(dir_kilosorted,sprintf('Mr%s-%s',monkeyName,sessionDate),driveID,'KS2_7_AC','ClusterInfo');

		cluster_info = importTSV_Y2018(fullfile(path_clustInfo,'cluster_info_new.tsv'));

		clust2proc = [];
		for idx_clust = 1:size(UberSTRF,2)
			clust2proc(end+1,1) = UberSTRF(idx_clust).ClustNum;
		end

		if clustType == 's'
			clust2proc = cluster_info.id(find(ismember(cluster_info.group,{'good'})));
			clustTypeName = 'Single Unit';
		elseif clustType == 'm'
			clust2proc = cluster_info.id(find(ismember(cluster_info.group,{'mua'})));
			clustTypeName = 'MUA';
		elseif clustType == 'b'
			clust2proc = cluster_info.id(find(ismember(cluster_info.group,{'good','mua'})));
			clustTypeName = 'Single Unit & MUA';
		elseif clustType == 'c'
			load(sprintf('~/STRF/clusters_goodSTRF_%s.mat',monkeyName));
			clust2proc = clusters_goodSTRF{find(contains(clusters_goodSTRF(:,1),sessionDate)),2};
			clustTypeName = 'Custom';
		else
			display('!!!Error: Wrong cluster type given!!!')
		end

		chanArr_par = [];
		Delay_par = [];
		Duration_par = [];
		BF_par = [];
		PLI_par = [];
		PLI2_par = [];
		DSI_par = [];
		PeakDelay_par = [];
		PeakEnvDelay_par = [];
		PeakBF_par = [];

		for idx_clust = 1:length(clust2proc)
			chanArr_par(1,idx_clust) = cluster_info.ch(find(cluster_info.id==clust2proc(idx_clust)));
			Delay_par(1,idx_clust) = UberSTRF(idx_clust).RFParam.Delay;
			Duration_par(1,idx_clust) = UberSTRF(idx_clust).RFParam.Duration;
			BF_par(1,idx_clust) = UberSTRF(idx_clust).RFParam.BFHz;
			PLI_par(1,idx_clust) = UberSTRF(idx_clust).RFParam.PLI;
			PLI2_par(1,idx_clust) = UberSTRF(idx_clust).RFParam.PLI2;
			DSI_par(1,idx_clust) = UberSTRF(idx_clust).RFParam.DSI;
			PeakDelay_par(1,idx_clust) = UberSTRF(idx_clust).RFParam.PeakDelay;
			PeakEnvDelay_par(1,idx_clust) = UberSTRF(idx_clust).RFParam.PeakEnvDelay;
			PeakBF_par(1,idx_clust) = UberSTRF(idx_clust).RFParam.PeakBF;
		end
		chanArr{idx_data,1} = chanArr_par;
		Delay{idx_data,1} = Delay_par;
		Duration{idx_data,1} = Duration_par;
		BF{idx_data,1} = BF_par;
		PLI{idx_data,1} = PLI_par;
		PLI2{idx_data,1} = PLI2_par;
		DSI{idx_data,1} = DSI_par;
		PeakDelay{idx_data,1} = PeakDelay_par;
		PeakEnvDelay{idx_data,1} = PeakEnvDelay_par;
		PeakBF{idx_data,1} = PeakBF_par;
	%catch
	%	fprintf('\n\n!!! Error occurred with data %s. Skipping...!!!\n',sessionDate)
	%end
end

if ifPlot == 'y'
	plot_strf_params_v3
end

