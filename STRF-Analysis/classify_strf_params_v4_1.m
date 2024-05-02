%function [STRFClass] = classify_strf_params_v4(monkeyName,clustType,p_crit,Tag)

	%if nargin < 4
		%Tag = sprintf('_%s',datestr(now,'mmddHHMM'));
		p_crit = .01;
	%end
	%if nargin < 4
	%	p_crit = .01;
	%end

	monkeyName = input('Enter monkey name: ','s');
	%monkeyName = 'Cassius';
	
	Tag = input('Enter a tag: ','s');

	%clustType = input('Which unit type to process?[(s)ingle, (m)ua, (b)oth, or (c)ustom]: ','s');
	clustType = 'b';
	
	%ifPlot = input('Plot the result? (y/n) ','s');
	ifPlot = 'y';

	dir_class = fullfile('~','STRF','STRFClass');

	dir_result = fullfile('~','STRF','STRFParams');
	dir_kilosorted = fullfile('~','kiloSorted_DMR');
	
	list_result = dir(fullfile(dir_result,sprintf('STRFParams_%s*8-100*.mat',monkeyName)));
	
	%load('targetInfo_lamProf.mat')
	load('areaColor.mat')
	sessionColor = [[0,0,0];[1,1,0];[1,0,1];[0,1,1];[1,0,0];[0,1,0];[0,0,1];[.5,0,0];[0,.5,0];[0,0,.5]];
	
	
	%if clustType == 's'
	%	clustTypeName = 'SingleUnit';
	%elseif clustType == 'm'
	%	clustTypeName = 'MUA';
	%elseif clustType == 'b'
	%	clustTypeName = 'SingleUnit & MUA';
	%else
	%	display('!!!Error: Wrong cluster type given!!!')
	%	return
	%end
	
	Area = {};
	idx_core = [];
	idx_belt = [];
	for idx_data = 1:length(list_result)
		%load(fullfile(dir_result,list_result(idx_data).name),'recArea');
		load(fullfile(dir_result,list_result(idx_data).name),'area');
		recArea = area;
		clearvars area;
		Area{idx_data,1} = recArea;
		if recArea == 'core'
			idx_core(end+1) = idx_data;
		else
			idx_belt(end+1) = idx_data;
		end
	end
	
	chanArr = {};
	
	
	for idx_data = 1:length(list_result)
		%try
			result = load(fullfile(dir_result,list_result(idx_data).name));
			%recArea = result.recArea;
			recArea = result.area;
			UberSTRF = result.UberSTRF;
	
			basename = split(list_result(idx_data).name,'.');
			basename = split(list_result(idx_data).name,'.');
			dataTok = split(basename{1},'_');
			sessionDate = dataTok{3};
			driveID = sprintf('%s_%s_%s',dataTok{4:6});
	
			fprintf('\n\nData for %s loaded. Area: %s\n',sessionDate,recArea)
	
			path_clustInfo = fullfile(dir_kilosorted,sprintf('Mr%s-%s',monkeyName,sessionDate),driveID,'KS2_7_AC','ClusterInfo');
	
			cluster_info = importTSV(fullfile(path_clustInfo,'cluster_info_new.tsv'));
	
			if clustType == 's'
				targClust = cluster_info.id(find(ismember(cluster_info.group,{'good'})));
				clustTypeName = 'SingleUnit';
			elseif clustType == 'm'
				targClust = cluster_info.id(find(ismember(cluster_info.group,{'mua'})));
				clustTypeName = 'MUA';
			elseif clustType == 'b'
				targClust = cluster_info.id(find(ismember(cluster_info.group,{'good','mua'})));
				clustTypeName = 'Both';
			elseif clustType == 'c'
				load(sprintf('~/STRF/clusters_goodSTRF_%s.mat',monkeyName));
				targClust = clusters_goodSTRF{find(contains(clusters_goodSTRF(:,1),sessionDate)),2};
				clustTypeName = 'Custom';
			else
				display('!!!Error: Wrong cluster type given!!!')
			end
	
			clust2proc = targClust(ismember(targClust,[UberSTRF.ClustNum]));
	
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
			T1_10_par = [];
			T2_10_par = [];
			T1_50_par = [];
			T2_50_par = [];
	
			parfor idx_clust = 1:length(clust2proc)
				for idx_UberSTRF = 1:length(UberSTRF)
					if UberSTRF(idx_UberSTRF).ClustNum == clust2proc(idx_clust)
						chanArr_par(1,idx_clust) = cluster_info.ch(find(cluster_info.id==clust2proc(idx_clust)));
						Delay_par(1,idx_clust) = UberSTRF(idx_UberSTRF).RFParam.Delay;
						Duration_par(1,idx_clust) = UberSTRF(idx_UberSTRF).RFParam.Duration;
						BF_par(1,idx_clust) = UberSTRF(idx_UberSTRF).RFParam.BFHz;
						PLI_par(1,idx_clust) = UberSTRF(idx_UberSTRF).RFParam.PLI;
						PLI2_par(1,idx_clust) = UberSTRF(idx_UberSTRF).RFParam.PLI2;
						DSI_par(1,idx_clust) = UberSTRF(idx_UberSTRF).RFParam.DSI;
						PeakDelay_par(1,idx_clust) = UberSTRF(idx_UberSTRF).RFParam.PeakDelay;
						PeakEnvDelay_par(1,idx_clust) = UberSTRF(idx_UberSTRF).RFParam.PeakEnvDelay;
						PeakBF_par(1,idx_clust) = UberSTRF(idx_UberSTRF).RFParam.PeakBF;
						T1_10_par(1,idx_clust) = UberSTRF(idx_UberSTRF).RFParam.t1_10;
						T2_10_par(1,idx_clust) = UberSTRF(idx_UberSTRF).RFParam.t2_10;
						T1_50_par(1,idx_clust) = UberSTRF(idx_UberSTRF).RFParam.t1_50;
						T2_50_par(1,idx_clust) = UberSTRF(idx_UberSTRF).RFParam.t2_50;
					end
				end
			end
			date{idx_data,1} = sessionDate;
			area{idx_data,1} = recArea;
			clusters{idx_data,1} = clust2proc;
			chanArr{idx_data,1} = chanArr_par';
			Delay{idx_data,1} = Delay_par';
			Duration{idx_data,1} = Duration_par';
			BF{idx_data,1} = BF_par';
			PLI{idx_data,1} = PLI_par';
			PLI2{idx_data,1} = PLI2_par';
			DSI{idx_data,1} = DSI_par';
			PeakDelay{idx_data,1} = PeakDelay_par';
			PeakEnvDelay{idx_data,1} = PeakEnvDelay_par';
			PeakBF{idx_data,1} = PeakBF_par';
			T1_10{idx_data,1} = T1_10_par';
			T2_10{idx_data,1} = T2_10_par';
			T1_50{idx_data,1} = T1_50_par';
			T2_50{idx_data,1} = T2_50_par';
		%catch
		%	fprintf('\n\n!!! Error occurred with data %s. Skipping...!!!\n',sessionDate)
		%end
	end

	totNumClust = 0;
	for k = 1:length(chanArr)
		totNumClust = totNumClust + numel(chanArr{k,1});
	end

	STRFClass.date = date;
	STRFClass.area = area;
	STRFClass.clusters = clusters;
	STRFClass.chanArr = chanArr;
	STRFClass.Delay = Delay;
	STRFClass.Duration = Duration;
	STRFClass.BF = BF;
	STRFClass.PLI = PLI;
	STRFClass.PLI2 = PLI2;
	STRFClass.DSI = DSI;
	STRFClass.PeakDelay = PeakDelay;
	STRFClass.PeakEnvDelay = PeakEnvDelay;
	STRFClass.PeakBF = PeakBF;
	STRFClass.T1_10 = T1_10;
	STRFClass.T2_10 = T2_10;
	STRFClass.T1_50 = T1_50;
	STRFClass.T2_50 = T2_50;
	STRFClass.totNumClust = totNumClust;
	fprintf('\n\nTotal number of clusters classified = %i.\n\n',totNumClust)

	save(sprintf('%s/STRFClass_%s_%s_%s.mat',dir_class,monkeyName,clustTypeName,Tag),'STRFClass')
	
	if ifPlot == 'y'
		plot_strf_params_v4_1
	end
%end % End of function definition	
