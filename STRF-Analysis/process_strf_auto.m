monkeyName = input('Enter monkey name [Cassius or Miyagi]: ','s');
clustType = input('Which unit type to process?[(s)ingle, (m)ua, (b)oth, or (c)ustom]: ','s');

warning('off')

dir_result = fullfile('~','STRF','Result');

load(sprintf('lamProf_%s.mat',monkeyName));
for idx_data = 2:length(lamProf)
%for idx_data = 2:5
	try
		fprintf('\n\nProcessing %s ...\n\n',lamProf(idx_data).dataName)
		% Parse data information
		dataTok = split(lamProf(idx_data).dataName,'_');
		monkeyID = dataTok{1};
		%if strcmp(monkeyID,'MrC')
		%	monkeyName = 'Cassius';
		%elseif strcmp(monkeyID,'MrM')
		%	monkeyName = 'Miyagi';
		%else
		%	display('!!! Error: Wrong monkey ID !!!')
		%	continue
		%end
		sessionDate = dataTok{2};
		driveID = sprintf('%s_AC_R1',dataTok{3});
		load(sprintf('coord_ac_%s',monkeyName));
		% Identify area
		for idx_area = 1:size(coord_ac,2)
			for idx_coord = 2:length(coord_ac(:,idx_area))
				if ~isempty(coord_ac{idx_coord,idx_area})
					if coord_ac{idx_coord,idx_area} == lamProf(idx_data).coord
						if contains(coord_ac{1,idx_area},{'A1','R'})
							area = 'core'
						elseif contains(coord_ac{1,idx_area},{'BeltM','BeltL'})
							area = 'belt'
						else
							display('!!! Error: Cannot identify area. !!!')
						end
					end
				end
			end
		end
		% Find clusters to calculate STRF
		clusterInfo = importTSV_Y2018(fullfile('~','kiloSorted_DMR',sprintf('Mr%s-%s',monkeyName,sessionDate),driveID,'KS2_7_AC','ClusterInfo','cluster_info_new.tsv'));

		clust2proc = [];

		if clustType == 's'
			clust2proc = cluster_info.id(find(ismember(cluster_info.group,{'good'})));
		elseif clustType == 'm'
			clust2proc = cluster_info.id(find(ismember(cluster_info.group,{'mua'})));
		elseif clustType == 'b'
			clust2proc = cluster_info.id(find(ismember(cluster_info.group,{'good','mua'})));
		elseif clustType == 'c'
			load(sprintf('~/STRF/clusters_goodSTRF_%s.mat',monkeyName));
			clust2proc = clusters_goodSTRF{find(contains(clusters_goodSTRF(:,1),sessionDate)),2};
		else
			display('!!!Error: Wrong cluster type given!!!')
			return
		end

		% Make triggers available to workers
		trigDir = dir(fullfile('~/STRF','Triggers',sprintf('%s/*%s*.mat',monkeyName,sessionDate)));
		TrigFile = fullfile(trigDir.folder,trigDir.name);
		load(TrigFile)
		
		calculate_strf_params

		%var2save = {'Fm','Fsd','Fss','Max','MaxFm','MaxRD','MdB','ModType','N','NBlocks','No1','No1A1','No1A2','No1B1','No1B2','PP','RD','RFParam','RTF','RTFVar','SModType','SPL','SPLN','STRF1','STRF1A_Ch1','STRF1A_Ch2','STRF1B_Ch1','STRF1B_Ch2','STRF1s','Sound','T1','T2','Tresh','TrigA','TrigB','TrigFile','UF','UberSTRF','Wo1','Wo1A1','Wo1A2','Wo1B1','Wo1B2','clust1','driveID','faxis','clust2proc','goodSTRFUnits','h_gabor','h_mod','h_strf','m','monkeyName','numClust','numCol','numRow','p','s_bin','sessionDate','spet1A','spet1B','spikeInd1','spikeTime1','spikeTimeFile','sprfile','sprtype','st_clu','taxis','theta','totNumClust','trigDir','area'};
		var2save = {'UberSTRF','area'};

		save(fullfile(dir_result,sprintf('STRFParams_%s_%s_%s.mat',monkeyName,sessionDate,driveID)),var2save{:});
	catch
		fprintf('\n\n!!! Error occurred while processing %s. Moving on!!!\n\n',lamProf(idx_data).dataName)
		continue
	end
	display('Clearing workspace...')
	close all
	clearvars -except monkeyName clustType dir_result lamProf idx_data
end


