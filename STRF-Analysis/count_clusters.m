function [clustCount] = count_clusters(monkeyName,p_crit)

	%monkeyName = 'Cassius';
	
	Const.NBoot = 8;
	Const.NB = 100;
	%p_crit = 1.e-32;
	
	
	load(sprintf('lamProf_%s.mat',monkeyName))
	
	for idx_data = 1:length(lamProf)
		nameRoot = split(lamProf(idx_data).dataName,'.');
		nameTok = split(nameRoot{1},'_');
		sessionDate = nameTok{2};
		driveIDRoot = sprintf('%s_%s',nameTok{3:4});
		dir_target = dir(sprintf('~/kiloSorted_DMR/Mr%s-%s/D*',monkeyName,sessionDate));
		numDrv = 0;
		for i = 1:length(dir_target)
			if contains(dir_target(i).name,driveIDRoot)
				driveID = dir_target(i).name;
				numDrv = numDrv+1;
			end
		end
		if numDrv > 1
			numDupDrv = numDupDrv +1;
		end
		fprintf('\nChecking %s-%s_%s\n',monkeyName,sessionDate,driveID)
		clustCount(idx_data).dataName = sprintf('%s-%s_%s.mat',monkeyName,sessionDate,driveID);
		try
			load(sprintf('~/kiloSorted_DMR/Mr%s-%s/%s/KS2_7_AC/ClusterInfo/spike_times_all_clust.mat',monkeyName,sessionDate,driveID))
			clustCount(idx_data).single = sum(arrayfun(@(x) contains(x,'good'), spikeTimesAllClust(:,3)));
			clustCount(idx_data).MUA = sum(arrayfun(@(x) contains(x,'mua'), spikeTimesAllClust(:,3)));
		catch
			clustCount(idx_data).status = sprintf('Cluster info is not available');
			fprintf('\nCluster info is not available for %s-%s_%s\n',monkeyName,sessionDate,driveID)
			continue
		end
	
		try
			load(sprintf('~/Research/Auditory/STRF/RI/RI_%s-%s_%s_%i-%i.mat',monkeyName,sessionDate,driveID,Const.NBoot,Const.NB))
			clustCount(idx_data).singleSRI = sum([STRFSig.p] < p_crit & arrayfun(@(x) strcmp(x.clustType,'good'),STRFSig));
			clustCount(idx_data).MUASRI = sum([STRFSig.p] < p_crit & arrayfun(@(x) strcmp(x.clustType,'mua'),STRFSig));
			clustCount(idx_data).status = 'Success';
		catch
			clustCount(idx_data).status = sprintf('RI is not available');
			fprintf('\nRI is not available for %s-%s_%s\n',monkeyName,sessionDate,driveID)
		end
	end
	totSingle = sum([clustCount.single])
	totMUA = sum([clustCount.MUA])
	totSingleSRI = sum([clustCount.singleSRI])
	totMUASRI = sum([clustCount.MUASRI])
end % End of function definition
