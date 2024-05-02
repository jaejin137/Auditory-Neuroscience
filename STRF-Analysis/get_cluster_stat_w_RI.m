function [ClustStat] = get_cluster_stat_w_RI(monkeyName,sessionDate,driveID,p_crit)

	load(sprintf('~/STRF/STRFParams/STRFParams_%s_%s_%s.mat',monkeyName,sessionDate,driveID))
	goodClust = find_good_clusters(monkeyName,sessionDate,driveID,p_crit);
	%ClustParam = nan(length(goodClust),4);
	idx_goodClust = [];
	for i = 1:size(goodClust,1)
		for j = 1:length(UberSTRF)
			if UberSTRF(j).ClustNum == goodClust{i,1}
				idx_goodClust(end+1,1) = j;
			end
		end
	end
	
	j = 1;
	for i = 1:length(idx_goodClust)
		try
			ClustParam(j,2) = UberSTRF(idx_goodClust(i)).RFParam.PeakDelay;
			ClustParam(j,3) = UberSTRF(idx_goodClust(i)).RFParam.PLI;
			ClustParam(j,4) = UberSTRF(idx_goodClust(i)).RFParam.PLI2;
			ClustParam(j,5) = UberSTRF(idx_goodClust(i)).RFParam.Duration;
			ClustParam(j,1) = UberSTRF(idx_goodClust(i)).ClustNum;	% assign lastly to prevent ClustParam from being in int32 format.
			j = j+1;
		end
	end
	
	corr_PLI_PeakDelay = corr(ClustParam(find(~isnan(ClustParam(:,3))),3),ClustParam(find(~isnan(ClustParam(:,2))),2));
	%corr_PLI2_PeakDelay = corr(ClustParam(find(~isnan(ClustParam(:,4))),4),ClustParam(find(~isnan(ClustParam(:,2))),2));
	corr_Duration_PeakDelay = corr(ClustParam(find(~isnan(ClustParam(:,5))),5),ClustParam(find(~isnan(ClustParam(:,2))),2));
	ClustStat.Area = recArea;
	ClustStat.Params = ClustParam;
	ClustStat.Corr_PLI_PeakDelay = corr_PLI_PeakDelay;
	%ClustStat.Corr_PLI2_PeakDelay = corr_PLI2_PeakDelay;
	ClustStat.Corr_Duration_PeakDelay = corr_Duration_PeakDelay;

end	% End of funciton definition
