function [ClustStat] = get_cluster_stat(monkeyName,sessionDate,driveID)

	load(sprintf('~/STRF/STRFParams/STRFParams_%s_%s_%s.mat',monkeyName,sessionDate,driveID))
	ClustParam = nan(length(UberSTRF),4);
	
	for i = 1:length(UberSTRF)
		try
			ClustParam(i,1) = UberSTRF(i).ClustNum;
			ClustParam(i,2) = UberSTRF(i).RFParam.PeakDelay;
			ClustParam(i,3) = UberSTRF(i).RFParam.PLI;
			ClustParam(i,4) = UberSTRF(i).RFParam.PLI2;
		end
	end
	
	corrPLI = corr(ClustParam(find(~isnan(ClustParam(:,3))),3),ClustParam(find(~isnan(ClustParam(:,2))),2));
	%corrPLI2 = corr(ClustParam(find(~isnan(ClustParam(:,4))),4),ClustParam(find(~isnan(ClustParam(:,2))),2));
	ClustStat.Area = recArea;
	ClustStat.Params = ClustParam;
	ClustStat.Corr_PLI_PeakDelay = corrPLI;
	%ClustStat.Corr_PLI2_PeakDelay = corrPLI2;

end	% End of funciton definition
