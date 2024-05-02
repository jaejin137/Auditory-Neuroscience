load('~/STRF/STRFParams_UberStat_w_RI_p001.mat')

UberStat_core = [];
UberStat_belt = [];
for idx_result = 1:size(UberStat,2)
	for idx_clust = 1:size(UberStat(idx_result).Params,1)
		if strcmp(UberStat(idx_result).Area,'core')
			UberStat_core(end+1:end+size(UberStat(idx_result).Params,1),1:2) = UberStat(idx_result).Params(:,2:3);
		elseif strcmp(UberStat(idx_result).Area,'belt')
			UberStat_belt(end+1:end+size(UberStat(idx_result).Params,1),1:2) = UberStat(idx_result).Params(:,2:3);
		end
	end
end

[h_PD, p_PD, ci_PD, stat_PD] = ttest2(UberStat_core(:,1),UberStat_belt(:,1));
[h_PLI, p_PLI, ci_PLI, stat_PLI] = ttest2(UberStat_core(:,2),UberStat_belt(:,2));
if p_PD < .05
	sig_PD = '*';
else
	sig_PD = '';
end
if p_PLI < .05
	sig_PLI = '*';
else
	sig_PLI = '';
end


Edges_PD = [0:10:max(max(UberStat_core(:,1)),max(UberStat_belt(:,1)))];
Edges_PLI = [0:0.02:max(max(UberStat_core(:,2)),max(UberStat_belt(:,2)))];
figure
subplot(1,2,1)
histogram(UberStat_core(:,1),Edges_PD,'FaceColor','r')
hold on
histogram(UberStat_belt(:,1),Edges_PD,'FaceColor','g')
xlabel('PeakDelay [ms]','FontSize',20)
ylabel('No. of Clusters','FontSize',20)
legend({'Core','Belt'})
%title('PeakDelay distribution across all clusters','FontSize',24)
title(sprintf('p%s = %.2d (two-sample t-test)',sig_PD,p_PD),'FontSize',20)
subplot(1,2,2)
histogram(UberStat_core(:,2),Edges_PLI,'FaceColor','r')
hold on
histogram(UberStat_belt(:,2),Edges_PLI,'FaceColor','g')
xlabel('Phase Locking Index','FontSize',20)
ylabel('No. of Clusters','FontSize',20)
legend({'Core','Belt'})
%title('PLI distribution across all clusters','FontSize',24)
title(sprintf('p%s = %.2d (two-sample t-test)',sig_PLI,p_PLI),'FontSize',20)
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1000,400]);
