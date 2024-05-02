function [] = plot_clust_stat(ClustStat)
	ClustParam = ClustStat.Params;
	corrPLI = ClustStat.corr_PLI_PeakDelay;

	%h_PLI = figure;
	plot(ClustParam(:,3),ClustParam(:,2),'.','MarkerSize',40)
	for i = 1:length(ClustParam)
		text(ClustParam(i,3),ClustParam(i,2),sprintf('%i',ClustParam(i,1)),'FontSize',16,'Color','red')
	end
	hold on
	[p_PLI,S_PLI] = polyfit(ClustParam(find(~isnan(ClustParam(:,3))),3),ClustParam(find(~isnan(ClustParam(:,2))),2),1);
	[y_PLI,delta_PLI] = polyval(p_PLI,ClustParam(find(~isnan(ClustParam(:,3))),3),S_PLI);
	[PLISort,idxSort] = sort(ClustParam(find(~isnan(ClustParam(:,3))),3));
	plot(ClustParam(find(~isnan(ClustParam(:,3))),3),y_PLI,'r-','LineWidth',2)
	xlabel('PLI'); ylabel('PeakDelay');
	title(sprintf('%s-%s (Corr=%.3f)',monkeyName,sessionDate,corrPLI))
	%drawnow
	
	%h_PLI2 = figure;
	%plot(ClustParam(:,4),ClustParam(:,2),'.','MarkerSize',40)
	%for i = 1:length(ClustParam)
	%	text(ClustParam(i,4),ClustParam(i,2),sprintf('%i',ClustParam(i,1)),'FontSize',16,'Color','red')
	%end
	%hold on
	%[p_PLI2,S_PLI2] = polyfit(ClustParam(find(~isnan(ClustParam(:,4))),4),ClustParam(find(~isnan(ClustParam(:,2))),2),1);
	%[y_PLI2,delta_PLI2] = polyval(p_PLI2,ClustParam(find(~isnan(ClustParam(:,4))),4),S_PLI2);
	%[PLI2Sort,idxSort2] = sort(ClustParam(find(~isnan(ClustParam(:,4))),4));
	%plot(ClustParam(find(~isnan(ClustParam(:,4))),4),y_PLI2,'r-','LineWidth',2)
	%xlabel('PLI2'); ylabel('PeakDelay');
	%title(sprintf('%s-%s (Corr=%.3f)',monkeyName,sessionDate,corrPLI2))
	%drawnow


	%close(h_PLI)
	%close(h_PLI2)

end	% End of function definition
