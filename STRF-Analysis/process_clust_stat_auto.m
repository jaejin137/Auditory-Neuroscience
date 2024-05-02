function process_clust_stat_auto(monkeyName,p_crit)
	% Use RI with specific p_crit value.
	%p_crit = 1.e-6;	% default = 0.01
	
	dir_result = fullfile('~','STRF','STRFParams');
	list_result = dir(sprintf('%s/*%s*.mat',dir_result,monkeyName));
	
	
	
	%numCol = 10;
	%numRow = ceil(length(list_result)/numCol);
	
	%h_stat = figure;
	%currPos = get(h_stat,'Position'); set(h_stat,'Position',[currPos(1),currPos(2),numCol*200,numRow*200]);
	
	UberStat_core = [];
	UberStat_belt = [];
	
	for idx_result = 1:length(list_result)
		nameRoot = split(list_result(idx_result).name,'.');
		nameTok = split(nameRoot{1},'_');
		%monkeyName = nameTok{2};
		sessionDate = nameTok{3};
		driveID = sprintf('%s_%s_%s',nameTok{4:6});
	
		fprintf('\nProcessing %s_%s_%s\n\n',monkeyName,sessionDate,driveID)
	
		try
			%UberStat(idx_result) = get_cluster_stat(monkeyName,sessionDate,driveID);
			UberStat(idx_result) = get_cluster_stat_w_RI(monkeyName,sessionDate,driveID,p_crit);
			if UberStat(idx_result).Area == 'core'
				numClustSess_core = size(UberStat(idx_result).Params,1);
				%UberStat_core(end+1:end+numClustSess_core,1:2) = UberStat(idx_result).Params(:,2:3);
				UberStat_core(end+1:end+numClustSess_core,1:3) = UberStat(idx_result).Params(:,[2 3 5]);
			elseif UberStat(idx_result).Area == 'belt'
				numClustSess_belt = size(UberStat(idx_result).Params,1);
				%UberStat_belt(end+1:end+numClustSess_belt,1:2) = UberStat(idx_result).Params(:,2:3);
				UberStat_belt(end+1:end+numClustSess_belt,1:3) = UberStat(idx_result).Params(:,[2 3 5]);
			end
	
			%ClustParam = UberStat(idx_result).Params;
			%corrPLI = UberStat(idx_result).Corr_PLI_PeakDelay;
			%
			%subplot(numRow,numCol,idx_result)
			%plot(ClustParam(:,3),ClustParam(:,2),'.','MarkerSize',20)
			%for i = 1:length(ClustParam)
			%	text(ClustParam(i,3),ClustParam(i,2),sprintf('%i',ClustParam(i,1)),'FontSize',10,'Color','black')
			%end
			%hold on
			%[p_PLI,S_PLI] = polyfit(ClustParam(find(~isnan(ClustParam(:,3))),3),ClustParam(find(~isnan(ClustParam(:,2))),2),1);
			%[y_PLI,delta_PLI] = polyval(p_PLI,ClustParam(find(~isnan(ClustParam(:,3))),3),S_PLI);
			%%[PLISort,idxSort] = sort(ClustParam(find(~isnan(ClustParam(:,3))),3));
			%% Red fit line if core area, Green line if belt area.
			%if UberStat(idx_result).Area == 'core'
			%	plot(ClustParam(find(~isnan(ClustParam(:,3))),3),y_PLI,'r-','LineWidth',4)
			%elseif UberStat(idx_result).Area == 'belt'
			%	plot(ClustParam(find(~isnan(ClustParam(:,3))),3),y_PLI,'g-','LineWidth',4)
			%else
			%	plot(ClustParam(find(~isnan(ClustParam(:,3))),3),y_PLI,'c-','LineWidth',4)
			%end
			%xlabel('PLI'); ylabel('PeakDelay');
			%title(sprintf('%s-%s\n(%s area, corr=%.3f)',monkeyName,sessionDate,UberStat(idx_result).Area,corrPLI))
			%drawnow
		catch
			fprintf('\n!!!Error occurred with %s_%s. Skipping ...\n\n',monkeyName,sessionDate)
		end
	end
	%save(sprintf('~/STRF/STRFParams_UberStat.mat'),'UberStat')
	save(sprintf('~/STRF/STRFParams_UberStat_w_RI_%s.mat',monkeyName),'monkeyName','UberStat*','p_crit')
	%saveas(h_stat,sprintf('~/STRF/STRFParams_PLI_vs_PeakDelay.png'))
	
	
	% Correlation between PLI and PeakDelay
	[p_pli_core,S_pli_core] = polyfit(UberStat_core(:,2),UberStat_core(:,1),1);
	[y_pli_core,delta_pli_core] = polyval(p_pli_core,UberStat_core(:,2),S_pli_core);
	[p_pli_belt,S_pli_belt] = polyfit(UberStat_belt(:,2),UberStat_belt(:,1),1);
	[y_pli_belt,delta_pli_belt] = polyval(p_pli_belt,UberStat_belt(:,2),S_pli_belt);
	% Linear correlation btw PeakDelay and PLI
	[rho_pli_core,pval_pli_core] = corr(UberStat_core(:,1),UberStat_core(:,2));
	[rho_pli_belt,pval_pli_belt] = corr(UberStat_belt(:,1),UberStat_belt(:,2));
	
	% Correlation between Duration and PeakDelay
	[p_dur_core,S_dur_core] = polyfit(UberStat_core(:,3),UberStat_core(:,1),1);
	[y_dur_core,delta_dur_core] = polyval(p_dur_core,UberStat_core(:,3),S_dur_core);
	[p_dur_belt,S_dur_belt] = polyfit(UberStat_belt(:,3),UberStat_belt(:,1),1);
	[y_dur_belt,delta_dur_belt] = polyval(p_dur_belt,UberStat_belt(:,3),S_dur_belt);
	% Linear correlation btw PeakDelay and PLI
	[rho_dur_core,pval_dur_core] = corr(UberStat_core(:,1),UberStat_core(:,3));
	[rho_dur_belt,pval_dur_belt] = corr(UberStat_belt(:,1),UberStat_belt(:,3));
	
	if pval_pli_core < 0.05
		corrSig_pli_core = '*';
	else
		corrSig_pli_core = '';
	end
	if pval_pli_belt < 0.05
		corrSig_pli_belt = '*';
	else
		corrSig_pli_belt = '';
	end
	
	if pval_dur_core < 0.05
		corrSig_dur_core = '*';
	else
		corrSig_dur_core = '';
	end
	if pval_dur_belt < 0.05
		corrSig_dur_belt = '*';
	else
		corrSig_dur_belt = '';
	end
	
	
	h_corr = figure;
	% PLI vs. PeakDelay
	subplot_tight(1,2,1,[.1 .1])
	% Core area
	h_1a = plot(UberStat_core(:,2),UberStat_core(:,1),'r.','MarkerSize',20);
	hold on
	h_1b = plot(UberStat_core(:,2),y_pli_core,'m-','LineWidth',4);
	% Belt area
	h_2a = plot(UberStat_belt(:,2),UberStat_belt(:,1),'g.','MarkerSize',20);
	h_2b = plot(UberStat_belt(:,2),y_pli_belt,'y-','LineWidth',4);
	xlabel('Phase Locking Index','FontSize',20)
	ylabel('PeakDelay [ms]','FontSize',20)
	legend([h_1a h_2a h_1b h_2b],{'Core','Belt','Lin.Fit_{Core}','Lin.Fit_{Belt}'},'FontSize',16)
	title(sprintf('Corr_{Core} = %.3f%s, Corr_{Belt} = %.3f%s    (p_{RI} = %.0d, N_{Core}=%i, N_{Belt}=%i)',rho_pli_core,corrSig_pli_core,rho_pli_belt,corrSig_pli_belt,p_crit,length(UberStat_core),length(UberStat_belt)),'FontSize',20)
	
	% Duration vs. PeakDelay
	subplot_tight(1,2,2,[.1 .1])
	% Core area
	h_1a = plot(UberStat_core(:,3),UberStat_core(:,1),'r.','MarkerSize',20);
	hold on
	h_1b = plot(UberStat_core(:,3),y_dur_core,'m-','LineWidth',4);
	% Belt area
	h_2a = plot(UberStat_belt(:,3),UberStat_belt(:,1),'g.','MarkerSize',20);
	h_2b = plot(UberStat_belt(:,3),y_dur_belt,'y-','LineWidth',4);
	xlabel('Duration [ms]','FontSize',20)
	ylabel('PeakDelay [ms]','FontSize',20)
	legend([h_1a h_2a h_1b h_2b],{'Core','Belt','Lin.Fit_{Core}','Lin.Fit_{Belt}'},'FontSize',16)
	title(sprintf('Corr_{Core} = %.3f%s, Corr_{Belt} = %.3f%s    (p_{RI} = %.0d, N_{Core}=%i, N_{Belt}=%i)',rho_dur_core,corrSig_dur_core,rho_dur_belt,corrSig_dur_belt,p_crit,length(UberStat_core),length(UberStat_belt)),'FontSize',20)
	
	currPos = get(h_corr,'Position'); set(h_corr,'Position',[currPos(1),currPos(2),1200,500]);
	
	
	
	%% Compare distribution of slope for core vs belt
	%Area_Corr = [];
	%i = 1; j = 1;
	%for idx_result = 1:length(UberStat)
	%	if ~isnan(UberStat(idx_result).Corr_PLI_PeakDelay);
	%		if UberStat(idx_result).Area == 'core'
	%			Area_Corr(i,1) = UberStat(idx_result).Corr_PLI_PeakDelay;
	%			i = i+1;
	%		elseif UberStat(idx_result).Area == 'belt'
	%			Area_Corr(j,2) = UberStat(idx_result).Corr_PLI_PeakDelay;
	%			j = j+1;
	%		end
	%	end
	%end
	%numCore = i-1;
	%numBelt = j-1;
	%
	%h_corr_bar = figure;
	%h_1a = bar(1,mean(Area_Corr(1:numCore,1)),'r');
	%hold on
	%h_1b = errorbar(1,mean(Area_Corr(1:numCore,1)),std(Area_Corr(1:numCore,1)),'k','LineWidth',4);
	%h_2a = bar(2,mean(Area_Corr(1:numBelt,2)),'g');
	%h_2b = errorbar(2,mean(Area_Corr(1:numBelt,2)),std(Area_Corr(1:numBelt,2)),'k','LineWidth',4);
	%legend([h_1a h_2a],{'Core','Belt'},'FontSize',16)
	%set(gca,'XTick',[])
	%set(gca,'XTickLabel',[])
	%title('Corr. PeakDelay vs. PLI','FontSize',20)
	%
	%saveas(h_corr_bar,sprintf('~/STRF/STRFParams_PLI_vs_PeakDelay_Bar.png'))

end	% End of function definition
