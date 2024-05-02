function plot_strf_class(monkeyName,clustType,Tag,p_crit)

	if nargin < 4
		p_crit = 1.;
	end

	%monkeyName = 'Cassius';
	%clustType = 'Both';
	%Tag = '06292007';
	%p_crit = 1.e-8;

	dir_class = fullfile('~','STRF','STRFClass');
	list_class = dir(sprintf('%s/STRFClass_%s_%s_%s.mat',dir_class,monkeyName,clustType,Tag));
	load(fullfile(dir_class,list_class.name));
	
	% Basic information for STRF results
	dir_result = fullfile('~','STRF','STRFParams');
	list_result = dir(sprintf('%s/*%s*.mat',dir_result,monkeyName));
	
	
	UberDelay_PeakDelay_core = [];
	UberDelay_PeakDelay_belt = [];
	UberDelay_T1_10_core = [];
	UberDelay_T1_10_belt = [];
	UberDelay_T1_50_core = [];
	UberDelay_T1_50_belt = [];
	%UberRI_core = [];
	%UberRI_belt = [];

	idx_core = find(contains(STRFClass.area,'core'));
	idx_belt = find(contains(STRFClass.area,'belt'));

	% T1_10 for core area
	for i = 1:length(idx_core)
		goodClust = find_good_clusters_simple(monkeyName,STRFClass.date{idx_core(i)},p_crit);
		if ~isempty(goodClust)
			targClust = intersect(STRFClass.clusters{idx_core(i)},cell2mat(goodClust(:,1)));
			idx_match = find(ismember(STRFClass.clusters{idx_core(i)},targClust));
			if ~isempty(idx_match)
				UberDelay_PeakDelay_core(end+1:end+length(idx_match),1) = STRFClass.PeakDelay{idx_core(i)}(idx_match);
				UberDelay_T1_10_core(end+1:end+length(idx_match),1) = STRFClass.T1_10{idx_core(i)}(idx_match);
				UberDelay_T1_50_core(end+1:end+length(idx_match),1) = STRFClass.T1_50{idx_core(i)}(idx_match);
			end
		end
	end
	% T1_10 for belt area
	for i = 1:length(idx_belt)
		goodClust = find_good_clusters_simple(monkeyName,STRFClass.date{idx_belt(i)},p_crit);
		if ~isempty(goodClust)
			targClust = intersect(STRFClass.clusters{idx_belt(i)},cell2mat(goodClust(:,1)));
			idx_match = find(ismember(STRFClass.clusters{idx_belt(i)},targClust));
			if ~isempty(idx_match)
				UberDelay_PeakDelay_belt(end+1:end+length(idx_match),1) = STRFClass.PeakDelay{idx_belt(i)}(idx_match);
				UberDelay_T1_10_belt(end+1:end+length(idx_match),1) = STRFClass.T1_10{idx_belt(i)}(idx_match);
				UberDelay_T1_50_belt(end+1:end+length(idx_match),1) = STRFClass.T1_50{idx_belt(i)}(idx_match);
			end
		end
	end

	h_delay = figure;

	Edges = [0:2:max(max(UberDelay_PeakDelay_core(:,1)),max(UberDelay_PeakDelay_belt(:,1)))];
	
	subplot(2,3,1)
	histogram(UberDelay_PeakDelay_core,Edges,'Normalization','probability','FaceColor','r','FaceAlpha',.5)
	xline(median(UberDelay_PeakDelay_core,'omitnan'),'b-.',sprintf('Median = %.1fms',median(UberDelay_PeakDelay_core,'omitnan')),'LineWidth',4,'FontSize',16)
	ylim([0 .2])
	xlabel('Onset Latency [ms]','FontSize',20)
	ylabel('Proportion of clusters','FontSize',20)
	%ylabel('No. of Clusters','FontSize',20)
	%legend({'Core','Belt'},'FontSize',16)
	%title(sprintf('p%s = %.2d (two-sample t-test)',sigLat,p_lat),'FontSize',20)
	title('PeakDelay for core','FontSize',20)
	
	subplot(2,3,4)
	histogram(UberDelay_PeakDelay_belt,Edges,'Normalization','probability','FaceColor','g','FaceAlpha',.5)
	xline(median(UberDelay_PeakDelay_belt,'omitnan'),'b-.',sprintf('Median = %.1fms',median(UberDelay_PeakDelay_belt,'omitnan')),'LineWidth',4,'FontSize',16)
	ylim([0 .2])
	xlabel('Onset Latency [ms]','FontSize',20)
	ylabel('Proportion of clusters','FontSize',20)
	%ylabel('No. of Clusters','FontSize',20)
	%legend({'Core','Belt'},'FontSize',16)
	%title(sprintf('p%s = %.2d (two-sample t-test)',sigLat,p_lat),'FontSize',20)
	title('PeakDelay for belt','FontSize',20)
	

	Edges = [0:2:max(max(UberDelay_T1_10_core(:,1)),max(UberDelay_T1_10_belt(:,1)))];
	
	subplot(2,3,2)
	histogram(UberDelay_T1_10_core,Edges,'Normalization','probability','FaceColor','r','FaceAlpha',.5)
	xline(median(UberDelay_T1_10_core,'omitnan'),'b-.',sprintf('Median = %.1fms',median(UberDelay_T1_10_core,'omitnan')),'LineWidth',4,'FontSize',16)
	ylim([0 .2])
	xlabel('Onset Latency [ms]','FontSize',20)
	ylabel('Proportion of clusters','FontSize',20)
	%ylabel('No. of Clusters','FontSize',20)
	%legend({'Core','Belt'},'FontSize',16)
	%title(sprintf('p%s = %.2d (two-sample t-test)',sigLat,p_lat),'FontSize',20)
	title('T1\_10 for core','FontSize',20)
	
	subplot(2,3,5)
	histogram(UberDelay_T1_10_belt,Edges,'Normalization','probability','FaceColor','g','FaceAlpha',.5)
	xline(median(UberDelay_T1_10_belt,'omitnan'),'b-.',sprintf('Median = %.1fms',median(UberDelay_T1_10_belt,'omitnan')),'LineWidth',4,'FontSize',16)
	ylim([0 .2])
	xlabel('Onset Latency [ms]','FontSize',20)
	ylabel('Proportion of clusters','FontSize',20)
	%ylabel('No. of Clusters','FontSize',20)
	%legend({'Core','Belt'},'FontSize',16)
	%title(sprintf('p%s = %.2d (two-sample t-test)',sigLat,p_lat),'FontSize',20)
	title('T1\_10 for belt','FontSize',20)
	
	Edges = [0:2:max(max(UberDelay_T1_50_core(:,1)),max(UberDelay_T1_50_belt(:,1)))];
	
	subplot(2,3,3)
	histogram(UberDelay_T1_50_core,Edges,'Normalization','probability','FaceColor','r','FaceAlpha',.5)
	xline(median(UberDelay_T1_50_core,'omitnan'),'b-.',sprintf('Median = %.1fms',median(UberDelay_T1_50_core,'omitnan')),'LineWidth',4,'FontSize',16)
	ylim([0 .2])
	xlabel('Onset Latency [ms]','FontSize',20)
	ylabel('Proportion of clusters','FontSize',20)
	%ylabel('No. of Clusters','FontSize',20)
	%legend({'Core','Belt'},'FontSize',16)
	%title(sprintf('p%s = %.2d (two-sample t-test)',sigLat,p_lat),'FontSize',20)
	title('T1\_50 for core','FontSize',20)
	
	subplot(2,3,6)
	histogram(UberDelay_T1_50_belt,Edges,'Normalization','probability','FaceColor','g','FaceAlpha',.5)
	xline(median(UberDelay_T1_50_belt,'omitnan'),'b-.',sprintf('Median = %.1fms',median(UberDelay_T1_50_belt,'omitnan')),'LineWidth',4,'FontSize',20)
	ylim([0 .2])
	xlabel('Onset Latency [ms]','FontSize',20)
	ylabel('Proportion of clusters','FontSize',20)
	%ylabel('No. of Clusters','FontSize',20)
	%legend({'Core','Belt'},'FontSize',16)
	%title(sprintf('p%s = %.2d (two-sample t-test)',sigLat,p_lat),'FontSize',20)
	title('T1\_50 for belt','FontSize',20)
	sgtitle(sprintf('%s (p_{RI,crit}=%.1e)',monkeyName,p_crit),'FontSize',24)

	% Resize the figure
	currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1200,720]);

	return




	for idx_sess = 1:length(STRFClass.date)
		%for idx_result = 1:length(list_result)
		%	if contains(list_result(idx_result).name,STRFClass.date{idx_sess})
		%		idx_match = idx_result;
		%	end
		%end
		%load(fullfile(dir_result,list_result(idx_match).name),'recArea');
		%nameRoot = split(list_result(idx_match).name,'.');
		%nameTok = split(nameRoot{1},'_');
		%monkeyName = nameTok{2};
		%sessionDate = nameTok{3};
		%driveID = sprintf('%s_%s_%s',nameTok{4:6});
		%fprintf('\n%s-%s_%s, area = %s\n\n',monkeyName,sessionDate,driveID,recArea)

		%try
			load(fullfile(dir_class,list_class(idx_sess).name),'meanDelay','latency_clust');
			load(sprintf('~/Research/Auditory/STRF/RI/RI_%s-%s_%s_%i-%i.mat',monkeyName,sessionDate,driveID,NBoot,NB),'STRFSig')
			% Find common clusters both in RI result and FR result.
			idx_clust_class = [];
			idx_clust_ri = [];
			for i = 1:size(latency_clust,1)
				for j = 1:size(STRFSig,2)
					if STRFSig(j).clustNum == latency_clust(i,1) & ~isnan(latency_clust(i,6)) & ~isnan(STRFSig(j).p)
						idx_clust_class(end+1,1) = i;
						idx_clust_ri(end+1,1) = j;
					end
				end
			end
			UberDelay(idx_sess).Area = recArea;
			if recArea == 'core'
				UberDelay(idx_sess).AreaIdx = 1;
				UberDelay_core(end+1:end+length(idx_clust_class),1) = latency_clust(idx_clust_class,6);
				UberRI_core(end+1:end+length(idx_clust_ri),1) = [STRFSig(idx_clust_ri).p]';
			elseif recArea == 'belt'
				UberDelay(idx_sess).AreaIdx = 2;
				UberDelay_belt(end+1:end+length(idx_clust_class),1) = latency_clust(idx_clust_class,6);
				UberRI_belt(end+1:end+length(idx_clust_ri),1) = [STRFSig(idx_clust_ri).p]';
			end
			UberDelay(idx_sess).meanLatency = meanDelay;
		%end
	end
	
	Area_Latency = [];
	
	i = 1; j = 1;
	for idx_sess = 1:length(UberDelay)
		if ~isnan(UberDelay(idx_sess).meanLatency);
			if UberDelay(idx_sess).Area == 'core'
				Area_Latency(i,1) = UberDelay(idx_sess).meanLatency;
				i = i+1;
			elseif UberDelay(idx_sess).Area == 'belt'
				Area_Latency(j,2) = UberDelay(idx_sess).meanLatency;
				j = j+1;
			end
		end
	end
	numCore = i-1;
	numBelt = j-1;
	
	
	[h_lat, p_lat, ci_lat, stat_lat] = ttest2(UberDelay_core(:,1),UberDelay_belt(:,1));
	if p_lat <.05
		sigLat = '*';
	else
		sigLat = '';
	end
	
	h_latency = figure;
	Edges = [0:1:max(max(UberDelay_core(:,1)),max(UberDelay_belt(:,1)))];

	nRow = 1;
	nCol = 3;
	marg = [.16 .08];

	% Plot distribution across all clusters, core vs. belt.
	subplot_tight(nRow,nCol,1,marg)
	%histogram(UberDelay_core,Edges,'FaceColor','r','FaceAlpha',.5)
	histogram(UberDelay_core,Edges,'Normalization','probability','FaceColor','r','FaceAlpha',.5)
	hold on
	%histogram(UberDelay_belt,Edges,'FaceColor','g','FaceAlpha',.5)
	histogram(UberDelay_belt,Edges,'Normalization','probability','FaceColor','g','FaceAlpha',.5)
	xlabel('Onset Latency [ms]','FontSize',20)
	%ylabel('No. of Clusters','FontSize',20)
	ylabel('Proportion of clusters','FontSize',20)
	legend({'Core','Belt'},'FontSize',16)
	title(sprintf('p%s = %.2d (two-sample t-test)',sigLat,p_lat),'FontSize',20)
	
	% Mean and STD, core vs. belt.
	subplot_tight(nRow,nCol,2,marg)
	%h_1a = bar(1,mean(Area_Latency(1:numCore,1)),'r');
	%hold on
	%h_1b = errorbar(1,mean(Area_Latency(1:numCore,1)),std(Area_Latency(1:numCore,1)),'k','LineWidth',4);
	%h_2a = bar(2,mean(Area_Latency(1:numBelt,2)),'g');
	%h_2b = errorbar(2,mean(Area_Latency(1:numBelt,2)),std(Area_Latency(1:numBelt,2)),'k','LineWidth',4);
	h_1a = bar(1,nanmean(UberDelay_core),'r');
	hold on
	h_1b = errorbar(1,nanmean(UberDelay_core),nanstd(UberDelay_core)/sqrt(length(UberDelay_core)),'k','LineWidth',4);
	h_2a = bar(2,mean(UberDelay_belt),'g');
	h_2b = errorbar(2,mean(UberDelay_belt),std(UberDelay_belt)/sqrt(length(UberDelay_belt)),'k','LineWidth',4);
	set(gca,'XTick',[])
	set(gca,'XTickLabel',[])
	ylabel('Onset Latency [ms]','FontSize',20)
	legend([h_1a h_2a],{'Core','Belt'},'FontSize',16)
	
	%% Onset latency vs RI
	subplot_tight(nRow,nCol,3,marg)
	plot(UberRI_core,UberDelay_core,'r.','MarkerSize',20)
	hold on
	plot(UberRI_belt,UberDelay_belt,'g.','MarkerSize',20)
	ylabel('Onset Latency [ms]','FontSize',20)
	xlabel('STRF Reliability Index','FontSize',20)
	legend({'Core','Belt'},'FontSize',16)
	
	currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),nCol*400,nRow*400]);
	saveas(h_latency,fullfile(dir_class,'OnsetLatency_Core_vs_Belt.png'))


end	% End of function definition
