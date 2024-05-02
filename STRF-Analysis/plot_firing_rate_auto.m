function plot_firing_rate_auto(monkeyName,binSize,stepSize)

	%monkeyName = 'Miyagi';
	% Basic information for firing rate results
	% bin size in ms.
	%binSize = 50;
	% step size in ms.
	%stepSize = 2;

	NBoot = 8;
	NB = 100;
	
	dir_fr = fullfile('~','STRF',sprintf('FR_%i_%i',binSize,stepSize));
	list_fr = dir(sprintf('%s/*%s*.mat',dir_fr,monkeyName));
	
	% Basic information for STRF results
	dir_result = fullfile('~','STRF','STRFParams');
	list_result = dir(sprintf('%s/*%s*.mat',dir_result,monkeyName));
	
	
	UberDelay = struct();
	UberDelay_core = [];
	UberDelay_belt = [];
	UberRI_core = [];
	UberRI_belt = [];


	for idx_fr = 1:length(list_fr)
	%for idx_fr = 4
		load(fullfile(dir_result,list_result(idx_fr).name),'recArea');
		nameRoot = split(list_result(idx_fr).name,'.');
		nameTok = split(nameRoot{1},'_');
		monkeyName = nameTok{2};
		sessionDate = nameTok{3};
		driveID = sprintf('%s_%s_%s',nameTok{4:6});
		fprintf('\n%s-%s_%s, area = %s\n\n',monkeyName,sessionDate,driveID,recArea)
		try
			load(fullfile(dir_fr,list_fr(idx_fr).name),'meanDelay','onset_clust');
			load(sprintf('~/Research/Auditory/STRF/RI/RI_%s-%s_%s_%i-%i.mat',monkeyName,sessionDate,driveID,NBoot,NB),'STRFSig')
			% Find common clusters both in RI result and FR result.
			idx_clust_fr = [];
			idx_clust_ri = [];
			for i = 1:size(onset_clust,1)
				for j = 1:size(STRFSig,2)
					if STRFSig(j).clustNum == onset_clust(i,1) & ~isnan(onset_clust(i,4)) & ~isnan(STRFSig(j).p)
						idx_clust_fr(end+1,1) = i;
						idx_clust_ri(end+1,1) = j;
					end
				end
			end
			UberDelay(idx_fr).Area = recArea;
			if recArea == 'core'
				UberDelay(idx_fr).AreaIdx = 1;
				UberDelay_core(end+1:end+length(idx_clust_fr),1) = onset_clust(idx_clust_fr,4);
				UberRI_core(end+1:end+length(idx_clust_ri),1) = [STRFSig(idx_clust_ri).p]';
			elseif recArea == 'belt'
				UberDelay(idx_fr).AreaIdx = 2;
				UberDelay_belt(end+1:end+length(idx_clust_fr),1) = onset_clust(idx_clust_fr,4);
				UberRI_belt(end+1:end+length(idx_clust_ri),1) = [STRFSig(idx_clust_ri).p]';
			end
			UberDelay(idx_fr).meanLatency = meanDelay;
		end
	end
	
	Area_Latency = [];
	
	i = 1; j = 1;
	for idx_fr = 1:length(UberDelay)
		if ~isnan(UberDelay(idx_fr).meanLatency);
			if UberDelay(idx_fr).Area == 'core'
				Area_Latency(i,1) = UberDelay(idx_fr).meanLatency;
				i = i+1;
			elseif UberDelay(idx_fr).Area == 'belt'
				Area_Latency(j,2) = UberDelay(idx_fr).meanLatency;
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
	Edges = [0:5:max(max(UberDelay_core(:,1)),max(UberDelay_belt(:,1)))];
	
	% Plot distribution across all clusters, core vs. belt.
	subplot_tight(1,2,1,[0.1,0.1])
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
	subplot_tight(1,2,2,[0.1,0.1])
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
	
	%% Onset latency vs 
	%subplot_tight(1,3,3,[0.1,0.04])
	%plot(UberRI_core,UberDelay_core,'r.','MarkerSize',20)
	%hold on
	%plot(UberRI_belt,UberDelay_belt,'g.','MarkerSize',20)
	%ylabel('Onset Latency [ms]','FontSize',20)
	%xlabel('STRF Reliability Index','FontSize',20)
	%legend({'Core','Belt'},'FontSize',16)
	%
	currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),900,400]);
	saveas(h_latency,fullfile(dir_fr,'OnsetLatency_Core_vs_Belt.png'))


end	% End of function definition
