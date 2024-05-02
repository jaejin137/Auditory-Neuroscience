function plot_firing_rate_auto_simple(monkeyName,binSize,stepSize,Tag,p_crit)

	if nargin < 5
		p_crit = 1.;
	end


	%monkeyName = 'Cassius';
	%Tag = '210608';
	%% Basic information for firing rate results
	%% bin size in ms.
	%binSize = 5;
	%% step size in ms.
	%stepSize = 1;
	%% critical p value for RI
	%p_crit = .01;


	NBoot = 8;
	NB = 100;
	
	dir_fr = fullfile('~','STRF','FiringRate',sprintf('FR_%i_%i_%s',binSize,stepSize,Tag));
	list_fr = dir(sprintf('%s/*%s*.mat',dir_fr,monkeyName));
	
	% Basic information for STRF results
	dir_result = fullfile('~','STRF','STRFParams');
	list_result = dir(sprintf('%s/*%s*.mat',dir_result,monkeyName));
	
	
	UberDelay = struct();
	UberDelay_FR_core = [];
	UberDelay_FR_belt = [];
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
		%try
			load(fullfile(dir_fr,list_fr(idx_fr).name),'meanDelay','latency_clust');
			load(sprintf('~/Research/Auditory/STRF/RI/RI_%s-%s_%s_%i-%i.mat',monkeyName,sessionDate,driveID,NBoot,NB),'STRFSig')
			% Find common clusters both in RI result and FR result.
			idx_clust_fr = [];
			idx_clust_ri = [];
			for i = 1:size(latency_clust,1)
				for j = 1:size(STRFSig,2)
					%if STRFSig(j).clustNum == latency_clust(i,1) & ~isnan(latency_clust(i,6)) & ~isnan(STRFSig(j).p)
					if STRFSig(j).clustNum == latency_clust(i,1) & ~isnan(latency_clust(i,6)) & STRFSig(j).p < p_crit
						idx_clust_fr(end+1,1) = i;
						idx_clust_ri(end+1,1) = j;
					end
				end
			end
			UberDelay(idx_fr).Area = recArea;
			if recArea == 'core'
				UberDelay(idx_fr).AreaIdx = 1;
				UberDelay_FR_core(end+1:end+length(idx_clust_fr),1) = latency_clust(idx_clust_fr,6);
				UberRI_core(end+1:end+length(idx_clust_ri),1) = [STRFSig(idx_clust_ri).p]';
			elseif recArea == 'belt'
				UberDelay(idx_fr).AreaIdx = 2;
				UberDelay_FR_belt(end+1:end+length(idx_clust_fr),1) = latency_clust(idx_clust_fr,6);
				UberRI_belt(end+1:end+length(idx_clust_ri),1) = [STRFSig(idx_clust_ri).p]';
			end
			UberDelay(idx_fr).meanLatency = meanDelay;
		%end
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
	
	
	[h_lat, p_lat, ci_lat, stat_lat] = ttest2(UberDelay_FR_core(:,1),UberDelay_FR_belt(:,1));
	if p_lat <.05
		sigLat = '*';
	else
		sigLat = '';
	end
	
	h_latency = figure;
	Edges = [0:2:max(max(UberDelay_FR_core(:,1)),max(UberDelay_FR_belt(:,1)))];

	nRow = 1;
	nCol = 3;
	marg = [.16 .08];

	% Plot distribution across all clusters, core vs. belt.
	%subplot_tight(nRow,nCol,1,marg)
	subplot(2,1,1)
	%histogram(UberDelay_FR_core,Edges,'FaceColor','r','FaceAlpha',.5)
	histogram(UberDelay_FR_core,Edges,'Normalization','probability','FaceColor','r','FaceAlpha',.5)
	xline(median(UberDelay_FR_core,'omitnan'),'b-.',sprintf('Median = %.1fms',median(UberDelay_FR_core,'omitnan')),'LineWidth',4,'FontSize',16)
	ylim([0 .2])
	xlabel('Onset Latency [ms]','FontSize',20)
	%ylabel('No. of Clusters','FontSize',20)
	ylabel('Proportion of clusters','FontSize',20)
	%legend({'Core','Belt'},'FontSize',16)
	%title(sprintf('p%s = %.2d (two-sample t-test)',sigLat,p_lat),'FontSize',20)
	title('Delay from FR for core','FontSize',20)
	%hold on
	subplot(2,1,2)
	%histogram(UberDelay_FR_belt,Edges,'FaceColor','g','FaceAlpha',.5)
	histogram(UberDelay_FR_belt,Edges,'Normalization','probability','FaceColor','g','FaceAlpha',.5)
	xline(median(UberDelay_FR_belt,'omitnan'),'b-.',sprintf('Median = %.1fms',median(UberDelay_FR_belt,'omitnan')),'LineWidth',4,'FontSize',16)
	ylim([0 .2])
	xlabel('Onset Latency [ms]','FontSize',20)
	%ylabel('No. of Clusters','FontSize',20)
	ylabel('Proportion of clusters','FontSize',20)
	%legend({'Core','Belt'},'FontSize',16)
	%title(sprintf('p%s = %.2d (two-sample t-test)',sigLat,p_lat),'FontSize',20)
	title('Delay from FR for belt','FontSize',20)
	sgtitle(sprintf('%s (bin=%ims,step=%ims,p_{RI,crit}=%.1e)',monkeyName,binSize,stepSize,p_crit),'FontSize',24)
	currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),330,720]);


end	% End of function definition
