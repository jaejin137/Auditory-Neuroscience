function [STRFSig] = get_RI_par(monkeyName,sessionDate,driveID,clustTypeName,custTag)

	%monkeyName = 'Cassius';
	%sessionDate = '190326';
	%driveID = 'D2_AC_R1';
	%clustTypeName = {'good','mua','unsorted'};

	%% Computational parameters
	% Power level (Needs to be automated.)
	PP = 112.5;
	% Bin size
	s_bin = 0.15;
	% Significance level
	alpha = .05;
	% Number of STRF Bootstrap samples
	NBoot = 8;
	% Number of bootstraps
	NB = 100;
	% Sampling rate (neural recordings)
	fs = 24414.0625;

	if nargin < 5
		custTag = sprintf('%i-%i',NBoot,NB);
	end
	
	% Turn off warning.
	warning('off')
	%close all
	
	% Load stimulus, trigger, and spike times
	sprfile = '~/STRF/Stimuli/Moving_ripple/DMR_50HZ/DNR_Cortex_96k5min_4_50.spr';
	
	%Trig = '~/STRF/Data/Cassius-190324/AudiResp_24_24-190324-175634_triggers.mat';
	trigDir = dir(fullfile('~/STRF','Triggers',sprintf('%s/*%s*.mat',monkeyName,sessionDate)));
	TrigFile = fullfile(trigDir.folder,trigDir.name);
	spikeTimeFile = fullfile('~/kiloSorted_DMR',sprintf('Mr%s-%s/%s/KS2_7_AC/ClusterInfo/spike_times_all_clust.mat',monkeyName,sessionDate,driveID));
	
	% Load trigger file and spike times (in seconds) for each cluster
	load(TrigFile)
	load(spikeTimeFile)
	st_clu = spikeTimesAllClust;
	clearvars spikeTimesAllClust;
	
	%% Define target data.
	%% Directory to session data
	%DIRSORTED = 'Cassius-190324';
	%% Area name
	%AREA = 'D2_AC_R1';
	%% Block name
	%BLOCK = 'AudiResp_24_24-190324-175634-ripplesF';
	%% Trigger file
	%FILETRIGGER = 'AudiResp_24_24-190324-175634_triggers.mat';
	%% Cluster numbers to process (Handpicked for now.)
	%clustNum = [69,116,168,214];
	%clustNum = cell2mat(st_clu(:,1));
	clustNum = cell2mat(st_clu(find(contains(st_clu(:,3),clustTypeName)),1));


	%% Load all the data needed.
	%% Load the sorted data
	%sp = loadKSdir(fullfile('.', DIRSORTED),AREA,BLOCK);
	%% Load trigger file.
	%load(fullfile('/Users/jaejin/Research/Auditory/STRF', DIRSORTED, FILETRIGGER))
	%% Define cluster indices
	%index_clusters = unique(sp.clu);

	tic
	%% Calculate RI for target clusters
	parfor idx_clust = 1:numel(clustNum)
	%for idx_clust = 1:numel(clustNum)
		fprintf('\nProcessing cluster %i, %i out of %i\n',clustNum(idx_clust),idx_clust,numel(clustNum))
		% Find the index of the good cluster.
		%index_temp = find(sp.clu==index_clusters(idx_clust));
		%index_temp = find(sp.clu==clustNum(idx_clust));
		 
		Trig = load(TrigFile);
		
		%spet = sp.st(index_temp)*fs;
		spet = st_clu{find(cell2mat(st_clu(:,1))==clustNum(idx_clust)),2};
		% Spike events for DMR trial A and B
		spetA = spet(spet<=(max(Trig.TrigA)));
		%spetB = spet(spet<=(max(Trig.TrigB)));
		% Cut and offset spetB
		spetB = spet(spet>=(min(Trig.TrigB)) & spet<=(max(Trig.TrigB)))-min(Trig.TrigB);
		% Offset Trig.TrigB
		Trig.TrigB = Trig.TrigB - min(Trig.TrigB);
		try
			% Bootstrap STRF
			fprintf('\nBootstrapping cluster %i, %i out of %i\n',clustNum(idx_clust),idx_clust,numel(clustNum))
			[taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN] = rtwstrfdbintboot(sprfile,0,s_bin,spet',uint64(Trig.TrigA),fs,80,30,'dB','MR',100,10,'float',NBoot);
			[taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN] = rtwstrfdbintboot(sprfile,0,s_bin,spet',uint64(Trig.TrigB),fs,80,30,'dB','MR',100,10,'float',NBoot);
			
			% Shuffle spike events for trial A and B
			spetAshuf = shufflespet(spetA);
			spetBshuf = shufflespet(spetB);
			
			% Bootstrap shuffled STRF
			fprintf('\nBootstrapping shuffled cluster %i, %i out of %i\n',clustNum(idx_clust),idx_clust,numel(clustNum))
			[taxis,faxis,STRF1Ash,STRF2Ash,PP,Wo1A,Wo2A,No1A,No2A,SPLN] = rtwstrfdbintboot(sprfile,0,s_bin,spetAshuf',uint64(Trig.TrigA),fs,80,30,'dB','MR',100,10,'float',NBoot);
			[taxis,faxis,STRF1Bsh,STRF2Bsh,PP,Wo1B,Wo2B,No1B,No2B,SPLN] = rtwstrfdbintboot(sprfile,0,s_bin,spetBshuf',uint64(Trig.TrigB),fs,80,30,'dB','MR',100,10,'float',NBoot);
			
			STRFBootData = struct;
			% Define STRFBootData
			STRFBootData.taxis = taxis;
			STRFBootData.faxis = faxis;
			STRFBootData.STRF1A = STRF1A;
			STRFBootData.STRF2A = STRF2A;
			STRFBootData.STRF1B = STRF1B;
			STRFBootData.STRF2B = STRF2B;
			STRFBootData.STRF1Ash = STRF1Ash;
			STRFBootData.STRF2Ash = STRF2Ash;
			STRFBootData.STRF1Bsh = STRF1Bsh;
			STRFBootData.STRF2Bsh = STRF2Bsh;
			STRFBootData.SPLN = sprfile;
			STRFBootData.No1A = No1A;
			STRFBootData.No1B = No1B;
			STRFBootData.No2A = No2A;
			STRFBootData.No2B = No2B;
			STRFBootData.Wo1A = Wo1A;
			STRFBootData.Wo1B = Wo1B;
			STRFBootData.Wo2A = Wo2A;
			STRFBootData.Wo2B = Wo2B;
			STRFBootData.PP = PP;
	    	STRFBootData.clustID = clustNum(idx_clust);
			
			% Calculate the reliability of the STRFs
			%STRFSigData = struct;
			fprintf('\nCalculating RI for cluster %i, %i out of %i\n',clustNum(idx_clust),idx_clust,numel(clustNum))
			[STRFSigData] = wstrfreliability_par(STRFBootData,alpha,NB);
			% Test significance of RI
			fprintf('\nTesting significance for cluster %i, %i out of %i\n',clustNum(idx_clust),idx_clust,numel(clustNum))
			% Use Whitney-Mann test
			[p,h,stat] = ranksum(STRFSigData.R12,STRFSigData.R12sh);
			% Use two sample t-test
			%[h,p,ci,stat] = ttest2(STRFSigData.R12,STRFSigData.R12sh);
	
			% Save the result
			STRFSig(idx_clust).clustNum = clustNum(idx_clust);
			STRFSig(idx_clust).clustType = strtrim(st_clu{find(cell2mat(st_clu(:,1))==clustNum(idx_clust)),3});
	    	%STRFSig(idx_clust).STRFBootData = STRFBootData;
			STRFSig(idx_clust).STRFSigData = STRFSigData;
			STRFSig(idx_clust).p = p;
			STRFSig(idx_clust).h = h;
			STRFSig(idx_clust).stat = stat;
			clearvars -except alpha clustNum clustTypeName custTag driveID fs monkeyName NB NBoot PP s_bin sessionDate STRFSig idx_clust
		catch
			fprintf('\n!!!Error occurred with cluster #%i. Skipping...!!!\n',clustNum(idx_clust))
		end
	end
	toc
	clearvars -except alpha clustNum clustTypeName custTag driveID fs monkeyName NB NBoot PP s_bin sessionDate STRFSig
	save(sprintf('~/STRF/RI/RI_%s-%s_%s_%s.mat',monkeyName,sessionDate,driveID,custTag))
	plot_RI
	print(sprintf('~/STRF/RI/RI_%s-%s_%s_%s.png',monkeyName,sessionDate,driveID,custTag),'-dpng') 


end	% End of function definition get_RI
