% Turn off warning.
warning off
close all

monkeyName = input('Enter monkey name (Cassius/Miyagi): ','s');
sessionDate = input('Enter session date (yymmdd): ','s');
driveID = input('Enter drive ID (e.g. D2_AC_R1): ','s');

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
% Power level (Needs to be automated.)
PP = 112.5;
%% Cluster numbers to process (Handpicked for now.)
%clustNum = [69,116,168,214];
%clustNum = cell2mat(st_clu(:,1));
clustNum = cell2mat(st_clu(find(contains(st_clu(:,3),{'good','mua'})),1));


%% Computational parameters
% Bin size
s_bin = 0.15;
% Significance level
alpha = .05;
% Number of bootstraps
NB = 10;
% Sampling rate (neural recordings)
fs = 24414.0625;


%% Load all the data needed.
%% Load the sorted data
%sp = loadKSdir(fullfile('.', DIRSORTED),AREA,BLOCK);
%% Load trigger file.
%load(fullfile('/Users/jaejin/Research/Auditory/STRF', DIRSORTED, FILETRIGGER))
%% Define cluster indices
%index_clusters = unique(sp.clu);


%% Calculate RI for target clusters
for s = 1:numel(clustNum)
	fprintf('\nProcessing cluster %i, %i out of %i\n',clustNum(s),s,numel(clustNum))
	% Find the index of the good cluster.
	%index_temp = find(sp.clu==index_clusters(s));
	%index_temp = find(sp.clu==clustNum(s));
	 
	%spet = sp.st(index_temp)*fs;
	spet = st_clu{find(cell2mat(st_clu(:,1))==clustNum(s)),2};

	% Spike events for DMR trial A and B
	spetA = spet(spet<=(max(TrigA)));
	%spetB = spet(spet<=(max(TrigB)));
	% Cut and offset spetB
	spetB = spet(spet>=(min(TrigB)) & spet<=(max(TrigB)))-min(TrigB);
	% Offset TrigB
	TrigB = TrigB - min(TrigB);
	
	try
		% Bootstrap STRF
		[taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN] = rtwstrfdbintboot(sprfile,0,s_bin,spet',uint64(TrigA),fs,80,30,'dB','MR',100,10,'float',8);
		[taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN] = rtwstrfdbintboot(sprfile,0,s_bin,spet',uint64(TrigB),fs,80,30,'dB','MR',100,10,'float',8);
		
		% Shuffle spike events for trial A and B
		spetAshuf = shufflespet(spetA);
		spetBshuf = shufflespet(spetB);
		
		% Bootstrap shuffled STRF
		[taxis,faxis,STRF1Ash,STRF2Ash,PP,Wo1A,Wo2A,No1A,No2A,SPLN] = rtwstrfdbintboot(sprfile,0,s_bin,spetAshuf',uint64(TrigA),fs,80,30,'dB','MR',100,10,'float',8);
		[taxis,faxis,STRF1Bsh,STRF2Bsh,PP,Wo1B,Wo2B,No1B,No2B,SPLN] = rtwstrfdbintboot(sprfile,0,s_bin,spetBshuf',uint64(TrigB),fs,80,30,'dB','MR',100,10,'float',8);
		
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
    	STRFBootData.CLUSTID = clustNum(s);
		
		% Calculate the reliability of the STRFs
		[STRFSig] = wstrfreliability(STRFBootData,alpha,NB);

		% Save the result
		STRFSigData{s,1} = clustNum(s);
		STRFSigData{s,2} = st_clu{s,3};
    	STRFSigData{s,3} = STRFBootData;
		STRFSigData{s,4} = STRFSig;

	catch
		fprintf('\n!!!Error occurred with cluster #%i. Skipping...!!!\n',clustNum(s))
	end
end


