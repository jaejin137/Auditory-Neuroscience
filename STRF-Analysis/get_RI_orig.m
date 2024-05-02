% Turn off warning.
warning off


%% Define target data.
% Directory to session data
DIRSORTED = 'Cassius-190324';
% Area name
AREA = 'D2_AC_R1';
% Block name
BLOCK = 'AudiResp_24_24-190324-175634-ripplesF';
% Trigger file
FILETRIGGER = 'AudiResp_24_24-190324-175634_triggers.mat';
% Power level (Needs to be automated.)
PP = 112.5;
% Cluster numbers to process (Handpicked for now.)
clusNum = [69,116,168,214];


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
% Load the sorted data
sp = loadKSdir(fullfile('.', DIRSORTED),AREA,BLOCK);
% Load trigger file.
load(fullfile('/Users/jaejin/Research/Auditory/STRF', DIRSORTED, FILETRIGGER))
% Define cluster indices
index_clusters = unique(sp.clu);


%% Calculate RI for target clusters
for s = 1:length(clusNum)
	fprintf('\nProcessing cluster %i, %i out of %i\n',clusNum(s),s,length(clusNum))
	% Find the index of the good cluster.
	%index_temp = find(sp.clu==index_clusters(s));
	index_temp = find(sp.clu==clusNum(s));
	 
	spet = sp.st(index_temp)*fs;
	
	% Define DMR's spr file.
	sprfile = './Moving_ripple/DMR_50HZ/DNR_Cortex_96k5min_4_50.spr';
	
	
	
	% Spike events for DMR trial A and B
	SpetA = spet(spet<=(max(TrigA)));
	SpetB = spet(spet<=(max(TrigB)));
	
	% Bootstrap STRF
	[taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN] = rtwstrfdbintboot(sprfile,0,s_bin,SpetA',TrigA,fs,80,30,'dB','MR',100,10,'float',8);
	[taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN] = rtwstrfdbintboot(sprfile,0,s_bin,SpetB',TrigB,fs,80,30,'dB','MR',100,10,'float',8);
	
	% Shuffle spike events for trial A and B
	spetAshuf = shufflespet(SpetA);
	spetBshuf = shufflespet(SpetB);
	
	% Bootstrap shuffled STRF
	[taxis,faxis,STRF1Ash,STRF2Ash,PP,Wo1A,Wo2A,No1A,No2A,SPLN] = rtwstrfdbintboot(sprfile,0,s_bin,spetAshuf',TrigA,fs,80,30,'dB','MR',100,10,'float',8);
	[taxis,faxis,STRF1Bsh,STRF2Bsh,PP,Wo1B,Wo2B,No1B,No2B,SPLN] = rtwstrfdbintboot(sprfile,0,s_bin,spetBshuf',TrigB,fs,80,30,'dB','MR',100,10,'float',8);
	
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
    STRFBootData.CLUSTID = clusNum(s);
	
	% Calculate the reliability of the STRFs
	[STRFSig] = wstrfreliability(STRFBootData,alpha,NB);

	% Save the result
	STRFSigData{s,1} = clusNum(s);
    STRFSigData{s,2} = STRFBootData;
	STRFSigData{s,3} = STRFSig;
end


%% Plot the result
for s = 1:length(clusNum)
	figure
	for i = 1:3
		if i == 1
			STRFSigTemp = STRFSigData{s,3};
			[N,X] = hist(STRFSigTemp.R1,[-1:.02:1]);
			[Nshuf,X] = hist(STRFSig.R1sh,[-1:.02:1]);
			subplot(2,2,i)
			bar(X,N/10);
			hold on
			bar(X,Nshuf/10,'r')
			xlim([-.2 .2])
			ylim([0 1])
			xlabel('Reliability Index')
			ylabel('Probability')
			legend('non-Shuffled','Shuffled')
			legend('boxoff')
			legend('Location','NorthWest')
			title('R1, R1shuf')
		elseif i == 2
			STRFSigTemp = STRFSigData{s,3};
			[N,X] = hist(STRFSigTemp.R2,[-1:.02:1]);
			[Nshuf,X] = hist(STRFSig.R2sh,[-1:.02:1]);
			subplot(2,2,i)
			bar(X,N/10);
			hold on
			bar(X,Nshuf/10,'r')
			xlim([-.2 .2])
			ylim([0 1])
			xlabel('Reliability Index')
			ylabel('Probability')
			legend('non-Shuffled','Shuffled')
			legend('boxoff')
			legend('Location','NorthWest')
			title('R2, R2shuf')
		elseif i == 3
			STRFSigTemp = STRFSigData{s,3};
			[N,X] = hist(STRFSigTemp.R12,[-1:.02:1]);
			[Nshuf,X] = hist(STRFSig.R12sh,[-1:.02:1]);
			subplot(2,2,i)
			bar(X,N/10);
			hold on
			bar(X,Nshuf/10,'r')
			xlim([-.2 .2])
			ylim([0 1])
			xlabel('Reliability Index')
			ylabel('Probability')
			legend('non-Shuffled','Shuffled')
			legend('boxoff')
			legend('Location','NorthWest')
			title('R12, R12shuf')
		end
	end
	suptitle(sprintf('Cluster %i', clusNum(s)))
end

