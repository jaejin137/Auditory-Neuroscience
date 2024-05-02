%function [UberSTRF] = calculate_strf_params(monkeyName,sessionDate,driveID,goodSTRFUnits,TrigA,TrigB)

ifPlot = 'y';

% Turn off warnings
warning('off')

monkeyName = input('Enter monkey name (Cassius/Miyagi): ','s');
sessionDate = input('Enter session date (yymmdd): ','s');
driveID = input('Enter drive ID (e.g. D2_AC_R1): ','s');
goodSTRFUnits = input('Enter good STRF units as a vector: ');

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

% Some constants
Fss = 24414.0625;
Fsd = 1000;
%load(fullfile('~/STRF/Data/Cassius-190324/channel_depth.dat'))

% Some more constants
T1 = 0;
s_bin = 0.15;
T2 = s_bin;
SPL = 80;
MdB = 30;
ModType = 'dB';
Sound = 'MR';
NBlocks = 100;
UF = 10;
sprtype='float';

% Even more constants
MdB=30;
p=0.001;
ModType='dB';
Sound='MR';
SModType='dB';
MaxFm=10;
MaxRD=4;
	
tic

% Sort clusters by depth.
%[depth,idx] = sort(channel_depth(:,3));
%st_clu_depth = st_clu(idx,:);
% Number of clusters
%numClust = length(st_clu_depth(:,1));

totNumClust = length(st_clu);


% Pick a cluster.
%m = 16;
%clust2proc = [69,116,163,168,198,201,214,229];
%clust2proc = cell2mat(spikeTimesGoodClust(:,1));
clust2proc = goodSTRFUnits;
numClust = length(clust2proc);

h_strf = figure;
%parfor m = 1:numClust
for m = 1:numClust
	warning('off');
	fprintf('\nProcessing clsuter #%i out of %i ...\n',m,numClust)
	%clust1 = st_clu_depth{m,1};
	clust1 = find(cell2mat(st_clu(:,1))==clust2proc(m));
	% Split spike times into trial A and B.
	%spikeTime1 = cell2mat(st_clu_depth(find(cell2mat(st_clu_depth(:,1))==clust1),2));
	spikeTime1 = cast(st_clu{clust1,2},'double');
	%spikeInd1 = spikeTime1*Fss;
	spikeInd1 = spikeTime1;
	spet1A = round(spikeInd1(spikeInd1<=max(TrigA)));
	spet1B = round(spikeInd1(spikeInd1>=min(TrigB) & spikeInd1<=max(TrigB)));
	
	% Get STRF from spike times.
	figure(h_strf)
	[taxis,faxis,STRF1A_Ch1,STRF1A_Ch2,PP,Wo1A1,Wo1A2,No1A1,No1A2,SPLN]=rtwstrfdbint(sprfile,T1,T2,spet1A',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
	[taxis,faxis,STRF1B_Ch1,STRF1B_Ch2,PP,Wo1B1,Wo1B2,No1B1,No1B2,SPLN]=rtwstrfdbint(sprfile,T1,T2,spet1B',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
	% Average STRF over trial A and B.
	STRF1 = (STRF1A_Ch1+STRF1B_Ch1)/2;
	No1 = round((No1A1 + No1B1)/2);
	Wo1 = round((Wo1A1 + Wo1B1)/2);
	
	% Get significant STRF.
	[STRF1s,Tresh]=wstrfstat(STRF1,p,No1,Wo1,PP,MdB,ModType,Sound,SModType);
	% Get STRF parameters.
	[RFParam]=strfparam(taxis,faxis,STRF1,Wo1,PP,Sound,MaxFm,MaxRD);
	
	% Gabor fit (Skipped for now due to convergence issue.)
	N=10;
	theta=5;
	%[GModel]=strfgaborfit(STRF1s,STRF1,taxis,faxis,N,theta);

	%This routine converts the STRF into a ripple transfer function (RTF) which
	%we can then plot - you can use the ‘y’ option to plot it
	[Fm,RD,RTF,RTFVar]=strf2rtf(taxis,faxis,STRF1s,MaxFm,MaxRD,'y');
	
	% Save the results
	UberSTRF(m).ClustNum = clust2proc(m);
	UberSTRF(m).taxis = taxis;
	UberSTRF(m).faxis = faxis;
	UberSTRF(m).PP = PP;
	UberSTRF(m).STRF1 = STRF1;
	UberSTRF(m).STRFs = STRF1s;
	UberSTRF(m).RFParam = RFParam;
	%UberSTRF(m).GModel = GModel;
	UberSTRF(m).Fm = Fm;
	UberSTRF(m).RD = RD;
	UberSTRF(m).RTF = RTF;
	UberSTRF(m).RTFVar = RTFVar;
end
	
toc

%% Plot the model fit. (Excerpted from Monty's code)
if ifPlot == 'y'
	h_strf = figure;
	%h_gabor = figure;
	h_mod = figure;
	numRow = ceil(sqrt(numClust));
	numCol = ceil(sqrt(numClust));
	for m = 1:numClust
		figure(h_strf)
		Max=max(max(abs(UberSTRF(m).STRF1)))*sqrt(UberSTRF(m).PP);
		subplot(numRow,numCol,m)
		imagesc(1000*UberSTRF(m).taxis,log2(UberSTRF(m).faxis/UberSTRF(m).faxis(1)),UberSTRF(m).STRF1*sqrt(UberSTRF(m).PP))
		caxis([-Max Max])
		set(gca,'YDir','normal')
		title('Original STRF')
		%figure(h_gabor)
		%subplot(numRow,numCol,m)
		%imagesc(1000*UberSTRF(m).taxis,log2(UberSTRF(m).faxis/UberSTRF(m).faxis(1)),UberSTRF(m).GModel.STRFm2*sqrt(UberSTRF(m).PP)), %Plot second order model for this case since it improves the fit quite a bit
		%caxis([-Max Max])
		%set(gca,'YDir','normal')
		%title('Gabor Model')
		figure(h_mod)
		subplot(numRow,numCol,m)
		imagesc(1000*UberSTRF(m).taxis,log2(UberSTRF(m).faxis/UberSTRF(m).faxis(1)),(UberSTRF(m).GModel.STRFm2-UberSTRF(m).STRF1)*sqrt(UberSTRF(m).PP)), %Plot Residuals
		xlabel('Temporal Mod (Hz)')
		ylabel('Spectral Mod (cycles/oct)')
		caxis([-Max Max])
		title('Error Residuals')
	end
end
