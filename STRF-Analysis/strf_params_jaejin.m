% Some constants
Fss = 24414.0625;
Fsd = 1000;

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
%MdB=30;
p=0.001;
ModType='dB';
Sound='MR';
SModType='dB';
%MaxFm=350;
MaxFm=50;
MaxRD=4;
	
tic


% "st_clu" is the spike times for all clusters.

% Find 'good' and 'mua' clusters.
clust_match = cellfun(@(x) contains(x,{'good','mua'}), st_clu(:,3), 'UniformOutput', 0);
clust2proc = cell2mat(st_clu(find(cell2mat(clust_match)),1));

numClust = length(clust2proc);

parfor m = 1:numClust
	clust1 = find(cell2mat(st_clu(:,1))==clust2proc(m));
	spikeTime1 = cast(st_clu{clust1,2},'double');
	spikeInd1 = spikeTime1;
	spet1A = round(spikeInd1(spikeInd1<=max(TrigA)));
	spet1B = round(spikeInd1(spikeInd1>=min(TrigB) & spikeInd1<=max(TrigB)));
	
	% Get STRF from spike times.
	[taxis,faxis,STRF1A_Ch1,STRF1A_Ch2,PP,Wo1A1,Wo1A2,No1A1,No1A2,SPLN]=rtwstrfdbint(sprfile,T1,T2,spet1A',TrigA,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
	[taxis,faxis,STRF1B_Ch1,STRF1B_Ch2,PP,Wo1B1,Wo1B2,No1B1,No1B2,SPLN]=rtwstrfdbint(sprfile,T1,T2,spet1B',TrigB,Fss,SPL,MdB,ModType,Sound,NBlocks,UF,sprtype);
	RTWSTRFDBINT(m).cluster = clust2proc(m);
	RTWSTRFDBINT(m).taxis = taxis;
	RTWSTRFDBINT(m).faxis = faxis;
	RTWSTRFDBINT(m).PP = PP;
	RTWSTRFDBINT(m).SPLN = SPLN;
	RTWSTRFDBINT(m).STRF1A1 = STRF1A_Ch1;
	RTWSTRFDBINT(m).STRF1A2 = STRF1A_Ch2;
	RTWSTRFDBINT(m).Wo1A1 = Wo1A1;
	RTWSTRFDBINT(m).Wo1A2 = Wo1A2;
	RTWSTRFDBINT(m).No1A1 = No1A1;
	RTWSTRFDBINT(m).No1A2 = No1A2;
	RTWSTRFDBINT(m).STRF1B1 = STRF1B_Ch1;
	RTWSTRFDBINT(m).STRF1B2 = STRF1B_Ch2;
	RTWSTRFDBINT(m).Wo1B1 = Wo1B1;
	RTWSTRFDBINT(m).Wo1B2 = Wo1B2;
	RTWSTRFDBINT(m).No1B1 = No1B1;
	RTWSTRFDBINT(m).No1B2 = No1B2;

	% Average STRF over trial A and B.
	STRF1 = (STRF1A_Ch1+STRF1B_Ch1)/2;
	No1 = round((No1A1 + No1B1)/2);
	Wo1 = round((Wo1A1 + Wo1B1)/2);
	
	% Get significant STRF.
	[STRF1s,Thresh]=wstrfstat(STRF1,p,No1,Wo1,PP,MdB,ModType,Sound,SModType);
	WSTRFSTAT(m).cluster = clust2proc(m);
	WSTRFSTAT(m).STRF1s = STRF1s;
	WSTRFSTAT(m).Thresh = Thresh;

	STRFPARAM(m) = strfparam(taxis,faxis,STRF1s,Wo1,PP,Sound,MaxFm,MaxRD);
	STRFPARAM(m).cluster = clust2proc(m);
	
end
	
toc

