%% This script calculates reliability indices of STRFs for (sorted) spike data
%% with a given trigger file.
%% Assumes a standard kilosort2 output in the sorted data directory.

% Directory path to the sorted spike data
DIR_SORTED = fullfile('.','Cassius-190324');
% File path to the trigger file
FILE_TRIGGER = fullfile(DIR_SORTED,'AudiResp_24_24-190324-175634_triggers.mat');

% Load the sorted spike data
sp = loadKSdir(DIR_SORTED);

% Load trigger file.
load(FILE_TRIGGER)

% Define cluster indices
index_clusters = unique(sp.clu);

% Find the index of the good cluster.
s = 69;	% <--- Good cluster
index_temp = find(sp.clu==s);

% Sampling rate (neural recordings)
fs = 24414.0625;
 
% Spking events
spet = sp.st(index_temp)*fs;

% Define DMR's spr file.
sprfile = './Moving_ripple/DMR_50HZ/DNR_Cortex_96k5min_4_50.spr';

% Bin size
s_bin = 0.15;

% Power level
PP = 80;

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
STRFBootData.CLUSTID = s;


% Calculate the reliability of the STRFs
alpha = 0.05;
NB = 10;

[STRFSig] = wstrfreliability(STRFBootData,alpha,NB);


