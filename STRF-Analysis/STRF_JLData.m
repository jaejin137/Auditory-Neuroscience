%calculating STRF with tdt output
%Trig times and spet in samples
%load('C:\work\Sorted\meta_171018_B2.mat')
%temp=SEV2mat('F:\SAM-171013\BLock-5','Channel',66,'EVENTNAME','xpz5');
% triggers=temp.xpz5.data;
% % fs=temp.xpz5.fs;
% triggers=meta.triggers;
% fs=meta.fs;
%
% [b,a]=butter(6,1000/(fs/2));
% trf=filter(b,a,triggers);
% [pks,locs]=findpeaks(trf/max(trf),'MinPeakHeight',.9);
% %%
% %there are 917 triggers on the file
%  [TrigA,TrigB]=trigfixstrf2(locs,400,917);

%%
data_file_name = 'sample_data_clust22';
dir_spr = 'Stimuli/Moving_ripple/DMR_50HZ/DNR_Cortex_96k5min_4_50.spr';

% load data
load(data_file_name);
spet = spet';
ch = 1;
fs = 24414.0625;

% calculate STRF
% TrigA
[taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN]=rtwstrfdbint(dir_spr,0,0.15,spet,TrigA,fs,80,30,'dB','MR',1700,5,'float');
% TrigB
[taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN]=rtwstrfdbint(dir_spr,0,0.15,spet,TrigB,fs,80,30,'dB','MR',1700,5,'float');

STRF1=(STRF1A+STRF1B)/2;
STRF2=(STRF2A+STRF2B)/2;
No1=No1A+No1B;
Wo1=(Wo1A+Wo1B)/2;
No2=No2A+No2B;
Wo2=(Wo2A+Wo2B)/2;
threshold = max(max(STRF1)) * 0.15;
i_exc = find(STRF1<=threshold); % excitatory part
i_inh = find(STRF1>=threshold); % inhibitory part
STRF1e = STRF1; STRF1i = STRF1;
STRF1e(i_exc) = threshold; STRF1i(i_inh) = threshold;

% get STRFData and STRF parameter
STRFData(ch) = struct('No1',No1,'Wo1',Wo1,'No2',No2,'Wo2',Wo2,'No1A',No1,'Wo1A',Wo1,'No2A',No2,'Wo2A',Wo2,'No1B',No1,'Wo1B',Wo1,'No2B',No2,'Wo2B',Wo2, ...
    'STRF1',STRF1,'STRF2',STRF2,'STRF1A',STRF1,'STRF2A',STRF2,'STRF1B',STRF1,'STRF2B',STRF2, ...
    'taxis',taxis,'faxis',faxis,'PP',PP,'SPLN',SPLN,'STRF1e',STRF1e,'STRF1i',STRF1i);
[STRF1s,Tresh1] = wstrfstat(STRFData(ch).STRF1,0.001,STRFData(ch).No1,STRFData(ch).Wo1,STRFData(ch).PP,30,'dB','MR','dB');
RF1P(ch) = strfparam(STRFData(ch).taxis,STRFData(ch).faxis,STRFData(ch).STRF1e,STRFData(ch).Wo1,STRFData(ch).PP,'MR',500,4,0.05,0.1,'y');

% display STRF
figure; %subplot(1,2,1)
taxis=(taxis)*1e3;
% faxis=(faxis)*1e3;
pcolor(taxis,log2(faxis/faxis(1)),(STRF1A+STRF1B)/2);
colormap jet;set(gca,'YDir','normal'); shading flat;
pkbf = RF1P(ch).PeakBF; % peak BF
best_freq = faxis(1) * 2^pkbf; % use peakBF insted of BFHz
latency = RF1P(ch).PeakDelay;
text(100,7.5,['BF = ' num2str(best_freq) ' Hz'],'FontWeight','bold');
text(100,7,['latency = ' num2str(latency) ' ms'],'FontWeight','bold');
set(gca,'YTickLabel',{'100','200','400','800','1.6k','3.2k','6.4k','12.8k','25.6k'});
xlabel('time [ms]'); ylabel('frequency [Hz]');
hold on;
plot(latency,pkbf,'+w','LineWidth',1.5);
title(data_file_name,'Interpreter','none');
