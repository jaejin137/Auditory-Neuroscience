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
%SnipTimeStamp (goes from 1-64)+SnipTimeStamp1(1-32)=>ML
%SnipTimeStamp1 (goes from 33-128)=>AL
addpath(genpath('matlab\keck'));
cs=1;
contaML=52;
contaML2=93;
contaAL=9;
% Block = '~20180516_Ripple_u03'; % single electrode
Block = '~20200114_RippleF2_spk2_d02'; % v-probe
% for ch=7 % single electrode
for ch=1:16 % v-probe
    %     [Data] = readtank_mwa_input_tb('C:\TDT\Synapse\Tanks\Domo-180227-115025','~20180516_Ripple_u03',s,'local');
        [Data] = readtank_mwa_input_tb('F:\TDT\Synapse\Tanks\Domo-200114',Block,ch,'local');
%     [Data] = readtank_mwa_input_tb('D:\TDT\Synapse\Tanks\Domo-180702-112023',Block,ch,'local');
    TrigTimes=round(Data.Fs*Data.Trig);
    %      [TrigA,TrigB]=trigfixstrf2(TrigTimes,400,917);
%          [TrigA,TrigB]=trigfixstrf2(TrigTimes,400,902/2); % single stim presentation
    [TrigA,TrigB]=trigfixstrf2(TrigTimes,400,899); % double stim presentation
    spet=(Data.SnipTimeStamp*Data.Fs);
    spet1=(Data.SnipTimeStamp1*Data.Fs1);
    fs=Data.Fs;
    %ML
    if ch<64
        %[taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN]=rtwstrfdb('C:\work\Penn\DNR_updated\DNR_Cortex_100k5min.spr',0,0.3,spet,TrigA,fs,80,30,'dB','MR',1300,'float');
        %         [taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN]=rtwstrfdb('C:\work\Penn\Moving_Ripple_generation\DNR_Cortex_96k5min.spr',0,0.3,spet,TrigA,fs,80,30,'dB','MR',300,'float');
%         [taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN]=rtwstrfdb('C:\STRF_files\ripple_96k\DNR_Cortex_96k5min.spr',0,0.15,spet,TrigA,fs,80,30,'dB','MR',300,'float');
        [taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN]=rtwstrfdbint('F:\STRF_files\ripple_96k\DNR_Cortex_96k5min_4_50.spr',0,0.15,spet,TrigA,fs,80,30,'dB','MR',1700,2,'float');
        
        %[taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN]=rtwstrfdb('C:\work\Penn\DNR_updated\DNR_Cortex_100k5min.spr',0,0.3,spet,TrigB,fs,80,30,'dB','MR',1300,'float');
        %         [taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN]=rtwstrfdb('C:\work\Penn\Moving_Ripple_generation\DNR_Cortex_96k5min.spr',0,0.3,spet,TrigB,fs,80,30,'dB','MR',300,'float');
%         [taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN]=rtwstrfdb('C:\STRF_files\ripple_96k\DNR_Cortex_96k5min.spr',0,0.15,spet,TrigB,fs,80,30,'dB','MR',300,'float');
        [taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN]=rtwstrfdbint('F:\STRF_files\ripple_96k\DNR_Cortex_96k5min_4_50.spr',0,0.15,spet,TrigB,fs,80,30,'dB','MR',1700,2,'float');
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
        %         STRFDataML(contraML) = struct('No1',No1,'Wo1',Wo1,'No2',No2,'Wo2',Wo2,'STRF1',STRF1,'STRF2',STRF2,'taxis',taxis,'faxis',faxis,'PP',PP,'SPLN',SPLN);
%         STRFData(ch) = struct('No1',No1,'Wo1',Wo1,'No2',No2,'Wo2',Wo2,'STRF1',STRF1,'STRF2',STRF2,'taxis',taxis,'faxis',faxis,'PP',PP,'SPLN',SPLN);
        STRFData(ch) = struct('No1',No1,'Wo1',Wo1,'No2',No2,'Wo2',Wo2,'STRF1',STRF1,'STRF2',STRF2,'taxis',taxis,'faxis',faxis,'PP',PP,'SPLN',SPLN,'STRF1e',STRF1e,'STRF1i',STRF1i);
        [STRF1s,Tresh1] = wstrfstat(STRFData(ch).STRF1,0.001,STRFData(ch).No1,STRFData(ch).Wo1,STRFData(ch).PP,30,'dB','MR','dB');
%         RF1P(ch) = strfparam(STRFData(ch).taxis,STRFData(ch).faxis,STRFData(ch).STRF1,STRFData(ch).Wo1,STRFData(ch).PP,'MR',500,4,0.05,0.1,'y');
        RF1P(ch) = strfparam(STRFData(ch).taxis,STRFData(ch).faxis,STRFData(ch).STRF1e,STRFData(ch).Wo1,STRFData(ch).PP,'MR',500,4,0.05,0.1,'y');
        %         BF(ch) = faxis(1) * 2^RF1P(ch).BF;
        
        figure; %subplot(1,2,1)
        taxis=(taxis)*1e3;
        % faxis=(faxis)*1e3;
        pcolor(taxis,log2(faxis/faxis(1)),(STRF1A+STRF1B)/2);
%         pcolor(taxis,log2(faxis/faxis(1)),STRF1A); % test
        colormap jet;set(gca,'YDir','normal'); shading flat;
%         text(100,7,['BF = ' num2str(RF1P(ch).BFHz) ' Hz'],'FontWeight','bold');
        best_freq = faxis(1) * 2^RF1P(ch).PeakBF; % use peakBF insted of BFHz
        text(100,7,['BF = ' num2str(best_freq) ' Hz'],'FontWeight','bold');
        saveFileName = [Block(2:10),'depth',Block(end-1:end),'_STRF_ch',num2str(ch)];
        print(['F:\Taku\03_STRF\' saveFileName],'-djpeg');
        %   contaML=contaML+1;
        close all
        %      if ch<32
        %      [taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN]=rtwstrfdb('C:\work\Penn\DNR_updated\DNR_Cortex_100k5min.spr',0,0.3,spet1,TrigA,fs,80,1300,'dB','MR',300,'float');
        %         [taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN]=rtwstrfdb('C:\work\Penn\DNR_updated\DNR_Cortex_100k5min.spr',0,0.3,spet1,TrigB,fs,80,1300,'dB','MR',300,'float');
        %         STRF1=(STRF1A+STRF1B)/2;
        %         STRF2=(STRF2A+STRF2B)/2;
        %         No1=No1A+No1B;
        %         Wo1=(Wo1A+Wo1B)/2;
        %         No2=No2A+No2B;
        %         Wo2=(Wo2A+Wo2B)/2;
        %         STRFDataML(contaML2) = struct('No1',No1,'Wo1',Wo1,'No2',No2,'Wo2',Wo2,'STRF1',STRF1,'STRF2',STRF2,'taxis',taxis,'faxis',faxis,'PP',PP,'SPLN',SPLN);
        %         figure;%subplot(1,2,1)
        % taxis=(taxis)*1e3;
        % % faxis=(faxis)*1e3;
        % pcolor(taxis,log2(faxis/faxis(1)),(STRF1A+STRF1B)/2);
        % colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
        %   print(['C:\work\STRF_tdt\ML\ML_' num2str(contaML2)],'-djpeg');
        %      close all
        %      contaML2=contaML2+1;
        %      end
    end
    %AL
    %      if ch>32
    %      [taxis,faxis,STRF1A,STRF2A,PP,Wo1A,Wo2A,No1A,No2A,SPLN]=rtwstrfdb('C:\work\Penn\DNR_updated\DNR_Cortex_100k5min.spr',0,0.3,spet1,TrigA,fs,80,1300,'dB','MR',1300,'float');
    %         [taxis,faxis,STRF1B,STRF2B,PP,Wo1B,Wo2B,No1B,No2B,SPLN]=rtwstrfdb('C:\work\Penn\DNR_updated\DNR_Cortex_100k5min.spr',0,0.3,spet1,TrigB,fs,80,1300,'dB','MR',1300,'float');
    %         STRF1=(STRF1A+STRF1B)/2;
    %         STRF2=(STRF2A+STRF2B)/2;
    %         No1=No1A+No1B;
    %         Wo1=(Wo1A+Wo1B)/2;
    %         No2=No2A+No2B;
    %         Wo2=(Wo2A+Wo2B)/2;
    %         STRFDataAL(contaAL) = struct('No1',No1,'Wo1',Wo1,'No2',No2,'Wo2',Wo2,'STRF1',STRF1,'STRF2',STRF2,'taxis',taxis,'faxis',faxis,'PP',PP,'SPLN',SPLN);
    %         figure;%subplot(1,2,1)
    % taxis=(taxis)*1e3;
    % % faxis=(faxis)*1e3;
    % pcolor(taxis,log2(faxis/faxis(1)),(STRF1A+STRF1B)/2);
    % colormap jet;set(gca,'YDir','normal'); shading flat;colormap jet;
    %   print(['C:\work\STRF_tdt\AL\AL_' num2str(contaAL)],'-djpeg');
    %      close all
    %        contaAL=contaAL+1;
    %      end
    
end




