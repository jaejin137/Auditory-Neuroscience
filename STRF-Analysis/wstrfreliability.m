%
%function [STRFSigData]=wstrfreliability(STRFBootData,alpha,NB)
%
%   FILE NAME   : STRF RELIABILITY
%   DESCRIPTION : Determines whether a STRF is reliable and significant
%                 above chance when compared against a spike train shuffled
%                 STRF. The program works with either one DMR trial (A) or
%                 two (A & B).
%
% STRFBootData  : Data Structure containing the following elements
%                 .taxis    - Time Axis
%                 .faxis    - Frequency Axis (Hz)
%                 .STRF1A   - STRF for channel 1 on trial A
%                 .STRF2A   - STRF for channel 2 on trial A
%                 .STRF1B   - STRF for channel 1 on trial B
%                 .STRF2B   - STRF for channel 2 on trial B
%                 .STRF1Ash - Shuffled STRF for channel 1 on trial A
%                 .STRF2Ash - Shuffled STRF for channel 2 on trial A
%                 .STRF1Bsh - Shuffled STRF for channel 1 on trial B
%                  (Optional)
%                 .STRF2Bsh - Shuffled STRF for channel 2 on trial B
%                  (Optional)
%                 .SPLN     - Sound Pressure Level per Frequency Band
%                 .No1A     - Number of spikes for trial A
%                 .No1B     - Number of spikes for trial B (Optional)
%                 .Wo1A     - Spike train for trial A
%                 .Wo1B     - Spike train for trial B (Optional)
%
%   alpha       : Significance level
%   NB          : Number of bootstraps for reliability index
%
%RETURNED VARIABLES
%   STRFSigData : Data structure containng results for significance testing
%		
%		.STRF1s - Significant STRF for channel 1
%		.STRF2s - Significant STRF for channel 2
%       .taxis	- Time Axis
%       .faxis	- Frequency Axis (Hz)
%		.Tresh1 - Threshold for channel 1
%		.Tresh2 - Threshold for channel 2
%		.sigma1 - Standard devation of noise STRF1
%		.sigma2 - Standard devation of niose STRF2
%		.P1     - Amplitude distribution for noise STRF1
%		.P2     - Amplitude distribition for noise STRF2
%		.X      - Amplitude for P1 and P2
%		.R1     - Reliability index for STRF1
%		.R2     - Reliability index for STRF2
%		.R12    - Combined reliability index for STRF1 and STRF2
%		.R1sh	- Reliability index for shuffled STRF1
%		.R2sh 	- Reliability index for shuffled STRF2
%		.R12sh 	- Combined reliability index for shuffled STRF1 and shuffled STRF2
%
% (C) Monty A. Escabi, May 2019
%
% Adopted and modified by Jaejin Lee
%

function [STRFSigData]=wstrfreliability(STRFBootData,alpha,NB)

%Concatenating Trials A and B
if isfield(STRFBootData,'STRF1B')
    %Concatenating Trials A and B STRFs
    STRF1=cat(3,STRFBootData.STRF1A,STRFBootData.STRF1B);
    STRF2=cat(3,STRFBootData.STRF2A,STRFBootData.STRF2B);
    STRF1sh=cat(3,STRFBootData.STRF1Ash,STRFBootData.STRF1Bsh);
    STRF2sh=cat(3,STRFBootData.STRF2Ash,STRFBootData.STRF2Bsh);
    
    %Concatenating Trials A and B Spike Counts and Spike Rates
    No1=cat(3,STRFBootData.No1A,STRFBootData.No1B);
    No2=cat(3,STRFBootData.No2A,STRFBootData.No2B);
    Wo1=cat(3,STRFBootData.Wo1A,STRFBootData.Wo1B);
    Wo2=cat(3,STRFBootData.Wo2A,STRFBootData.Wo2B);
else
    %If there is no Trial B, select only trial A STRF
    STRF1=STRFBootData.STRF1A;
    STRF2=STRFBootData.STRF2A;
    STRF1sh=STRFBootData.STRF1Ash;     %Shuffled STRF
    STRF2sh=STRFBootData.STRF2Ash;     %Shuffled STRF
    
    %Selecting Spike Counts and Spike Rates
    No1=STRFBootData.No1A;
    No2=STRFBootData.No2A;
    Wo1=STRFBootData.Wo1A;
    Wo2=STRFBootData.Wo2A;
end
PP=STRFBootData.PP;

%Bootstrapping Reliability Index for Shuffled and Unshuffled STRF
for i=1:NB
     
     %Display Progress
     clc
     display(['Bootstrap iteration: ' int2str(i)])
     
     %Number of STRF blocks used
     L=length(No1);     %Should be an interger multiple of 2
     
     %Estimating significant strfs for randomly choosen halfs of the data
     k=randperm(L);
     k1=k(1:L/2);
     k2=k(5:L/2+1);
     STRF1a=mean(STRF1(:,:,k1),3);
     STRF1b=mean(STRF1(:,:,k2),3);
     STRF1ash=mean(STRF1sh(:,:,k1),3);
     STRF1bsh=mean(STRF1sh(:,:,k2),3);
     Wo1a=mean(Wo1(k1));
     Wo1b=mean(Wo1(k2));
     No1a=mean(No1(k1));
     No1b=mean(No1(k2));
     k=randperm(L);
     k1=k(1:L/2);
     k2=k(5:L/2+1);
     STRF2a=mean(STRF2(:,:,k1),3);
     STRF2b=mean(STRF2(:,:,k2),3);
     STRF2ash=mean(STRF2sh(:,:,k1),3);
     STRF2bsh=mean(STRF2sh(:,:,k2),3);
     Wo2a=mean(Wo1(k1));
     Wo2b=mean(Wo1(k2));
     No2a=mean(No1(k1));
     No2b=mean(No1(k2));
     [STRF1a_sig]=wstrfstat(STRF1a,alpha,No1a,Wo1a,PP,30,'dB','MR','dB');
     [STRF1b_sig]=wstrfstat(STRF1b,alpha,No1b,Wo1b,PP,30,'dB','MR','dB');
     [STRF1ash_sig]=wstrfstat(STRF1ash,alpha,No1a,Wo1a,PP,30,'dB','MR','dB');
     [STRF1bsh_sig]=wstrfstat(STRF1bsh,alpha,No1b,Wo1b,PP,30,'dB','MR','dB');
     [STRF2a_sig]=wstrfstat(STRF2a,alpha,No2a,Wo2a,PP,30,'dB','MR','dB');
     [STRF2b_sig]=wstrfstat(STRF2b,alpha,No2b,Wo2b,PP,30,'dB','MR','dB');
     [STRF2ash_sig]=wstrfstat(STRF2ash,alpha,No2a,Wo2a,PP,30,'dB','MR','dB');
     [STRF2bsh_sig]=wstrfstat(STRF2bsh,alpha,No2b,Wo2b,PP,30,'dB','MR','dB');
    
     %Estimating Relability for original STRF data
     r=corrcoef(STRF1a_sig,STRF1b_sig);
     R1(i)=r(1,2);
     r=corrcoef(STRF2a_sig,STRF2b_sig);
     R2(i)=r(1,2);
     r=corrcoef([STRF1a_sig STRF2a_sig],[STRF1b_sig STRF2b_sig]);
     R12(i)=r(1,2);
        
     %Estimating Relability forshuffled STRF data
     r=corrcoef(STRF1ash_sig,STRF1bsh_sig);
     R1sh(i)=r(1,2);
     r=corrcoef(STRF2ash_sig,STRF2bsh_sig);
     R2sh(i)=r(1,2);
     r=corrcoef([STRF1ash_sig STRF2ash_sig],[STRF1bsh_sig STRF2bsh_sig]);
     R12sh(i)=r(1,2);
end

%Output Data Structure
STRFSigData.R1 = R1;
STRFSigData.R2 = R2;
STRFSigData.R12 = R12;
STRFSigData.R1sh = R1sh;
STRFSigData.R2sh = R2sh;
STRFSigData.R12sh = R12sh;
