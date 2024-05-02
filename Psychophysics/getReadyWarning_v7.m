function getReadyWarning_v7(freqA, stepSize, A_loc, RP, calibrationFile)

%% Load Calibration File
% [file, path] = uigetfile('*.mat', 'Select the Calibration File');
% fileName = strcat(path,file);
% if ~exist('calibrationFile','var')
%     fileName = 'tone_calib_current_2015051102.mat';
%     calibrationFile = load(fileName);
% end

%% SET PARAMETERS
% RP = TDTRP('C:\Users\arl\Desktop\Stim Parsing\Continuous_Play_v5.rcx', 'RX8');
toneDuration = 150; % In ms
toneSilence = 50; % In ms
rampLength = 5; % In ms
fs = RP.GetSFreq(); % Sampling Rate
% Hardcoded Base dB
ampA = 62;
% Frequency for B tone
freqB = freqA*(2^(stepSize/12));
% Warning sound at the center speaker
B_loc = A_loc;
% Number of warning tones
numTones = 2;
toneBuffer = numTones*(ceil(fs*(toneDuration+toneSilence)/1000)+1);
% Length of single tone vector
toneLeng = ceil(fs*(toneDuration+toneSilence)/1000)+1;


%% SET UP CALIBRATION ARRAY
% Allocate 3-dimensional array for calibration map (Freqs, dBs, Speakers)
[~,sCount] = size(calibrationFile.cTDT);
[~,fCount] = size(calibrationFile.cTDT(1).freqMappings);
[dCount,~] = size(calibrationFile.cTDT(1).freqMappings(1).dBMappings);

calibmap = zeros(sCount,fCount,dCount);

% Assign voltage values for each speaker, frequency, and dB value.
for i=1:sCount
    for j=1:fCount
        for k=1:dCount
            calibmap(i,j,k) = calibrationFile.cTDT(i).freqMappings(j).dBMappings(k,2);
        end
    end
end


%% Time vectors.
tt = 0:1/fs:toneDuration/1000; % Time Matrix for Sound Onset
ts = 0:1/fs:toneSilence/1000; % Time Matrix for Silence
cosRamp = (rampLength/1000)*fs; % Cos Ramp Duration


%% Construct full stream of tones and upload it ot RX8
fprintf('\nConstructing warning... ')
%% Construct full streams
% Arrays for full stream
streamA = zeros(1,numTones*toneLeng);
% streamT = zeros(1,numTones*toneLeng);
streamB = zeros(1,numTones*toneLeng);

% Initialize tone number
toneNum = 1;

% First tone
A_voltage = setCalibrationByIndex(ampA, freqA, A_loc, calibrationFile, calibmap);
toneA = A_voltage*sin(2*pi*freqA*tt);
toneA = pa_ramp(toneA, cosRamp, fs); % Add Ramp
toneA = toneA'; % Rotate Matrix
toneA = [toneA 0*ts]; % Add Silence
% toneBuffer
% toneLeng
% size(toneA)
streamA(1:toneLeng) = toneA;
toneNum = toneNum+1;

% Second tone
ampB = ampA;
B_voltage = setCalibrationByIndex(ampB, freqB, B_loc, calibrationFile, calibmap);
toneB = B_voltage*sin(2*pi*freqB*tt); % Generate Waveform
toneB = pa_ramp(toneB,cosRamp, fs); % Add Ramp
toneB = toneB'; % Rotate Matrix
toneB = [toneB 0*ts]; % Add Silence
streamB((toneNum-1)*toneLeng+1:toneNum*toneLeng) = toneB;

% Merge stream A and stream B
if A_loc == B_loc
    streamA = streamA+streamB;
%     streamB = 0*streamB;
end
fprintf('done!')


%% Upload to RX8
fprintf('  -->  Uploading Warning... ')
% Full warning stream A to speaker location A
RP.SetTagVal('WBuffer', toneBuffer); % Set buffer size
RP.SetTagVal('DAC_WChan', A_loc); % Speaker Switching Code
RP.WriteTagVEX('Wdatain', 0, 'F32', streamA);

fprintf('done!\n')

% Check sizes
% size(tt)
% size(ts)
% size(toneA)
% toneLeng
% numTones
% size(streamA)
% toneBuffer
% return

% Plot stream
% figure(98)
% plot(streamA)
% hold on
% plot(streamB)


    

%% SCRIPT FOR FINDING INDEXING VOLTAGE FROM CALIBRATION FILE
% function V = setCalibrationBySearch(dB, frequency, speaker, calibrationFile)
% freqMaps = calibrationFile(speaker).freqMappings;
% for i = 1:numel(freqMappings)
%     if freqMaps(i).value == frequency
%         decibelMaps = freqMaps(i).dBMappings;
%         for j = 1:numel(decibelMaps)
%             if decibelMaps(j,1) == dB
%                 V = decibelMaps(j,2);
%                 break
%             end
%         end
%     break    
%     end
% end

function V = setCalibrationByIndex(dB, frequency, speaker, calibrationFile, calibmap)
if speaker == 23
    i = 9;
else
    i = speaker-8; % Index 1: Speaker Number
end
fMaps = calibrationFile.cTDT(i).freqMappings;
fList = [];
for x = 1:numel(fMaps)
    fList = [fList fMaps(x).value];
end
fDiff = frequency - fList;
[~, j] = min(abs(fDiff)); % Index 2: Frequency

dMaps = calibrationFile.cTDT(i).freqMappings(j).dBMappings;
dList = [];
for y = 1:size(dMaps,1)
    dList = [dList dMaps(y,1)];
end
dDiff = dB - dList;
[~, k] = min(abs(dDiff)); % Index 3: Decibel

V = calibmap(i,j,k);

% Interpolate for Frequency
if fDiff(j) > 0
    fAmps = calibmap(i,:,k);
    V = interp1(fList, fAmps, frequency);
end

% Interpolate for dB
if dDiff(k) > 0
    V = V * (10^(dDiff(k)/20));
end
