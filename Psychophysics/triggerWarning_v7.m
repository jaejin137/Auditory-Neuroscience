function triggerWarning_v7(RP)

%% Basic parameters
toneDuration = 150; % In ms
toneSilence = 50; % In ms
fs = RP.GetSFreq(); % Sampling Rate
% % Warning sound at the center speaker
% B_loc = A_loc;
% Null channel
N_loc = 18;
% Number of warning tones
numTones = 2;
toneBuffer = numTones*(ceil(fs*(toneDuration+toneSilence)/1000)+1);

%% Trigger the stimulus
fprintf('\nTriggering warning... ')
RP.SoftTrg(1);
fprintf('done!\n')


% Ensure tone is played in full
nIndex = 0;
while nIndex < toneBuffer
    curIndex = RP.GetTagVal('Windex');
    if (curIndex > nIndex) nIndex = curIndex;
    else break;
    end
end

% Stop Warning
RP.SoftTrg(2);
% RP.Halt;
% Avoid coflict in channel number with one of stimulus tones
RP.SetTagVal('DAC_WChan',N_loc);

%% Reset
% RP.ZeroTag('Wdatain');
% RP.Halt;
    

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
