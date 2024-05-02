function generateTonePattern(freqA, stepSize, duration, targetOnset, A_loc, B_loc, calibrationFile)

%% CHECK FOR EXCEPTIONS
if (duration > 7600) error('Trial Duration Exceeds 7600ms');
end

if (targetOnset > duration) error('Target Time Out of Trial Duration Bounds');
end

%% SET PARAMETERS
RP = TDTRP('C:\Users\arl\Desktop\Stim Parsing\Continuous_Play.rcx', 'RX8');
toneDuration = 150; % In ms
toneSilence = 50; % In ms
rampLength = 20; % In ms
fs = RP.GetSFreq(); % Sampling Rate
numberTones = floor(duration/(2*(toneDuration+toneSilence)));
ampA = 62; % Hardcoded Base dB
deltaAmps = [4 -4 8 -8 14]; % Amplitude Variations in B Tone
target = 0;

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

%% CALCULATE TARGET TONE
if (targetOnset > 0)
    target = round(targetOnset/(2*(toneDuration + toneSilence)))+1;
    if (target > numberTones) target = numberTones; % Check for End Case
    end
end

%% CALCULATE TONE STREAMS
tonesPlayed = 0;
streamA = [];
streamB = [];
while (tonesPlayed < numberTones)
    t = 0:1/fs:toneDuration / 1000; % Time Matrix for Sound Onset
    st = 0:1/fs:toneSilence / 1000; % Time Matrix for Silence
    ft = 0:1/fs:(toneSilence+toneDuration) / 1000; % Time Matrix for Full Trial
    cosRamp = (rampLength/1000)*fs; % Cos Ramp Duration
    
    freqB = freqA*(2^(stepSize/12));
    ampB = ampA + deltaAmps(randi(numel(deltaAmps))); % Add Amplitude Variability in B Tone

    soundA = sin(2*pi*freqA*t);
    soundA = pa_ramp(soundA, cosRamp, fs); % Add Ramp
    soundA = soundA'; % Rotate Matrix
    soundA = [soundA 0*st 0*ft]; % Add Silence
    
    soundB = (ampB*0.01)*sin(2*pi*freqB*t); % Generate Waveform
    soundB = pa_ramp(soundB,cosRamp, fs); % Add Ramp
    soundB = soundB'; % Rotate Matrix
    soundB = [0*ft soundB 0*st]; % Add Silence
    
%     streamA = [streamA soundA];
%     streamB = [streamB soundB];
    
    %% SEND TO RX-8
    toneBuffer = fs*(2*(toneDuration+toneSilence))/1000;
    if (tonesPlayed == target-1)A_voltage = setCalibrationByIndex(ampA+12, freqA, A_loc, calibrationFile, calibmap); % Generate Target Tone
    else A_voltage = setCalibrationByIndex(ampA, freqA, A_loc, calibrationFile, calibmap);
    end
    B_voltage = setCalibrationByIndex(ampB, freqB, B_loc, calibrationFile, calibmap);
    
    RP.WriteTagVEX('Rdatain', 0, 'F32', soundA);
    RP.SetTagVal('R_Voltage', A_voltage);
    RP.SetTagVal('RBuffer', toneBuffer); % Set buffer size
    RP.SetTagVal('DAC_RChan', A_loc); % Speaker Switching Code

    RP.WriteTagVEX('Ldatain', 0, 'F32', soundB);
    RP.SetTagVal('L_Voltage', B_voltage);
    RP.SetTagVal('LBuffer', toneBuffer); % Set buffer size
    RP.SetTagVal('DAC_LChan', B_loc); % Speaker Switching Code
    RP.SoftTrg(1);
    
    % Timer to Terminate Each Tone Pair
    nIndex = 0;
    while nIndex < toneBuffer
        curIndex = RP.GetTagVal('index');
        if (curIndex > nIndex) nIndex = curIndex;
        else break;
        end
    end
    tonesPlayed = tonesPlayed + 1;
end
RP.ZeroTag('Rdatain');
RP.ZeroTag('Ldatain');
RP.Halt;

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
if (speaker == 23) i = 9;
else i = speaker-8; % Index 1: Speaker Number
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
if (fDiff(j) > 0)
    fAmps = calibmap(i,:,k);
    V = interp1(fList, fAmps, frequency);
end

% Interpolate for dB
if (dDiff(k) > 0)
    V = V * (10^(dDiff(k)/20));
end

%% PLAY TONES AND PLOT
% disp(streamA);
% sound(streamA,fs);
% sound(streamB,fs);
% tt = 2*linspace(0,stimDuration,length(streamA));
% plot(tt, streamA, tt, streamB);
% xlabel('Duration (ms)');
% ylabel('Volume');
% ylim([-1 1]);
% legend('A Tones','B Tones')


%% OLD CODE TO LOAD AND PLAY ENTIRE STREAM RATHER THAN INDIVIDUAL TONES
%         % Load Up Buffer with Tones A and B and Play
%         RP.WriteTagVEX('Rdatain', 0, 'F32', streamA);
%         RP.SetTagVal('RBuffer', fs); % Set buffer size
%         RP.SetTagVal('DAC_RChan', A_loc) % Speaker Switching Code
%         
%         RP.WriteTagVEX('Ldatain', 0, 'F32', streamB);
%         RP.SetTagVal('LBuffer', fs); % Set buffer size
%         RP.SetTagVal('DAC_LChan', B_loc) % Speaker Switching Code
%         RP.SoftTrg(1);
%         
%         % Timer to terminate
%         lastIndex = RP.GetTagVal('index');
%         nIndex = 0;
%         
%         while nIndex < length(streamA)
%             curIndex = RP.GetTagVal('index');
%             if (curIndex > lastIndex) diff = curIndex - lastIndex;
%             else diff = (curIndex+fs) - lastIndex;
%             end
%             nIndex = nIndex + diff;
%             lastIndex = curIndex;
%         end
%         
%         % Clear Buffer
%         RP.ZeroTag('Rdatain');
%         RP.ZeroTag('Ldatain');
%         RP.ZeroTag('DAC_RChan');
%         RP.ZeroTag('DAC_LChan');
