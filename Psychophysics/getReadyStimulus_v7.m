function [td, streamCurr] = getReadyStimulus_v7(freqA, stepSize, targetOnset, A_loc, B_loc, RP, calibrationFile)
tic

%% SET PARAMETERS
% RP = TDTRP('C:\Users\arl\Desktop\Stim Parsing\Continuous_Play_v5.rcx', 'RX8');
toneDuration = 50; % In ms
toneSilence = 50; % In ms
rampLength = 7.5; % In ms
fs = RP.GetSFreq(); % Sampling Rate
% Hardcoded Base dB
ampA = 62;
% Increment in intensity in dB for target tone.
incA = 15;
% Frequency for B tone
freqB = freqA*(2^(stepSize/12));
% Amplitude Variations in B Tone.
deltaAmps = [-5 5 10 20];         
% Total durtion of stimulus.
duration = 7200; % Changed from 8000ms to 9000ms.
% Time for Catch Trial.
timeCatch = 4200;
% Calculate tone buffer size and number of tones.
numTones = floor(duration/(toneDuration+toneSilence));
toneBuffer = numTones*(ceil(fs*(toneDuration+toneSilence)/1000+1));
% Length of single tone vector
toneLeng = ceil(fs*(toneDuration+toneSilence)/1000);
% number of B tones in A-B-B-... sequence.
numBTones = 2;
% number of repetitions of A-B-B-... sequence.
totABBRep = floor(duration/((1+numBTones)*(toneDuration+toneSilence)));


% Calculate target tone location in terms of A-B-B repetition number.
if targetOnset > 0
    targetLoc = floor(targetOnset/((1+numBTones)*(toneDuration+toneSilence)))+1;
    if targetLoc > totABBRep
        display('*** Warning! Target is out of range! ***')
        targetLoc = totABBRep; % Check for End Case
    elseif targetOnset >= timeCatch
        fprintf('*** Catch Trial ***')
    end
end


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

% Time vectors.
tt = 0:1/fs:toneDuration/1000; % Time vector for a single tone
ts = 0:1/fs:toneSilence/1000; % Time vector for a silence
cosRamp = (rampLength/1000)*fs; % Cos Ramp Duration


%% Construct full stream of tones
fprintf('\nConstructing stimulus... ') 
% Arrays for full stream of each tone
streamA = zeros(1,totABBRep*(1+numBTones)*toneLeng);
streamT = zeros(1,totABBRep*(1+numBTones)*toneLeng);
streamB = zeros(1,totABBRep*(1+numBTones)*toneLeng);

% A-B-B Repetition number
ABBRepNum = 1;

while ABBRepNum < totABBRep    
    % Target tone only if targetOnset < timeCatch
    if ABBRepNum == targetLoc && targetOnset < timeCatch
        T_voltage = setCalibrationByIndex(ampA+incA, freqA, A_loc, calibrationFile, calibmap);
        toneT = T_voltage*sin(2*pi*freqA*tt);
        toneT = pa_ramp(toneT, cosRamp, fs); % Add Ramp
        toneT = toneT'; % Rotate Matrix
        toneT = [toneT 0*ts]; % Add Silence
        streamT((ABBRepNum-1)*(1+numBTones)*toneLeng+1:((ABBRepNum-1)*(1+numBTones)+1)*toneLeng) = toneT;
    else
        A_voltage = setCalibrationByIndex(ampA, freqA, A_loc, calibrationFile, calibmap);       
        toneA = A_voltage*sin(2*pi*freqA*tt);
        toneA = pa_ramp(toneA, cosRamp, fs); % Add Ramp
        toneA = toneA'; % Rotate Matrix
        toneA = [toneA 0*ts]; % Add Silence
        streamA((ABBRepNum-1)*(1+numBTones)*toneLeng+1:((ABBRepNum-1)*(1+numBTones)+1)*toneLeng) = toneA;
    end
    
    % Regular B tones before target
    BToneNum = 1;
    for k=1:numBTones
        ampB = ampA + deltaAmps(randi(numel(deltaAmps))); % Add Amplitude Variability in B Tone
        B_voltage = setCalibrationByIndex(ampB, freqB, B_loc, calibrationFile, calibmap);
        toneB = B_voltage*sin(2*pi*freqB*tt); % Generate Waveform
        toneB = pa_ramp(toneB,cosRamp, fs); % Add Ramp
        toneB = toneB'; % Rotate Matrix
        toneB = [toneB 0*ts]; % Add Silence
        streamB(((ABBRepNum-1)*(1+numBTones)+BToneNum)*toneLeng+1:(((ABBRepNum-1)*(1+numBTones)+BToneNum+1)*toneLeng)) = toneB;
        BToneNum = BToneNum+1;
    end
    ABBRepNum = ABBRepNum+1;
end

  
% When both tones are played from the same speaker
if A_loc == B_loc
    streamA = streamA+streamB;
    streamB = 0*streamB;
end
fprintf('done!')

%% Upload to RX8
fprintf('  -->  Uploading Stimulus... ')
% common stuffs
RP.SetTagVal('Buffer', toneBuffer); % Set buffer size
% full stream A to speaker location A
RP.SetTagVal('DAC_RChan', A_loc); % Speaker Switching Code
RP.WriteTagVEX('Rdatain', 0, 'F32', streamA);
% full stream T to speaker location A
RP.WriteTagVEX('Tdatain', 0, 'F32', streamT);
% full stream B to speaker location B
if A_loc ~= B_loc
    RP.SetTagVal('DAC_LChan', B_loc); % Speaker Switching Code
    RP.WriteTagVEX('Ldatain', 0, 'F32', streamB);
end
fprintf('done!\n')

% Record currnet stimulus.
streamCurr{1} = streamA;
streamCurr{2} = streamB;
streamCurr{3} = streamT;
streamCurr{4} = [A_loc B_loc];
td = 0:1/fs:(length(streamA)-1)/fs; % Time vector for total duration

toc

% %% Plot stream
% figure(99)
% clf(99)
% h1 = plot(streamA);
% hold on
% h2 = plot(streamT);
% h3 = plot(streamB);
% drawnow



% %% Set RX8
% % RP = TDTRP('C:\Users\arl\Desktop\Stim Parsing\Continuous_Play_v4.rcx', 'RX8');
% %% Basic parameters
% toneDuration = 150; % In ms
% toneSilence = 50; % In ms
% fs = RP.GetSFreq(); % Sampling Rate
% % Total durtion of stimulus.
% duration = 8000;
% % Calculate tone buffer size and number of tones.
% numTones = floor(duration/(toneDuration+toneSilence));
% toneBuffer = ceil(fs*numTones*(toneDuration+toneSilence)/1000)+1;
% %% Trigger the stimulus
% RP.SoftTrg(1);
% display('*** Stimulus triggered! ***')
% % Timer to Terminate Each Tone Pair
% nIndex = 0;
% while nIndex < toneBuffer
%     curIndex = RP.GetTagVal('index');
%     if (curIndex > nIndex) nIndex = curIndex;
%     else break;
%     end
% end
% % Reset
% RP.ZeroTag('Rdatain');
% RP.ZeroTag('Ldatain');
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
