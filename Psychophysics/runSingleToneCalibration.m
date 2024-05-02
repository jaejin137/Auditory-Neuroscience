function runSingleToneCalibration()
% Runs calibrations for 9-speaker configuration in 2 stages:
%
%   1) Calibrates the decibel values for a specific frequency as
%   decibelMappings
%   2) Calibrates the frequency values for a specific speaker as
%   freqMappings
%   3) Loops through all 9 speakers and saves data as speakerMappings
%
% Data structure is [Target Speaker #, [Frequency Mappings]] --> [Target Frequency,
% [Decibel Mappings]] --> [Target Decibel, Voltage]
% Access as speakerMapping(i).value, speakerMapping(i).freqMappings -->
% freqMappings(i).value, freqMappings(i).dBMappings

%% SET PARAMETERS
RP = TDTRP('C:\Users\arl\Desktop\Stim Parsing\calib_tone.rcx', 'RX8');
toneDuration = 50; % In ms
rampLength = 5; % In ms
fs = RP.GetSFreq(); % Sampling Rate
t = 0:1/fs:toneDuration/1000;
cosRamp = (rampLength/1000)*fs;

VUL = 4.25; % Voltage Upper Limit
ld = 0.1; % low delta for incrementing voltage
hd = 0.3; % high delta for incrementing voltage
initV = 0.01;
acc = 0.1; % dB SPL threshold

db_list = 42:5:72; % Set dB list here

% Set frequency list here
start_freq = 300;
end_freq = 15000;
freq_list = [];

curr_freq = start_freq;
while curr_freq < end_freq
    freq_list(end+1) = curr_freq;
    curr_freq = curr_freq * (10^0.075);
end

speaker_list = 9:16; % Assign speaker numbers
speaker_list = [speaker_list 23];

speakerMappings = {}; % Populate with frequency mappings
RP.SetTagVal('SerialSize', fs*(toneDuration/1000));
micBuffer = fs*(toneDuration/1000);
RP.SetTagVal('InBufferSize', micBuffer);

%% LOOP THROUGH SPEAKERS
for x = 1:numel(speaker_list)
    RP.SetTagVal('SpeakerNum', speaker_list(x));
    disp(strcat('Calibrating Speaker_', num2str(x)));
    %% LOOP THROUGH TARGET FREQUENCIES
    freqMappings = {}; % Populate with decibel mappings
    for i = 1:numel(freq_list)
        %% TONE GENERATION
        decibelMappings = (ones(numel(db_list),2))*nan; % Populate with voltages
        sound = sin(2*pi*freq_list(i)*t);
        sound = pa_ramp(sound,cosRamp, fs);
        sound = sound';
        t0 = length(sound)/fs;
        RP.WriteTagVEX('Voc', 0, 'F32', sound);
        %% LOOP THROUGH TARGET DECIBELS
        for j = 1:numel(db_list)
            RP.SetTagVal('VocSc', initV);
            RP.SoftTrg(1); % Play the tone
            pause(t0) % Turn off sound and stop recording
            RP.SoftTrg(2);
            pause(0.2) % Stop recording
            t_outA = RP.ReadTagVEX('In1', 0, micBuffer, 'F32', 'F64', 1);
            outA = sqrt(mean(t_outA.^2)); % RMS pascal;
            curr_dB = 20.*(log10((outA./0.00002)));
            disp(strcat(num2str(initV),'V'));
            disp(strcat(num2str(curr_dB),'dB'));
            diff = db_list(j) - curr_dB;
            disp(diff);
            
            V=initV;
            
            %% CALIBRATION AND CALCULATIONS FOR VOLTAGE-DB MATCHING
            while abs(diff) >= acc
                if diff > 1
                    V = V*(1+(diff/5)*hd);
                elseif diff > 0 && diff <= 1
                    V = V*(1+(diff/5)*ld);
                elseif diff < -1
                    V = V/(1+abs(diff/5)*hd);
                elseif diff < 0 && diff >= -1
                    V = V/(1+(diff/5)*ld);
                end
                if V > VUL
                    disp('Voltage Set is Greater Than Limit');
                    diff = 0; % Leave that as NaN
                else
                    RP.SetTagVal('VocSc', V);
                    RP.SoftTrg(1);
                    pause(t0) % Turn off sound and stop recording
                    RP.SoftTrg(2);
                    pause(0.2) % Stop recording
                    t_outA = RP.ReadTagVEX('In1', 0, micBuffer, 'F32', 'F64', 1);
                    outA = sqrt(mean(t_outA.^2)); % RMS pascal;
                    curr_dB = 20.*(log10((outA./0.00002)));
                    disp(strcat(num2str(V),'V'));
                    disp(strcat(num2str(curr_dB),'dB'));
                    diff = db_list(j) - curr_dB;
                end
            end
            decibelMappings(j,1) = db_list(j);
            decibelMappings(j,2) = V;
            decibelMappings(j,3) = curr_dB;
            disp('----------CALIBRATED----------');
        end
        disp(decibelMappings);
        freqMappings(i).value = freq_list(i);
        freqMappings(i).dBMappings = decibelMappings;
    end
    speakerMappings(x).value = speaker_list(x);
    speakerMappings(x).freqMappings = freqMappings;
    disp(speakerMappings(x).value);
end

%% SAVE THE CALIBRATION DATA
str = strcat('tone_calib_current');
cTDT = speakerMappings;
cd 'C:\Users\arl\Documents\MATLAB\Calibration'
uisave('cTDT',str);
RP.Halt;
