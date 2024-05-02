function triggerStimulus_v7(RP,sPort)

%%% This function only triggers already uploaded stimulus in RX8 and
%%% monitor a stop signal


%% Prepare serial port communication
% delete(instrfindall)
% global sPort;
% sPort = serial('COM1', 'BaudRate', 9600, 'DataBits', 8);
% set(sPort, 'Terminator', 'CR/LF');
% fopen(sPort);

%% Basic parameters
toneDuration = 50; % In ms
toneSilence = 50; % In ms
fs = RP.GetSFreq(); % Sampling Rate
% Total durtion of stimulus.
duration = 7200;
% Calculate tone buffer size and number of tones.
numTones = floor(duration/(toneDuration+toneSilence));
toneBuffer = numTones*(ceil(fs*(toneDuration+toneSilence)/1000+1));

%% Trigger the stimulus
fprintf('\nTriggering stimulus... ')
RP.SoftTrg(3);
fprintf('done!\n')


%% Monitor stop signal and abort if necessary
fprintf('\nMonitoring stop signal...\n')
% Wait for Stop Signal from Serial
stopFlag = findTrialStop(sPort);
if stopFlag == 1
    RP.SoftTrg(4);
%     RP.Halt;
%     RP.ClearCOF;
    return;
end

%% Timer to Terminate Each Tone Pair
nIndex = 0;
while nIndex < toneBuffer
    curIndex = RP.GetTagVal('index');
    if (curIndex > nIndex) nIndex = curIndex;
    else break;
    end
end

%% Reset
% RP.SoftTrg(4);
% RP.Halt;
% RP.ZeroTag('Rdatain');
% RP.ZeroTag('Ldatain');




function [stopFlag] = findTrialStop(sPort)
    checkStopFlag = 1;
    while checkStopFlag == 1
        stopFlag = [];
        dat = fscanf(sPort);
        if strfind(dat,'9')
            stopFlag = 1;
            checkStopFlag = 0;
            disp('Aborting Trial');
        end
    end

    
