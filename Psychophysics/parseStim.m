function parseStim()
%% INITIAL SETUP
% Find TDT Circuit
% OPEN LabVIEW 2009
% File --> Open BigGreenUberLoa --> Core --> ProjectExplorer -->
% BigGreenProject
% Right click on RT PXI Target and Deploy All
% File --> Open --> BigreenUberLoa --> Diagnostics --> SerialTest
% Adjust parameters and hit arrow button to run
circuitLoc = fullfile('C:','Users','arl','Desktop','Stim Parsing','Continuous_Play.rcx');
RP = TDTRP('C:\Users\arl\Desktop\Stim Parsing\Continuous_Play.rcx','RX8');

fs = RP.GetSFreq(); % Sampling Rate

% Load Calibration File
[file, path] = uigetfile('*.mat', 'Select the Calibration File');
fileName = strcat(path,file);
calibrationFile = load(fileName);

% Initialize Serial Port with COM Channel
global sPort;
sPort = serial('COM1', 'BaudRate', 9600, 'DataBits', 8);
set(sPort, 'Terminator', 'CR/LF');
fopen(sPort);

trialEnd = 0;
while (trialEnd ~= 1)
    %% LOAD AND TEST CIRCUIT
    disp('Loading Circuit...');
    if RP.LoadCOF(circuitLoc)
        disp(strcat('Loaded...',circuitLoc));
    else
        disp('Failed to Load Circuit');
    end
    
    %% WAIT FOR DATA FROM SERIAL
    % Wait for Stimulus String from Serial
    fprintf('Waiting for Stimulus String...')
    checkingForStimFlag = 1;
    while checkingForStimFlag
        dat = fscanf(sPort);
        fprintf('%s\n',dat);
        
        if length(dat) >= 6
            try
                [A_freq, semiSteps, duration, targetOnset, A_loc, B_loc] = parse(dat);
                if semiSteps == 99 
                    semiSteps=0.5;
                end
            catch
                disp('Error: Stimulus String in Incorrect Format')
                continue
            end
            if paramFlag==1
                checkingForStimFlag = 0;
            end
        end
    end
    
    % Wait for Start Signal from Serial
    abortFlag = findTrialStart;
    if abortFlag == 1
        disp('Aborting Trial');
        continue;
    end
    
    % Run Circuit
    if RP.Run
        disp('Circuit Running');
    else
        disp('Error in Running Circuit');
    end
    
    %% GENERATE AND PLAY TONES
    if all(bitget(RP.GetStatus,1:3))
        generateTonePattern(A_freq, semiSteps, duration, targetOnset, A_loc, B_loc, calibrationFile);
    end
    
    % Wait for Stop Signal from Serial
    stopFlag = findTrialStop;
    if stopFlag == 1
        RP.Halt;
        RP.ClearCOF;
        continue;
    end
end

%% SCRIPTS FOR TRIAL START STOP
    function [stopFlag] = findTrialStop
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
    end

    function abortFlag = findTrialStart
        checkAbortFlag = 1;
        while checkAbortFlag == 1
            abortFlag = [];
            
            dat2 = fscanf(sPort);
            if strfind(dat2, '9')
                abortFlag = 1;
                checkAbortFlag = 0;
            end
            if strfind(dat2, '5')
                checkAbortFlag = 0;
            end
        end
    end

%% SCRIPT FOR STRING PARSING
    function [freq, steps, duration, target, loc1, loc2] = parse(data)
        % Initialize Values
        freq='';
        steps='';
        duration = '';
        target = '';
        loc1='';
        loc2='';
        paramFlag=0;
        
        % Suppose Stimulus String has form 'hi0000.00.0000.0000.0.0bye'
        hi_index = strfind(data,'hi');
        bye_index = strfind(data,'bye');
        
        % Check for Proper String Format
        if (bye_index - hi_index) == 23
            paramFlag=1;
        else
            return
        end
        freq = str2double(data(hi_index+2:hi_index+5));
        steps = str2double(data(hi_index+7:hi_index+8)); % As Semitonal Value
        duration = str2double(data(hi_index+10:hi_index+13));
        target = str2double(data(hi_index+15:hi_index+18));
        loc1 = str2double(data(hi_index+20));
        loc2 = str2double(data(hi_index+22));
    end

fclose(sPort);
delete(sPort);
clear sPort
close all

end
