function parseStim_v7()

%%% This version is a possible fix for the delay problem due to the
%%% uploading time of stimulus to RX8 by constructing the stimulus and
%%% uploading it to RX8 right after it finishes parsing and simply
%%% triggering already uploaded stimulus when it gets stimulus onset signal
%%% from serial port by Labview.

% filename for mat file to save all stimuli for all trials.
ansSaveStim = input('Save the stimli(y/n)? ','s');
if strcmp(ansSaveStim,'y')
    userinStim = input('Enter filename to save stimuli: ','s');
    fnameStim = sprintf('%s_stim.mat',userinStim);
end


%% INITIAL SETUP
% Find TDT Circuit
circuitLoc = fullfile('C:','Tasks','Human_V2','Continuous_Play_FF.rcx');
% Set RX8
RP = TDTRP('C:\Tasks\Human_V2\Continuous_Play_FF.rcx','RX8');


%% Initialize Serial Port with COM Channel
delete(instrfindall)
global sPort;
sPort = serial('COM1', 'BaudRate', 9600, 'DataBits', 8);
set(sPort, 'Terminator', 'CR/LF');
fopen(sPort);




%% Load Calibration File
% [file, path] = uigetfile('*.mat', 'Select the Calibration File');
% fileName = strcat(path,file);
if ~exist('calibrationFile','var')
    fileName = 'tone_calib_current_2015051102.mat';
    calibrationFile = load(fileName);
end

% Clear
RP.Halt;
RP.ClearCOF;

%% Run trials
trialEnd = 0;
% Trial number
trialNum = 1;
while (trialEnd ~= 1)
    fprintf('\n\n===================== Start of Trial ==========================\n')
    fprintf('Trial number %i.\n',trialNum)
    
    %% LOAD AND TEST CIRCUIT
    fprintf('\nLoading Circuit... ');
    if RP.LoadCOF(circuitLoc)
%         disp(strcat('Loaded...',circuitLoc));
        fprintf('done!\n')
    else
        disp('Failed to Load Circuit');
    end
    
    % Run circuit
    if RP.Run
%         disp('Circuit Running');
    else
        disp('Error in Running Circuit');
    end    
    
    %% Get warning ready early on
    getReadyWarning_v7(4000,-1,13,RP,calibrationFile);
    
    %% WAIT FOR DATA FROM SERIAL
    % Wait for Stimulus String from Serial
    fprintf('\nWaiting for Stimulus String...')
    checkingForStimFlag = 1;
    fprintf('  --> Parsing... ')
    while checkingForStimFlag
        dat = fscanf(sPort);
%        fprintf('%s\n',dat);

        if length(dat) >= 6
            try
                [A_freq, semiSteps, targetOnset, A_loc, B_loc] = parse(dat);
                if semiSteps == 91
                    semiSteps=0.25;
                elseif semiSteps == 92
                    semiSteps=0.5;
                elseif semiSteps == 93
                    semiSteps=0.75;
                elseif semiSteps == 95
                    semiSteps=1.25;
                elseif semiSteps == 96
                    semiSteps=1.5;
                elseif semiSteps == 97
                    semiSteps=1.75;
                elseif semiSteps == 99
                    semiSteps=2.5;
                end
                fprintf('done!\n')
            catch
                disp('Error: Stimulus String in Incorrect Format')
                continue
            end
            if paramFlag==1
                checkingForStimFlag = 0;
            end
        end
    end
    % Redefine speaker locations just for convenience.
    if A_loc==23 && B_loc==23
        SpkA = 9;
        SpkB = 9;
    elseif A_loc==23
        SpkA = 9;
        SpkB = B_loc-8;
    elseif B_loc==23
        SpkA = A_loc-8;
        SpkB = 9;
    else
        SpkA = A_loc-8;
        SpkB = B_loc-8;
    end  
    
    
    %% Construct and upload stimulus to RX8
    % Get stimulus ready
    [td, streamCurr] = getReadyStimulus_v7(A_freq,semiSteps,targetOnset,A_loc,B_loc,RP,calibrationFile);
    streamAll(trialNum,:) = streamCurr;
    % Plot stimulus
    figure(91)
    clf
    hs1 = plot(td,streamCurr{1});
    hold on
    hs2 = plot(td,streamCurr{2});
    hs3 = plot(td,streamCurr{3});
    xlim([0 max(td)]);
    ylim([-2 2]);
    title(sprintf('Trial #%i: dF=%0.2f semitones, Speaker #%i #%i',trialNum,semiSteps,SpkA,SpkB))
    xlabel('time [seconds]')
    ylabel('amplitude [volts]')
    legend('tone A','tone B','target')
    drawnow

    % Save stimulus using built-in "matfile" command
    if strcmp(ansSaveStim,'y')
        tic
        fprintf('\n\nSaving stimulus... ')
        if trialNum==1
            save(fnameStim,'streamAll','-v7.3')
        else
            m = matfile(fnameStim,'Writable',true);
            m.streamAll(trialNum,:) = streamCurr;
        end
        fprintf('done!\n')
        toc
    end
    
   
    %% Wait for trial start signal from serial
    abortFlag = findTrialStart;
    if abortFlag == 1
        disp('Aborting Trial');
%         RP.Halt
        continue;
    end

    
    %% Trigger warning tone
    if all(bitget(RP.GetStatus,1:3))
        triggerWarning(RP);
    end

    %% Wait for stimulus start signal from serial
    abortFlag = findTrialStart;
    if abortFlag == 1
        disp('Aborting Trial');
        continue;
    end
    
    %% Run circuit
%     if RP.Run
% %         disp('Circuit Running');
%     else
%         disp('Error in Running Circuit');
%     end
    
    %% Trigger stimulus and monitor stop signal
    if all(bitget(RP.GetStatus,1:3))
        triggerStimulus(RP,sPort);
    end
    trialNum = trialNum+1;
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
    function [freq, steps, target, loc1, loc2] = parse(data)
        % Initialize Values
        freq='';
        steps='';
%         duration = '';
        target = '';
        loc1='';
        loc2='';
        paramFlag=0;
        
        % Suppose Stimulus String has form 'hi0000.00.0000.0.0bye'

        hi_index = strfind(data,'hi');
        bye_index = strfind(data,'bye');
        
        % Check for Proper String Format
%         if (bye_index - hi_index) == 23
        if (bye_index - hi_index) == 18
            paramFlag=1;
        else
            return
        end
        freq = str2double(data(hi_index+2:hi_index+5));
        steps = str2double(data(hi_index+7:hi_index+8)); % As Semitonal Value
%         duration = str2double(data(hi_index+10:hi_index+13));
        target = str2double(data(hi_index+10:hi_index+13));
        loc1 = str2double(data(hi_index+15));
        loc2 = str2double(data(hi_index+17));
        % Shift speaker numbers.
        if loc1 == 9
            loc1 = 23;
        else
            loc1 = loc1+8;
        end
        if loc2 == 9
            loc2 = 23;
        else
            loc2 = loc2+8;
        end
        
    end

fclose(sPort);
delete(sPort);
clear sPort
close all

%% Reset RX8
RP.ZeroTag('Wdatain');
RP.ZeroTag('Rdatain');
RP.ZeroTag('Ldatain');
RP.Halt;

end
