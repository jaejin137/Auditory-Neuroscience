freqA=300;
freqStepMin=0;
freqStepStep=.25;
freqStepMax=12;
ampFact=0;
dur=3000;
chA=9; % Center-1
chB=16; % Center-1

% Load calibration file.
calibrationFile=load('tone_calib_current_2015051102.mat')

% Loop generateTonePattern
for freqStep = freqStepMin:freqStepStep:freqStepMax
	fprintf('\nfreqStep = %d\n',freqStep)
	generateTonePattern(freqA,freqStep,ampFact,dur,0,chA,chB,calibrationFile)
end
