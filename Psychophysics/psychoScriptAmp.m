freqA=1500;
freqStep=0;
ampFactMin=0;
ampFactStep=1;
ampFactMax=18;
dur=2000;
chA=9; % Center-1
chB=16; % Center-1

% Load calibration file.
calibrationFile=load('tone_calib_current_2015051102.mat')

% Loop generateTonePattern
for ampFact = ampFactMin:ampFactStep:ampFactMax
	fprintf('\nampFact = %i\n',ampFact)
	generateTonePattern(freqA,freqStep,ampFact,dur,0,chA,chB,calibrationFile)
end
