%% This script calculate ITPC(inter-trial phase coherence) for a given
%% Hilbert transform of LFP, which is in a cell structure, and plot the
%% result.


% Get the number of trials.
NumTrial = max(size(NeuralRawTr));
% Determine the size of longest trials.
MaxCellSize = 1;
for i=1:NumTrial
    CellSize = max(size(NeuralRawTr{1,i}));
    if CellSize > MaxCellSize
        MaxCellSize = CellSize;
    end
end
% Determine the size of shortest trials.
MinCellSize = 1000000000000;
for i=1:NumTrial
    CellSize = max(size(NeuralRawTr{1,i}));
    if CellSize < MinCellSize
        MinCellSize = CellSize;
    end
end

% Extract LFP(<600Hz) from raw waveform and parse different bands.
[b_LFP,a_LFP] = butter(4,600/(NeuralFs/2),'low');
[b_thet,a_thet] = butter(4,[4 8]/(NeuralFs/2),'bandpass');
[b_alph,a_alph] = butter(4,[9 12]/(NeuralFs/2),'bandpass');
[b_beta,a_beta] = butter(4,[13 30]/(NeuralFs/2),'bandpass');
[b_gamm,a_gamm] = butter(4,[30 90]/(NeuralFs/2),'bandpass');
[b_epsil,a_epsil] = butter(4,[90 150]/(NeuralFs/2),'bandpass');
% Suppress all warnings.
warning('off')
for i=1:NumTrial
	NeuralLFPTr{i} = filtfilt(b_LFP,a_LFP,double(NeuralRawTr{i}));
	% Get Hilbert transform of LFP.
	NeuralLFPHTTr{i} = hilbert(NeuralLFPTr{i});
	NeuralThetTr{i} = filtfilt(b_thet,a_thet,double(NeuralLFPTr{i}));
	NeuralAlphTr{i} = filtfilt(b_alph,a_alph,double(NeuralLFPTr{i}));
	NeuralBetaTr{i} = filtfilt(b_beta,a_beta,double(NeuralLFPTr{i}));
	NeuralGammTr{i} = filtfilt(b_gamm,a_gamm,double(NeuralLFPTr{i}));
	NeuralEpsilTr{i} = filtfilt(b_epsil,a_epsil,double(NeuralLFPTr{i}));
end
warning('on')


% Convert LFP from cell structure to matrix.
NeuralLFPHTMat = zeros(NumTrial,MinCellSize);
for i=1:NumTrial
    NeuralLFPHTMat(i,:) = NeuralLFPHTTr{1,i}(1,1:MinCellSize);
end

% Get the phase.
Phi = atan(imag(NeuralLFPHTMat)./real(NeuralLFPHTMat));
% Get the cosine value. 
CosPhi = cos(Phi);
% Get the sine value.
SinPhi = sin(Phi);
% Sum of cosine over trials at each time.
SumCosPhi = sum(cos(Phi),1);
% Sum of sine over trials at each time. 
SumSinPhi = sum(sin(Phi),1);
% Calculate inter-trial phase coherence.
ITPC = sqrt(SumCosPhi.^2+SumSinPhi.^2)/NumTrial;
% Plot the resulting ITPC.
figure
plot(ITPC)
ylim([0 1])
title('ITPC');

