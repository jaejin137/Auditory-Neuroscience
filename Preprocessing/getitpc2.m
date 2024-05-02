%% This script calculate ITPC(inter-trial phase coherence) for a given
%% Hilbert transform of MUA, which is in a cell structure, and plot the
%% result.


% Get the number of trials.
NumTrial = max(size(NeuralMUAHTTr));
% Determine the size of longest trials.
MaxCellSize = 1;
for i=1:NumTrial
    CellSize = max(size(NeuralMUAHTTr{1,i}));
    if CellSize > MaxCellSize
        MaxCellSize = CellSize;
    end
end
% Determine the size of shortest trials.
MinCellSize = 1000000000000;
for i=1:NumTrial
    CellSize = max(size(NeuralMUAHTTr{1,i}));
    if CellSize < MinCellSize
        MinCellSize = CellSize;
    end
end

% Convert MUA from cell structure to matrix.
NeuralMUAHTMat = zeros(NumTrial,MinCellSize);
for i=1:NumTrial
    NeuralMUAHTMat(i,:) = NeuralMUAHTTr{1,i}(1,1:MinCellSize);
end

% Get the phase.
Phi = atan(imag(NeuralMUAHTMat)./real(NeuralMUAHTMat));
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
figure; plot(ITPC); title('ITPC');
