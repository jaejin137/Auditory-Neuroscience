if ~exist('BlockName')
	display('Please specify BlockName.')
	break;
end

BlockName = strrep(BlockName,'.mat','');

load(sprintf('%s.mat',BlockName))
% Assign data with new names.
fprintf('\n========== Processing %s ==========\n',BlockName)
NeuralRaw = RawData.streams.NRaw.data;
TTLv = RawData.streams.TTLv.data(1,:);
NeuralFs = RawData.streams.NRaw.fs;
TTLvFs = RawData.streams.TTLv.fs;
FsRatio = NeuralFs/TTLvFs;
% Threshold by subtracting mean of TTLv assigning only 0 or 1.
TTLvOn = TTLv-mean(TTLv)>0;
% Take the differences to get the direction of jumps.
TTLvOnDif = diff(TTLvOn);
% Find Onset and Offset points(indices) based on the jump direction.
TTLvOnset = find(TTLvOnDif == 1);
TTLvOffset = find(TTLvOnDif == -1);
% Parse the neural data based on the onset and offset points from TTLv.
% Iterate over all trials.
% When TTLv starts with on state and only one pulse exists.
if isempty(TTLvOnset) & ~isempty(TTLvOffset)
    NeuralRawTr{1} = NeuralRaw(1:TTLvOffset(1)*FsRatio);
% When TTLv starts with off state and only one pulse exists.
elseif ~isempty(TTLvOnset) & isempty(TTLvOffset)
    NeuralRawTr{1} = NeuralRaw(TTLvOnset(1)*FsRatio:end);           
% When TTLv starts with off state:
elseif TTLvOnset(1)<TTLvOffset(1)
    for i=1:nnz(TTLvOffset)
        NeuralRawTr{i} = NeuralRaw(TTLvOnset(i)*FsRatio:TTLvOffset(i)*FsRatio);
    end
% When TTLv starts with on state:
elseif TTLvOnset(1)>TTLvOffset(1)
    NeuralRawTr{1} = NeuralRaw(1:TTLvOffset(1)*FsRatio);
    for i=2:nnz(TTLvOffset)
        NeuralRawTr{i} = NeuralRaw(TTLvOnset(i-1)*FsRatio:TTLvOffset(i)*FsRatio);
    end
end       	
% Clear NeuralRaw, RawData, and TTLv stuffs to save space.
clear NeuralRaw RawData TTLv
% Process to get MUA for each trials
for i=1:length(NeuralRawTr)
    fprintf('Processing Trial #%i ...\n',i)
	[b1,a1] = butter(4,[500 3000]/(NeuralFs/2),'bandpass');
	[b2,a2] = butter(4,600/(NeuralFs/2),'low');
	NeuralBPFTr{i} = filtfilt(b1,a1,double(NeuralRawTr{i}));
	NeuralRecTr{i} = abs(NeuralBPFTr{i});
	NeuralMUATr{i} = filtfilt(b2,a2,double(NeuralRecTr{i}));
	NeuralMUAHTTr{i} = hilbert(NeuralMUATr{i});
end
% Save all the MUAs for all the trials as a structure in one mat file for each block.
save(sprintf('%s_MUA.mat',BlockName),'NeuralRawTr','NeuralMUATr','NeuralMUAHTTr')
fprintf('\n========== MUA extraction completed for %s ==========\n',BlockName)
