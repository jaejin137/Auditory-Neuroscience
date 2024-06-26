% Get the job type.
ans_jobtype = 0;
while ~ans_jobtype
	jobtype = input('Please specify job type, (c)onvert only, (e)xtract MUA only, or (b)oth: ','s');
	if ~(strcmp(jobtype,'c') || strcmp(jobtype,'e') || strcmp(jobtype,'b'))
		fprintf('Wrong answer!\n')
	else
		ans_jobtype = 1;
	end
end

% Simple warning and confirming.
if strcmp(jobtype,'c') || strcmp(jobtype,'b')
	ans_activex = 0;
	while ~ans_activex
		display('Activex Control is required for the job.\n');
		activex = input('Is activex installed on your system?(y/n): ','s')
		if ~(strcmp(activex,'y') || strcmp(activex,'n'))
			fprintf('Wrong answer!\n')
		else
			ans_activex = 1;
		end
	end
	if (strcmp(activex,'n'))
		display('Requested job cannot proceed without ActiveX Control. Sorry.\n')
	end
end

% Switch to a different task.
switch jobtype

case 'c'
	% Get the tank name.
	TankName = input('Enter the tank name: ','s');
	MatDirName = input('Enter a new directory name to store converted mat files: ','s');
	% Get into the tank to look up all blocks.
	cd(sprintf('%s',TankName))
	BlockList = dir;
	% Exclude anything that are NOT data blocks.
	for j=length(BlockList):-1:1
	    % Check if regular files are included as blocks.
	    if ~BlockList(j).isdir
	        BlockList(j) = [];
	    end
	    % Check if anything starting with a dot are included.
	    DirName = BlockList(j).name;
	    if DirName=='.'
	        BlockList(j) = [];
        end
        % Check if the TankName itself is included.
        %if strcmp(DirName,TankName)
        %    BlockList(j) = [];
        %end
	end
	% Get out of the tank and create a directory to store converted raw data.
	cd ..
	mkdir(sprintf('%s',MatDirName))
	cd(sprintf('%s',MatDirName))
	
	% Start to convert each block to mat files and extract MUA.
	display('Converting tank data into mat files ...\n')
	for k=1:length(BlockList)
		BlockName = BlockList(k).name;
		fprintf('\n========== TankName = %s, BlockName = %s ==========\n',TankName,BlockName);
		% Excute TDT2mat script provided by TDT.
		TDT2mat(TankName,BlockName);
        % Change name
        RawData = ans;
		% Save the result in mat format.
		save(sprintf('%s.mat',BlockName),'RawData')
		% Clear ans generated by TDT2mat script.
		clear ans;
	end
	display('Conversion completed!')

	% Go back to the parent directory.
	cd ..

case 'e'
	MatDirName = input('Enter the directory name containing raw waveforms in mat format: ','s');
	cd(sprintf('%s',MatDirName))
	BlockList = dir;
	% Exclude anything that are NOT data blocks.
	for j=length(BlockList):-1:1
	    DirName = BlockList(j).name;
        % Check if anything starting with a dot are included.
	    if DirName=='.'
	        BlockList(j) = [];
        % Check if anything already processed are included.
        elseif ~isempty(strfind(DirName,'_MUA'))
            BlockList(j) = [];
        end
	end

	% Extract MUA from raw waveform
	display('Extracting MUA from raw waveform ...')
	for k=1:length(BlockList)
		BlockName = strrep(BlockList(k).name,'.mat','');
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
        clear Neural*
	end
	display('MUA Extraction completed!')

case 'b'
    % Get the tank name.
	TankName = input('Enter the tank name: ','s');
	MatDirName = input('Enter a new directory name to store converted mat files: ','s');
	% Get into the tank to look up all blocks.
	cd(sprintf('%s',TankName))
	BlockList = dir;
	% Exclude anything that are NOT data blocks.
	for j=length(BlockList):-1:1
	    % Check if regular files are included as blocks.
	    if ~BlockList(j).isdir
	        BlockList(j) = [];
	    end
	    % Check if anything starting with a dot are included.
	    DirName = BlockList(j).name;
	    if DirName=='.'
	        BlockList(j) = [];
        end
        % Check if the TankName itself is included.
        %if strcmp(DirName,TankName)
        %    BlockList(j) = [];
        %end
	end
	% Get out of the tank and create a directory to store converted raw data.
	cd ..
	mkdir(sprintf('%s',MatDirName))
	cd(sprintf('%s',MatDirName))
	
	% Start to convert each block to mat files and extract MUA.
	display('Converting tank data into mat files ...\n')
	for k=1:length(BlockList)
		BlockName = BlockList(k).name;
		fprintf('\n========== TankName = %s, BlockName = %s ==========\n',TankName,BlockName);
		% Excute TDT2mat script provided by TDT.
		TDT2mat(TankName,BlockName);
        % Change name
        RawData = ans;
		% Save the result in mat format.
		save(sprintf('%s.mat',BlockName),'RawData')
		% Clear ans generated by TDT2mat script.
		clear ans;
	end
	display('Conversion completed!')

	% Extract MUA from raw waveform
	display('Extracting MUA from raw waveform ...')
	for k=1:length(BlockList)
		BlockName = strrep(BlockList(k).name,'.mat','');
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
			NeuralBPFTr{i} = filtfilt(b1,a1,NeuralRawTr{i});
			NeuralRecTr{i} = abs(NeuralBPFTr{i});
			NeuralMUATr{i} = filtfilt(b2,a2,NeuralRecTr{i});
			NeuralMUAHTTr{i} = hilbert(NeuralMUATr{i});
		end
		% Save all the MUAs for all the trials as a structure in one mat file for each block.
		save(sprintf('%s_MUA.mat',BlockName),'NeuralRawTr','NeuralMUATr','NeuralMUAHTTr')
	end
	display('MUA Extraction completed!')
	
end

