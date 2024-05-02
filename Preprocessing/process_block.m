% This script takes a block name, search the corresponding raw waveform
% and the behavioral data, get band specific LFPs and MUA, group them into
% two decision groups, segregation and intergration, and plot figures. 
%
% Required packages: b2r.m (colormap for ribbon plot)
%


% Clean start.
clear
close all

% Record processing time.
t_tic = tic;
t_run = clock;

% Suppress all warnings.
%warning('query','all');
warning('off','all');

%% Define some global constants, etc.
TANKDIR = 'Test1';
RAWDIR = 'Test1_Raw';
BEHAVDIR = 'BehaviorReplicated';
RESULTDIR = 'Test1_Proc';
DF = [.5 3 5 12];
AUDTTHRESH = 0.015; % threshold for meaningful auditory stimulus.
JOYTHRESH = 1000; % threshold for meaningful joystick movements.
TMWNDWTH = 20; % width of time window(in msec) for similarity index.
TMWNDMOV = 10; % size of moving of time window(in msec) for similarity index.
MINDUR = 500; % minimum post-stimulus duration considered as long enough.
LEFTOFFSET = 116.9; % offset for left-cut(in msec) to include last two tones.


%% Load the file and get the basic info about raw waveform
% Get the block name to process.
fprintf('\nWarning! Please make sure that the current directory contains raw waveforms. Thank you!\n\n')
%fprintf('Please hit any key to continue.\n\n')
%pause
blockname = input('Enter the block name (without .mat) to process: ','s');

% Check if the target neural data exists.
currdir = strsplit(pwd,'/');
if ~strcmp(currdir(end),RAWDIR)
	fprintf('\nThe current directory is not the RAWDIR you defined!\n')
	fprintf('Please change to the correct directory. Thank you!\n\n')
	return
else
	filelist = strsplit(ls);
	if ~ismember(sprintf('%s.mat',blockname),filelist)
	    fprintf('\nTarget block does NOT exist in the current directory.\n')
	    fprintf('Please change to the correct directory.\n')
	    return
	else
		fprintf('\nGood! Target neural data exists.\n')
	end
end
% Check if the target behavioral data exists.
cd ..
cd(BEHAVDIR)
currdirlist = strsplit(ls);
if ~ismember(sprintf('%s.mat',blockname),currdirlist)
	fprintf('\nThe corresponding behavior data could NOT be found!\n')
	fprintf('Sorry. Try another block!\n\n')
	% Go back to the directory that has raw waveforms.
	cd ..
	cd(RAWDIR)
	return
else
	fprintf('Better! Target behavioral data also exists.\n\n')
	cd ..
	cd(RAWDIR)
end

ans_saveworkspace = 0;
while ~ans_saveworkspace
    saveworkspace = input('Save the workspace?(y/n) ','s');
    if strcmp(saveworkspace,'y') | strcmp(saveworkspace,'n')
        ans_saveworkspace = 1;
	else
        fprintf('Wrong Answer!\n')
	end
end
ans_savefig = 0;
while ~ans_savefig
    savefig = input('Save the figures?(y/n) ','s');
    if strcmp(savefig,'y') | strcmp(savefig,'n')
        ans_savefig = 1;
	else
        fprintf('Wrong Answer!\n')
	end
end
%if strcmp(saveworkspace,'y') | strcmp(savefig,'y')
%	RESULTDIR = input('Enter the directory name to save the result into: ','s');
%end

% Load neural data.
fprintf('\nLoading neural data ...\n',blockname)
load(sprintf('%s.mat',blockname))
fprintf('Neural data successfully loaded!\n',blockname)

% Rename the neural data.
fprintf('\nStart processing block %s ...\n',blockname)
nraw = RawData.streams.NRaw.data;
audt = RawData.streams.Audt.data;
ttlv = RawData.streams.TTLv.data(1,:);
fs_neural = RawData.streams.NRaw.fs;
fs_ttlv = RawData.streams.TTLv.fs;
fs_ratio = fs_neural/fs_ttlv;
clear RawData

%% Split the block into individual trials based on ttlv pulses.
% Threshold by subtracting mean of ttlv assigning only 0 or 1.
ttlv_on = ttlv-mean(ttlv)>0;
% Take the differences to get the direction of jumps.
ttlv_on_dif = diff(ttlv_on);
% Find Onset and Offset points(indices) based on the jump direction.
ttlv_onset = find(ttlv_on_dif == 1);
ttlv_offset = find(ttlv_on_dif == -1);

% Parse the neural data based on the onset and offset points from ttlv.
% Iterate over all trials.
% When ttlv starts with on state and only one pulse exists.
if isempty(ttlv_onset) & ~isempty(ttlv_offset)
    audt_sec{1} = audt(1:ttlv_offset(1)*fs_ratio);
    nraw_sec{1} = nraw(1:ttlv_offset(1)*fs_ratio);
% When ttlv starts with off state and only one pulse exists.
elseif ~isempty(ttlv_onset) & isempty(ttlv_offset)
    audt_sec{1} = audt(ttlv_onset(1)*fs_ratio:end);           
    nraw_sec{1} = nraw(ttlv_onset(1)*fs_ratio:end);           
% When ttlv starts with off state:
elseif ttlv_onset(1)<ttlv_offset(1)
    for i=1:nnz(ttlv_offset)
        audt_sec{i} = audt(ttlv_onset(i)*fs_ratio:ttlv_offset(i)*fs_ratio);
        nraw_sec{i} = nraw(ttlv_onset(i)*fs_ratio:ttlv_offset(i)*fs_ratio);
    end
% When ttlv starts with on state:
elseif ttlv_onset(1)>ttlv_offset(1)
    audt_sec{1} = audt(1:ttlv_offset(1)*fs_ratio);
    nraw_sec{1} = nraw(1:ttlv_offset(1)*fs_ratio);
    for i=2:nnz(ttlv_offset)
        audt_sec{i} = audt(ttlv_onset(i-1)*fs_ratio:ttlv_offset(i)*fs_ratio);
        nraw_sec{i} = nraw(ttlv_onset(i-1)*fs_ratio:ttlv_offset(i)*fs_ratio);
    end
end

% Delete some workspaces to free memory
clear audt nraw

% Reshape audt_sec and nraw_sec.
audt_sec = audt_sec';
nraw_sec = nraw_sec';

% Take out null trials where there is no auditory stimulus.
truetrials = [];
nulltrials = [];
for i=1:length(audt_sec)
	if nnz(audt_sec{i}>AUDTTHRESH)
		truetrials(end+1) = i;
	else
		nulltrials(end+1) = i;
	end
end
%audt_sec = audt_sec(~cellfun('isempty',audt_sec));
%nraw_sec = nraw_sec(~cellfun('isempty',nraw_sec));

% Get the number of trials and the length of longest and shortest trial.
num_tottrials = length(audt_sec);
num_truetrials = length(truetrials);
fprintf('\nTotal number of trials = %i.\n',num_tottrials)
secleng = zeros(num_tottrials,1);
for i=1:num_tottrials
	secleng(i) = length(nraw_sec{i});
end
seclengmax = max(secleng);
seclengmin = min(secleng);

%% Select good trials based on the behavioral data.
% Load behavioral data.
cd ..
cd(BEHAVDIR)
load(sprintf('%s.mat',blockname))
% Rename the behavioral data. 
behavdata = eval(sprintf('%s',blockname))';
% Parse the behavioral data.
behavdec_tr = zeros(num_tottrials,6);
if length(behavdata)==num_tottrials
	for i=1:num_tottrials
	    behavdec_tr(i,:) = ambigdec_merge_neuralandbehav(blockname,i);
	end
	% Go back to the directory that has raw waveforms.
	cd ..
	cd(RAWDIR)
else
	fprintf('\nOops! The number of trials between neural and behavior data does NOT match!\n')
	fprintf('Please check out the number of trials on both data.\n')
	fprintf('Sorry!\n\n')
	% Go back to the directory that has raw waveforms.
	cd ..
	cd(RAWDIR)
	return
end

% Find the correspondence between behavioral data and neural data.
joyleng = zeros(num_tottrials,1);
fs_joy = zeros(num_tottrials,1);
for i=1:num_tottrials
	joyleng(i) = length(behavdata(i).JoystickStatus);
	fs_joy(i) = fs_neural*joyleng(i)/secleng(i);
end
joylengmax = max(joyleng);

% Define "time" based on fs_neural, fs_joy, and fs_ttlv.
t_neural = zeros(num_tottrials,seclengmax);
t_joy = zeros(num_tottrials,joylengmax);
for i=1:num_tottrials
	t_neural(i,1:secleng(i)) = (1:secleng(i))/fs_neural;
	t_joy(i,1:joyleng(i)) = (1:joyleng(i))/fs_joy(i);
end

% Find the onset time of the first stimulus and the offset time of the last stimulus in tone sequence for each trial.
n_stim = zeros(num_tottrials,2);
for i=1:num_tottrials
	if ~ismember(i,nulltrials)
		n_stim(i,1) = min(find(audt_sec{i}>AUDTTHRESH));
		n_stim(i,2) = max(find(audt_sec{i}>AUDTTHRESH));
	end
end

% Find good trials to process.
goodtrials = [];
errtrials = [];
for i=1:num_tottrials
%	fprintf('Checking trial #%i ...\n',i)
%	if i==1
%		figure
%	end
	if ismember(i,truetrials)
		% Test if joystick started at 0 and ended at 0.
		if behavdata(i).JoystickStatus(1)<JOYTHRESH & behavdata(i).JoystickStatus(end)<JOYTHRESH
			% Test if there was only one movement regardless of direction.
			if nnz(diff(behavdata(i).JoystickStatus)>1000)==1 & nnz(diff(behavdata(i).JoystickStatus)<-1000)==1
				% Test if it happened after the end of tone sequence.
%				if min(find(diff(behavdata(i).JoystickStatus)>1000))/fs_joy(i)>(n_stim(i,2)/fs_neural+TMWNDWTH/1000)
				if min(find(diff(behavdata(i).JoystickStatus)>1000))/fs_joy(i)-(n_stim(i,2)/fs_neural) > MINDUR/1000
					goodtrials(end+1) = i;
					%fprintf('    --> Good trial!\n')
				else
					%fprintf('    --> Too early!\n')
					errtrials(end+1) = i;
				end
			else
				%fprintf('    --> Multiple movements!\n')
				errtrials(end+1) = i;
			end
		else
			%fprintf('    --> Wrong start or end!\n')
			errtrials(end+1) = i;
		end
		% Plot the auditory stimuli and the joystick status.
		%hold off
		%plot(t_joy(1,1:joyleng(i)),behavdata(i).JoystickStatus)
		%hold on
		%plot(t_neural(1,1:secleng(i)),audt_sec{i}/max(audt_sec{i})*max(behavdata(i).JoystickStatus),'g')
		%legend(sprintf('trial #%i',i))
		%hold off
		%pause
	end
end
% number of good trials.
num_goodtrials = length(goodtrials);
if num_goodtrials
	fprintf('Number of good trials = %i.\n',num_goodtrials)
	fprintf('Good trials are %s.\n\n',num2str(goodtrials))
	pause(2)
	%ans_ifcontinue = 0;
	%while ~ans_ifcontinue
	%	ifcontinue = input('Continue?(y/n): ','s');
	%	if strcmp(ifcontinue,'y')
	%		fprintf('\nContinue to process ...')
	%		ans_ifcontinue = 1;
	%	elseif strcmp(ifcontinue,'n')
	%		fprintf('\nProcess aborted!\n\n')
	%		return
	%	else
	%		fprintf('\nWrong answer!')
	%	end
	%end
else
	fprintf('\nSorry. No good trial was found. Please adjust the criteria.\n')
	fprintf('Thank you!\n\n');
	return
end

%% Remove all the elements other than good trials.
badtrials = [nulltrials'; errtrials']';
joyleng(badtrials) = [];
fs_joy(badtrials) = [];
t_joy(badtrials,:) = [];
t_neural(badtrials,:) = [];
n_stim(badtrials,:) = [];
behavdec_tr(badtrials,:) = [];
behavdata(badtrials) = [];
audt_sec(badtrials) = [];
nraw_sec(badtrials) = [];

% Determine left-cut and right-cut time for neural data.
n_cut = zeros(num_goodtrials,2);
n_dur = zeros(num_goodtrials,1);
n_leftoffset = round(LEFTOFFSET/1000*fs_neural);
for i=1:num_goodtrials
	% left-cut time
	n_cut(i,1) = n_stim(i,2)-n_leftoffset+1;
	% right-cut time
	n_cut(i,2) = round(min(find(diff(behavdata(i).JoystickStatus)>1000))/fs_joy(i)*fs_neural)-1;
	n_dur(i) = n_cut(i,2)-n_cut(i,1)+1;
end
n_durmax = max(n_dur);
n_durmin = min(n_dur);

% Assign each section with cut-out to each row of nraw matrix.
nraw_tr = zeros(num_goodtrials,n_durmax);
for i=1:num_goodtrials
	nraw_tr(i,1:n_dur(i)) = nraw_sec{i}(n_cut(i,1):n_cut(i,2));
end

% Delete some workspaces to free memory
clear audt_sec nraw_sec


%% Extract LFP trial by trial.
fprintf('\nStart processing LFP ...\n\n')

% Extract LFP(<600Hz) from raw waveform and split different bands.
fc_lfp = 600; % cuttoff frequency for LFP
pb_theta = [4 8]; % passband for theta.
pb_alpha = [8 12]; % passband for alpha.
pb_beta = [13 30]; % passband for beta.
pb_gamma = [30 90]; % passband for gamma.
pb_epsilon = [90 150]; % passband for epsilon.

% Get coefficients for nth order Butterworth filter for each band.
order = 4; % 4th order
[b_lfp,a_lfp] = butter(order,fc_lfp/(fs_neural/2),'low');
[b_theta,a_theta] = butter(order,pb_theta/(fs_neural/2),'bandpass');
[b_alpha,a_alpha] = butter(order,pb_alpha/(fs_neural/2),'bandpass');
[b_beta,a_beta] = butter(order,pb_beta/(fs_neural/2),'bandpass');
[b_gamma,a_gamma] = butter(order,pb_gamma/(fs_neural/2),'bandpass');
[b_epsilon,a_epsilon] = butter(order,pb_epsilon/(fs_neural/2),'bandpass');

% Preallocate workspace for lfps.
lfp_tr = zeros(num_goodtrials,n_durmax);
lfp_ht_tr = zeros(num_goodtrials,n_durmax);
lfp_theta_tr = zeros(num_goodtrials,n_durmax);
lfp_alpha_tr = zeros(num_goodtrials,n_durmax);
lfp_beta_tr = zeros(num_goodtrials,n_durmax);
lfp_gamma_tr = zeros(num_goodtrials,n_durmax);
lfp_epsilon_tr = zeros(num_goodtrials,n_durmax);

for i=1:num_goodtrials
    fprintf('Processing LFP for Trial #%i ...\n',goodtrials(i))
	lfp_tr(i,1:n_dur(i)) = filtfilt(b_lfp,a_lfp,nraw_tr(i,1:n_dur(i)));
	% Get Hilbert transform of LFP.
	lfp_ht_tr(i,1:n_dur(i)) = hilbert(lfp_tr(i,1:n_dur(i)));
	lfp_theta_tr(i,1:n_dur(i)) = filtfilt(b_theta,a_theta,lfp_tr(i,1:n_dur(i)));
	lfp_alpha_tr(i,1:n_dur(i)) = filtfilt(b_alpha,a_alpha,lfp_tr(i,1:n_dur(i)));
	lfp_beta_tr(i,1:n_dur(i)) = filtfilt(b_beta,a_beta,lfp_tr(i,1:n_dur(i)));
	lfp_gamma_tr(i,1:n_dur(i)) = filtfilt(b_gamma,a_gamma,lfp_tr(i,1:n_dur(i)));
	lfp_epsilon_tr(i,1:n_dur(i)) = filtfilt(b_epsilon,a_epsilon,lfp_tr(i,1:n_dur(i)));
end


%% Extract MUA trial by trial.
fprintf('\nStart processing MUA ...\n\n')
nraw_bpf_tr = zeros(num_goodtrials,n_durmax);
nraw_rec_tr = zeros(num_goodtrials,n_durmax);
mua_tr = zeros(num_goodtrials,n_durmax);
mua_ht_tr = zeros(num_goodtrials,n_durmax);
for i=1:num_goodtrials
    fprintf('Processing MUA for Trial #%i ...\n',goodtrials(i))
	pb_mua = [500 3000]; % Define passband for MUA.
	fc_mua = 600; % Define cuffoff frequency for low pass filter.
	[b_bp,a_bp] = butter(order,pb_mua/(fs_neural/2),'bandpass');
	[b_low,a_low] = butter(order,fc_mua/(fs_neural/2),'low');
	nraw_bpf_tr(i,1:n_dur(i)) = filtfilt(b_bp,a_bp,nraw_tr(i,1:n_dur(i)));
	nraw_rec_tr(i,1:n_dur(i)) = abs(nraw_bpf_tr(i,1:n_dur(i)));
	mua_tr(i,1:n_dur(i)) = filtfilt(b_low,a_low,nraw_rec_tr(i,1:n_dur(i)));
	mua_ht_tr(i,1:n_dur(i)) = hilbert(mua_tr(i,1:n_dur(i)));
end
clear nraw_bpf_tr nraw_rec_tr


%% Classify all trials into two groups, integration(0) and segregation(1).
% Get the number of occurences for each case
num_int_p5 = nnz(behavdec_tr(:,2)==0 & behavdec_tr(:,1)==.5);
num_int_3 = nnz(behavdec_tr(:,2)==0 & behavdec_tr(:,1)==3);
num_int_5 = nnz(behavdec_tr(:,2)==0 & behavdec_tr(:,1)==5);
num_int_12 = nnz(behavdec_tr(:,2)==0 & behavdec_tr(:,1)==12);
num_seg_p5 = nnz(behavdec_tr(:,2)==1 & behavdec_tr(:,1)==.5);
num_seg_3 = nnz(behavdec_tr(:,2)==1 & behavdec_tr(:,1)==3);
num_seg_5 = nnz(behavdec_tr(:,2)==1 & behavdec_tr(:,1)==5);
num_seg_12 = nnz(behavdec_tr(:,2)==1 & behavdec_tr(:,1)==12);
num_err = nnz(behavdec_tr(:,2)==66);
num_unk = nnz(behavdec_tr(:,2)==66);
% Split lfp into integration and segregation.
lfp_int_p5 = lfp_tr(find(behavdec_tr(:,2)==0 & behavdec_tr(:,1)==.5),:);
lfp_int_3 = lfp_tr(find(behavdec_tr(:,2)==0 & behavdec_tr(:,1)==3),:);
lfp_int_5 = lfp_tr(find(behavdec_tr(:,2)==0 & behavdec_tr(:,1)==5),:);
lfp_int_12 = lfp_tr(find(behavdec_tr(:,2)==0 & behavdec_tr(:,1)==12),:);
lfp_seg_p5 = lfp_tr(find(behavdec_tr(:,2)==1 & behavdec_tr(:,1)==.5),:);
lfp_seg_3 = lfp_tr(find(behavdec_tr(:,2)==1 & behavdec_tr(:,1)==3),:);
lfp_seg_5 = lfp_tr(find(behavdec_tr(:,2)==1 & behavdec_tr(:,1)==5),:);
lfp_seg_12 = lfp_tr(find(behavdec_tr(:,2)==1 & behavdec_tr(:,1)==12),:);
% Split lfp_ht into integration and segregation.
lfp_ht_int_p5 = lfp_ht_tr(find(behavdec_tr(:,2)==0 & behavdec_tr(:,1)==.5),:);
lfp_ht_int_3 = lfp_ht_tr(find(behavdec_tr(:,2)==0 & behavdec_tr(:,1)==3),:);
lfp_ht_int_5 = lfp_ht_tr(find(behavdec_tr(:,2)==0 & behavdec_tr(:,1)==5),:);
lfp_ht_int_12 = lfp_ht_tr(find(behavdec_tr(:,2)==0 & behavdec_tr(:,1)==12),:);
lfp_ht_seg_p5 = lfp_ht_tr(find(behavdec_tr(:,2)==1 & behavdec_tr(:,1)==.5),:);
lfp_ht_seg_3 = lfp_ht_tr(find(behavdec_tr(:,2)==1 & behavdec_tr(:,1)==3),:);
lfp_ht_seg_5 = lfp_ht_tr(find(behavdec_tr(:,2)==1 & behavdec_tr(:,1)==5),:);
lfp_ht_seg_12 = lfp_ht_tr(find(behavdec_tr(:,2)==1 & behavdec_tr(:,1)==12),:);

%% Get SI(similarity indices) using peak amplitudes for a given df using moving time window.
fprintf('\nCalculating similarity indices. Please wait ...\n')
% Convert time to the number of data points.
n_twwidth = round(TMWNDWTH/1000*fs_neural);
n_twmov = round(TMWNDMOV/1000*fs_neural);
num_bin = ceil((n_durmin-n_twwidth)/n_twmov);
num_df = length(DF);
% For moving forward in time
peak_int_p5_fw = zeros(num_int_p5,num_bin);
peak_int_3_fw = zeros(num_int_3,num_bin);
peak_int_5_fw = zeros(num_int_5,num_bin);
peak_int_12_fw = zeros(num_int_12,num_bin);
peak_seg_p5_fw = zeros(num_seg_p5,num_bin);
peak_seg_3_fw = zeros(num_seg_3,num_bin);
peak_seg_5_fw = zeros(num_seg_5,num_bin);
peak_seg_12_fw = zeros(num_seg_12,num_bin);
peak_int_mean_fw = zeros(num_df,num_bin);
peak_seg_mean_fw = zeros(num_df,num_bin);
% For moving backward in time
peak_int_p5_bw = zeros(num_int_p5,num_bin);
peak_int_3_bw = zeros(num_int_3,num_bin);
peak_int_5_bw = zeros(num_int_5,num_bin);
peak_int_12_bw = zeros(num_int_12,num_bin);
peak_seg_p5_bw = zeros(num_seg_p5,num_bin);
peak_seg_3_bw = zeros(num_seg_3,num_bin);
peak_seg_5_bw = zeros(num_seg_5,num_bin);
peak_seg_12_bw = zeros(num_seg_12,num_bin);
peak_int_mean_bw = zeros(num_df,num_bin);
peak_seg_mean_bw = zeros(num_df,num_bin);

% Get peak-based SI by moving forward in time.
si_peak_fw = zeros(num_df,num_bin);
for j=1:num_bin
	for i=1:num_int_p5
		peak_int_p5_fw(i,j) = max(abs(lfp_int_p5(i,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	end
	for i=1:num_int_3
		peak_int_3_fw(i,j) = max(abs(lfp_int_3(i,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	end
	for i=1:num_int_5
		peak_int_5_fw(i,j) = max(abs(lfp_int_5(i,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	end
	for i=1:num_int_12
		peak_int_12_fw(i,j) = max(abs(lfp_int_12(i,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	end
	for i=1:num_seg_p5
		peak_seg_p5_fw(i,j) = max(abs(lfp_seg_p5(i,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	end
	for i=1:num_seg_3
		peak_seg_3_fw(i,j) = max(abs(lfp_seg_3(i,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	end
	for i=1:num_seg_5
		peak_seg_5_fw(i,j) = max(abs(lfp_seg_5(i,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	end
	for i=1:num_seg_12
		peak_seg_12_fw(i,j) = max(abs(lfp_seg_12(i,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	end
	peak_int_mean_fw(1,j) = sum(peak_int_p5_fw(:,j),1)/num_int_p5;
	peak_int_mean_fw(2,j) = sum(peak_int_3_fw(:,j),1)/num_int_3;
	peak_int_mean_fw(3,j) = sum(peak_int_5_fw(:,j),1)/num_int_5;
	peak_int_mean_fw(4,j) = sum(peak_int_12_fw(:,j),1)/num_int_12;
	peak_seg_mean_fw(1,j) = sum(peak_seg_p5_fw(:,j),1)/num_seg_p5;
	peak_seg_mean_fw(2,j) = sum(peak_seg_3_fw(:,j),1)/num_seg_3;
	peak_seg_mean_fw(3,j) = sum(peak_seg_5_fw(:,j),1)/num_seg_5;
	peak_seg_mean_fw(4,j) = sum(peak_seg_12_fw(:,j),1)/num_seg_12;
	for k=1:num_df
		si_peak_fw(k,j) = (peak_int_mean_fw(k,j)-peak_seg_mean_fw(k,j))/(peak_int_mean_fw(k,j)+peak_seg_mean_fw(k,j));
	end
end

% Get peak-based SI by moving backward in time.
si_peak_bw = zeros(num_df,num_bin);
for j=1:num_bin
	for i=1:num_int_p5
		lfp_int_p5_bw(i,:) = circshift(lfp_int_p5(i,:),nnz(~lfp_int_p5(i,:)),2);
		peak_int_p5_bw(i,j) = max(abs(lfp_int_p5_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	for i=1:num_int_3
		lfp_int_3_bw(i,:) = circshift(lfp_int_3(i,:),nnz(~lfp_int_3(i,:)),2);
		peak_int_3_bw(i,j) = max(abs(lfp_int_3_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	for i=1:num_int_5
		lfp_int_5_bw(i,:) = circshift(lfp_int_5(i,:),nnz(~lfp_int_5(i,:)),2);
		peak_int_5_bw(i,j) = max(abs(lfp_int_5_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	for i=1:num_int_12
		lfp_int_12_bw(i,:) = circshift(lfp_int_12(i,:),nnz(~lfp_int_12(i,:)),2);
		peak_int_12_bw(i,j) = max(abs(lfp_int_12_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	for i=1:num_seg_p5
		lfp_seg_p5_bw(i,:) = circshift(lfp_seg_p5(i,:),nnz(~lfp_seg_p5(i,:)),2);
		peak_seg_p5_bw(i,j) = max(abs(lfp_seg_p5_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	for i=1:num_seg_3
		lfp_seg_3_bw(i,:) = circshift(lfp_seg_3(i,:),nnz(~lfp_seg_3(i,:)),2);
		peak_seg_3_bw(i,j) = max(abs(lfp_seg_3_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	for i=1:num_seg_5
		lfp_seg_5_bw(i,:) = circshift(lfp_seg_5(i,:),nnz(~lfp_seg_5(i,:)),2);
		peak_seg_5_bw(i,j) = max(abs(lfp_seg_5_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	for i=1:num_seg_12
		lfp_seg_12_bw(i,:) = circshift(lfp_seg_12(i,:),nnz(~lfp_seg_12(i,:)),2);
		peak_seg_12_bw(i,j) = max(abs(lfp_seg_12_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	peak_int_mean_bw(1,j) = sum(peak_int_p5_bw(:,j),1)/num_int_p5;
	peak_int_mean_bw(2,j) = sum(peak_int_3_bw(:,j),1)/num_int_3;
	peak_int_mean_bw(3,j) = sum(peak_int_5_bw(:,j),1)/num_int_5;
	peak_int_mean_bw(4,j) = sum(peak_int_12_bw(:,j),1)/num_int_12;
	peak_seg_mean_bw(1,j) = sum(peak_seg_p5_bw(:,j),1)/num_seg_p5;
	peak_seg_mean_bw(2,j) = sum(peak_seg_3_bw(:,j),1)/num_seg_3;
	peak_seg_mean_bw(3,j) = sum(peak_seg_5_bw(:,j),1)/num_seg_5;
	peak_seg_mean_bw(4,j) = sum(peak_seg_12_bw(:,j),1)/num_seg_12;
	for k=1:num_df
		si_peak_bw(k,j) = (peak_int_mean_bw(k,j)-peak_seg_mean_bw(k,j))/(peak_int_mean_bw(k,j)+peak_seg_mean_bw(k,j));
	end
end

% Mean and s.d. of SI over all time bins
% Forward in time.
%si_peak_fw(:,ceil(find(isnan(si_peak_fw))/num_df)) = [];
si_peak_fw_mean = mean(si_peak_fw,2);
si_peak_fw_std = std(si_peak_fw,1,2);
% Backward in time.
%si_peak_bw(:,ceil(find(isnan(si_peak_bw))/num_df)) = [];
si_peak_bw_mean = mean(si_peak_bw,2);
si_peak_bw_std = std(si_peak_bw,1,2);

%% Get ITPC(inter-trial phase coherence) using phase.
% For moving forward in time.
phase_int_p5_fw = zeros(num_int_p5,n_twwidth,num_bin);
phase_int_3_fw = zeros(num_int_3,n_twwidth,num_bin);
phase_int_5_fw = zeros(num_int_5,n_twwidth,num_bin);
phase_int_12_fw = zeros(num_int_12,n_twwidth,num_bin);
phase_seg_p5_fw = zeros(num_seg_p5,n_twwidth,num_bin);
phase_seg_3_fw = zeros(num_seg_3,n_twwidth,num_bin);
phase_seg_5_fw = zeros(num_seg_5,n_twwidth,num_bin);
phase_seg_12_fw = zeros(num_seg_12,n_twwidth,num_bin);
itpc_int_mean_fw = zeros(num_df,num_bin);
itpc_seg_mean_fw = zeros(num_df,num_bin);
% For moving backward in time.
phase_int_p5_bw = zeros(num_int_p5,n_twwidth,num_bin);
phase_int_3_bw = zeros(num_int_3,n_twwidth,num_bin);
phase_int_5_bw = zeros(num_int_5,n_twwidth,num_bin);
phase_int_12_bw = zeros(num_int_12,n_twwidth,num_bin);
phase_seg_p5_bw = zeros(num_seg_p5,n_twwidth,num_bin);
phase_seg_3_bw = zeros(num_seg_3,n_twwidth,num_bin);
phase_seg_5_bw = zeros(num_seg_5,n_twwidth,num_bin);
phase_seg_12_bw = zeros(num_seg_12,n_twwidth,num_bin);
itpc_int_mean_bw = zeros(num_df,num_bin);
itpc_seg_mean_bw = zeros(num_df,num_bin);

% Get ITPC-based SI by moving forward in time.
si_itpc_fw = zeros(num_df,num_bin);
for j=1:num_bin
	% For .5 semitone intergration trials.
	phase_int_p5_fw(:,:,j) = atan(imag(lfp_ht_int_p5(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth))./real(lfp_ht_int_p5(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	sum_cos_int_p5_fw = sum(cos(phase_int_p5_fw(:,:,j)),1);
	sum_sin_int_p5_fw = sum(sin(phase_int_p5_fw(:,:,j)),1);
	itpc_int_p5_fw = sqrt(sum_cos_int_p5_fw.^2+sum_sin_int_p5_fw.^2)/num_int_p5;
	itpc_int_mean_fw(1,j) = sum(itpc_int_p5_fw,2)/n_twwidth;
	% For 3 semitone intergration trials.
	phase_int_3_fw(:,:,j) = atan(imag(lfp_ht_int_3(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth))./real(lfp_ht_int_3(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	sum_cos_int_3_fw = sum(cos(phase_int_3_fw(:,:,j)),1);
	sum_sin_int_3_fw = sum(sin(phase_int_3_fw(:,:,j)),1);
	itpc_int_3_fw = sqrt(sum_cos_int_3_fw.^2+sum_sin_int_3_fw.^2)/num_int_3;
	itpc_int_mean_fw(2,j) = sum(itpc_int_3_fw,2)/n_twwidth;
	% For 5 semitone intergration trials.
	phase_int_5_fw(:,:,j) = atan(imag(lfp_ht_int_5(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth))./real(lfp_ht_int_5(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	sum_cos_int_5_fw = sum(cos(phase_int_5_fw(:,:,j)),1);
	sum_sin_int_5_fw = sum(sin(phase_int_5_fw(:,:,j)),1);
	itpc_int_5_fw = sqrt(sum_cos_int_5_fw.^2+sum_sin_int_5_fw.^2)/num_int_5;
	itpc_int_mean_fw(3,j) = sum(itpc_int_5_fw,2)/n_twwidth;
	% For 12 semitone intergration trials.
	phase_int_12_fw(:,:,j) = atan(imag(lfp_ht_int_12(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth))./real(lfp_ht_int_12(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	sum_cos_int_12_fw = sum(cos(phase_int_12_fw(:,:,j)),1);
	sum_sin_int_12_fw = sum(sin(phase_int_12_fw(:,:,j)),1);
	itpc_int_12_fw = sqrt(sum_cos_int_12_fw.^2+sum_sin_int_12_fw.^2)/num_int_12;
	itpc_int_mean_fw(4,j) = sum(itpc_int_12_fw,2)/n_twwidth;
	% For .5 semitone segregation trials.
	phase_seg_p5_fw(:,:,j) = atan(imag(lfp_ht_seg_p5(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth))./real(lfp_ht_seg_p5(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	sum_cos_seg_p5_fw = sum(cos(phase_seg_p5_fw(:,:,j)),1);
	sum_sin_seg_p5_fw = sum(sin(phase_seg_p5_fw(:,:,j)),1);
	itpc_seg_p5_fw = sqrt(sum_cos_seg_p5_fw.^2+sum_sin_seg_p5_fw.^2)/num_seg_p5;
	itpc_seg_mean_fw(1,j) = sum(itpc_seg_p5_fw,2)/n_twwidth;
	% For 3 semitone segregation trials.
	phase_seg_3_fw(:,:,j) = atan(imag(lfp_ht_seg_3(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth))./real(lfp_ht_seg_3(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	sum_cos_seg_3_fw = sum(cos(phase_seg_3_fw(:,:,j)),1);
	sum_sin_seg_3_fw = sum(sin(phase_seg_3_fw(:,:,j)),1);
	itpc_seg_3_fw = sqrt(sum_cos_seg_3_fw.^2+sum_sin_seg_3_fw.^2)/num_seg_3;
	itpc_seg_mean_fw(2,j) = sum(itpc_seg_3_fw,2)/n_twwidth;
	% For 5 semitone segregation trials.
	phase_seg_5_fw(:,:,j) = atan(imag(lfp_ht_seg_5(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth))./real(lfp_ht_seg_5(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	sum_cos_seg_5_fw = sum(cos(phase_seg_5_fw(:,:,j)),1);
	sum_sin_seg_5_fw = sum(sin(phase_seg_5_fw(:,:,j)),1);
	itpc_seg_5_fw = sqrt(sum_cos_seg_5_fw.^2+sum_sin_seg_5_fw.^2)/num_seg_5;
	itpc_seg_mean_fw(3,j) = sum(itpc_seg_5_fw,2)/n_twwidth;
	% For 12 semitone segregation trials.
	phase_seg_12_fw(:,:,j) = atan(imag(lfp_ht_seg_12(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth))./real(lfp_ht_seg_12(:,(j-1)*n_twmov+1:(j-1)*n_twmov+n_twwidth)));
	sum_cos_seg_12_fw = sum(cos(phase_seg_12_fw(:,:,j)),1);
	sum_sin_seg_12_fw = sum(sin(phase_seg_12_fw(:,:,j)),1);
	itpc_seg_12_fw = sqrt(sum_cos_seg_12_fw.^2+sum_sin_seg_12_fw.^2)/num_seg_12;
	itpc_seg_mean_fw(4,j) = sum(itpc_seg_12_fw,2)/n_twwidth;
	% Get SI using temporal means of itpc for integration and segretation.
	for k=1:4
		si_itpc_fw(k,j) = itpc_int_mean_fw(k,j)-itpc_seg_mean_fw(k,j);
	end
end

% Get ITPC-based SI by moving backward in time.
si_itpc_bw = zeros(num_df,num_bin);
for j=1:num_bin
	% For .5 semitone intergration trials.
	for i=1:num_int_p5
		lfp_ht_int_p5_bw(i,:) = circshift(lfp_ht_int_p5(i,:),sum(lfp_ht_int_p5(i,:)==0,2),2);
		phase_int_p5_bw(i,:,j) = atan(imag(lfp_ht_int_p5_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov))./real(lfp_ht_int_p5_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	sum_cos_int_p5_bw = sum(cos(phase_int_p5_bw(:,:,j)),1);
	sum_sin_int_p5_bw = sum(sin(phase_int_p5_bw(:,:,j)),1);
	itpc_int_p5_bw = sqrt(sum_cos_int_p5_bw.^2+sum_sin_int_p5_bw.^2)/num_int_p5;
	itpc_int_mean_bw(1,j) = sum(itpc_int_p5_bw,2)/n_twwidth;
	% For 3 semitone intergration trials.
	for i=1:num_int_3
		lfp_ht_int_3_bw(i,:) = circshift(lfp_ht_int_3(i,:),sum(lfp_ht_int_3(i,:)==0,2),2);
		phase_int_3_bw(i,:,j) = atan(imag(lfp_ht_int_3_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov))./real(lfp_ht_int_3_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	sum_cos_int_3_bw = sum(cos(phase_int_3_bw(:,:,j)),1);
	sum_sin_int_3_bw = sum(sin(phase_int_3_bw(:,:,j)),1);
	itpc_int_3_bw = sqrt(sum_cos_int_3_bw.^2+sum_sin_int_3_bw.^2)/num_int_3;
	itpc_int_mean_bw(2,j) = sum(itpc_int_3_bw,2)/n_twwidth;
	% For 5 semitone intergration trials.
	for i=1:num_int_5
		lfp_ht_int_5_bw(i,:) = circshift(lfp_ht_int_5(i,:),sum(lfp_ht_int_5(i,:)==0,2),2);
		phase_int_5_bw(i,:,j) = atan(imag(lfp_ht_int_5_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov))./real(lfp_ht_int_5_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	sum_cos_int_5_bw = sum(cos(phase_int_5_bw(:,:,j)),1);
	sum_sin_int_5_bw = sum(sin(phase_int_5_bw(:,:,j)),1);
	itpc_int_5_bw = sqrt(sum_cos_int_5_bw.^2+sum_sin_int_5_bw.^2)/num_int_5;
	itpc_int_mean_bw(3,j) = sum(itpc_int_5_bw,2)/n_twwidth;
	% For 12 semitone intergration trials.
	for i=1:num_int_12
		lfp_ht_int_12_bw(i,:) = circshift(lfp_ht_int_12(i,:),sum(lfp_ht_int_12(i,:)==0,2),2);
		phase_int_12_bw(i,:,j) = atan(imag(lfp_ht_int_12_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov))./real(lfp_ht_int_12_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	sum_cos_int_12_bw = sum(cos(phase_int_12_bw(:,:,j)),1);
	sum_sin_int_12_bw = sum(sin(phase_int_12_bw(:,:,j)),1);
	itpc_int_12_bw = sqrt(sum_cos_int_12_bw.^2+sum_sin_int_12_bw.^2)/num_int_12;
	itpc_int_mean_bw(4,j) = sum(itpc_int_12_bw,2)/n_twwidth;
	% For .5 semitone segregation trials.
	for i=1:num_seg_p5
		lfp_ht_seg_p5_bw(i,:) = circshift(lfp_ht_seg_p5(i,:),sum(lfp_ht_seg_p5(i,:)==0,2),2);
		phase_seg_p5_bw(i,:,j) = atan(imag(lfp_ht_seg_p5_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov))./real(lfp_ht_seg_p5_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	sum_cos_seg_p5_bw = sum(cos(phase_seg_p5_bw(:,:,j)),1);
	sum_sin_seg_p5_bw = sum(sin(phase_seg_p5_bw(:,:,j)),1);
	itpc_seg_p5_bw = sqrt(sum_cos_seg_p5_bw.^2+sum_sin_seg_p5_bw.^2)/num_seg_p5;
	itpc_seg_mean_bw(1,j) = sum(itpc_seg_p5_bw,2)/n_twwidth;
	% For 3 semitone segregation trials.
	for i=1:num_seg_3
		lfp_ht_seg_3_bw(i,:) = circshift(lfp_ht_seg_3(i,:),sum(lfp_ht_seg_3(i,:)==0,2),2);
		phase_seg_3_bw(i,:,j) = atan(imag(lfp_ht_seg_3_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov))./real(lfp_ht_seg_3_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	sum_cos_seg_3_bw = sum(cos(phase_seg_3_bw(:,:,j)),1);
	sum_sin_seg_3_bw = sum(sin(phase_seg_3_bw(:,:,j)),1);
	itpc_seg_3_bw = sqrt(sum_cos_seg_3_bw.^2+sum_sin_seg_3_bw.^2)/num_seg_3;
	itpc_seg_mean_bw(2,j) = sum(itpc_seg_3_bw,2)/n_twwidth;
	% For 5 semitone segregation trials.
	for i=1:num_seg_5
		lfp_ht_seg_5_bw(i,:) = circshift(lfp_ht_seg_5(i,:),sum(lfp_ht_seg_5(i,:)==0,2),2);
		phase_seg_5_bw(i,:,j) = atan(imag(lfp_ht_seg_5_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov))./real(lfp_ht_seg_5_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	sum_cos_seg_5_bw = sum(cos(phase_seg_5_bw(:,:,j)),1);
	sum_sin_seg_5_bw = sum(sin(phase_seg_5_bw(:,:,j)),1);
	itpc_seg_5_bw = sqrt(sum_cos_seg_5_bw.^2+sum_sin_seg_5_bw.^2)/num_seg_5;
	itpc_seg_mean_bw(3,j) = sum(itpc_seg_5_bw,2)/n_twwidth;
	% For 12 semitone segregation trials.
	for i=1:num_seg_12
		lfp_ht_seg_12_bw(i,:) = circshift(lfp_ht_seg_12(i,:),sum(lfp_ht_seg_12(i,:)==0,2),2);
		phase_seg_12_bw(i,:,j) = atan(imag(lfp_ht_seg_12_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov))./real(lfp_ht_seg_12_bw(i,end-((j-1)*n_twmov+n_twwidth)+1:end-(j-1)*n_twmov)));
	end
	sum_cos_seg_12_bw = sum(cos(phase_seg_12_bw(:,:,j)),1);
	sum_sin_seg_12_bw = sum(sin(phase_seg_12_bw(:,:,j)),1);
	itpc_seg_12_bw = sqrt(sum_cos_seg_12_bw.^2+sum_sin_seg_12_bw.^2)/num_seg_12;
	itpc_seg_mean_bw(4,j) = sum(itpc_seg_12_bw,2)/n_twwidth;
	% Get SI using temporal means of itpc for integration and segretation.
	for k=1:4
		si_itpc_bw(k,j) = itpc_int_mean_bw(k,j)-itpc_seg_mean_bw(k,j);
	end
end

% Mean and s.d. of SI over all time bins
% Forward in time
%si_itpc_fw(:,ceil(find(isnan(si_itpc_fw))/num_df)) = [];
si_itpc_fw_mean = mean(si_itpc_fw,2);
si_itpc_fw_std = std(si_itpc_fw,1,2);
% Backward in time
%si_itpc_bw(:,ceil(find(isnan(si_itpc_bw))/num_df)) = [];
si_itpc_bw_mean = mean(si_itpc_bw,2);
si_itpc_bw_std = std(si_itpc_bw,1,2);

fprintf('Done successfully!\n\n')
t_proc = toc(t_tic);

%============================ Plot ================================

%% Plot some basic figures for the first 10 trials.
fprintf('Generating plots ...\n\n')
if ~exist('plotrange')
	fprintf('No "plotrange" given. Default plot range (1~10) automatically selected.\n')
	fprintf('Help: plotrange = [min max]\n')
	plotrange = [1 10];
end
% Determine the number of row of subplots.
numrow = ceil((max(plotrange)-min(plotrange)+1)/2);

% Plot LFP, gamma, and alpha waves in an array.

% Define time sequence for conversion of unit of x axis to ms.
t_neural = (1:n_durmax)/fs_neural;
t_neural_ms = (1:n_durmax)/fs_neural*1000;
t_neural_max = max(t_neural);
t_neural_ms_max = max(t_neural_ms);

% Open a new figure window.
f1 = figure;
for i=min(plotrange):min(max(plotrange),num_goodtrials)
    subplot(numrow,6,mod(3*(i-min(plotrange)),numrow*6)+1);
    plot(t_neural_ms,lfp_tr(i,:));
	xlim([0 t_neural_ms_max])
	ylim([-50 50])
	set(gca,'fontsize',6)
	xlabel('[ms]','fontsize',6)
    title(sprintf('LFP of Tr%i',i),'fontsize',8);

	subplot(numrow,6,mod(3*(i-min(plotrange)),numrow*6)+2);
    plot(t_neural_ms,lfp_gamma_tr(i,:));
	xlim([0 t_neural_ms_max])
	ylim([-50 50])
	set(gca,'fontsize',6)
	xlabel('[ms]','fontsize',6)
    title('Gamma Wave','fontsize',8);

	subplot(numrow,6,mod(3*(i-min(plotrange)),numrow*6)+3);
    plot(t_neural_ms,lfp_alpha_tr(i,:));
	xlim([0 t_neural_ms_max])
	ylim([-50 50])
	set(gca,'fontsize',6)
	xlabel('[ms]','fontsize',6)
    title('Alpha Wave','fontsize',8);
end
position_f1 = get(f1,'position');

% Plot Raw, MUA, and Hilbert Transform of MUA in an array.
f2 = figure;
for i=min(plotrange):min(max(plotrange),num_goodtrials)
    subplot(numrow,6,mod(3*(i-min(plotrange)),numrow*6)+1);
    plot(t_neural_ms,nraw_tr(i,:));
    hold on;
    plot(t_neural_ms,mua_tr(i,:),'g');
	xlim([0 t_neural_ms_max])
	ylim([-150 150])
	set(gca,'fontsize',6)
	xlabel('[ms]','fontsize',6)
    legend({'Raw','MUA'},'location','northeast','fontsize',6)
	legend('boxoff')
    title(sprintf('Raw and MUA of Tr%i',i),'fontsize',8);

	subplot(numrow,6,mod(3*(i-min(plotrange)),numrow*6)+2);
    plot(t_neural_ms,abs(mua_ht_tr(i,:)));
	xlim([0 t_neural_ms_max])
	ylim([0 50])
	set(gca,'fontsize',6)
	xlabel('[ms]','fontsize',6)
    title('Hilbert Transform','fontsize',8);

	subplot(numrow,6,mod(3*(i-min(plotrange)),numrow*6)+3);
    periodogram(mua_tr(i,:));
	set(gca,'fontsize',6)
	xlabel('Normalized Frequency','fontsize',6)
	ylabel('Power/Frequency','fontsize',6)
    title('Periodogram','fontsize',8);
end
position_f2 = get(f2,'position');


% Plot the resulting SI.
% Forward in time
TICKS = [1,2,3,4];
TICKLABELS = {'.5','3','5','12'};
f3 = figure;
subplot(2,3,1)
h1 = ribbon2(si_peak_fw');
zlim([-.5 .5])
set(gca,'xtick',TICKS,'xticklabel',TICKLABELS,'fontsize',6)
xlabel('\Delta f','fontsize',6)
ylabel('Bin Number','fontsize',6)
zlabel('SI_{peak,fw}','fontsize',6)
set(h1,'facecolor','interp')
colormap(b2r(-.5,.5))

subplot(2,3,2)
h2 = ribbon2(si_itpc_fw');
zlim([-.5 .5])
set(gca,'xtick',TICKS,'xticklabel',TICKLABELS,'fontsize',6)
xlabel('\Delta f','fontsize',6)
ylabel('Bin Number','fontsize',6)
zlabel('SI_{ITPC,fw}','fontsize',6)
set(h2,'facecolor','interp')
colormap(b2r(-.5,.5))

subplot(2,3,3)
h31 = errorbar(si_peak_fw_mean,si_peak_fw_std);
hold on;
h32 = errorbar(si_itpc_fw_mean,si_itpc_fw_std,'r');
% Reference horizontal line at 0
plot(get(gca,'xlim'),[0 0],'g')
set(gca,'xtick',TICKS,'xticklabel',TICKLABELS,'fontsize',6)
ylim([-.5 .5])
ylabel('SI','fontsize',6)
xlabel('\Delta f','fontsize',6)
legend({'SI_{peak,fw}','SI_{ITPC,fw}'},'location','northeast','fontsize',6)
legend('boxoff')
title('Similarity Indices vs \Delta f','fontsize',8);

% Backward in time
subplot(2,3,4)
h4 = ribbon2(si_peak_bw');
zlim([-.5 .5])
set(gca,'xtick',TICKS,'xticklabel',TICKLABELS,'fontsize',6)
xlabel('\Delta f','fontsize',6)
ylabel('Bin Number','fontsize',6)
zlabel('SI_{peak,bw}','fontsize',6);
set(h4,'facecolor','interp')
colormap(b2r(-.5,.5))

subplot(2,3,5)
h5 = ribbon2(si_itpc_bw');
zlim([-.5 .5])
set(gca,'xtick',TICKS,'xticklabel',TICKLABELS,'fontsize',6)
xlabel('\Delta f','fontsize',6)
ylabel('Bin Number','fontsize',6)
zlabel('SI_{ITPC,bw}','fontsize',6);
set(h5,'facecolor','interp')
colormap(b2r(-.5,.5))

subplot(2,3,6)
h61 = errorbar(si_peak_bw_mean,si_peak_bw_std);
hold on;
h62 = errorbar(si_itpc_bw_mean,si_itpc_bw_std,'r');
plot(get(gca,'xlim'),[0 0],'g') % Reference horizontal line at 0
ylim([-.5 .5])
set(gca,'xtick',TICKS,'xticklabel',TICKLABELS,'fontsize',6)
ylabel('SI','fontsize',6)
xlabel('\Delta f','fontsize',6)
legend({'SI_{peak,bw}','SI_{ITPC,bw}'},'location','northeast','fontsize',6)
legend('boxoff')
title('Similarity Indices vs \Delta f','fontsize',8);
position_f3 = get(f3,'position');

%% Save the result.
% Check if the directory to save the result exists.
if strcmp(saveworkspace,'y') | strcmp(savefig,'y')
    cd ..
    dirlist = strsplit(ls);
    if ~ismember(RESULTDIR,dirlist)
        mkdir(RESULTDIR)
    end
    cd(RESULTDIR)
	% Save the workspace
	if strcmp(saveworkspace,'y')
		fprintf('Saving the result ...\n')
		PROCFILE = sprintf('%s_proc_%i%02i%02i%02i%02i.mat',blockname,t_run(1:5));
		save(PROCFILE)
	end
	% Save figures
	if strcmp(savefig,'y')
		fprintf('\nSaving the figures ...\n')
		FIGLFP = sprintf('%s_LFP_%i%02i%02i%02i%02i.eps',blockname,t_run(1:5));
		FIGMUA = sprintf('%s_MUA_%i%02i%02i%02i%02i.eps',blockname,t_run(1:5));
		FIGSI = sprintf('%s_SI_%i%02i%02i%02i%02i.eps',blockname,t_run(1:5));
		for i=1:3
			if i==1
				set(f1,'position',[1 1 1200 1000],'renderer','painters','paperorientation','landscape','paperpositionmode','auto')
				print(f1,'-depsc','-tiff','-r600',FIGLFP)
				set(f1,'position',position_f1)
			elseif i==2
				set(f2,'position',[1 1 1200 1000],'renderer','painters','paperorientation','landscape','paperpositionmode','auto')
				print(f2,'-depsc','-tiff','-r600',FIGMUA)
				set(f2,'position',position_f2)
			elseif i==3
				set(f3,'position',[1 1 600 400],'renderer','painters','paperorientation','landscape','paperpositionmode','auto')
				print(f3,'-depsc','-tiff','-r600',FIGSI)
				set(f3,'position',position_f3)
			end
		end
	end
	cd ..
	cd(RAWDIR)
end


%% Completion.
fprintf('\nProcessing block %s completed!\n',blockname)
fprintf('\nThank you for your business with us!\n\n',blockname)
% Turn on warning.
warning('on','all');
%warning('query','all');
