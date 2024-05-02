%%============ Prepare to load multiple session data ==========================
% Get the path to the data folder
%datPath = sprintf('./Abu/Behavior/%s',input('Enter the BEHAVIORALdata file name to process: ','s'));
%datPath = sprintf('./Abu/Neural/%s',input('Enter the NEURAL data folder name to process: ','s'));
%datPath = 'abu_chht_20171116d';
%datPath = '2017-11-16_16-12-13' ;

tag_continue = 0;
while tag_continue == 0
	ifContinue = input('Continue previous process?(y/n) ','s');
	if ifContinue == 'n'
		clear
		datPath = './Abu/Data';
		key_common = input('Enter session keyword: ','s');
		key_except = {};
		tag_except = 1;
		ind_except = 1;
		while tag_except==1
			key_except{ind_except} = input('Any exceptions? (hit Enter if none) ','s');
			if isempty(key_except{ind_except})
				key_except = key_except(~cellfun('isempty',key_except));
				tag_except = 0;
			else
				ind_except = ind_except+1;
			end
		end
		
		ifFilter = input('Filter the data over different rhythms?(y/n) ','s');
		%ifFilter = 'n';
		
		% Get the data file list
		cmd_list = sprintf('ls -d %s/%s*',datPath,key_common);
		[status,list] = system(cmd_list);
		dname_tmp = strsplit(list);
		for i=1:length(dname_tmp)
			dname_tmp_split = strsplit(dname_tmp{i},'/');
			dname{i,1} = dname_tmp_split{end};
		end
		dname = dname(~cellfun(@isempty,dname));
		dname = sort(dname(~contains(dname,key_except)));
		numsess = length(dname);
		% Define global workspaces
		grandNumTr = 0;
		grandChProc = [];
		grandBehavSummary = cell([],10);
		grandChDataBroadDns = cell([],256);
		grandChStatBroad = nan([],2,256);
		grandChDataThetaDns = cell([],256);
		grandChStatTheta = nan([],2,256);
		grandChDataAlphaDns = cell([],256);
		grandChStatAlpha = nan([],2,256);
		grandChDataBetaDns = cell([],256);
		grandChStatBeta = nan([],2,256);
		grandChDataGammaDns = cell([],256);
		grandChStatGamma = nan([],2,256);
		grandChDataHFODns = cell([],256);
		grandChStatHFO = nan([],2,256);
		l_init = 1;
		tag_continue = 1;
	
	elseif ifContinue == 'y'
		tag_skip = 0;
		while tag_skip == 0
			ifSkip = input('Skip to the next session data?(y/n) ','s');
			if ifSkip == 'y'
				l_init = l+1;
				tag_skip = 1;
			elseif ifSkip == 'n'
				l_init = l;
				tag_skip = 1;
			else
				return
			end
		end
		tag_continue = 1;
	else
		display('Wrong answer!')
	end
end

tic

%%=============== Load all the target session data set one by one ==============
for l=l_init:numel(dname)
	[status,list] = system(sprintf('ls %s/%s | grep .mat',datPath,dname{l}));
	list = strsplit(list);
	targetBeh{l} = list{1};
	%%========================== Load the behavioral data.======================
	display('Loading behavioral data ...')
	T = load(sprintf('%s/%s/%s',datPath,dname{l},targetBeh{l}));
	varName = fieldnames(T);
	BehavDatRaw = T.(varName{1});
	clear T
	numtrbeh = length(BehavDatRaw);
	% Select only valid trials.
	validtr = [];
	for i=1:length(BehavDatRaw)
		if ~isempty(BehavDatRaw(i).Error) && BehavDatRaw(i).StimulusOn(2)>0
			validtr(end+1) = i;
		end
	end
	numtr = length(validtr);
	BehavDat = BehavDatRaw(validtr);
	
	clear BehavDatRaw
	
	% Get the parameters for each trial.
	behavSummary = cell(numtr,10);
	for i=1:numtr
		behavSummary(i,2:4) = strread(BehavDat(i).CurrentParam,'%s','delimiter','&');
		behavSummary(i,1:2) = strsplit(behavSummary{i,2},'T');
		behavSummary{i,1} = strcat(behavSummary{i,1},'T');
		if contains(BehavDat(i).CurrentParam,'HT') && BehavDat(i).Error(1) == 0
			behavSummary{i,5} = 'h';
			behavSummary{i,6} = BehavDat(i).LeverAcquired(2) - BehavDat(i).StimulusOn(3);
		elseif contains(BehavDat(i).CurrentParam,'HT') && BehavDat(i).Error(1) == 1
			behavSummary{i,5} = 'm';
		elseif contains(BehavDat(i).CurrentParam,'CRT') && BehavDat(i).Error(1) == 0
			behavSummary{i,5} = 'c';
		elseif BehavDat(i).Error(1) == 2 || BehavDat(i).Error(1) == 999
			behavSummary{i,5} = 'f';
			behavSummary{i,6} = BehavDat(i).Error(2) - BehavDat(i).StimulusOn(2);	% Time to false lever release
		else
			fprintf('*** Cannot identify error type for trial #%i! Exiting! ***\n',i)
			return
		end
	end
	%clear BehavDat

	% indices of each response type(hit,mis,fa,cr) and TNR values
	respLabel = {'h' 'm' 'f' 'c'};
	tnrLabel = num2cell(unique(str2num(cell2mat(behavSummary(:,4)))));
	trind = cell(numel(respLabel),1);
	tnrind = cell(numel(tnrLabel),1);
	for i=1:numel(respLabel)
		trind{i,1} = find(cell2mat(behavSummary(:,5))==respLabel{i});
	end
	for i=1:numel(tnrLabel)
		tnrind{i,1} = find(str2num(cell2mat(behavSummary(:,4)))==tnrLabel{i});
	end
	
	display('Behavioral data successfully loaded.')
	%%======================== End of loading behavioral data =====================
	
	
	%%================================== Loading ADC data =========================
	
	ROI = 1000;	% ROI window in msec.
	ROIEnh = 400;
	ROIOFFSET = 360; % ROI offset in msec.
	fs = 20000.0;
	ups = 1;
	dns = 5;
	splt_pad = 0.;	% Recommend 2.0 sec.
	tStart = 0;
	
	display ('Loading ADC data ...')

	%% Optimize TTL threshold to match the number of trials.
	% Load TTL
	ttl = load_open_ephys_data(sprintf('%s/%s/100_ADC1.continuous',datPath,dname{l}));
	tt = [0:1/fs:(length(ttl)-1)/fs]';
	% initial values
	ttlthr = .06;
	numtrttl = 0;
	delta_ttlthr = .001;
	cnt_ttlopt = 0;
	lim_ttlopt = 100;
	while numtrttl ~= numtrbeh
		cnt_ttlopt = cnt_ttlopt+1;
		if numtrttl>numtrbeh
			ttlthr = ttlthr - delta_ttlthr;
		else
			ttlthr = ttlthr + delta_ttlthr;
		end
		ttlOn = ttl>ttlthr;
		ttlUp = find(diff(ttlOn)>.5);
		ttlDn = find(diff(ttlOn)<-.5);
		ttl_minleng = min(length(ttlUp),length(ttlDn));
		ind_trcut = [ttlUp(1:ttl_minleng) ttlDn(1:ttl_minleng)];
		numtrttl = length(ind_trcut);
		fprintf('%i %.3f %i %i\n',cnt_ttlopt, ttlthr, numtrttl, numtrbeh)
		if cnt_ttlopt > lim_ttlopt
			display('*** Error: Disparity detected in the number of trials. Aborted! ***')
			fprintf('*** No. of trials from Behavior:%i, No. of trials from TTL:%i. ***\n',numtrbeh,numtrttl)
			break
		end
	end
	if numtrttl ~= numtrbeh
		continue
	end
	clear ttlOn
	display('$$$ First sanity check passed! $$$')
	

	% Load lever
	lever = load_open_ephys_data(sprintf('%s/%s/100_ADC2.continuous',datPath,dname{l}));
	leverTrRaw = cell(numtrttl,1);
	for i=1:numtrttl
		leverTrRaw{i} = lever(ind_trcut(i,1):ind_trcut(i,2));
	end
	clear lever
	leverTr = leverTrRaw(validtr);
	clear leverTrRaw
	% indices of lever releases
	leverRel = nan(numtr,1);
	% hit
	for i=1:length(trind{1})
		leverRel(trind{1}(i)) = BehavDat(trind{1}(i)).LeverAcquired(2);
	end
	% FA
	for i=1:length(trind{3})
		leverRel(trind{3}(i)) = BehavDat(trind{3}(i)).Error(2);
	end


	% Load sound stimuli
	stimtmp = load_open_ephys_data(sprintf('%s/%s/100_ADC3.continuous',datPath,dname{l}));

	stimTrRaw = cell(numtrttl,1);
	for i=1:numtrttl
		stimTrRaw{i} = stimtmp(ind_trcut(i,1):ind_trcut(i,2));
	end
	stimTr = stimTrRaw(validtr);
	clear stimTrRaw
	
	% Define absolute indices of sound-on interval for each valid trials.
	%ind_sndon = cell(numtr,1);
	%ind_snd = nan(numtr,2);
	%stimthr = .01;
	%for i = 1:numtr
	%	ind_sndon{i} = find(abs(stimTr{i}-mean(stimTr{i}))>stimthr);
	%	ind_snd(i,:) = [ind_sndon{i}(1) ind_sndon{i}(end) ind_sndon{i}(end)-ind_sndon{i}(1)+1];
	%end
	clear stimTr %ind_sndon
	
	% Mean response time in msec.
	meanRespTime = round(mean(cell2mat(behavSummary(trind{1},6))));

	% Define indices of ROI within each trial for each responses.	
	ind_roi = nan(numtr,3);
	if ~isempty(trind{1})	% Hit
		for i=1:length(trind{1})
			ind_roi(trind{1}(i),1:2) = [(BehavDat(trind{1}(i)).StimulusOn(3)-ROI/2+ROIOFFSET)/1000*fs (BehavDat(trind{1}(i)).StimulusOn(3)+ROI/2+ROIOFFSET)/1000*fs-1];
		end
	end
	if ~isempty(trind{2})	% Miss
		for i=1:length(trind{2})
			ind_roi(trind{2}(i),1:2) = [(BehavDat(trind{2}(i)).StimulusOn(3)-ROI/2+ROIOFFSET)/1000*fs (BehavDat(trind{2}(i)).StimulusOn(3)+ROI/2+ROIOFFSET)/1000*fs-1];
		end
	end
	if ~isempty(trind{3})	% FA
		for i=1:length(trind{3})
			ind_roi(trind{3}(i),1:2) = [(leverRel(trind{3}(i))-meanRespTime-ROI/2+ROIOFFSET)/1000*fs (leverRel(trind{3}(i))-meanRespTime+ROI/2+ROIOFFSET)/1000*fs-1];
		end
	end
	if ~isempty(trind{4})	% CR
		for i=1:length(trind{4})
			ind_roi(trind{4}(i),1:2) = [(BehavDat(trind{4}(i)).TimeStamp(2)-ROI/2+ROIOFFSET)/1000*fs (BehavDat(trind{4}(i)).TimeStamp(2)+ROI/2+ROIOFFSET)/1000*fs-1];
		end
	end
	%if ~isempty(trind{1})	% Hit
	%	ind_roi(trind{1},1:2) = [ind_snd(trind{1},1)+(str2num(cell2mat(behavSummary(trind{1},2)))-ROI/2)/1000*fs ind_snd(trind{1},1)+(str2num(cell2mat(behavSummary(trind{1},2)))+ROI/2)/1000*fs-1];
	%end
	%if ~isempty(trind{2})	% Miss
	%	ind_roi(trind{2},1:2) = [ind_snd(trind{2},1)+(str2num(cell2mat(behavSummary(trind{2},2)))-ROI/2)/1000*fs ind_snd(trind{2},1)+(str2num(cell2mat(behavSummary(trind{2},2)))+ROI/2)/1000*fs-1];
	%end
	%if ~isempty(trind{3})	% FA
	%	ind_roi(trind{3},1:2) = [ind_snd(trind{3},1)+(cell2mat(behavSummary(trind{3},6))-meanRespTime-ROI/2)/1000*fs ind_snd(trind{3},1)+(cell2mat(behavSummary(trind{3},6))-meanRespTime+ROI/2)/1000*fs-1];
	%end
	%if ~isempty(trind{4})	% CR
	%	ind_roi(trind{4},1:2) = [ind_snd(trind{4},1)+(str2num(cell2mat(behavSummary(trind{4},2)))-ROI/2)/1000*fs ind_snd(trind{4},1)+(str2num(cell2mat(behavSummary(trind{4},2)))+ROI/2)/1000*fs-1];
	%end
	ind_roi(:,3) = ind_roi(:,2)-ind_roi(:,1)+1;

	% Define absolute indices of ROI
	absind_roi = nan(numtr,2);
	for i=1:numtr
		absind_roi(i,1) = ind_trcut(validtr(i),1)+ind_roi(i,1);
		absind_roi(i,2) = ind_trcut(validtr(i),1)+ind_roi(i,2);
	end



	% Mark ROIs and plot them
	trMark = zeros(length(stimtmp),1);
	for i=1:numtrttl
		trMark(ind_trcut(i,1):ind_trcut(i,2)) = 0.1;
	end
	%sndMark = zeros(length(stimtmp),1);
	%for i=1:numtr
	%	sndMark(ind_trcut(validtr(i),1)+ind_snd(i,1):ind_trcut(validtr(i),1)+ind_snd(i,2)) = 0.08;
	%end
	roiMark = zeros(length(stimtmp),1);
	for i=1:numtr
		roiMark(ind_trcut(validtr(i),1)+ind_roi(i,1):ind_trcut(validtr(i),1)+ind_roi(i,2)) = 0.06;
	end

	% Plot ROIs
	figure
	plot(tt,stimtmp-mean(stimtmp))
	hold on
	plot(tt,trMark)
	plot(tt,roiMark)
	xlim([60 120])
	title('Raw Stimuli and ROIs')
	drawnow
	clear stimtmp trMark sndMark roiMark tt ttl

	display('ADC data successfully loaded.')
	%%============================= End of loading ADC data =======================
	
	%%======================== Loading neural data ================================
	datPathImp = sprintf('%s/%s',datPath,dname{l});
	rmBadCh = 'y';
	impThrMax = 1.0e5;
	numchGiven = 256;
	
	display('Loading neural data ...')
	%% Get good channels.
	if rmBadCh == 'y'
		goodChNum = OpenE_Ztest(datPathImp,0,impThrMax);
	else
		goodChNum = [1:256];
	end
	numGoodCh = length(goodChNum);
	chproc = goodChNum(1:min(numGoodCh,numchGiven));
	numch = length(chproc);
	
	% Create number of channels by power of 2 padding the number of good channels.
	PN = 0:log2(256);
	numchpow = [2.^PN(2.^PN<numch) numch];
	
	%% Load and preprocess broad data channel by channel
	if ifFilter == 'n'
		% Some numbers for filter
		BUTTORDER = 2; % order of Butterworth filter. Minimum: 2nd order.
		PB_60 = [59 61]; % passband for theta.
		[b_60,a_60] = butter(BUTTORDER/2,PB_60/(fs/2),'stop');
		% Define array for downsampled signal.
		chdata_broad_dns = cell(numtr,256);
		% Define stat arrays.
		chstat_broad = nan(numtr,2,256);
		parfor j=1:256
			if ismember(j,chproc)
				chdata_broad_par = cell(numtr,1);
				chdata_broad_dns_par = cell(numtr,1);
				% Load the good channel data, cut it, and filter 60Hz noise.
				chfname = sprintf('100_CH%i.continuous',j);
				chrawtot = load_open_ephys_data(sprintf('%s/%s/%s',datPath,dname{l},chfname));
				chraw2proc = filtfilt(b_60,a_60,chrawtot);
				chrawtot = [];
	
				% Process over broadband.
				chbroad2proc = zscore(chraw2proc);
				chraw2proc = [];
				% Rectify the broad data to process.
				%chbroad2proc = abs(chbroad2proc);
				% Normalize the broad data to process using mean power.
				%chbroad2proc = chbroad2proc/mean(chbroad2proc.^2);
				% Split individual trials and downsample the splitted broad band data.
				for i=1:numtr
					% Cut ROI
					chdata_broad_par{i,1} = chbroad2proc(absind_roi(i,1):absind_roi(i,2));
					% Find optimal minPeakProminence and select peak
					optPP = 3;
					delPP = .5;
					locs = [];
					while isempty(locs)
						[pks locs] = findpeaks(chdata_broad_par{i,1},'minPeakProminence',optPP);
						optPP = optPP-delPP;
					end
					pkloc = locs(find(pks==max(pks)));
					% Sanity check
					if pkloc-ROIEnh/2/1000*fs < 1 || pkloc+ROIEnh/2/1000*fs > length(chdata_broad_par{i,1})
						pkloc = round(length(chdata_broad_par{i,1})/2);
					end

					chdata_broad_par{i,1} = chdata_broad_par{i,1}(pkloc-ROIEnh/2/1000*fs:pkloc+ROIEnh/2/1000*fs);
					chdata_broad_dns_par{i,1} = resample(chdata_broad_par{i,1},ups,dns);
				end
				chbroad2proc = [];
				chdata_broad_dns(:,j) = chdata_broad_dns_par;
				chdata_broad_dns_par = [];
				% Get mean and variance over each trial for each channel
				trmean_broad = nan(numtr,1);
				trvar_broad = nan(numtr,1);
				for i=1:numtr
					trmean_broad(i) = mean(chdata_broad_par{i});
					trvar_broad(i) = var(chdata_broad_par{i});
				end
				chdata_broad_par = [];
				% Merge mean and variance
				chstat_broad(:,:,j) = [trmean_broad trvar_broad];
			end
		end
	else
		% Some numbers for filter
		BUTTORDER = 2; % order of Butterworth filter. Minimum: 2nd order.
		FC_LFP = 600; % cuttoff frequency for LFP
		PB_60 = [59 61]; % passband for theta.
		PB_THETA = [4 8]; % passband for theta.
		PB_ALPHA = [8 12]; % passband for alpha.
		PB_BETA = [13 30]; % passband for beta.
		PB_GAMMA = [30 90]; % passband for gamma.
		PB_HFO = [90 600]; % passband for HFO.
		%PB_MUA = [500 3000]; % Define passband for MUA.
		%FC_MUA = 600; % Define cuffoff frequency for low pass filter.
		% Filter lfp, theta, alpha, beta, gamma after removing 60Hz noise.
		[b_60,a_60] = butter(BUTTORDER/2,PB_60/(fs/2),'stop');
		%[b_lfp,a_lfp] = butter(BUTTORDER/2,FC_LFP/(fs/2),'low');
		[b_theta,a_theta] = butter(BUTTORDER/2,PB_THETA/(fs/2),'bandpass');
		[b_alpha,a_alpha] = butter(BUTTORDER/2,PB_ALPHA/(fs/2),'bandpass');
		[b_beta,a_beta] = butter(BUTTORDER/2,PB_BETA/(fs/2),'bandpass');
		[b_gamma,a_gamma] = butter(BUTTORDER/2,PB_GAMMA/(fs/2),'bandpass');
		[b_hfo,a_hfo] = butter(BUTTORDER/2,PB_HFO/(fs/2),'bandpass');
		% Define array for downsampled signal.
		chdata_broad_dns = cell(numtr,256);
		chdata_theta_dns = cell(numtr,256);
		chdata_alpha_dns = cell(numtr,256);
		chdata_beta_dns = cell(numtr,256);
		chdata_gamma_dns = cell(numtr,256);
		chdata_hfo_dns = cell(numtr,256);
		% Define stat arrays.
		chstat_broad = nan(numtr,2,256);
		chstat_theta = nan(numtr,2,256);
		chstat_alpha = nan(numtr,2,256);
		chstat_beta = nan(numtr,2,256);
		chstat_gamma = nan(numtr,2,256);
		chstat_hfo = nan(numtr,2,256);
	
		parfor j=1:253
			if ismember(j,chproc)
				% Load the good channel data, cut it, and filter 60Hz noise.
				chfname = sprintf('100_CH%i.continuous',j);
				chrawtot = load_open_ephys_data(sprintf('%s/%s/%s',datPath,dname{l},chfname));
				chraw2proc = filtfilt(b_60,a_60,chrawtot);
				chrawtot = [];

				chdata_broad_par = cell(numtr,1);
				chdata_broad_dns_par = cell(numtr,1);
				chdata_theta_par = cell(numtr,1);
				chdata_theta_dns_par = cell(numtr,1);
				chdata_alpha_par = cell(numtr,1);
				chdata_alpha_dns_par = cell(numtr,1);
				chdata_beta_par = cell(numtr,1);
				chdata_beta_dns_par = cell(numtr,1);
				chdata_gamma_par = cell(numtr,1);
				chdata_gamma_dns_par = cell(numtr,1);
				chdata_hfo_par = cell(numtr,1);
				chdata_hfo_dns_par = cell(numtr,1);
	
				% Process over broadband.
				chbroad2proc = zscore(chraw2proc);
				% Split individual trials and downsample the splitted broad band data.
				for i=1:numtr
					chdata_broad_par{i,1} = chbroad2proc(absind_roi(i,1):absind_roi(i,2));
					chdata_broad_dns_par{i,1} = resample(chdata_broad_par{i,1},ups,dns);
				end
				chbroad2proc = [];
				chdata_broad_dns(:,j) = chdata_broad_dns_par;
				chdata_broad_dns_par = [];
				% Get mean and variance over each trial for each channel
				trmean_broad = nan(numtr,1);
				trvar_broad = nan(numtr,1);
				for i=1:numtr
					trmean_broad(i) = mean(chdata_broad_par{i});
					trvar_broad(i) = var(chdata_broad_par{i});
				end
				chdata_broad_par = [];
				% Merge mean and variance
				chstat_broad(:,:,j) = [trmean_broad trvar_broad];
				
				% Process over theta band.
				chtheta2proc = zscore(filtfilt(b_theta,a_theta,chraw2proc));
				% Split individual trials and downsample the splitted theta band data.
				for i=1:numtr
					chdata_theta_par{i,1} = chtheta2proc(absind_roi(i,1):absind_roi(i,2));
					chdata_theta_dns_par{i,1} = resample(chdata_theta_par{i,1},ups,dns);
				end
				chtheta2proc = [];
				chdata_theta_dns(:,j) = chdata_theta_dns_par;
				chdata_theta_dns_par = [];
				% Get mean and variance over each trial for each channel
				trmean_theta = nan(numtr,1);
				trvar_theta = nan(numtr,1);
				for i=1:numtr
					trmean_theta(i) = mean(chdata_theta_par{i});
					trvar_theta(i) = var(chdata_theta_par{i});
				end
				chdata_theta_par = [];
				% Merge mean and variance
				chstat_theta(:,:,j) = [trmean_theta trvar_theta];
	
				% Process over alpha band.
				chalpha2proc = zscore(filtfilt(b_alpha,a_alpha,chraw2proc));
				% Split individual trials and downsample the splitted alpha band data.
				for i=1:numtr
					chdata_alpha_par{i,1} = chalpha2proc(absind_roi(i,1):absind_roi(i,2));
					chdata_alpha_dns_par{i,1} = resample(chdata_alpha_par{i,1},ups,dns);
				end
				chalpha2proc = [];
				chdata_alpha_dns(:,j) = chdata_alpha_dns_par;
				chdata_alpha_dns_par = [];
				% Get mean and variance over each trial for each channel
				trmean_alpha = nan(numtr,1);
				trvar_alpha = nan(numtr,1);
				for i=1:numtr
					trmean_alpha(i) = mean(chdata_alpha_par{i});
					trvar_alpha(i) = var(chdata_alpha_par{i});
				end
				chdata_alpha_par = [];
				% Merge mean and variance
				chstat_alpha(:,:,j) = [trmean_alpha trvar_alpha];
	
				% Process over beta band.
				chbeta2proc = zscore(filtfilt(b_beta,a_beta,chraw2proc));
				% Split individual trials and downsample the splitted beta band data.
				for i=1:numtr
					chdata_beta_par{i,1} = chbeta2proc(absind_roi(i,1):absind_roi(i,2));
					chdata_beta_dns_par{i,1} = resample(chdata_beta_par{i,1},ups,dns);
				end
				chbeta2proc = [];
				chdata_beta_dns(:,j) = chdata_beta_dns_par;
				chdata_beta_dns_par = [];
				% Get mean and variance over each trial for each channel
				trmean_beta = nan(numtr,1);
				trvar_beta = nan(numtr,1);
				for i=1:numtr
					trmean_beta(i) = mean(chdata_beta_par{i});
					trvar_beta(i) = var(chdata_beta_par{i});
				end
				chdata_beta_par = [];
				% Merge mean and variance
				chstat_beta(:,:,j) = [trmean_beta trvar_beta];
	
				% Process over gamma band.
				chgamma2proc = zscore(filtfilt(b_gamma,a_gamma,chraw2proc));
				% Split individual trials and downsample the splitted gamma band data.
				for i=1:numtr
					chdata_gamma_par{i,1} = chgamma2proc(absind_roi(i,1):absind_roi(i,2));
					chdata_gamma_dns_par{i,1} = resample(chdata_gamma_par{i,1},ups,dns);
				end
				chgamma2proc = [];
				chdata_gamma_dns(:,j) = chdata_gamma_dns_par;
				chdata_gamma_dns_par = [];
				% Get mean and variance over each trial for each channel
				trmean_gamma = nan(numtr,1);
				trvar_gamma = nan(numtr,1);
				for i=1:numtr
					trmean_gamma(i) = mean(chdata_gamma_par{i});
					trvar_gamma(i) = var(chdata_gamma_par{i});
				end
				chdata_gamma_par = [];
				% Merge mean and variance
				chstat_gamma(:,:,j) = [trmean_gamma trvar_gamma];
	
				% Process over hfo band.
				chhfo2proc = zscore(filtfilt(b_hfo,a_hfo,chraw2proc));
				% Split individual trials and downsample the splitted hfo band data.
				for i=1:numtr
					chdata_hfo_par{i,1} = chhfo2proc(absind_roi(i,1):absind_roi(i,2));
					chdata_hfo_dns_par{i,1} = resample(chdata_hfo_par{i,1},ups,dns);
				end
				chhfo2proc = [];
				chdata_hfo_dns(:,j) = chdata_hfo_dns_par;
				chdata_hfo_dns_par = [];
				% Get mean and variance over each trial for each channel
				trmean_hfo = nan(numtr,1);
				trvar_hfo = nan(numtr,1);
				for i=1:numtr
					trmean_hfo(i) = mean(chdata_hfo_par{i});
					trvar_hfo(i) = var(chdata_hfo_par{i});
				end
				chdata_hfo_par = [];
				% Merge mean and variance
				chstat_hfo(:,:,j) = [trmean_hfo trvar_hfo];
				chraw2proc = [];
				
			end
		end
	end
	grandChProc = [grandChProc; chproc'];
	clear chproc
	grandNumTr = grandNumTr+numtr;
	clear numtr
	grandBehavSummary = [grandBehavSummary; behavSummary];
	clear behavSummary
	grandChDataBroadDns = [grandChDataBroadDns; chdata_broad_dns];
	clear chdata_broad_dns
	grandChStatBroad = [grandChStatBroad; chstat_broad];
	clear chstat_broad
	if ifFilter == 'y'
		grandChDataThetaDns = [grandChDataThetaDns; chdata_theta_dns];
		clear chdata_theta_dns
		grandChStatTheta = [grandChStatTheta; chstat_theta];
		clear chstat_theta
		grandChDataAlphaDns = [grandChDataAlphaDns; chdata_alpha_dns];
		clear chdata_alpha_dns
		grandChStatAlpha = [grandChStatAlpha; chstat_alpha];
		clear chstat_alpha
		grandChDataBetaDns = [grandChDataBetaDns; chdata_beta_dns];
		clear chdata_beta_dns
		grandChStatBeta = [grandChStatBeta; chstat_beta];
		clear chstat_beta
		grandChDataGammaDns = [grandChDataGammaDns; chdata_gamma_dns];
		clear chdata_gamma_dns
		grandChStatGamma = [grandChStatGamma; chstat_gamma];
		clear chstat_gamma
		grandChDataHFODns = [grandChDataHFODns; chdata_hfo_dns];
		clear chdata_hfo_dns
		grandChStatHFO = [grandChStatHFO; chstat_hfo];
		clear chstat_hfo
	end
end

chproc = unique(grandChProc);
clear grandChProc
numtr = grandNumTr;
clear grandNumTr
behavSummary = grandBehavSummary;
clear grandBehavSummary
chdata_broad_dns = grandChDataBroadDns;
clear grandChDataBroadDns
chstat_broad = grandChStatBroad;
clear grandStatBroad
if ifFilter == 'y'
	chdata_theta_dns = grandChDataThetaDns;
	clear grandChDataThetaDns
	chstat_theta = grandChStatTheta;
	clear grandStatTheta
	chdata_alpha_dns = grandChDataAlphaDns;
	clear grandChDataAlphaDns
	chstat_alpha = grandChStatAlpha;
	clear grandStatAlpha
	chdata_beta_dns = grandChDataBetaDns;
	clear grandChDataBetaDns
	chstat_beta = grandChStatBeta;
	clear grandStatBeta
	chdata_gamma_dns = grandChDataGammaDns;
	clear grandChDataGammaDns
	chstat_gamma = grandChStatGamma;
	clear grandStatGamma
	chdata_hfo_dns = grandChDataHFODns;
	clear grandChDataHFODns
	chstat_hfo = grandChStatHFO;
	clear grandStatHFO
end



%% Calculate evoked potential
% Average over all channels
for i=1:numtr
	if ifFilter == 'n'
		chdata_erp_broad_dns(i,1) = num2cell(mean(cell2mat(chdata_broad_dns(i,:)),2),1);
	elseif ifFilter == 'y'
		chdata_erp_broad_dns(i,1) = num2cell(mean(cell2mat(chdata_broad_dns(i,:)),2),1);
		chdata_erp_theta_dns(i,1) = num2cell(mean(cell2mat(chdata_theta_dns(i,:)),2),1);
		chdata_erp_alpha_dns(i,1) = num2cell(mean(cell2mat(chdata_alpha_dns(i,:)),2),1);
		chdata_erp_beta_dns(i,1) = num2cell(mean(cell2mat(chdata_beta_dns(i,:)),2),1);
		chdata_erp_gamma_dns(i,1) = num2cell(mean(cell2mat(chdata_gamma_dns(i,:)),2),1);
		chdata_erp_hfo_dns(i,1) = num2cell(mean(cell2mat(chdata_hfo_dns(i,:)),2),1);
	end	
end

display('Neural data successfully loaded and processed.')
%%========================= End of loading neural data ========================



return

%%============================== Decoding neural signal =======================
% Build a table for all classes.
classInd = cell(length(tnrLabel),length(respLabel));
className = cell(length(tnrLabel),length(respLabel));
classNum = 1;
for j=1:length(respLabel)
	for i=1:length(tnrLabel)
		classInd{i,j} = find(strcmp(behavSummary(:,5),respLabel(j)).*strcmp(behavSummary(:,4),num2str(tnrLabel{i})));
		className{i,j} = {strcat(num2str(tnrLabel{i}),respLabel{j})};
		behavSummary(classInd{i,j},7) = {strcat(num2str(tnrLabel{i}),respLabel{j})};
		behavSummary(classInd{i,j},8) = {classNum};
		classNum = classNum+1;
	end
end

% Build a table for all conditions over which to run SVM
% Create unbiased sets of condition indices.
% condition index #1: hit trials across all TNR values.
% condition index #2: miss trials across all TNR values.
% condition index #3: FA trials across all TNR values.
% condition index #4: CR trials across all TNR values.
% condition index #5-(end-1): different TNR values across hit, FA, miss.
% condition index end: collapsed across all rows and columns(ignore both TNR and response type)
condIndUB = cell(length(respLabel)+length(tnrLabel)+1,1);
for i=1:numel(respLabel)
	if sum(~cellfun(@isempty,classInd(:,i))) >= 2
		minoccur = min(cellfun(@length,classInd(~cellfun(@isempty,classInd(:,i)),i)));
		if minoccur >= 2
			condIndUB{i} = reshape(cell2mat(cellfun(@(x) x(1:minoccur),classInd(~cellfun(@isempty,classInd(:,i)),i),'UniformOutput',false)),1,[]);
		end
	end
end
for i=1:length(tnrLabel)
	if sum(~cellfun(@isempty,classInd(i,:))) >= 2
		minoccur = min(cellfun(@length,classInd(i,~cellfun(@isempty,classInd(i,:)))));
		if minoccur >= 2
			condIndUB{length(respLabel)+i} = reshape(cell2mat(cellfun(@(x) x(1:minoccur),classInd(i,~cellfun(@isempty,classInd(i,:))),'UniformOutput',false)),1,[]);
		end
	end
end
% Collapse class indices across all TNR and create condition indices.
for j=1:length(respLabel)
	classIndAllTNR{1,j} = vertcat(classInd{:,j});
end
minoccur = min(cellfun(@length,classIndAllTNR(~cellfun(@isempty,classIndAllTNR))));
if minoccur >= 2
	condIndUB{end} = reshape(cell2mat(cellfun(@(x) x(1:minoccur),classIndAllTNR(~cellfun(@isempty,classIndAllTNR)),'UniformOutput',false)),1,[]);
end
% Remove null cell elements
condIndUB = condIndUB(find(~cellfun(@isempty,condIndUB)));

% Find all classes in each condition.
classUniq = cell(length(condIndUB),1);
for k=1:length(condIndUB)
	if k == length(condIndUB)
		classLabel = 5;
	else
		classLabel = 8;
	end
	classUniq{k} = unique(cell2mat(behavSummary(condIndUB{k},classLabel)));
end


%% For classification using SVM model
%%------------------------------- Bootstrapping -------------------------------

% Get SVM models using increasing number of channels.
optLibSVM = 0;
itermax = 100;
% Redefine channels to bootstrap over.
numch = length(chproc);
numchpow = [2.^PN(2.^PN<numch) numch];
numchbstr = numchpow;


% Generate a table of shuffled channel and trial indices for bootstrap for all conditions.
chpick = cell(itermax,numel(numchbstr),length(condIndUB));
rndTrInd = cell(itermax,numel(numchbstr),length(condIndUB));
for k=1:length(condIndUB)
	for j=1:numel(numchbstr)
		for i=1:itermax
			chpick(i,j,k) = {chproc(randperm(numch,numchbstr(j)))};
			tmpTrInd = condIndUB{k};
			rndTrInd(i,j,k) = {tmpTrInd(randperm(length(condIndUB{k})))};
		end
	end
end


% Define arrays for accuracy
accCV_broad = nan(itermax,numel(numchbstr),length(condIndUB));
accCV_broad_rnd = nan(itermax,numel(numchbstr),length(condIndUB));
if ifFilter == 'y'
	accCV_theta = nan(itermax,numel(numchbstr),length(condIndUB));
	accCV_theta_rnd = nan(itermax,numel(numchbstr),length(condIndUB));
	accCV_alpha = nan(itermax,numel(numchbstr),length(condIndUB));
	accCV_alpha_rnd = nan(itermax,numel(numchbstr),length(condIndUB));
	accCV_beta = nan(itermax,numel(numchbstr),length(condIndUB));
	accCV_beta_rnd = nan(itermax,numel(numchbstr),length(condIndUB));
	accCV_gamma = nan(itermax,numel(numchbstr),length(condIndUB));
	accCV_gamma_rnd = nan(itermax,numel(numchbstr),length(condIndUB));
	accCV_hfo = nan(itermax,numel(numchbstr),length(condIndUB));
	accCV_hfo_rnd = nan(itermax,numel(numchbstr),length(condIndUB));
end

% Set SVM optimization options.
if optLibSVM == 1;
	% Find the best C and gamma for broad band.
	[C_best_broad gamma_best_broad] = parfindcgam(trPar(binInd,classLabel),chstat_broad(binInd,:,chpick{iternum,k}));
	libsvmopt_broad = sprintf('-s 0 -t 0 -c %d -g %d -v 10 -q',C_best_broad,gamma_best_broad);
	if ifFilter == 'y'
		% Find the best C and gamma for theta band.
		[C_best_theta gamma_best_theta] = parfindcgam(trPar(binInd,classLabel),chstat_theta(binInd,:,chpick{iternum,k}));
		libsvmopt_theta = sprintf('-s 0 -t 0 -c %d -v 10 -g %d -q',C_best_theta,gamma_best_theta);
		% Find the best C and gamma for alpha band.
		[C_best_alpha gamma_best_alpha] = parfindcgam(trPar(binInd,classLabel),chstat_alpha(binInd,:,chpick{iternum,k}));
		libsvmopt_alpha = sprintf('-s 0 -t 0 -c %d -v 10 -g %d -q',C_best_alpha,gamma_best_alpha);
		% Find the best C and gamma for beta band.
		[C_best_beta gamma_best_beta] = parfindcgam(trPar(binInd,classLabel),chstat_beta(binInd,:,chpick{iternum,k}));
		%libsvmopt_beta = sprintf('-s 0 -t 0 -c %d -v 10 -g %d -q',C_best_beta,gamma_best_beta);
		% Find the best C and gamma for gamma band.
		[C_best_gamma gamma_best_gamma] = parfindcgam(trPar(binInd,classLabel),chstat_gamma(binInd,:,chpick{iternum,k}));
		libsvmopt_gamma = sprintf('-s 0 -t 0 -c %d -v 10 -g %d -q',C_best_gamma,gamma_best_gamma);
		% Find the best C and gamma for HFO band.
		[C_best_gamma gamma_best_gamma] = parfindcgam(trPar(binInd,classLabel),chstat_hfo(binInd,:,chpick{iternum,k}));
		libsvmopt_hfo = sprintf('-s 0 -t 0 -c %d -v 10 -g %d -q',C_best_gamma,gamma_best_gamma);
	end
else
	% (Use default options for linear kernel.)
	libsvmopt_broad = '-s 0 -t 0 -c 1 -v 10 -q';
	if ifFilter == 'y'
		libsvmopt_theta = '-s 0 -t 0 -c 1 -v 10 -q';
		libsvmopt_alpha = '-s 0 -t 0 -c 1 -v 10 -q';
		libsvmopt_beta = '-s 0 -t 0 -c 1 -v 10 -q';
		libsvmopt_gamma = '-s 0 -t 0 -c 1 -v 10 -q';
		libsvmopt_hfo = '-s 0 -t 0 -c 1 -v 10 -q';
	end
end


% Bootstrap over k channels
display('Starting bootstrap...')
for k=1:length(condIndUB)
	if numel(classUniq{k}) >= 2
		if k == length(condIndUB)
			classLabel = 5;
		else
			classLabel = 8;
		end
		parfor j = 1:numel(numchbstr)
			for iternum = 1:itermax
				accCV_broad(iternum,j,k) = libsvmtrain(double(cell2mat(behavSummary(condIndUB{k},classLabel))),chstat_broad(condIndUB{k},1:2,chpick{iternum,j,k}),libsvmopt_broad);
				accCV_broad_rnd(iternum,j,k) = libsvmtrain(double(cell2mat(behavSummary(rndTrInd{iternum,j,k},classLabel))),chstat_broad(condIndUB{k},1:2,chpick{iternum,j,k}),libsvmopt_broad);
				if ifFilter == 'y'
					accCV_theta(iternum,j,k) = libsvmtrain(double(cell2mat(behavSummary(condIndUB{k},classLabel))),chstat_theta(condIndUB{k},1:2,chpick{iternum,j,k}),libsvmopt_theta);
					accCV_theta_rnd(iternum,j,k) = libsvmtrain(double(cell2mat(behavSummary(rndTrInd{iternum,j,k},classLabel))),chstat_theta(condIndUB{k},1:2,chpick{iternum,j,k}),libsvmopt_theta);
					accCV_alpha(iternum,j,k) = libsvmtrain(double(cell2mat(behavSummary(condIndUB{k},classLabel))),chstat_alpha(condIndUB{k},1:2,chpick{iternum,j,k}),libsvmopt_alpha);
					accCV_alpha_rnd(iternum,j,k) = libsvmtrain(double(cell2mat(behavSummary(rndTrInd{iternum,j,k},classLabel))),chstat_alpha(condIndUB{k},1:2,chpick{iternum,j,k}),libsvmopt_alpha);
					accCV_beta(iternum,j,k) = libsvmtrain(double(cell2mat(behavSummary(condIndUB{k},classLabel))),chstat_beta(condIndUB{k},1:2,chpick{iternum,j,k}),libsvmopt_beta);
					accCV_beta_rnd(iternum,j,k) = libsvmtrain(double(cell2mat(behavSummary(rndTrInd{iternum,j,k},classLabel))),chstat_beta(condIndUB{k},1:2,chpick{iternum,j,k}),libsvmopt_beta);
					accCV_gamma(iternum,j,k) = libsvmtrain(double(cell2mat(behavSummary(condIndUB{k},classLabel))),chstat_gamma(condIndUB{k},1:2,chpick{iternum,j,k}),libsvmopt_gamma);
					accCV_gamma_rnd(iternum,j,k) = libsvmtrain(double(cell2mat(behavSummary(rndTrInd{iternum,j,k},classLabel))),chstat_gamma(condIndUB{k},1:2,chpick{iternum,j,k}),libsvmopt_gamma);
					accCV_hfo(iternum,j,k) = libsvmtrain(double(cell2mat(behavSummary(condIndUB{k},classLabel))),chstat_hfo(condIndUB{k},1:2,chpick{iternum,j,k}),libsvmopt_hfo);
					accCV_hfo_rnd(iternum,j,k) = libsvmtrain(double(cell2mat(behavSummary(rndTrInd{iternum,j,k},classLabel))),chstat_hfo(condIndUB{k},1:2,chpick{iternum,j,k}),libsvmopt_hfo);
				end
			end
		end
	end
end

%% Calculate decoding accuracy for different number of channels.
for k = 1:length(condIndUB)
	if numel(classUniq{k}) >= 2
		for j=1:numel(numchbstr)
			accCV_broad_mean(j,k) = mean(accCV_broad(:,j,k),1);
			accCV_broad_std(j,k) = std(accCV_broad(:,j,k),1);
			accCV_broad_sem(j,k) = std(accCV_broad(:,j,k),1)/sqrt(itermax);
			accCV_broad_rnd_mean(j,k) = mean(accCV_broad_rnd(:,j,k),1);
			accCV_broad_rnd_std(j,k) = std(accCV_broad_rnd(:,j,k),1);
			accCV_broad_rnd_sem(j,k) = std(accCV_broad_rnd(:,j,k),1)/sqrt(itermax);
			if ifFilter == 'y'
				accCV_theta_mean(j,k) = mean(accCV_theta(:,j,k),1);
				accCV_theta_std(j,k) = std(accCV_theta(:,j,k),1);
				accCV_theta_sem(j,k) = std(accCV_theta(:,j,k),1)/sqrt(itermax);
				accCV_theta_rnd_mean(j,k) = mean(accCV_theta_rnd(:,j,k),1);
				accCV_theta_rnd_std(j,k) = std(accCV_theta_rnd(:,j,k),1);
				accCV_theta_rnd_sem(j,k) = std(accCV_theta_rnd(:,j,k),1)/sqrt(itermax);
				accCV_alpha_mean(j,k) = mean(accCV_alpha(:,j,k),1);
				accCV_alpha_std(j,k) = std(accCV_alpha(:,j,k),1);
				accCV_alpha_sem(j,k) = std(accCV_alpha(:,j,k),1)/sqrt(itermax);
				accCV_alpha_rnd_mean(j,k) = mean(accCV_alpha_rnd(:,j,k),1);
				accCV_alpha_rnd_std(j,k) = std(accCV_alpha_rnd(:,j,k),1);
				accCV_alpha_rnd_sem(j,k) = std(accCV_alpha_rnd(:,j,k),1)/sqrt(itermax);
				accCV_beta_mean(j,k) = mean(accCV_beta(:,j,k),1);
				accCV_beta_std(j,k) = std(accCV_beta(:,j,k),1);
				accCV_beta_sem(j,k) = std(accCV_beta(:,j,k),1)/sqrt(itermax);
				accCV_beta_rnd_mean(j,k) = mean(accCV_beta_rnd(:,j,k),1);
				accCV_beta_rnd_std(j,k) = std(accCV_beta_rnd(:,j,k),1);
				accCV_beta_rnd_sem(j,k) = std(accCV_beta_rnd(:,j,k),1)/sqrt(itermax);
				accCV_gamma_mean(j,k) = mean(accCV_gamma(:,j,k),1);
				accCV_gamma_std(j,k) = std(accCV_gamma(:,j,k),1);
				accCV_gamma_sem(j,k) = std(accCV_gamma(:,j,k),1)/sqrt(itermax);
				accCV_gamma_rnd_mean(j,k) = mean(accCV_gamma_rnd(:,j,k),1);
				accCV_gamma_rnd_std(j,k) = std(accCV_gamma_rnd(:,j,k),1);
				accCV_gamma_rnd_sem(j,k) = std(accCV_gamma_rnd(:,j,k),1)/sqrt(itermax);
				accCV_hfo_mean(j,k) = mean(accCV_hfo(:,j,k),1);
				accCV_hfo_std(j,k) = std(accCV_hfo(:,j,k),1);
				accCV_hfo_sem(j,k) = std(accCV_hfo(:,j,k),1)/sqrt(itermax);
				accCV_hfo_rnd_mean(j,k) = mean(accCV_hfo_rnd(:,j,k),1);
				accCV_hfo_rnd_std(j,k) = std(accCV_hfo_rnd(:,j,k),1);
				accCV_hfo_rnd_sem(j,k) = std(accCV_hfo_rnd(:,j,k),1)/sqrt(itermax);
			end
		end
	end
end

display('Bootstrap finished.')


%% Classify each group.
libsvmopt = '-s 0 -t 0 -c 1 -q';
for k=1:length(condIndUB)
	if numel(classUniq{k}) >= 2
		svmModel_broad(k) = libsvmtrain(double(cell2mat(behavSummary(condIndUB{k},classLabel))),chstat_broad(condIndUB{k},1:2,chproc),libsvmopt);
		if ifFilter == 'y'
			svmModel_theta(k) = libsvmtrain(double(cell2mat(behavSummary(condIndUB{k},classLabel))),chstat_theta(condIndUB{k},1:2,chproc),libsvmopt);
			svmModel_alpha(k) = libsvmtrain(double(cell2mat(behavSummary(condIndUB{k},classLabel))),chstat_alpha(condIndUB{k},1:2,chproc),libsvmopt);
			svmModel_beta(k) = libsvmtrain(double(cell2mat(behavSummary(condIndUB{k},classLabel))),chstat_beta(condIndUB{k},1:2,chproc),libsvmopt);
			svmModel_gamma(k) = libsvmtrain(double(cell2mat(behavSummary(condIndUB{k},classLabel))),chstat_gamma(condIndUB{k},1:2,chproc),libsvmopt);
			svmModel_hfo(k) = libsvmtrain(double(cell2mat(behavSummary(condIndUB{k},classLabel))),chstat_hfo(condIndUB{k},1:2,chproc),libsvmopt);
		end
	end
end


if numch == 256
	% Redefine the data consistent with spatical channel distribution.
	LUT = [78 250 249 229 230 186 177 165 119 111 100 56 55 39 38 82;
	       164 245 241 234 238 180 185 171 117 121 107 48 53 47 42 120;
	       174 247 243 232 236 178 175 169 122 113 104 108 59 51 40 112;
	       242 244 239 235 176 181 168 172 114 115 106 110 52 49 43 44;
	       240 251 231 228 184 187 173 170 116 118 101 102 58 57 41 46;
	       70 246 237 233 182 179 183 167 123 103 109 105 50 54 45 204;
	       33 158 66 124 31 36 25 37 199 205 197 202 196 206 200 203;
	       35 1 7 63 3 61 5 29 207 254 217 252 215 209 211 201;
	       62 2 60 6 94 88 98 0 213 140 194 210 198 208 192 212;
	       96 90 92 91 95 93 97 193 20 143 139 138 142 132 144 136;
	       99 89 65 67 125 71 127 248 27 147 148 151 149 141 145 133;
	       64 69 126 160 155 163 157 255 21 75 68 76 72 150 146 154;
	       153 161 159 189 135 191 137 162 26 83 87 81 77 79 73 74;
	       131 128 129 188 130 190 134 156 16 32 86 30 84 24 85 28;
	       216 224 220 226 214 222 218 166 22 11 10 14 8 80 12 34;
	       223 219 227 221 225 195 253 152 18 23 19 9 17 13 15 4];
	LUT = LUT + 1; %change to 1 base instead of 0 base
	% For broad data
	%for k=1:numch
	%	chstat_broad_rearr(:,:,k) = chstat_broad(:,:,LUT(k));
	%end
	chstat_broad_newarr = permute(chstat_broad(:,:,LUT),[3 1 2]);
	clear chstat_broad_rearr
	for j=1:numtr
		chstat_broad_mean_grid(:,:,j) = reshape(chstat_broad_newarr(:,j,1),16,16);
		chstat_broad_var_grid(:,:,j) = reshape(chstat_broad_newarr(:,j,2),16,16);
	end
	clear chstat_broad_newarr
	if ifFilter == 'y'
		% For theta band
		%for k=1:numch
		%	chstat_theta_rearr(:,:,k) = chstat_theta(:,:,LUT(k));
		%end
		chstat_theta_newarr = permute(chstat_theta(:,:,LUT),[3 1 2]);
		clear chstat_theta_rearr
		for j=1:numtr
			chstat_theta_mean_grid(:,:,j) = reshape(chstat_theta_newarr(:,j,1),16,16);
			chstat_theta_var_grid(:,:,j) = reshape(chstat_theta_newarr(:,j,2),16,16);
		end
		clear chstat_theta_newarr
		% For alpha band
		%for k=1:numch
		%	chstat_alpha_rearr(:,:,k) = chstat_alpha(:,:,LUT(k));
		%end
		chstat_alpha_newarr = permute(chstat_alpha(:,:,LUT),[3 1 2]);
		clear chstat_alpha_rearr
		for j=1:numtr
			chstat_alpha_mean_grid(:,:,j) = reshape(chstat_alpha_newarr(:,j,1),16,16);
			chstat_alpha_var_grid(:,:,j) = reshape(chstat_alpha_newarr(:,j,2),16,16);
		end
		clear chstat_alpha_newarr
		% For beta band
		%for k=1:numch
		%	chstat_beta_rearr(:,:,k) = chstat_beta(:,:,LUT(k));
		%end
		chstat_beta_newarr = permute(chstat_beta(:,:,LUT),[3 1 2]);
		clear chstat_beta_rearr
		for j=1:numtr
			chstat_beta_mean_grid(:,:,j) = reshape(chstat_beta_newarr(:,j,1),16,16);
			chstat_beta_var_grid(:,:,j) = reshape(chstat_beta_newarr(:,j,2),16,16);
		end
		clear chstat_beta_newarr
		% For gamma band
		%for k=1:numch
		%	chstat_gamma_rearr(:,:,k) = chstat_gamma(:,:,LUT(k));
		%end
		chstat_gamma_newarr = permute(chstat_gamma(:,:,LUT),[3 1 2]);
		clear chstat_gamma_rearr
		for j=1:numtr
			chstat_gamma_mean_grid(:,:,j) = reshape(chstat_gamma_newarr(:,j,1),16,16);
			chstat_gamma_var_grid(:,:,j) = reshape(chstat_gamma_newarr(:,j,2),16,16);
		end
		clear chstat_gamma_newarr
		% For hfo band
		%for k=1:numch
		%	chstat_hfo_rearr(:,:,k) = chstat_hfo(:,:,LUT(k));
		%end
		chstat_hfo_newarr = permute(chstat_hfo(:,:,LUT),[3 1 2]);
		clear chstat_hfo_rearr
		for j=1:numtr
			chstat_hfo_mean_grid(:,:,j) = reshape(chstat_hfo_newarr(:,j,1),16,16);
			chstat_hfo_var_grid(:,:,j) = reshape(chstat_hfo_newarr(:,j,2),16,16);
		end
		clear chstat_hfo_newarr
	end
end

display('Decoding done.')
%%============================== End of bootstrapping =========================

toc



%%============================== Plotting =====================================
display('Generating figures... Please wait.')

%% Plot impedance measure.
%OpenE_Ztest(datPathImp,1,impThrMax);

%% Plot evoked potential for each trial.
xrange = [-ROI/2000:1/(fs/(dns/ups)):ROI/2000];
xrange = xrange(1:length(chdata_erp_broad_dns{trind{1}(1),1}));
num2plot_erp = 16;
if ifFilter == 'n'
	figure
	suptitle('Broad Band E.R.P.')
	for i=1:num2plot_erp
		subplot(ceil(sqrt(num2plot_erp)),ceil(sqrt(num2plot_erp)),i)
		plot(xrange,chdata_erp_broad_dns{i,1},'linewidth',2)
		xlim([-ROI/2000 ROI/2000])
		ylim([-10 10]);
		title(behavSummary{i,5})
	end
else
	% Broad band
	figure
	suptitle('Broad Band E.R.P.')
	for i=1:num2plot_erp
		subplot(ceil(sqrt(num2plot_erp)),ceil(sqrt(num2plot_erp)),i)
		plot(xrange,chdata_erp_broad_dns{i,1},'linewidth',2)
		hold on
		xlim([-ROI/2000 ROI/2000])
		ylim([-10 10]);
		title(behavSummary{i,5})
	end
	% Theta band
	figure
	suptitle('Theta Band E.R.P.')
	for i=1:num2plot_erp
		subplot(ceil(sqrt(num2plot_erp)),ceil(sqrt(num2plot_erp)),i)
		plot(xrange,chdata_erp_theta_dns{i,1},'linewidth',2)
		hold on
		xlim([-ROI/2000 ROI/2000])
		ylim([-10 10]);
		title(behavSummary{i,5})
	end
	% Alpha band
	figure
	suptitle('Alpha Band E.R.P.')
	for i=1:num2plot_erp
		subplot(ceil(sqrt(num2plot_erp)),ceil(sqrt(num2plot_erp)),i)
		plot(xrange,chdata_erp_alpha_dns{i,1},'linewidth',2)
		hold on
		xlim([-ROI/2000 ROI/2000])
		ylim([-10 10]);
		title(behavSummary{i,5})
	end
	% Beta band
	figure
	suptitle('Beta Band E.R.P.')
	for i=1:num2plot_erp
		subplot(ceil(sqrt(num2plot_erp)),ceil(sqrt(num2plot_erp)),i)
		plot(xrange,chdata_erp_beta_dns{i,1},'linewidth',2)
		hold on
		xlim([-ROI/2000 ROI/2000])
		ylim([-10 10]);
		title(behavSummary{i,5})
	end
	% Gamma band
	figure
	suptitle('Gamma Band E.R.P.')
	for i=1:num2plot_erp
		subplot(ceil(sqrt(num2plot_erp)),ceil(sqrt(num2plot_erp)),i)
		plot(xrange,chdata_erp_gamma_dns{i,1},'linewidth',2)
		hold on
		xlim([-ROI/2000 ROI/2000])
		ylim([-10 10]);
		title(behavSummary{i,5})
	end
	% HFO band
	figure
	suptitle('HFO Band E.R.P.')
	for i=1:num2plot_erp
		subplot(ceil(sqrt(num2plot_erp)),ceil(sqrt(num2plot_erp)),i)
		plot(xrange,chdata_erp_hfo_dns{i,1},'linewidth',2)
		hold on
		xlim([-ROI/2000 ROI/2000])
		ylim([-10 10]);
		title(behavSummary{i,5})
	end
end
%if ifFilter == 'n'
%	for i=1:length(trind{1})
%		plot(xrange,chdata_erp_broad_dns{trind{1}(i),1},'linewidth',2)
%		hold on
%	end
%	xlim([-ROI/2000 ROI/2000])
%	ylim([-10 10]);
%	title('Broad Band E.R.P.','fontsize',8)
%else
%	suptitle('Event Related Potential')
%	% broad band
%	subplot(2,3,1)
%	for i=1:length(trind{1})
%		plot(xrange,chdata_erp_broad_dns{trind{1}(i),1},'linewidth',2)
%		hold on
%	end
%	xlim([-ROI/2000 ROI/2000])
%	%ylim([-10 10]);
%	title('Broad Band','fontsize',8)
%	% theta band
%	subplot(2,3,2)
%	for i=1:length(trind{1})
%		plot(xrange,chdata_erp_theta_dns{trind{1}(i),1},'linewidth',2)
%		hold on
%	end
%	xlim([-ROI/2000 ROI/2000])
%	%ylim([-10 10]);
%	title('Theta Band','fontsize',8)
%	% alpha band
%	subplot(2,3,3)
%	for i=1:length(trind{1})
%		plot(xrange,chdata_erp_alpha_dns{trind{1}(i),1},'linewidth',2)
%		hold on
%	end
%	xlim([-ROI/2000 ROI/2000])
%	%ylim([-10 10]);
%	title('Alpha Band','fontsize',8)
%	% beta band
%	subplot(2,3,4)
%	for i=1:length(trind{1})
%		plot(xrange,chdata_erp_beta_dns{trind{1}(i),1},'linewidth',2)
%		hold on
%	end
%	xlim([-ROI/2000 ROI/2000])
%	%ylim([-10 10]);
%	title('Beta Band','fontsize',8)
%	% gamma band
%	subplot(2,3,5)
%	for i=1:length(trind{1})
%		plot(xrange,chdata_erp_gamma_dns{trind{1}(i),1},'linewidth',2)
%		hold on
%	end
%	xlim([-ROI/2000 ROI/2000])
%	%ylim([-10 10]);
%	title('Gamma Band','fontsize',8)
%	% hfo band
%	subplot(2,3,6)
%	for i=1:length(trind{1})
%		plot(xrange,chdata_erp_hfo_dns{trind{1}(i),1},'linewidth',2)
%		hold on
%	end
%	xlim([-ROI/2000 ROI/2000])
%	%ylim([-10 10]);
%	title('HFO Band','fontsize',8)
%end

% Plot % correct vs numch.
if ifFilter == 'n'
	for j=1:length(condIndUB)
		if numel(classUniq{j}) >= 2
			figure
			% Broad
			%subplot(3,3,j)
			errorbar(accCV_broad_mean(:,j),accCV_broad_sem(:,j),'ko-','LineWidth',2)
			hold on
			errorbar(accCV_broad_rnd_mean(:,j),accCV_broad_rnd_sem(:,j),'g--','LineWidth',2)
			ylim([0 100])
			xlabel('Number of channels')
			ylabel('% Correct')
			title(sprintf('Condition %i',j))
			Xticks = 1:length(numchpow);
			XtickLabels = numchpow;
			set(gca,'Xtick',Xticks,'XtickLabel',XtickLabels)
			legend('Broad','Broad_{btstr}')
			legend boxoff
		end
	end
else
	for j=1:length(condIndUB)
		if numel(classUniq{j}) >= 2
			figure
			suptitle(sprintf('Condition %i',j))
			% Broad
			subplot(2,3,1)
			errorbar(accCV_broad_mean(:,j),accCV_broad_sem(:,j),'k-','LineWidth',2)
			hold on
			errorbar(accCV_broad_rnd_mean(:,j),accCV_broad_rnd_sem(:,j),'g--','LineWidth',2)
			ylim([0 100])
			xlabel('Number of channels')
			ylabel('% Correct')
			Xticks = 1:length(numchpow);
			XtickLabels = numchpow;
			set(gca,'Xtick',Xticks,'XtickLabel',XtickLabels)
			legend('Broad','Broad_{btstr}')
			legend boxoff
			% Theta
			subplot(2,3,2)
			errorbar(accCV_theta_mean(:,j),accCV_theta_sem(:,j),'k-','LineWidth',2)
			hold on
			errorbar(accCV_theta_rnd_mean(:,j),accCV_theta_rnd_sem(:,j),'g--','LineWidth',2)
			ylim([0 100])
			xlabel('Number of channels')
			ylabel('% Correct')
			Xticks = 1:length(numchpow);
			XtickLabels = numchpow;
			set(gca,'Xtick',Xticks,'XtickLabel',XtickLabels)
			legend('Theta','Theta_{btstr}')
			legend boxoff
			% Alpha
			subplot(2,3,3)
			errorbar(accCV_alpha_mean(:,j),accCV_alpha_sem(:,j),'b-','LineWidth',2)
			hold on
			errorbar(accCV_alpha_rnd_mean(:,j),accCV_alpha_rnd_sem(:,j),'c--','LineWidth',2)
			ylim([0 100])
			xlabel('Number of channels')
			ylabel('% Correct')
			Xticks = 1:length(numchpow);
			XtickLabels = numchpow;
			set(gca,'Xtick',Xticks,'XtickLabel',XtickLabels)
			legend('Alpha','Alpha_{btstr}')
			legend boxoff
			% Beta
			subplot(2,3,4)
			errorbar(accCV_beta_mean(:,j),accCV_beta_sem(:,j),'g-','LineWidth',2)
			hold on
			errorbar(accCV_beta_rnd_mean(:,j),accCV_beta_rnd_sem(:,j),'y--','LineWidth',2)
			ylim([0 100])
			xlabel('Number of channels')
			ylabel('% Correct')
			Xticks = 1:length(numchpow);
			XtickLabels = numchpow;
			set(gca,'Xtick',Xticks,'XtickLabel',XtickLabels)
			legend('Beta','Beta_{btstr}')
			legend boxoff
			% Gamma
			subplot(2,3,5)
			errorbar(accCV_gamma_mean(:,j),accCV_gamma_sem(:,j),'r-','LineWidth',2)
			hold on
			errorbar(accCV_gamma_rnd_mean(:,j),accCV_gamma_rnd_sem(:,j),'m--','LineWidth',2)
			ylim([0 100])
			xlabel('Number of channels')
			ylabel('% Correct')
			Xticks = 1:length(numchpow);
			XtickLabels = numchpow;
			set(gca,'Xtick',Xticks,'XtickLabel',XtickLabels)
			legend('Gamma','Gamma_{btstr}')
			legend boxoff
			% HFO
			subplot(2,3,6)
			errorbar(accCV_hfo_mean(:,j),accCV_hfo_sem(:,j),'r-','LineWidth',2)
			hold on
			errorbar(accCV_hfo_rnd_mean(:,j),accCV_hfo_rnd_sem(:,j),'m--','LineWidth',2)
			ylim([0 100])
			xlabel('Number of channels')
			ylabel('% Correct')
			Xticks = 1:length(numchpow);
			XtickLabels = numchpow;
			set(gca,'Xtick',Xticks,'XtickLabel',XtickLabels)
			% Resize figure
			%currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1400,800]);
		end
	end
end

return

% Plot variance vs mean.
if ifFilter == 'n'
	for i=1:length(condIndUB)
		if numel(classUniq{i}) >= 2
			figure
			for k = 1:length(chproc)
				gscatter(chstat_broad(condIndUB{i},1,chproc(k)),chstat_broad(condIndUB{i},2,chproc(k)),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_broad(i).SVs(:,1),svmModel_broad(i).SVs(:,2),'ko')
			%for k = chproc
			%	scatter3(chstat_broad(condIndUB{i},1,k),chstat_broad(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)),str2num(cell2mat(behavSummary(condIndUB{i},4))),cell2mat(behavSummary(condIndUB{i},8)),'filled')
			%	drawnow
			%	hold on
			%end
%			xlim([-3 3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			%zlabel('TNR')
			%ZtickLabels = tnrLabel;
			%set(gca,'ZtickLabel',ZtickLabels)
			title(sprintf('Broad, Condition %i',i))
			legend([className{classUniq{i}}])
			legend boxoff
		end
	end

else

	for i=1:length(condIndUB)
		if numel(classUniq{i}) >= 2
			figure
			suptitle(sprintf('Condition %i',i))
			% For broad data
			subplot(2,3,1)
			for k = chproc
				gscatter(chstat_broad(condIndUB{i},1,k),chstat_broad(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_broad(i).SVs(:,1),svmModel_broad(i).SVs(:,2),'ko')
%			xlim([-3 3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			title('Broad')
			legend([className{classUniq{i}}])
			legend boxoff

			% For theta band
			subplot(2,3,2)
			for k = chproc
				gscatter(chstat_theta(condIndUB{i},1,k),chstat_theta(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_theta(i).SVs(:,1),svmModel_theta(i).SVs(:,2),'ko')
%			xlim([-.3 .3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			title('Theta')
			legend([className{classUniq{i}}])
			legend boxoff

			% For alpha band
			subplot(2,3,3)
			for k = chproc
				gscatter(chstat_alpha(condIndUB{i},1,k),chstat_alpha(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_alpha(i).SVs(:,1),svmModel_alpha(i).SVs(:,2),'ko')
%			xlim([-.3 .3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			title('Alpha')
			legend([className{classUniq{i}}])
			legend boxoff

			% For beta band
			subplot(2,3,4)
			for k = chproc
				gscatter(chstat_beta(condIndUB{i},1,k),chstat_beta(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_beta(i).SVs(:,1),svmModel_beta(i).SVs(:,2),'ko')
%			xlim([-.3 .3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			title('Beta')
			legend([className{classUniq{i}}])
			legend boxoff

			% For gamma band
			subplot(2,3,5)
			for k = chproc
				gscatter(chstat_gamma(condIndUB{i},1,k),chstat_gamma(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_gamma(i).SVs(:,1),svmModel_gamma(i).SVs(:,2),'ko')
%			xlim([-.3 .3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			title('Gamma')
			legend([className{classUniq{i}}])
			legend boxoff

			% For hfo band
			subplot(2,3,6)
			for k = chproc
				gscatter(chstat_hfo(condIndUB{i},1,k),chstat_hfo(condIndUB{i},2,k),cell2mat(behavSummary(condIndUB{i},classLabel)))
				drawnow
				hold on
			end
			plot(svmModel_hfo(i).SVs(:,1),svmModel_hfo(i).SVs(:,2),'ko')
%			xlim([-.3 .3])
%			ylim([0 80])
			xlabel('Mean')
			ylabel('Var')
			title('HFO')
			legend([className{classUniq{i}}])
			legend boxoff
			% Resize figure
			%currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1400,800]);
		end
	end
end


%% Plot mean trial by trial activities over all channel grid
%figure
%subplot(1,2,1)
%h1 = imagesc(chstat_broad_mean_grid(:,:,1));
%set(h1,'Alphadata',~isnan(chstat_broad_mean_grid(:,:,1)))
%caxis([0 1])
%colorbar
%title(sprintf('Mean, Tr#%03i',1),'fontsize',8)
%subplot(1,2,2)
%h2 = imagesc(chstat_broad_var_grid(:,:,2));
%set(h2,'Alphadata',~isnan(chstat_broad_var_grid(:,:,1)))
%caxis([0 1])
%colorbar
%title(sprintf('Var, Tr#%03i',1),'fontsize',8)
%%currPos = get(gcf,'Position');
%set(gcf,'Position',[%currPos(1),currPos(2),700,200]);
%for k=2:2*numtr
%	set(h1,'cdata',chstat_broad_mean_grid(:,:,k));
%	set(h2,'cdata',chstat_broad_var_grid(:,:,k));
%	if mod(k,2)
%		subplot(1,2,1)
%		title(sprintf('Mean, Tr#%i',k),'fontsize',8)
%		subplot(1,2,2)
%		title(sprintf('Var, Tr#%i',k),'fontsize',8)
%	end
%	pause(.5)
%	dbroadnow;
%end


