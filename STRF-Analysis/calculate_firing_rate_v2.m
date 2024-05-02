function [] = calculate_firing_rate_v2(monkeyName,sessionDate,driveID,binTime_ms,stepTime_ms)

	warning('off')

	%monkeyName='Cassius';
	%sessionDate='190326';
	%driveID='D2_AC_R1';
	%NBoot=8;
	%NB=100;
	%p_crit=.01;
	
	%clearvars -except Trig*

	% Bin size in time (in milliseconds)
	%binTime_ms = 20;
	% Moving step size in time (in milliseconds)
	%stepTime_ms = 1;
	% Time after each trigger in ms to calculate the firing rate for
	TaT_ms = 2000;

	% Some constants
	Fss = 24414.0625;
	sigFact = 2.0;	% Significance factor
	numHiBin = 3;	% number of consecutive "hi" bins to define onset (firing rate higher than sigFact*SD of mean firing rate).
	
	% Define baseline duration in ms (default=200ms)
	blDur_ms = 2000;

	% Hardwired forbidden response window after onset of trigger in ms.
	forbWindow = 0;

	% Bin array size
	binSize = floor(binTime_ms/1000*Fss);
	% Moving step array size
	stepSize = floor(stepTime_ms/1000*Fss);
	
	Tag = datestr(now,'yymmdd');
	dir_out = sprintf('~/STRF/FR_%i_%i_%s',binTime_ms,stepTime_ms,Tag);

	if ~exist(dir_out,'dir')
		mkdir(dir_out)
	end
	
	
	% Triggers
	trigDir = dir(fullfile('~/STRF','Triggers',sprintf('%s/*%s*.mat',monkeyName,sessionDate)));
	TrigFile = fullfile(trigDir.folder,trigDir.name);
	% Load trigger file and spike times (in seconds) for each cluster
	load(TrigFile)
	TrigA = TrigA;  % For parfor to be albe to access trigger file
	TrigB = TrigB;
	% Convert triggers in seconds.
	TrigA_sec = TrigA'/Fss;
	TrigB_sec = TrigB'/Fss;
	
	% Load spike time data for all clusters.
	spikeTimeFile = fullfile('~/kiloSorted_DMR',sprintf('Mr%s-%s/%s/KS2_7_AC/ClusterInfo/spike_times_all_clust.mat',monkeyName,sessionDate,driveID));
	%load(sprintf('%s_%s/ClusterInfo/spike_times_all_clust.mat',monkeyName,sessionDate))
	load(spikeTimeFile)
	spikeTimesAllClust = spikeTimesAllClust;
	
	
	% Find baseline on and off points
	blOn(1) = TrigA_sec(1) - blDur_ms/1000;	% Baseline onset time in seconds
	blOff(1) = TrigA_sec(1);	% Baseline offset time in seconds
	blOn(2) = TrigB_sec(1) - blDur_ms/1000;	% Baseline onset time in seconds
	blOff(2) = TrigB_sec(1);	% Baseline offset time in seconds
	%blCenter = TrigA_sec(end)+(TrigB_sec(1) - TrigA_sec(end))/2;
	%blOn = blCenter-blDur_ms/1000;
	%blOff = blCenter+blDur_ms/1000;
	bidx_blOn = floor(blOn/(stepTime_ms/1000));
	bidx_blOff = floor(blOff/(stepTime_ms/1000));
	% Define real time and baseline offset time
	tt_bin = [0:stepTime_ms/1000:700]';
	tt_bin_offset_ms = (tt_bin-blOff(1))*1000+.5*stepTime_ms;
	% Define real time and baseline offset time for TaT A
	tt_bin_TaT_ms = [0:stepTime_ms:TaT_ms]';
	
	
	% Find 'good' and 'mua' clusters
	clust_match = cellfun(@(x) contains(x,{'good','mua'}), spikeTimesAllClust(:,3), 'UniformOutput', 0);
	idx_matchingClust = find(cell2mat(clust_match));
	clust2proc = spikeTimesAllClust(find(cell2mat(clust_match)),1);
	
	% Firing rate of each cluster for entire time course
	fr_clust = cell(length(clust2proc),4);
	fr_clust(:,1) = clust2proc;
	% Number of spikes for each cluster for entire time course
	ns_clust = cell(length(clust2proc),4);
	ns_clust(:,1) = clust2proc;
	
	% Baseline firing rate for each cluster
	%fr_bl_clust = nan(length(clust2proc),2);
	%fr_bl_clust(:,1) = cell2mat(clust2proc);
	% Baseline firing rate for each cluster
	fr_bl_stat = cell(length(clust2proc),2);
	fr_bl_stat(:,1) = clust2proc;
	% Significant firing rate for each cluster
	bidx_hi_clust = cell(length(clust2proc),2);
	bidx_hi_clust(:,1) = clust2proc;
	fr_hi_clust = cell(length(clust2proc),2);
	fr_hi_clust(:,1) = clust2proc;
	fr_lo_clust = cell(length(clust2proc),2);
	fr_lo_clust(:,1) = clust2proc;
	
	numClust = length(clust2proc);
	parfor idx_clust = 1:numClust
	%for idx_clust = 1:numClust
		% Spike time in seconds (for cluster #3)
		%fprintf('\nCluster #%i\n\n',spikeTimesAllClust{idx_matchingClust(idx_clust),1});
		st_sec = double(spikeTimesAllClust{idx_matchingClust(idx_clust),2})/Fss;
		
		% Get firing rate for individual bins
		ns_bin = [];
		fr_bin = [];
		ns_bin_A = [];
		fr_bin_A = [];
		ns_bin_B = [];
		fr_bin_B = [];

		% Firing rate for entire time course	
		for bidx_step = 1:floor(max(st_sec)/(stepTime_ms/1000))
			ns_bin(bidx_step,1) = numel(find(st_sec > stepTime_ms*(bidx_step-1)/1000 & st_sec < (stepTime_ms*(bidx_step-1)+binTime_ms)/1000));
			fr_bin(bidx_step,1) = ns_bin(bidx_step,1)/(binTime_ms/1000);
		end
		% Firing rate for TaT
		for bidx_step = 1:floor(TaT_ms/stepTime_ms)
			% For trigger A
			ns_bin_A(bidx_step,1) = numel(find(st_sec > TrigA_sec(1)+stepTime_ms*(bidx_step-1)/1000 & st_sec < TrigA_sec(1)+(stepTime_ms*(bidx_step-1)+binTime_ms)/1000));
			fr_bin_A(bidx_step,1) = ns_bin_A(bidx_step,1)/(binTime_ms/1000);
			% For trigger B
			ns_bin_B(bidx_step,1) = numel(find(st_sec > TrigB_sec(1)+stepTime_ms*(bidx_step-1)/1000 & st_sec < TrigB_sec(1)+(stepTime_ms*(bidx_step-1)+binTime_ms)/1000));
			fr_bin_B(bidx_step,1) = ns_bin_B(bidx_step,1)/(binTime_ms/1000);
		end

		% Average firing rate for baseline period for each trigger
		fr_bl_mean = [];
		fr_bl_std = [];
		fr_bl_mean_A = mean(fr_bin(bidx_blOn(1):bidx_blOff(1)));
		fr_bl_std_A = std(fr_bin(bidx_blOn(1):bidx_blOff(1)));
		fr_bl_mean_B = mean(fr_bin(bidx_blOn(2):bidx_blOff(2)));
		fr_bl_std_B = std(fr_bin(bidx_blOn(2):bidx_blOff(2)));
		% Significant firing rate
		% For entire time course
		fr_hi = nan(length(fr_bin),1);
		fr_lo = nan(length(fr_bin),1);
		bidx_fr_hi = find(fr_bin > sigFact*fr_bl_std_A);
		bidx_fr_lo = find(fr_bin <= sigFact*fr_bl_std_A);
		fr_hi(bidx_fr_hi) = fr_bin(bidx_fr_hi);
		fr_lo(bidx_fr_lo) = fr_bin(bidx_fr_lo);
		% For each trigger
		fr_hi_A = nan(length(fr_bin_A),1);
		fr_lo_A = nan(length(fr_bin_A),1);
		fr_hi_B = nan(length(fr_bin_B),1);
		fr_lo_B = nan(length(fr_bin_B),1);
		bidx_fr_hi_A = find(fr_bin_A > sigFact*fr_bl_std_A);
		bidx_fr_lo_A = find(fr_bin_A <= sigFact*fr_bl_std_A);
		bidx_fr_hi_B = find(fr_bin_B > sigFact*fr_bl_std_B);
		bidx_fr_lo_B = find(fr_bin_B <= sigFact*fr_bl_std_B);
		fr_hi_A(bidx_fr_hi_A) = fr_bin_A(bidx_fr_hi_A);
		fr_lo_A(bidx_fr_lo_A) = fr_bin_A(bidx_fr_lo_A);
		fr_hi_B(bidx_fr_hi_B) = fr_bin_B(bidx_fr_hi_B);
		fr_lo_B(bidx_fr_lo_B) = fr_bin_B(bidx_fr_lo_B);
		% Assign significant firing rate for entire time
		bidx_hi_clust_tmp(idx_clust,1) = {bidx_fr_hi};
		fr_hi_clust_tmp(idx_clust,1) = {fr_hi};
		fr_lo_clust_tmp(idx_clust,1) = {fr_lo};
		% Assign significant firing rate for TaT A
		bidx_hi_clust_A_tmp(idx_clust,1) = {bidx_fr_hi_A};
		fr_hi_clust_A_tmp(idx_clust,1) = {fr_hi_A};
		fr_lo_clust_A_tmp(idx_clust,1) = {fr_lo_A};
		% Assign significant firing rate for TaT B
		bidx_hi_clust_B_tmp(idx_clust,1) = {bidx_fr_hi_B};
		fr_hi_clust_B_tmp(idx_clust,1) = {fr_hi_B};
		fr_lo_clust_B_tmp(idx_clust,1) = {fr_lo_B};

		ns_clust_tmp(idx_clust,1) = {ns_bin};
		fr_clust_tmp(idx_clust,1) = {fr_bin};
		ns_clust_A_tmp(idx_clust,1) = {ns_bin_A};
		fr_clust_A_tmp(idx_clust,1) = {fr_bin_A};
		ns_clust_B_tmp(idx_clust,1) = {ns_bin_B};
		fr_clust_B_tmp(idx_clust,1) = {fr_bin_B};
		fr_bl_stat_A_tmp(idx_clust,1) = {[fr_bl_mean_A fr_bl_std_A]};	
		fr_bl_stat_B_tmp(idx_clust,1) = {[fr_bl_mean_B fr_bl_std_B]};
	end

	ns_clust(:,2:4) = [ns_clust_tmp ns_clust_A_tmp ns_clust_B_tmp];
	fr_clust(:,2:4) = [fr_clust_tmp fr_clust_A_tmp fr_clust_B_tmp];
	fr_bl_stat(:,2:3) = [fr_bl_stat_A_tmp fr_bl_stat_B_tmp];
	bidx_hi_clust(:,2:4)= [bidx_hi_clust_tmp bidx_hi_clust_A_tmp bidx_hi_clust_B_tmp];
	fr_hi_clust(:,2:4)= [fr_hi_clust_tmp fr_hi_clust_A_tmp fr_hi_clust_B_tmp];
	fr_lo_clust(:,2:4)= [fr_lo_clust_tmp fr_lo_clust_A_tmp fr_lo_clust_B_tmp];

	
	% Find idices of onset (i.e. the first bin of first N consecutive hi bins) for all clusters.
	lim_onset = 0.1;	% Limit on onset time in seconds after the onset of trigger
	latency_clust = nan(length(clust2proc),8);
	latency_clust(:,1) = cell2mat(clust2proc);
	for idx_clust = 1:numClust
		% For trigger A (using entire firing rate)
		diff_bidx_hi = diff(bidx_hi_clust{idx_clust,2});
		for i = 1:length(diff_bidx_hi)-(numHiBin-1)
			if isequal(diff_bidx_hi(i:i+numHiBin-1),[1;1;1]) && bidx_hi_clust{idx_clust,2}(i)>bidx_blOff(1)+floor(forbWindow/stepTime_ms)+1 && bidx_hi_clust{idx_clust,2}(i)<bidx_blOff(1)+lim_onset*1000/stepTime_ms
				% Index of delayed onset
				latency_clust(idx_clust,2) = bidx_hi_clust{idx_clust,2}(i);
				% Time of delayed onset in seconds
				latency_clust(idx_clust,3) = (bidx_hi_clust{idx_clust,2}(i)-1+.5)*stepTime_ms/1000;
				% Onset delay in ms
				latency_clust(idx_clust,4) = (bidx_hi_clust{idx_clust,2}(i)-1+.5)*stepTime_ms-blOff(1)*1000;
				break
			end
		end
		% For trigger A (using trimmed firing rate for trigger A)
		diff_bidx_hi_A = diff(bidx_hi_clust{idx_clust,3});
		for i = 1:length(diff_bidx_hi_A)-(numHiBin-1)
			if isequal(diff_bidx_hi_A(i:i+numHiBin-1),[1;1;1]) && bidx_hi_clust{idx_clust,3}(i)>floor(forbWindow/stepTime_ms)+1 && bidx_hi_clust{idx_clust,3}(i)<lim_onset*1000/stepTime_ms
				% Index of delayed onset
				latency_clust(idx_clust,5) = bidx_hi_clust{idx_clust,3}(i);
				% Onset delay in ms
				latency_clust(idx_clust,6) = (bidx_hi_clust{idx_clust,3}(i)-1+.5)*stepTime_ms;
				break
			end
		end
		% For trigger B (using trimmed firing rate for trigger B)
		diff_bidx_hi_B = diff(bidx_hi_clust{idx_clust,4});
		for i = 1:length(diff_bidx_hi_B)-(numHiBin-1)
			if isequal(diff_bidx_hi_B(i:i+numHiBin-1),[1;1;1]) && bidx_hi_clust{idx_clust,4}(i)>floor(forbWindow/stepTime_ms)+1 && bidx_hi_clust{idx_clust,4}(i)<lim_onset*1000/stepTime_ms
				% Index of delayed onset
				latency_clust(idx_clust,7) = bidx_hi_clust{idx_clust,4}(i);
				% Onset delay in ms
				latency_clust(idx_clust,8) = (bidx_hi_clust{idx_clust,4}(i)-1+.5)*stepTime_ms;
				break
			end
		end
	end

	meanDelay = nanmean(latency_clust(:,4));
	meanDelay_A = nanmean(latency_clust(:,6));
	meanDelay_B = nanmean(latency_clust(:,8));
	fprintf('\nMean delay = %.1fms\n\n',meanDelay);
	fprintf('\nMean delay = %.1fms for trigger A\n\n',meanDelay_A);
	fprintf('\nMean delay = %.1fms for trigger B\n\n',meanDelay_B);

	save(fullfile(dir_out,sprintf('FR_%s_%s_%i_%i.mat',monkeyName,sessionDate,binTime_ms,stepTime_ms)))

	plot_latency_v2

end	% End of function definition
