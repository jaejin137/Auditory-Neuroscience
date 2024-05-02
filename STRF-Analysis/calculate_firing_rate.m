function [] = calculate_firing_rate(monkeyName,sessionDate,driveID,binTime_ms,stepTime_ms)

	warning('off')

	%monkeyName='Cassius';
	%sessionDate='190326';
	%driveID='D2_AC_R1';
	%NBoot=8;
	%NB=100;
	%p_crit=.01;
	
	%clearvars -except Trig*

	% Bin size in time (in milliseconds)
	%binTime_ms = 5;
	% Moving step size in time (in milliseconds)
	%stepTime_ms = 1;

	% Some constants
	Fss = 24414.0625;
	sigFact = 2.0;	% Significance factor
	numHiBin = 3;	% number of consecutive "hi" bins to define onset (firing rate higher than sigFact*SD of mean firing rate).
	
	% Define baseline duration in ms
	blDur_ms = 200;

	% Hardwired forbidden response window after onset of trigger in ms.
	forbWindow = 0;

	% Bin array size
	binSize = floor(binTime_ms/1000*Fss);
	% Moving step array size
	stepSize = floor(stepTime_ms/1000*Fss);
	
	dir_out = sprintf('~/STRF/FR_%i_%i',binTime_ms,stepTime_ms);

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
	blOn = TrigA_sec(1) - blDur_ms/1000;	% Baseline onset time in seconds
	blOff = TrigA_sec(1);	% Baseline offset time in seconds
	%blCenter = TrigA_sec(end)+(TrigB_sec(1) - TrigA_sec(end))/2;
	%blOn = blCenter-blDur_ms/1000;
	%blOff = blCenter+blDur_ms/1000;
	bidx_blOn = floor(blOn/(stepTime_ms/1000));
	bidx_blOff = floor(blOff/(stepTime_ms/1000));
	% Define real time and baseline offset time
	tt_bin = [0:stepTime_ms/1000:700]';
	tt_bin_offset_ms = (tt_bin-blOff)*1000+.5*stepTime_ms;
	
	
	% Find 'good' and 'mua' clusters
	clust_match = cellfun(@(x) contains(x,{'good','mua'}), spikeTimesAllClust(:,3), 'UniformOutput', 0);
	idx_matchingClust = find(cell2mat(clust_match));
	clust2proc = spikeTimesAllClust(find(cell2mat(clust_match)),1);
	
	% Firing rate of each cluster
	fr_clust = cell(length(clust2proc),2);
	fr_clust(:,1) = clust2proc;
	% Number of spikes for each cluster
	ns_clust = cell(length(clust2proc),2);
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
		% Spike time in seconds (for cluster #3)
		st_sec = double(spikeTimesAllClust{idx_matchingClust(idx_clust),2})/Fss;
		
		% Get firing rate for individual bins
		fr_bin = [];
		ns_bin = [];
		
		%for bidx_step = 1:floor(max(st_sec{idx_clust,2})/(stepTime_ms/1000))
		for bidx_step = 1:floor(max(st_sec)/(stepTime_ms/1000))
			ns_bin(bidx_step,1) = numel(find(st_sec > stepTime_ms*(bidx_step-1)/1000 & st_sec < (stepTime_ms*(bidx_step-1)+binTime_ms)/1000));
			fr_bin(bidx_step,1) = numel(find(st_sec > stepTime_ms*(bidx_step-1)/1000 & st_sec < (stepTime_ms*(bidx_step-1)+binTime_ms)/1000))/(binTime_ms/1000);
		end
		fr_clust(idx_clust,2) = {fr_bin};
		ns_clust(idx_clust,2) = {ns_bin};
		
		% Average firing rate
		fr_bl_mean = [];
		fr_bl_std = [];
		%fr_bl = numel(find(st_sec > blOn & st_sec < blOff))/(blDur_ms/1000);
		%fr_bl_clust(idx_clust,2) = fr_bl;
		fr_bl_mean = mean(fr_bin(bidx_blOn:bidx_blOff));
		fr_bl_std = std(fr_bin(bidx_blOn:bidx_blOff));
		fr_bl_stat(idx_clust,2) = {[fr_bl_mean fr_bl_std]};
		%fr_bl_stat(idx_clust,3) = fr_bl_std;
	
		% Significant firing rate
		fr_hi = nan(length(fr_bin),1);
		fr_lo = nan(length(fr_bin),1);
		bidx_fr_hi = find(fr_bin > sigFact*fr_bl_std);
		bidx_fr_lo = find(fr_bin <= sigFact*fr_bl_std);
		fr_hi(bidx_fr_hi) = fr_bin(bidx_fr_hi);
		fr_lo(bidx_fr_lo) = fr_bin(bidx_fr_lo);
		% Assign significant firing rate
		bidx_hi_clust(idx_clust,2) = {bidx_fr_hi};
		fr_hi_clust(idx_clust,2) = {fr_hi};
		fr_lo_clust(idx_clust,2) = {fr_lo};
	end
	
	% Find idices of onset (i.e. the first bin of first N consecutive hi bins) for all clusters.
	lim_onset = 0.1;	% Limit on onset time in seconds after the onset of trigger
	onset_clust = nan(length(clust2proc),4);
	onset_clust(:,1) = cell2mat(clust2proc);
	for idx_clust = 1:numClust
		diff_bidx_hi = diff(bidx_hi_clust{idx_clust,2});
		for i = 1:length(diff_bidx_hi)-(numHiBin-1)
			if isequal(diff_bidx_hi(i:i+numHiBin-1),[1;1;1]) && bidx_hi_clust{idx_clust,2}(i)>bidx_blOff+floor(forbWindow/stepTime_ms)+1 && bidx_hi_clust{idx_clust,2}(i)<bidx_blOff+lim_onset*1000/stepTime_ms
				onset_clust(idx_clust,2) = bidx_hi_clust{idx_clust,2}(i);
				onset_clust(idx_clust,3) = (bidx_hi_clust{idx_clust,2}(i)-1+.5)*stepTime_ms/1000;
				onset_clust(idx_clust,4) = (bidx_hi_clust{idx_clust,2}(i)-1+.5)*stepTime_ms-blOff*1000;
				break
			end
		end
	end

	meanDelay = nanmean(onset_clust(:,4));
	fprintf('\nMean delay = %.1fms\n\n',meanDelay);

	
	% Plot firing rate for each cluster.
	%numCol = ceil(sqrt(length(clust2proc)));
	%numRow = floor(sqrt(length(clust2proc)));
	%if numRow*numCol < length(clust2proc)
	%	numRow = ceil(sqrt(length(clust2proc)));
	%end
	numCol = 5;
	numRow = ceil(length(clust2proc)/numCol);
	t_lim_ms = 150;	% Time limit in ms.
	fr_max = 400;
	
	h_fr = figure;
	currPos = get(h_fr,'Position'); set(h_fr,'Position',[currPos(1),currPos(2),numCol*200,numRow*150])
	sgtitle(sprintf('%s-%s',monkeyName,sessionDate),'FontSize',20)
	for idx_clust = 1:length(clust2proc)
		subplot(numRow,numCol,idx_clust)
		ylim_up_def = 50;
		ylim_up = min(fr_max,max(ylim_up_def,max(fr_clust{idx_clust,2})));
		%if ylim_up == ylim_up_def
			plot(tt_bin_offset_ms(1:length(fr_clust{idx_clust,2})),fr_clust{idx_clust,2})
			%plot(tt_bin(1:length(fr_clust{idx_clust,2})),fr_clust{idx_clust,2},'.-','MarkerSize',20)
		%else
			%plot(tt_bin(1:length(fr_clust{idx_clust,2})),fr_clust{idx_clust,2},'r.-','MarkerSize',20)
		%end
		hold on
		xline(0,'y','LineWidth',2)
		try
			xline(onset_clust(idx_clust,4),'r:',{sprintf('%.1fms',onset_clust(idx_clust,4))},'LineWidth',2,'FontSize',16);
		end
		%yline(fr_bl_clust(idx_clust,2),'g','LineWidth',2);
		yline(fr_bl_stat{idx_clust,2}(1,1),'g','LineWidth',2);
		%xlabel('Bin No.')
		xlabel('Time [ms]')
		ylabel('Firing Rate')
		%xlim([tt_bin_offset_ms(1) tt_bin_offset_ms(length(fr_clust{idx_clust,2}))])
		xlim([-200 200])
		ylim([0 ylim_up])
		drawnow
	end
	
	% Plot significant firing rate
	h_fr_sig = figure;
	currPos = get(h_fr_sig,'Position'); set(h_fr_sig,'Position',[currPos(1),currPos(2),numCol*200,numRow*150])
	sgtitle(sprintf('%s-%s (meanDelay = %.1fms)',monkeyName,sessionDate,meanDelay),'FontSize',20)
	for idx_clust = 1:length(clust2proc)
		subplot(numRow,numCol,idx_clust)
		ylim_up_def = 50;
		ylim_up = min(fr_max,max(ylim_up_def,max(fr_hi_clust{idx_clust,2})));
		plot(tt_bin_offset_ms(1:length(fr_clust{idx_clust,2})),fr_hi_clust{idx_clust,2},'r.','MarkerSize',10)
		hold on
		plot(tt_bin_offset_ms(1:length(fr_clust{idx_clust,2})),fr_lo_clust{idx_clust,2},'b.','MarkerSize',10)
		xline(0,'y','LineWidth',2)
		try
			xline(onset_clust(idx_clust,4),'r:',{sprintf('%.1fms',onset_clust(idx_clust,4))},'LineWidth',2,'FontSize',16);
		end
		%yline(fr_bl_clust(idx_clust,2),'g','LineWidth',2);
		yline(fr_bl_stat{idx_clust,2}(1,1),'g','LineWidth',2);
		%xlabel('Bin No.')
		xlabel('Time [ms]')
		ylabel('Firing Rate')
		%xlim([blOn blOff+.1])
		%xlim([blOff-.100 blOff+.150])
		xlim([0 t_lim_ms])
		ylim([0 ylim_up])
		drawnow
	end
	
	%% Plot the number of spikes per bin
	%h_ns = figure;
	%for idx_clust = 1:length(clust2proc)
	%	subplot(numRow,numCol,idx_clust)
	%	plot(tt_bin_offset_ms(1:length(fr_clust{idx_clust,2})),ns_clust{idx_clust,2})
	%	%plot(tt_bin_offset_ms(1:length(fr_clust{idx_clust,2})),ns_clust{idx_clust,2},'.-','MarkerSize',20)
	%	%xlabel('Bin No.')
	%	xlabel('Time [ms]')
	%	ylabel('Spikes / Bin')
	%	xlim([tt_bin_offset_ms(1) tt_bin_offset_ms(length(fr_clust{idx_clust,2}))])
	%end
	%currPos = get(h_ns,'Position'); set(h_ns,'Position',[currPos(1),currPos(2),numCol*200,numRow*150])
	
	% Filename tag
	%Tag = datestr(now,'mmddHHMM');

	%save(sprintf('~/STRF/FiringRate/FR_%s_%s_%s.mat',monkeyName,sessionDate,Tag))
	%saveas(h_fr,sprintf('~/STRF/FiringRate/FR_%s_%s_%s.png',monkeyName,sessionDate,Tag),'png')
	%saveas(h_fr_sig,sprintf('~/STRF/FiringRate/FRSig_%s_%s_%s.png',monkeyName,sessionDate,Tag),'png')
	save(fullfile(dir_out,sprintf('FR_%s_%s.mat',monkeyName,sessionDate)))
	saveas(h_fr,fullfile(dir_out,sprintf('FR_%s_%s_Cnt.png',monkeyName,sessionDate)),'png')
	saveas(h_fr_sig,fullfile(dir_out,sprintf('FR_%s_%s.png',monkeyName,sessionDate)),'png')


end	% End of function definition
