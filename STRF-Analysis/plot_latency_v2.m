%Tag = datestr(now,'yymmdd');
dir_out = fullfile('~','STRF',sprintf('FR_%i_%i_%s',binTime_ms,stepTime_ms,Tag));
if ~exist(dir_out,'dir')
	mkdir(dir_out)
end


% Plot firing rate for each cluster.
%numCol = ceil(sqrt(length(clust2proc)));
%numRow = floor(sqrt(length(clust2proc)));
%if numRow*numCol < length(clust2proc)
%	numRow = ceil(sqrt(length(clust2proc)));
%end
numCol = 10;
numRow = ceil(length(clust2proc)/(numCol/2));
t_lim_ms = 100;	% Time limit in ms.
fr_max = 400;

h_fr = figure;
currPos = get(h_fr,'Position'); set(h_fr,'Position',[currPos(1),currPos(2),numCol*150,numRow*100])
for idx_clust = 1:length(clust2proc)
	% Trigger A
	subplot_tight(numRow,numCol,2*(idx_clust-1)+1,[.1,.04])
	ylim_up_def = 50;
	ylim_up = min(fr_max,max(ylim_up_def,max(fr_clust{idx_clust,3})));
	plot(tt_bin_TaT_ms(1:end-1), fr_clust{idx_clust,3})
	hold on
	try
		xline(latency_clust(idx_clust,6),'r:',{sprintf('%.1fms',latency_clust(idx_clust,6))},'LineWidth',2,'FontSize',16);
	end
	yline(fr_bl_stat{idx_clust,2}(1,1),'g','LineWidth',2);
	xlim([0 t_lim_ms])
	ylim([0 ylim_up])
	xlabel('Time [ms]')
	ylabel('Firing Rate')
	title(sprintf('#%i TrigA',clust2proc{idx_clust}))
	drawnow

	% Trigger B
	subplot_tight(numRow,numCol,2*idx_clust,[.1,.04])
	ylim_up_def = 50;
	ylim_up = min(fr_max,max(ylim_up_def,max(fr_clust{idx_clust,4})));
	plot(tt_bin_TaT_ms(1:end-1), fr_clust{idx_clust,4})
	hold on
	try
		xline(latency_clust(idx_clust,8),'r:',{sprintf('%.1fms',latency_clust(idx_clust,8))},'LineWidth',2,'FontSize',16);
	end
	yline(fr_bl_stat{idx_clust,3}(1,1),'g','LineWidth',2);
	xlim([0 t_lim_ms])
	ylim([0 ylim_up])
	xlabel('Time [ms]')
	ylabel('Firing Rate')
	title(sprintf('#%i TrigB',clust2proc{idx_clust}))
	drawnow
end
%sgtitle(sprintf('%s-%s',monkeyName,sessionDate),'FontSize',24)

% Plot significant firing rate
h_fr_sig = figure;
currPos = get(h_fr_sig,'Position'); set(h_fr_sig,'Position',[currPos(1),currPos(2),numCol*150,numRow*100])
for idx_clust = 1:length(clust2proc)
	% For trigger A
	subplot_tight(numRow,numCol,2*(idx_clust-1)+1,[.1,.04])
	ylim_up_def = 50;
	ylim_up = min(fr_max,max(ylim_up_def,max(fr_hi_clust{idx_clust,3})));
	plot(tt_bin_TaT_ms(1:end-1),fr_hi_clust{idx_clust,3},'r.','MarkerSize',10)
	hold on
	plot(tt_bin_TaT_ms(1:end-1),fr_lo_clust{idx_clust,3},'b.','MarkerSize',10)
	try
		xline(latency_clust(idx_clust,6),'r:',{sprintf('%.1fms',latency_clust(idx_clust,6))},'LineWidth',2,'FontSize',16);
	end
	yline(fr_bl_stat{idx_clust,2}(1,1),'g','LineWidth',2);
	xlim([0 t_lim_ms])
	ylim([0 ylim_up])
	xlabel('Time [ms]')
	ylabel('Firing Rate')
	title(sprintf('#%i TrigA',clust2proc{idx_clust}))
	drawnow
	
	% For trigger B
	subplot_tight(numRow,numCol,2*idx_clust,[.1,.04])
	ylim_up_def = 50;
	ylim_up = min(fr_max,max(ylim_up_def,max(fr_hi_clust{idx_clust,4})));
	plot(tt_bin_TaT_ms(1:end-1),fr_hi_clust{idx_clust,4},'r.','MarkerSize',10)
	hold on
	plot(tt_bin_TaT_ms(1:end-1),fr_lo_clust{idx_clust,4},'b.','MarkerSize',10)
	try
		xline(latency_clust(idx_clust,8),'r:',{sprintf('%.1fms',latency_clust(idx_clust,8))},'LineWidth',2,'FontSize',16);
	end
	yline(fr_bl_stat{idx_clust,3}(1,1),'g','LineWidth',2);
	xlim([0 t_lim_ms])
	ylim([0 ylim_up])
	xlabel('Time [ms]')
	ylabel('Firing Rate')
	title(sprintf('#%i TrigB',clust2proc{idx_clust}))
	drawnow
end
%sgtitle(sprintf('%s-%s (meanDelay = %.1fms)',monkeyName,sessionDate,meanDelay),'FontSize',24)

%save(fullfile(dir_out,sprintf('FR_%s_%s.mat',monkeyName,sessionDate)))
saveas(h_fr,fullfile(dir_out,sprintf('FR_%s_%s.png',monkeyName,sessionDate)),'png')
saveas(h_fr_sig,fullfile(dir_out,sprintf('FR_%s_%s_Latency.png',monkeyName,sessionDate)),'png')

