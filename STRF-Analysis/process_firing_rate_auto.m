function process_firing_rate_auto(binSize,stepSize)

	% Basic information for firing rate results
	% bin size in ms.
	%binSize = 5;
	% step size in ms.
	%stepSize = 1;
	
	Tag = datestr(now,'yymmdd');
	dir_fr = fullfile('~','STRF','FiringRate',sprintf('FR_%i_%i_%s',binSize,stepSize,Tag));
	
	% Basic information for STRF results
	dir_result = fullfile('~','STRF','STRFParams');
	list_result = dir(sprintf('%s/*.mat',dir_result));
	
	for idx_result = 1:length(list_result)
		nameRoot = split(list_result(idx_result).name,'.');
		nameTok = split(nameRoot{1},'_');
		monkeyName = nameTok{2};
		sessionDate = nameTok{3};
		driveID = sprintf('%s_%s_%s',nameTok{4:6});
	
		fprintf('\nProcessing %s_%s\n\n',monkeyName,sessionDate)
	
		close all
		try
			calculate_firing_rate_v3(monkeyName,sessionDate,driveID,binSize,stepSize,Tag)
		catch
			fprintf('\n!!!Error occurred with %s_%s!!! Skipping...\n',monkeyName,sessionDate)
		end
	end


end	% End of function definition

