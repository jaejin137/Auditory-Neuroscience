%clear
%close all
warning('off')

monkeyName = input('Enter monkey name [Cassius or Miyagi]: ','s');
%clustType = input('Which unit type to process?[(s)ingle, (m)ua, (b)oth, or (c)ustom]: ','s');

dir_result = fullfile('~','STRF','STRFParams');

load(sprintf('lamProf_%s.mat',monkeyName),'lamProf');

if ~exist('dataLog')
	dataLog = struct;
end

numDupDrv = 0;

%for idx_data = 1:length(lamProf)
for idx_data = 2
	fprintf('\n\nProcessing %s ...\n\n',lamProf(idx_data).dataName)
	% Parse data information
	nameRoot = split(lamProf(idx_data).dataName,'.');
	nameTok = split(nameRoot{1},'_');
	sessionDate = nameTok{2};
	driveIDRoot = sprintf('%s_%s',nameTok{3:4});
	dir_target = dir(sprintf('~/kiloSorted_DMR/Mr%s-%s/D*',monkeyName,sessionDate));
	% Sanity check
	numDrv = 0;
	for i = 1:length(dir_target)
		if contains(dir_target(i).name,driveIDRoot)
			driveID = dir_target(i).name;
			numDrv = numDrv+1;
		end
	end
	if numDrv == 0
		fprintf('\n!!!Error: No drives matched in %s!!!\n\n',sessionDate)
		continue
	elseif numDrv > 1
		fprintf('\n!!!Error: Multiple drives matched in %s!!!\n\n',sessionDate)
		numDupDrv = numDupDrv +1;
	end
	% Simple log data
	dataLog(idx_data).monkeyName = monkeyName;
	dataLog(idx_data).sessionDate = sessionDate;
	dataLog(idx_data).driveID = driveID;
	dataLog(idx_data).status = 'progress';
	load(sprintf('coord_ac_%s',monkeyName));

	% Identify area
	for idx_area = 1:size(coord_ac,2)
		for idx_coord = 2:length(coord_ac(:,idx_area))
			if ~isempty(coord_ac{idx_coord,idx_area})
				if coord_ac{idx_coord,idx_area} == lamProf(idx_data).coord
					if contains(coord_ac{1,idx_area},{'A1','R'})
						area = 'core'
					elseif contains(coord_ac{1,idx_area},{'BeltM','BeltL'})
						area = 'belt'
					else
						display('!!! Error: Cannot identify area. !!!')
					end
				end
			end
		end
	end

	%try	
		% Calculate STRF parameters
		calculate_strf_params_v4(monkeyName,sessionDate,driveID)
		% Save the result
		var2save = {'UberSTRF','area'};
		save(fullfile(dir_result,sprintf('STRFParams_%s_%s_%s.mat',monkeyName,sessionDate,driveID)),var2save{:});
		dataLog(idx_data).status = 'success';
	%catch
	%	dataLog(idx_data).status = 'fail';
	%	fprintf('\nError occurred with %s_%s_%s... Skipping...\n',monkeyName,sessionDate,driveID)
	%end
	display('Clearing workspace...')
	clearvars -except monkeyName dir_result lamProf idx_data dataLog
end


