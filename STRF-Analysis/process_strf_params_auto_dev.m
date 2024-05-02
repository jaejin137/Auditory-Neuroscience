%clear
%close all
warning('off')

monkeyName = input('Enter monkey name [Cassius or Miyagi]: ','s');
clustType = input('Which unit type to process?[(s)ingle, (m)ua, (b)oth, or (a)ll]: ','s');

dir_result = fullfile('~','STRF','RI');

load(sprintf('lamProf_%s.mat',monkeyName),'lamProf')

if ~exist('dataLog')
	dataLog = struct;
end

numDupDrv = 0;

for idx_data = 1:length(lamProf)
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


	if clustType == 's'
		clustTypeName = {'good'};
	elseif clustType == 'm'
		clustTypeName = {'mua'};
	elseif clustType == 'b'
		clustTypeName = {'good','mua'};
	elseif clustType == 'a'
		clustTypeName = {'good','mua','unsorted'};
	else
		fprintf('\n!!!Error: Wrong cluster type given!!!\n')
	end

	fprintf('\nProcessing %s_%s_%s...\n',monkeyName,sessionDate,driveID)
	%pause(3)

	try
		% get RI
		custTag = datestr(now,'_yymmdd');
		get_RI_par(monkeyName,sessionDate,driveID,clustTypeName,custTag);
		dataLog(idx_data).status = 'success';
	catch
		dataLog(idx_data).status = 'fail';
		fprintf('\nError occurred with %s_%s_%s... Skipping...\n',monkeyName,sessionDate,driveID)
	end
end
