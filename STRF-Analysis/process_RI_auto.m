clear
close all
warning('off')

monkeyName = input('Enter monkey name [Cassius or Miyagi]: ','s');
clustType = input('Which unit type to process?[(s)ingle, (m)ua, (b)oth, or (a)ll]: ','s');


dir_result = fullfile('~','STRF','RI');

dir_target = dir(sprintf('~/kiloSorted_DMR/Mr%s*/*AC*',monkeyName));

for idx_drv = 2:length(dir_target)
	dirTok = split(dir_target(idx_drv).folder,'/');
	nameTok = split(dirTok(end),'-');
	%monkeyName = nameTok{1};
	sessionDate = nameTok{2};
	driveID = dir_target(idx_drv).name;
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
		custTag = datestr(now,'_mmdd');
		get_RI_par(monkeyName,sessionDate,driveID,clustTypeName,custTag);
	catch
		fprintf('\nError occurred with %s_%s_%s... Skipping...\n',monkeyName,sessionDate,driveID)
	end
end
