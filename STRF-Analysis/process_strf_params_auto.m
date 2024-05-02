monkeyName = input('Enter monkey name [Cassius or Miyagi]: ','s');
clustType = input('Which unit type to process?[(s)ingle, (m)ua or (b)oth]: ','s');
Const.p_crit = input('Enter p_crit for RI thresholding: ');	% Default p_crit=0.01

Const.NBoot = 8;
Const.NB = 100;
%Const.p_crit = 1.e-33;
%Const.p_crit = .01;

warning('off')

dir_result = fullfile('~','STRF','STRFParams');

load(sprintf('~/STRF/lamProf_%s.mat',monkeyName));

tAllDataStart = tic;

for idx_data = 1:length(lamProf)
	fprintf('\n\nProcessing %s ...\n\n',lamProf(idx_data).dataName)
	dataLog(idx_data).status = 'progress';

	% Parsing
	try
		% New parser
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
			dataLog(idx_data).status = 'failed';
			dataLog(idx_data).message = 'no matching drive found';
			fprintf('\n\n!!!Error: No matching drive found for %s. Moving on!!!\n\n',lamProf(idx_data).dataName)
			clearvars -except monkeyName clustType dir_result lamProf idx_data Const dataLog tAllDataStart
			continue
		elseif numDrv > 1
			fprintf('\n!!!Error: Multiple drives matched in %s!!!\n\n',sessionDate)
			dataLog(idx_data).status = 'failed';
			dataLog(idx_data).message = 'multiple matching drive found';
			fprintf('\n\n!!!Error: Multiple matching drive found for %s. Moving on!!!\n\n',lamProf(idx_data).dataName)
			clearvars -except monkeyName clustType dir_result lamProf idx_data Const dataLog tAllDataStart
			continue
		end
		% Initial log
		dataLog(idx_data).monkeyName = monkeyName;
		dataLog(idx_data).sessionDate = sessionDate;
		dataLog(idx_data).driveID = driveID;
	catch
		dataLog(idx_data).status = 'failed';
		dataLog(idx_data).message = 'cannot parse recording info';
		fprintf('\n\n!!!Error: Cannot parse recording info for %s. Moving on!!!\n\n',lamProf(idx_data).dataName)
		clearvars -except monkeyName clustType dir_result lamProf idx_data Const dataLog tAllDataStart
		continue
	end
	

	% Find coordinates
	try
		load(sprintf('coord_ac_%s',monkeyName));
	catch
		dataLog(idx_data).status = 'failed';
		dataLog(idx_data).message = 'coordinate not found';
		fprintf('\n\n!!!Error: Cannot find coordinates for %s. Moving on!!!\n\n',lamProf(idx_data).dataName)
		clearvars -except monkeyName clustType dir_result lamProf idx_data Const dataLog tAllDataStart
		continue
	end

	% Identify area
	try
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
	catch
		dataLog(idx_data).status = 'failed';
		dataLog(idx_data).message = 'cannot identify area';
		fprintf('\n\n!!!Error: Cannot identify the area for %s. Moving on!!!\n\n',lamProf(idx_data).dataName)
		clearvars -except monkeyName clustType dir_result lamProf idx_data Const dataLog tAllDataStart
		continue
	end

	% Find cluster info file
	try
		clusterInfo = importTSV(fullfile('~','kiloSorted_DMR',sprintf('Mr%s-%s',monkeyName,sessionDate),driveID,'KS2_7_AC','ClusterInfo','cluster_info_new.tsv'));
	catch
		dataLog(idx_data).status = 'failed';
		dataLog(idx_data).message = 'cluster info not available';
		fprintf('\n\n!!!Error: Cannot cluster info not available for %s. Moving on!!!\n\n',lamProf(idx_data).dataName)
		clearvars -except monkeyName clustType dir_result lamProf idx_data Const dataLog tAllDataStart
		continue
	end

	% Find good clusters based on Reliability Indices
	try	
		goodSTRFUnits = find_good_clusters(monkeyName,sessionDate,driveID,Const.NBoot,Const.NB,Const.p_crit);
		clust2proc = cell2mat(goodSTRFUnits(:,1));
	catch
		dataLog(idx_data).status = 'failed';
		dataLog(idx_data).message = 'good clusters not found';
		fprintf('\n\n!!!Error: Cannot find good clusters for %s. Moving on!!!\n\n',lamProf(idx_data).dataName)
		clearvars -except monkeyName clustType dir_result lamProf idx_data Const dataLog tAllDataStart
		continue
	end

	% Calculate STRF parameters
	try
		UberSTRF = calculate_strf_params(monkeyName,sessionDate,driveID,clust2proc);
	catch
		dataLog(idx_data).status = 'failed';
		dataLog(idx_data).message = 'calculating STRF parameters failed';
		fprintf('\n\n!!!Error: Calculating STRF parameters failed for %s. Moving on!!!\n\n',lamProf(idx_data).dataName)
		clearvars -except monkeyName clustType dir_result lamProf idx_data Const dataLog tAllDataStart
		continue
	end

	% Save result
	try
		var2save = {'monkeyName','sessionDate','driveID','area','Const','UberSTRF','dataLog'};
		dataLog(idx_data).status = 'success';
		save(fullfile(dir_result,sprintf('STRFParams_%s_%s_%s.mat',monkeyName,sessionDate,driveID)),var2save{:});
	catch
		dataLog(idx_data).status = 'failed';
		dataLog(idx_data).message = 'cannot save result';
		fprintf('\n\n!!!Error: Cannot save the result for %s. Moving on!!!\n\n',lamProf(idx_data).dataName)
		clearvars -except monkeyName clustType dir_result lamProf idx_data Const dataLog tAllDataStart
		continue
	end
	
	% Clear
	display('Clearing workspace...')
	close all
	clearvars -except monkeyName clustType dir_result lamProf idx_data Const dataLog tAllDataStart
end

tAllDataEnd = toc(tAllDataStart)

