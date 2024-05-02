monkeyName = input('Enter monkey name [Cassius or Miyagi]: ','s');

warning('off')

dir_result = fullfile('~','STRF','STRFParams');
list_result = dir(fullfile(dir_result,sprintf('*STRFParams_%s*.mat',monkeyName)));

load(sprintf('~/STRF/lamProf_%s.mat',monkeyName));


for idx_result = 1:length(list_result)
	fprintf('\n\nChecking %s ...\n\n',list_result(idx_result).name)

	% Check if area exists.
	if isempty(who('-file',fullfile(dir_result,list_result(idx_result).name),'area'))
	fprintf('\n\nArea missing for %s. Identifying area ...\n\n',list_result(idx_result).name)

		% Parsing
		try
			% New parser
			nameRoot = split(list_result(idx_result).name,'.');
			nameTok = split(nameRoot{1},'_');
			monkeyName = nameTok{2};
			sessionDate = nameTok{3};
			driveID = sprintf('%s_%s_%s',nameTok{4:6});
			dir_target = dir(sprintf('~/kiloSorted_DMR/Mr%s-%s/D*',monkeyName,sessionDate));
		catch
			fprintf('\n\n!!!Error: Cannot parse recording info for %s. Moving on!!!\n\n',lamProf(idx_result).dataName)
		end
		

		% Find coordinate file.
		try
			load(sprintf('coord_ac_%s',monkeyName));
		catch
			fprintf('\n\n!!!Error: Cannot find coordinate file for %s. Moving on!!!\n\n',monkeyName)
		end
		
		% Find current coordinate
		lamProf_date = cellfun(@num2str,{lamProf.date},'UniformOutput',false)';
		currCoord = lamProf(find(strcmp(lamProf_date,sessionDate))).coord;
		
		% Identify area
		try
			for idx_area = 1:size(coord_ac,2)
				for idx_coord = 2:length(coord_ac(:,idx_area))
					if ~isempty(coord_ac{idx_coord,idx_area})
						if coord_ac{idx_coord,idx_area} == currCoord
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
			fprintf('\n\n!!!Error: Cannot identify the area for %s. Moving on!!!\n\n',sessionDate)
		end
		save(list_result(idx_result).name,'area','-append','-v7.3','-nocompression')
		fprintf('\nAppending area completed for %s!\n\n',list_result(idx_result).name)
	else
		fprintf('\nArea exists for %s. Move on!\n\n',list_result(idx_result).name)
	end

end

