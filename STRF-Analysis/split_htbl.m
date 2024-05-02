coordind_core = cellfun(@(x) [-x(2)+8,x(1)+8], coord_core, 'UniformOutput', false);
coordind_belt = cellfun(@(x) [-x(2)+8,x(1)+8], coord_belt, 'UniformOutput', false);

htbl_core_supra = zeros(15,15,14);
htbl_core_infra = zeros(15,15,14);
htbl_belt_supra = zeros(15,15,14);
htbl_belt_infra = zeros(15,15,14);

for j = 1:14
	for i = 1:length(coordind_core)
		tmpind = num2cell([coordind_core{i},j]);
		if ~isnan(htbl_supra(tmpind{:}))
			htbl_core_supra(tmpind{:}) = htbl_supra(tmpind{:});
		end
	end
	for i = 1:length(coordind_core)
		tmpind = num2cell([coordind_core{i},j]);
		if ~isnan(htbl_infra(tmpind{:}))
			htbl_core_infra(tmpind{:}) = htbl_infra(tmpind{:});
		end
	end
	for i = 1:length(coordind_belt)
		tmpind = num2cell([coordind_belt{i},j]);
		if ~isnan(htbl_supra(tmpind{:}))
			htbl_belt_supra(tmpind{:}) = htbl_supra(tmpind{:});
		end
	end
	for i = 1:length(coordind_belt)
		tmpind = num2cell([coordind_belt{i},j]);
		if ~isnan(htbl_infra(tmpind{:}))
			htbl_belt_infra(tmpind{:}) = htbl_infra(tmpind{:});
		end
	end
end
