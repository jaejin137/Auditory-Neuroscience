

for i=1:length(si_tot)
	if isempty(si_tot(i).si_peak_fw)
		badblockindx(end+1) = i;
	end
end
badblockindx
si_tot_cp = si_tot
si_tot_cp(badblockindx) = [];
si_tot_mean = zeros(4,4,length(si_tot_cp)

for i=1:length(si_tot_cp)
	si_tot_mean(:,1,i) = mean(si_tot_cp(i).si_peak_fw,2);
	si_tot_mean(:,2,i) = mean(si_tot_cp(i).si_peak_bw,2);
	si_tot_mean(:,3,i) = mean(si_tot_cp(i).si_itpc_fw,2);
	si_tot_mean(:,4,i) = mean(si_tot_cp(i).si_itpc_bw,2);
end


si_tot_mean_sum = zeros(4,4)
num_samples = zeros(4,4)
for i=1:4
	for j=1:4
		for k=1:length(si_tot_mean)
			if ~isnan(si_tot_mean(i,j,k))
				si_tot_mean_sum(i,j) = si_tot_mean_sum(i,j) + si_tot_mean(i,j,k);
				num_samples(i,j) = num_samples(i,j)+1;
			end
		end
	end
end


