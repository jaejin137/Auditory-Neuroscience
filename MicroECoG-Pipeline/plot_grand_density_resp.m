data = table2cell(datafilenames(2:end,1));
GrandDecAccDensity_Resp = nan(6,3,length(data));
for k=1:length(data)
	try
		if k==1
			GrandDecAccDensity_Resp(:,:,1) = struct2array(load(data{1,1},'decAccDensityResp'));
		end
		GrandDecAccDensity_Resp(:,:,k) = struct2array(load(data{k,1},'decAccDensityResp'));
	catch
		continue
	end
end
num_density_data = sum(~isnan(GrandDecAccDensity_Resp(1,1,:)));

GrandDecAccDensity_Resp_Mean = nanmean(GrandDecAccDensity_Resp,3);
GrandDecAccDensity_Resp_SEM = nanstd(GrandDecAccDensity_Resp,0,3)/sqrt(num_density_data);

figure
label_band = {'Broad','Theta','Alpha','Beta','Gamma','HFO'};
for i=1:6
	subplot(2,3,i)
	errorbar(GrandDecAccDensity_Resp_Mean(i,:)',GrandDecAccDensity_Resp_SEM(i,:)','o-','LineWidth',2)
	xlim([.9 3.1]); ylim([50 100])
	xticks([1 2 3])
	xticklabels({'Low','Mid','High'})
	xlabel('Channel Density'); ylabel('C.V. Accuracy');
	legend on; legend('boxoff')
	legend(label_band{i},'Location','southeast')
	set(gca,'FontSize',14)
	currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),700,420]);
end
