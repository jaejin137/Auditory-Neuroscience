monkeyName = input('Enter monkey name: ','s');
%monkeyName = 'Cassius';

%ifPlot = input('Plot the result? (y/n) ','s');
ifPlot = 'y';

dir_result = fullfile('~','STRF','Result');
dir_kilosorted = fullfile('~','kiloSorted_DMR');

list_result = dir(fullfile(dir_result,sprintf('STRFParams_%s*.mat',monkeyName)));

%load('targetInfo_lamProf.mat')
load('areaColor.mat')
sessionColor = [[0,0,0];[1,1,0];[1,0,1];[0,1,1];[1,0,0];[0,1,0];[0,0,1];[.5,0,0];[0,.5,0];[0,0,.5]];

%if ifPlot == 'y'
	h = figure;
	currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1000,600]);
%end


idx_core = [];
idx_belt = [];

for idx_data = 1:length(list_result)
%for idx_data = 1:4
	try

		load(fullfile(dir_result,list_result(idx_data).name));
		list_result(idx_data).area = area;

		if list_result(idx_data).area == 'core'
			idx_core(end+1) = idx_data;
		else
			idx_belt(end+1) = idx_data;
		end

		basename = split(list_result(idx_data).name,'.');
		dataTok = split(basename{1},'_');
		sessionDate = dataTok{3};
		driveID = sprintf('%s_%s_%s',dataTok{4:6});
		
		fprintf('\n\nData for %s loaded. Area: %s\n',sessionDate,area)

		path_clustInfo = fullfile(dir_kilosorted,sprintf('Mr%s-%s',monkeyName,sessionDate),driveID,'KS2_7_AC','ClusterInfo');

		cluster_info = importTSV(fullfile(path_clustInfo,'cluster_info_new.tsv'));

		goodClust = [];
		for idx_clust = 1:size(UberSTRF,2)
			goodClust(end+1,1) = UberSTRF(idx_clust).ClustNum;
		end

		%Channel = nan(length(list_result),256);
		for idx_clust = 1:length(goodClust)
			Channel(idx_data,idx_clust) = cluster_info.ch(find(cluster_info.id==goodClust(idx_clust)));
			Delay(idx_data,idx_clust) = UberSTRF(idx_clust).RFParam.Delay;
			Duration(idx_data,idx_clust) = UberSTRF(idx_clust).RFParam.Duration;
			BF(idx_data,idx_clust) = UberSTRF(idx_clust).RFParam.BF;
			PLI(idx_data,idx_clust) = UberSTRF(idx_clust).RFParam.PLI;
			PLI2(idx_data,idx_clust) = UberSTRF(idx_clust).RFParam.PLI2;
			DSI(idx_data,idx_clust) = UberSTRF(idx_clust).RFParam.DSI;
			PeakDelay(idx_data,idx_clust) = UberSTRF(idx_clust).RFParam.PeakDelay;
			PeakEnvDelay(idx_data,idx_clust) = UberSTRF(idx_clust).RFParam.PeakEnvDelay;
			PeakBF(idx_data,idx_clust) = UberSTRF(idx_clust).RFParam.PeakBF;

			if ifPlot == 'y'
				figure(h)
				subplot(3,4,1)
				scatter3(UberSTRF(idx_clust).RFParam.Delay,UberSTRF(idx_clust).RFParam.Duration,UberSTRF(idx_clust).RFParam.BF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
				xlabel('Delay'); ylabel('Duration'); zlabel('BF');
				%drawnow
				hold on
				subplot(3,4,2)
				scatter3(UberSTRF(idx_clust).RFParam.Delay,UberSTRF(idx_clust).RFParam.Duration,UberSTRF(idx_clust).RFParam.BF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
				xlabel('Delay'); ylabel('Duration'); zlabel('BF');
				view([90 0])
				%drawnow
				hold on
				subplot(3,4,3)
				scatter3(UberSTRF(idx_clust).RFParam.Delay,UberSTRF(idx_clust).RFParam.Duration,UberSTRF(idx_clust).RFParam.BF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
				xlabel('Delay'); ylabel('Duration'); zlabel('BF');
				view([0 0])
				%drawnow
				hold on
				subplot(3,4,4)
				scatter3(UberSTRF(idx_clust).RFParam.Delay,UberSTRF(idx_clust).RFParam.Duration,UberSTRF(idx_clust).RFParam.BF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
				xlabel('Delay'); ylabel('Duration'); zlabel('BF');
				view([0 90])
				%drawnow
				hold on

				figure(h)
				subplot(3,4,5)
				scatter3(UberSTRF(idx_clust).RFParam.PeakDelay,UberSTRF(idx_clust).RFParam.PeakEnvDelay,UberSTRF(idx_clust).RFParam.PeakBF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
				xlabel('PeakDelay'); ylabel('PeakEnvDelay'); zlabel('PeakBF');
				%drawnow
				hold on
				subplot(3,4,6)
				scatter3(UberSTRF(idx_clust).RFParam.PeakDelay,UberSTRF(idx_clust).RFParam.PeakEnvDelay,UberSTRF(idx_clust).RFParam.PeakBF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
				xlabel('PeakDelay'); ylabel('PeakEnvDelay'); zlabel('PeakBF');
				view([90 0])
				%drawnow
				hold on
				subplot(3,4,7)
				scatter3(UberSTRF(idx_clust).RFParam.PeakDelay,UberSTRF(idx_clust).RFParam.PeakEnvDelay,UberSTRF(idx_clust).RFParam.PeakBF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
				xlabel('PeakDelay'); ylabel('PeakEnvDelay'); zlabel('PeakBF');
				view([0 0])
				%drawnow
				hold on
				subplot(3,4,8)
				scatter3(UberSTRF(idx_clust).RFParam.PeakDelay,UberSTRF(idx_clust).RFParam.PeakEnvDelay,UberSTRF(idx_clust).RFParam.PeakBF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
				xlabel('PeakDelay'); ylabel('PeakEnvDelay'); zlabel('PeakBF');
				view([0 90])
				%drawnow
				hold on

				figure(h)
				subplot(3,4,9)
				scatter3(UberSTRF(idx_clust).RFParam.PLI,UberSTRF(idx_clust).RFParam.PLI2,UberSTRF(idx_clust).RFParam.DSI,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
				xlabel('PLI'); ylabel('PLI2'); zlabel('DSI');
				%drawnow
				hold on
				subplot(3,4,10)
				scatter3(UberSTRF(idx_clust).RFParam.PLI,UberSTRF(idx_clust).RFParam.PLI2,UberSTRF(idx_clust).RFParam.DSI,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
				xlabel('PLI'); ylabel('PLI2'); zlabel('DSI');
				view([90 0])
				%drawnow
				hold on
				subplot(3,4,11)
				scatter3(UberSTRF(idx_clust).RFParam.PLI,UberSTRF(idx_clust).RFParam.PLI2,UberSTRF(idx_clust).RFParam.DSI,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
				xlabel('PLI'); ylabel('PLI2'); zlabel('DSI');
				view([0 0])
				%drawnow
				hold on
				subplot(3,4,12)
				scatter3(UberSTRF(idx_clust).RFParam.PLI,UberSTRF(idx_clust).RFParam.PLI2,UberSTRF(idx_clust).RFParam.DSI,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
				xlabel('PLI'); ylabel('PLI2'); zlabel('DSI');
				view([0 90])
				%drawnow
				hold on
			end
		end
	catch
		fprintf('\n\n!!! Error occurred with data %s. Skipping...!!!\n',sessionDate)
	end
	drawnow
	all_figs = findobj(0, 'type', 'figure');
	delete(setdiff(all_figs, h.Number));
	clearvars -except monkeyName dir_kilosorted dir_result list_result sessionColor areaColor idx_data idx_core idx_belt h ifPlot Channel Delay Duration BF PLI PLI2 DSI PeakDelay PeakEnvDelay PeakBF
end

if ifPlot ~= 'y'
	close(h)
end


%% k-means clustering
% Define feature array.
strf_featarr(:,:,1) = BF;
strf_featarr(:,:,2) = Delay;
strf_featarr(:,:,3) = Duration;
strf_featarr(:,:,4) = PeakBF;
strf_featarr(:,:,5) = PeakDelay;
strf_featarr(:,:,6) = PeakEnvDelay;
strf_featarr(:,:,7) = DSI;
strf_featarr(:,:,8) = PLI;
strf_featarr(:,:,9) = PLI2;
% Reshape feature array into 2D array.
for i = 1:size(strf_featarr,1)
	for j = 1:size(strf_featarr,2)
		strf_featarr2d(size(strf_featarr,2)*(i-1)+j,:) = strf_featarr(i,j,:);
	end
end
% Clean up null data and redefine 2D feature array.
idx_nozeros = [];
for i = 1:size(strf_featarr,1)*size(strf_featarr,1)
	if ~strf_featarr2d(i,:) == 0
		idx_nozeros(end+1) = i;
	end
end
feat2d = strf_featarr2d(idx_nozeros,:);
% Get k-means clustering for first 3 features.
[idx_area,C_area] = kmeans(feat2d(:,1:3),2);

% Plot scatter plot in feature space.
h_k = figure;
for i = 1:length(feat2d)
	scatter3(feat2d(i,1),feat2d(i,2),feat2d(i,3),20,idx_area(i),'filled')
	hold on
	drawnow
end
xlabel('Delay');
ylabel('Duration');
zlabel('BF')
title('k-means clustering')


label_feat = {'BF','Delay','Duration','PeakBF','PeakDelay','PeakEnvDelay','DSI','PLI','PLI2'};


%% Plot distribution of STRF parameters for core and belt
h_dist = figure;

numBin = 20;
for l = 1:size(strf_featarr,3)
	figure(h_dist)
	subplot(3,3,l)
	histogram(nonzeros(strf_featarr(idx_core,:,l)),numBin,'FaceColor','r')
	hold on
	histogram(nonzeros(strf_featarr(idx_belt,:,l)),numBin,'FaceColor','g')
	title(label_feat{l})
	legend('core','belt')
	legend('location','northeast')
	drawnow
end
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1000,800]);


%% Plot STRF parameters vs channels
h_ch = figure;

for l = 1:size(strf_featarr,3)
	figure(h_ch)
	subplot(3,3,l)
	scatter(nonzeros(strf_featarr(idx_core,:,l)),nonzeros(Channel(idx_core,:)),20,'r','filled')
	hold on
	scatter(nonzeros(strf_featarr(idx_belt,:,l)),nonzeros(Channel(idx_belt,:)),20,'g','filled')
	set(gca,'YDir','reverse')
	title(label_feat{l})
	ylabel('Channels')
	legend('core','belt')
	%legend('box','off')
	legend('location','northeast')
	drawnow
end
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1000,800]);


%% Plot bar graph of STRF parameters
h_bar = figure;

for l = 1:size(strf_featarr,3)
	figure(h_bar)
	% For core
	subplot(3,3,l)
	data2plot= nonzeros(strf_featarr(idx_core,:,l));
	idx2plot = find(isfinite(data2plot));
	h_b(1) = bar(1,mean(data2plot(idx2plot)),'r');
	hold on
	errorbar(1,mean(data2plot(idx2plot)),std(data2plot(idx2plot))/sqrt(length(data2plot(idx2plot))),'k','LineWidth',2)
	% For belt
	data2plot= nonzeros(strf_featarr(idx_belt,:,l));
	idx2plot = find(isfinite(data2plot));
	h_b(2) = bar(2,mean(data2plot(idx2plot)),'g');
	errorbar(2,mean(data2plot(idx2plot)),std(data2plot(idx2plot))/sqrt(length(data2plot(idx2plot))),'k','LineWidth',2)
	set(gca,'xticklabels',[])
	legend(h_b,'core','belt')
	legend('location','southeast')
	title(label_feat{l})
	drawnow
end
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1000,800]);


%% Plot STRF parameters vs channels
ch_num = 1:24;
h_ch_bar = figure;
for l = 1:size(strf_featarr,3)
    figure(h_ch_bar)
    subplot(3,3,l)
    % For core
	data2plot_core = nonzeros(strf_featarr(idx_core,:,l));
	idx2plot_core = find(isfinite(data2plot_core));
	ch2plot_core = nonzeros(Channel(idx_core,:));
	ch2plot_core = ch2plot_core(idx2plot_core);
	data2plot_belt = nonzeros(strf_featarr(idx_belt,:,l));
	idx2plot_belt = find(isfinite(data2plot_belt));
	ch2plot_belt = nonzeros(Channel(idx_belt,:));
	ch2plot_belt = ch2plot_belt(idx2plot_belt);
	for i = ch_num
		idx_ch_core{i} = find(ch2plot_core == ch_num(i));
		idx_ch_belt{i} = find(ch2plot_belt == ch_num(i));

		bardata(i,:) = [mean(data2plot_core(idx_ch_core{i})),mean(data2plot_belt(idx_ch_belt{i}))];
	end
	h_bh = barh(bardata);
	h_bh(1).FaceColor = 'r';
	h_bh(1).EdgeColor = 'r';
	h_bh(2).FaceColor = 'g';
	h_bh(2).EdgeColor = 'g';
	set(gca,'YDir','reverse')
    legend(h_bh,'core','belt')
    legend('location','southeast')
    title(label_feat{l})
	ylabel('Channels')
	drawnow
end
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1000,800]);



