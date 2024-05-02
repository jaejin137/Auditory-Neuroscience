monkeyName = input('Enter monkey name: ','s');
%monkeyName = 'Cassius';

clustType = input('Cluster unit type to process: [(s)ingle, (m)ua, or (b)oth] ','s');

%ifPlot = input('Plot the result? (y/n) ','s');
ifPlot = 'n';

dir_result = fullfile('~','STRF','STRFParams');
dir_kilosorted = fullfile('~','kiloSorted_DMR');

list_result = dir(fullfile(dir_result,sprintf('STRFParams_%s*.mat',monkeyName)));
if isempty(list_result)
	fprintf("\n!!!Error: Could not retrieve STRF files!!!\n\n")
	return
end

%load('targetInfo_lamProf.mat')
load('areaColor.mat')
sessionColor = [[0,0,0];[1,1,0];[1,0,1];[0,1,1];[1,0,0];[0,1,0];[0,0,1];[.5,0,0];[0,.5,0];[0,0,.5]];


if clustType == 's'
	clustTypeName = 'Single Unit';
elseif clustType == 'm'
	clustTypeName = 'MUA';
elseif clustType == 'b'
	clustTypeName = 'Single Unit & MUA';
else
	display('!!!Error: Wrong cluster type given!!!')
	return
end


if ifPlot == 'y'
	h = figure;
else
	h = [];
end


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

		clust2proc = [];
		%for idx_clust = 1:size(UberSTRF,2)
		%	clust2proc(end+1,1) = UberSTRF(idx_clust).ClustNum;
		%end

		if clustType == 's'
			clust2proc = cluster_info.id(find(ismember(cluster_info.group,{'good'})));
			clustTypeName = 'Single Unit';
		elseif clustType == 'm'
			clust2proc = cluster_info.id(find(ismember(cluster_info.group,{'mua'})));
			clustTypeName = 'MUA';
		elseif clustType == 'b'
			clust2proc = cluster_info.id(find(ismember(cluster_info.group,{'good','mua'})));
			clustTypeName = 'Single Unit & MUA';
		else
			display('!!!Error: Wrong cluster type given!!!')
			return
		end

		chanArr = nan(length(list_result),256);
		for idx_clust = 1:length(clust2proc)
			idx_match = [];
			for i=1:length(UberSTRF)
				if UberSTRF(i).ClustNum == clust2proc(idx_clust)
					idx_match = i;
				end
			end
			if ~isempty(idx_match)
				chanArr(idx_data,idx_match) = cluster_info.ch(find(cluster_info.id==clust2proc(idx_match)));
				Delay(idx_data,idx_match) = UberSTRF(idx_match).RFParam.Delay;
				Duration(idx_data,idx_match) = UberSTRF(idx_match).RFParam.Duration;
				BF(idx_data,idx_match) = UberSTRF(idx_match).RFParam.BF;
				PLI(idx_data,idx_match) = UberSTRF(idx_match).RFParam.PLI;
				PLI2(idx_data,idx_match) = UberSTRF(idx_match).RFParam.PLI2;
				DSI(idx_data,idx_match) = UberSTRF(idx_match).RFParam.DSI;
				PeakDelay(idx_data,idx_match) = UberSTRF(idx_match).RFParam.PeakDelay;
				PeakEnvDelay(idx_data,idx_match) = UberSTRF(idx_match).RFParam.PeakEnvDelay;
				PeakBF(idx_data,idx_match) = UberSTRF(idx_match).RFParam.PeakBF;

				if ifPlot == 'y'
					figure(h)
					subplot(3,4,1)
					scatter3(UberSTRF(idx_match).RFParam.Delay,UberSTRF(idx_match).RFParam.Duration,UberSTRF(idx_match).RFParam.BF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
					xlabel('Delay'); ylabel('Duration'); zlabel('BF');
					%drawnow
					hold on
					subplot(3,4,2)
					scatter3(UberSTRF(idx_match).RFParam.Delay,UberSTRF(idx_match).RFParam.Duration,UberSTRF(idx_match).RFParam.BF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
					xlabel('Delay'); ylabel('Duration'); zlabel('BF');
					view([90 0])
					%drawnow
					hold on
					subplot(3,4,3)
					scatter3(UberSTRF(idx_match).RFParam.Delay,UberSTRF(idx_match).RFParam.Duration,UberSTRF(idx_match).RFParam.BF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
					xlabel('Delay'); ylabel('Duration'); zlabel('BF');
					view([0 0])
					%drawnow
					hold on
					subplot(3,4,4)
					scatter3(UberSTRF(idx_match).RFParam.Delay,UberSTRF(idx_match).RFParam.Duration,UberSTRF(idx_match).RFParam.BF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
					xlabel('Delay'); ylabel('Duration'); zlabel('BF');
					view([0 90])
					%drawnow
					hold on
					
					figure(h)
					subplot(3,4,5)
					scatter3(UberSTRF(idx_match).RFParam.PeakDelay,UberSTRF(idx_match).RFParam.PeakEnvDelay,UberSTRF(idx_match).RFParam.PeakBF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
					xlabel('PeakDelay'); ylabel('PeakEnvDelay'); zlabel('PeakBF');
					%drawnow
					hold on
					subplot(3,4,6)
					scatter3(UberSTRF(idx_match).RFParam.PeakDelay,UberSTRF(idx_match).RFParam.PeakEnvDelay,UberSTRF(idx_match).RFParam.PeakBF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
					xlabel('PeakDelay'); ylabel('PeakEnvDelay'); zlabel('PeakBF');
					view([90 0])
					%drawnow
					hold on
					subplot(3,4,7)
					scatter3(UberSTRF(idx_match).RFParam.PeakDelay,UberSTRF(idx_match).RFParam.PeakEnvDelay,UberSTRF(idx_match).RFParam.PeakBF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
					xlabel('PeakDelay'); ylabel('PeakEnvDelay'); zlabel('PeakBF');
					view([0 0])
					%drawnow
					hold on
					subplot(3,4,8)
					scatter3(UberSTRF(idx_match).RFParam.PeakDelay,UberSTRF(idx_match).RFParam.PeakEnvDelay,UberSTRF(idx_match).RFParam.PeakBF,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
					xlabel('PeakDelay'); ylabel('PeakEnvDelay'); zlabel('PeakBF');
					view([0 90])
					%drawnow
					hold on
					
					figure(h)
					subplot(3,4,9)
					scatter3(UberSTRF(idx_match).RFParam.PLI,UberSTRF(idx_match).RFParam.PLI2,UberSTRF(idx_match).RFParam.DSI,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
					xlabel('PLI'); ylabel('PLI2'); zlabel('DSI');
					%drawnow
					hold on
					subplot(3,4,10)
					scatter3(UberSTRF(idx_match).RFParam.PLI,UberSTRF(idx_match).RFParam.PLI2,UberSTRF(idx_match).RFParam.DSI,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
					xlabel('PLI'); ylabel('PLI2'); zlabel('DSI');
					view([90 0])
					%drawnow
					hold on
					subplot(3,4,11)
					scatter3(UberSTRF(idx_match).RFParam.PLI,UberSTRF(idx_match).RFParam.PLI2,UberSTRF(idx_match).RFParam.DSI,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
					xlabel('PLI'); ylabel('PLI2'); zlabel('DSI');
					view([0 0])
					%drawnow
					hold on
					subplot(3,4,12)
					scatter3(UberSTRF(idx_match).RFParam.PLI,UberSTRF(idx_match).RFParam.PLI2,UberSTRF(idx_match).RFParam.DSI,10,areaColor{strcmp(areaColor(:,1),list_result(idx_data).area),2},'filled')
					xlabel('PLI'); ylabel('PLI2'); zlabel('DSI');
					view([0 90])
					%drawnow
					hold on
				end
			end
		end
	catch
		fprintf('\n\n!!! Error occurred with data %s. Skipping...!!!\n',sessionDate)
	end
	%drawnow
	%hold on
	all_figs = findobj(0, 'type', 'figure');
	if ifPlot=='y'
		delete(setdiff(all_figs, h.Number));
	end
	clearvars -except monkeyName clustType clustTypeName dir_kilosorted dir_result list_result sessionColor areaColor idx_data idx_core idx_belt h ifPlot chanArr Delay Duration BF PLI PLI2 DSI PeakDelay PeakEnvDelay PeakBF
end

if ifPlot == 'y'
	sgtitle(clustTypeName)
	currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1000,600]);
else
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
for i = 1:size(strf_featarr,1)*size(strf_featarr,2)
	if ~strf_featarr2d(i,:) == 0
		idx_nozeros(end+1) = i;
	end
end
feat2d = strf_featarr2d(idx_nozeros,:);
% Get k-means clustering for first 3 features.


label_feat = {'BF','Delay','Duration','PeakBF','PeakDelay','PeakEnvDelay','DSI','PLI','PLI2'};
% Plot scatter plot in feature space.
h_k = figure;
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),700,300]);
for j = 1:2
	[idx_area,C_area] = kmeans(feat2d(:,3*(j-1)+1:3*j),2,'Replicates',5,'Display','final');
	subplot(1,2,j)
	for i = 1:length(feat2d)
		scatter3(feat2d(i,3*(j-1)+1),feat2d(i,3*(j-1)+2),feat2d(i,3*(j-1)+3),20,idx_area(i),'filled')
		hold on
	end
	plot3(C_area(:,1),C_area(:,2),C_area(:,3),'r.','MarkerSize',40)
	xlabel(label_feat{3*(j-1)+1});
	ylabel(label_feat{3*(j-1)+2});
	zlabel(label_feat{3*(j-1)+3});
	drawnow
end
sgtitle('k-means clustering','FontSize',24)



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
sgtitle(clustTypeName)
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
sgtitle(clustTypeName)
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1000,800]);


%% Plot STRF parameters vs channels
h_ch = figure;

for l = 1:size(strf_featarr,3)
	figure(h_ch)
	subplot(3,3,l)
	for j = 1:size(chanArr,2)
		scatter(strf_featarr(idx_core,j,1),chanArr(idx_core,j),20,'r','filled')
		hold on
		scatter(strf_featarr(idx_belt,j,1),chanArr(idx_belt,j),20,'g','filled')
	end
	set(gca,'YDir','reverse')
	title(label_feat{l})
	ylabel('Channels')
	legend('core','belt')
	%legend('box','off')
	legend('location','northeast')
	drawnow
end
sgtitle(clustTypeName)
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1000,800]);


%% Plot STRF parameters vs channels
Channels = 1:24;
h_ch_bar = figure;
for l = 1:size(strf_featarr,3)
    figure(h_ch_bar)
    subplot(3,3,l)
	featarr_core = strf_featarr(idx_core,:);
	chanArr_core = chanArr(idx_core,:);
	featarr_belt = strf_featarr(idx_belt,:);
	chanArr_belt = chanArr(idx_belt,:);
	for chan = Channels
		feat_chan_core = featarr_core(find(chanArr_core==chan));
		feat_chan_belt = featarr_belt(find(chanArr_belt==chan));
		bardata(chan,:) = [mean(feat_chan_core),mean(feat_chan_belt)];
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
sgtitle(clustTypeName)
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1000,800]);



