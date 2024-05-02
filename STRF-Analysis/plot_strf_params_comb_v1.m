% Just plot the results of STRF parameters
datatag = input('Enter the data tag: ','s');
savetag = input('Enter a save tag: ','s');


%dir_result = fullfile('~','STRF','STRFParams');
dir_strf = fullfile('~','STRF');
dir_strfclass = fullfile('~','STRF','STRFClass');
dir_kilosorted = fullfile('~','kiloSorted_DMR');

load(fullfile(dir_strf,'areaColor.mat'));

clustType = 'b';
clustTypeName = 'Both';
p_crit = .01;

monkName = {'Cassius','Miyagi'};

fprintf('\n\nLoading data ...\n\n')

for idx_monk = 1:2
	load(fullfile(dir_strfclass,sprintf('STRFClass_%s_%s_%s.mat',monkName{idx_monk},clustTypeName,datatag)));
	if idx_monk == 1
		STRFClassComb = STRFClass;
		STRFClassComb.name(1:length(STRFClass.date),1) = {'Cassius'};
		clearvars STRFClass;
	else
		STRFClassComb.date(end+1:end+length(STRFClass.date)) = STRFClass.date;
		STRFClassComb.area(end+1:end+length(STRFClass.date)) = STRFClass.area;
		STRFClassComb.clusters(end+1:end+length(STRFClass.date)) = STRFClass.clusters;
		STRFClassComb.chanArr(end+1:end+length(STRFClass.date)) = STRFClass.chanArr;
		STRFClassComb.Delay(end+1:end+length(STRFClass.date)) = STRFClass.Delay;
		STRFClassComb.Duration(end+1:end+length(STRFClass.date)) = STRFClass.Duration;
		STRFClassComb.BF(end+1:end+length(STRFClass.date)) = STRFClass.BF;
		STRFClassComb.PLI(end+1:end+length(STRFClass.date)) = STRFClass.PLI;
		STRFClassComb.PLI2(end+1:end+length(STRFClass.date)) = STRFClass.PLI2;
		STRFClassComb.DSI(end+1:end+length(STRFClass.date)) = STRFClass.DSI;
		STRFClassComb.PeakDelay(end+1:end+length(STRFClass.date)) = STRFClass.PeakDelay;
		STRFClassComb.PeakEnvDelay(end+1:end+length(STRFClass.date)) = STRFClass.PeakEnvDelay;
		STRFClassComb.PeakBF(end+1:end+length(STRFClass.date)) = STRFClass.PeakBF;
		STRFClassComb.T1_10(end+1:end+length(STRFClass.date)) = STRFClass.T1_10;
		STRFClassComb.T2_10(end+1:end+length(STRFClass.date)) = STRFClass.T2_10;
		STRFClassComb.T1_50(end+1:end+length(STRFClass.date)) = STRFClass.T1_50;
		STRFClassComb.T2_50(end+1:end+length(STRFClass.date)) = STRFClass.T2_50;
		STRFClassComb.totNumClust = STRFClassComb.totNumClust + STRFClass.totNumClust;
		STRFClassComb.name(end+1:end+length(STRFClass.date)) = {'Miyagi'};
		clearvars STRFClass;
	end
end

Area = STRFClassComb.area;
clusters = STRFClassComb.clusters;
chanArr = STRFClassComb.chanArr;
Delay = STRFClassComb.Delay;
Duration = STRFClassComb.Duration;
BF = STRFClassComb.BF;
PLI = STRFClassComb.PLI;
PLI2 = STRFClassComb.PLI2;
DSI = STRFClassComb.DSI;
PeakDelay = STRFClassComb.PeakDelay;
PeakEnvDelay = STRFClassComb.PeakEnvDelay;
PeakBF = STRFClassComb.PeakBF;
T1_10 = STRFClassComb.T1_10;
T2_10 = STRFClassComb.T2_10;
T1_50 = STRFClassComb.T1_50;
T2_50 = STRFClassComb.T2_50;

numSess = length(STRFClassComb.date);
idx_core = [];
idx_belt = [];
for k = 1:numSess
	if strcmp(Area(k),'core')
		idx_core(end+1) = k;
	elseif strcmp(Area(k),'belt')
		idx_belt(end+1) = k;
	else
		fprintf('\n\n!!!Error: Unidentified area detected!!!\n\n');
		return
	end
end


fprintf('\n\nPlotting results ...\n\n')

%% Feature spaces

h_feat = figure;
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1400,900]);
for idx_data = 1:numSess
	subplot(3,4,1)
	scatter3(Delay{idx_data},Duration{idx_data},BF{idx_data},20,areaColor{strcmp(areaColor(:,1),Area{idx_data}),2},'filled')
	xlabel('Delay (ms)'); ylabel('Duration (ms)'); zlabel('BF (Hz)');
	hold on
	
	subplot(3,4,2)
	scatter3(Delay{idx_data},Duration{idx_data},BF{idx_data},20,areaColor{strcmp(areaColor(:,1),Area{idx_data}),2},'filled')
	xlabel('Delay (ms)'); ylabel('Duration (ms)'); zlabel('BF (Hz)');
	view([90 0])
	hold on
	
	subplot(3,4,3)
	scatter3(Delay{idx_data},Duration{idx_data},BF{idx_data},20,areaColor{strcmp(areaColor(:,1),Area{idx_data}),2},'filled')
	xlabel('Delay (ms)'); ylabel('Duration (ms)'); zlabel('BF (Hz)');
	view([0 0])
	hold on
	
	subplot(3,4,4)
	scatter3(Delay{idx_data},Duration{idx_data},BF{idx_data},20,areaColor{strcmp(areaColor(:,1),Area{idx_data}),2},'filled')
	xlabel('Delay (ms)'); ylabel('Duration (ms)'); zlabel('BF (Hz)');
	view([0 90])
	hold on
	
	subplot(3,4,5)
	scatter3(PeakDelay{idx_data},PeakEnvDelay{idx_data},PeakBF{idx_data},20,areaColor{strcmp(areaColor(:,1),Area{idx_data}),2},'filled')
	xlabel('PeakDelay (ms)'); ylabel('PeakEnvDelay (ms)'); zlabel('PeakBF (oct)');
	hold on
	
	subplot(3,4,6)
	scatter3(PeakDelay{idx_data},PeakEnvDelay{idx_data},PeakBF{idx_data},20,areaColor{strcmp(areaColor(:,1),Area{idx_data}),2},'filled')
	xlabel('PeakDelay (ms)'); ylabel('PeakEnvDelay (ms)'); zlabel('PeakBF (oct)');
	view([90 0])
	hold on
	
	subplot(3,4,7)
	scatter3(PeakDelay{idx_data},PeakEnvDelay{idx_data},PeakBF{idx_data},20,areaColor{strcmp(areaColor(:,1),Area{idx_data}),2},'filled')
	xlabel('PeakDelay (ms)'); ylabel('PeakEnvDelay (ms)'); zlabel('PeakBF (oct)');
	view([0 0])
	hold on
	
	subplot(3,4,8)
	scatter3(PeakDelay{idx_data},PeakEnvDelay{idx_data},PeakBF{idx_data},20,areaColor{strcmp(areaColor(:,1),Area{idx_data}),2},'filled')
	xlabel('PeakDelay (ms)'); ylabel('PeakEnvDelay (ms)'); zlabel('PeakBF (oct)');
	view([0 90])
	hold on
	
	subplot(3,4,9)
	scatter3(PLI{idx_data},PLI2{idx_data},DSI{idx_data},20,areaColor{strcmp(areaColor(:,1),Area{idx_data}),2},'filled')
	xlabel('PLI'); ylabel('PLI2'); zlabel('DSI');
	hold on
	
	subplot(3,4,10)
	scatter3(PLI{idx_data},PLI2{idx_data},DSI{idx_data},20,areaColor{strcmp(areaColor(:,1),Area{idx_data}),2},'filled')
	xlabel('PLI'); ylabel('PLI2'); zlabel('DSI');
	view([90 0])
	hold on
	
	subplot(3,4,11)
	scatter3(PLI{idx_data},PLI2{idx_data},DSI{idx_data},20,areaColor{strcmp(areaColor(:,1),Area{idx_data}),2},'filled')
	xlabel('PLI'); ylabel('PLI2'); zlabel('DSI');
	view([0 0])
	hold on
	
	subplot(3,4,12)
	scatter3(PLI{idx_data},PLI2{idx_data},DSI{idx_data},20,areaColor{strcmp(areaColor(:,1),Area{idx_data}),2},'filled')
	xlabel('PLI'); ylabel('PLI2'); zlabel('DSI');
	view([0 90])
	hold on
end
drawnow
sgtitle(sprintf('%s (Red:Core, Green: Belt) (p_{sig}=%.2e)',clustTypeName,p_crit),'Interpreter','none','FontWeight','bold','FontSize',16)
fprintf('\nSaving STRFParams_Features_%s.png\n',clustTypeName)
saveas(h_feat,sprintf('~/STRF/Figures/STRFParams_Features_%s_%s.png',clustTypeName,savetag),'png')


%% Construct feature array
% Define feature cell array.
strf_feat = cell(numSess,9);
strf_feat(:,1) = BF;
strf_feat(:,2) = Delay;
strf_feat(:,3) = Duration;
strf_feat(:,4) = PeakBF;
strf_feat(:,5) = PeakDelay;
strf_feat(:,6) = PeakEnvDelay;
strf_feat(:,7) = DSI;
strf_feat(:,8) = PLI;
strf_feat(:,9) = PLI2;

% Change feature cell array to full array
maxNumClust = 0;
for idx_data = 1:length(strf_feat)
	if length(strf_feat{idx_data,1}) > maxNumClust
		maxNumClust = length(strf_feat{idx_data,1});
	end
end
strf_featArr = nan(numSess,maxNumClust,9);
for idx_feat = 1:9
	for idx_data = 1:length(strf_feat)
		strf_featArr(idx_data,1:length(strf_feat{idx_data,1}),idx_feat) = strf_feat{idx_data,idx_feat};
	end
end

% Reshape feature array into 2D array.
for i = 1:size(strf_featArr,1)
	for j = 1:size(strf_featArr,2)
		strf_featArr2d(size(strf_featArr,2)*(i-1)+j,:) = strf_featArr(i,j,:);
	end
end

% Clean up null data and redefine 2D feature array.
idx__nnan = find(~isnan(strf_featArr2d(:,1)));
feat2d = strf_featArr2d(idx__nnan,:);

feat_label = {'BF (Hz)','Delay (ms)','Duration (ms)','PeakBF (Octave)','PeakDelay (ms)','PeakEnvDelay (ms)','DSI','PLI','PLI2'};
feat_name = {'BF','Delay','Duration','PeakBF','PeakDelay','PeakEnvDelay','DSI','PLI','PLI2'};
feat_unit = {'(Hz)','(ms)','(ms)','(Octave)','(ms)','(ms)','','',''};


%% k-means clustering
% Get k-means clustering for first 3 features.
idx_kmean = {};
C_kmean = {};
feat2d_ninf = feat2d;
feat2d_ninf(find(isinf(feat2d))) = nan;
[idx_kmean{1,1},C_kmean{1,1}] = kmeans(feat2d_ninf(:,1:3),2);
[idx_kmean{2,1},C_kmean{2,1}] = kmeans(feat2d_ninf(:,4:6),2);
[idx_kmean{3,1},C_kmean{3,1}] = kmeans(feat2d_ninf(:,7:9),2);

% Plot scatter plot in feature space.
h_k = figure;
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1300,300]);
color_kmean = {'m','c'};
% feature space 1
figure(h_k)
subplot(1,3,1)
for i = 1:length(feat2d_ninf)
	scatter3(feat2d_ninf(i,2),feat2d_ninf(i,3),feat2d_ninf(i,1),20,color_kmean{idx_kmean{1}(i)},'filled')
	hold on
end
xlabel('Delay (ms)');
ylabel('Duration (ms)');
zlabel('BF (Hz)')
drawnow
% feature space 2
figure(h_k)
subplot(1,3,2)
for i = 1:length(feat2d_ninf)
	scatter3(feat2d_ninf(i,5),feat2d_ninf(i,6),feat2d_ninf(i,4),20,color_kmean{idx_kmean{2}(i)},'filled')
	hold on
end
xlabel('PeakDelay (ms)');
ylabel('PeakEnvDelay (ms)');
zlabel('PeakBF (oct)')
drawnow
% feature space 1 (Be carefull with inf and nan!)
figure(h_k)
subplot(1,3,3)
for i = 1:length(feat2d_ninf)
	if ~isnan(idx_kmean{3}(i))
		scatter3(feat2d_ninf(i,8),feat2d_ninf(i,9),feat2d_ninf(i,7),20,color_kmean{idx_kmean{3}(i)},'filled')
	end
	hold on
end
xlabel('PLI');
ylabel('PLI2');
zlabel('DSI')
drawnow
sgtitle(sprintf('%s, k-means clustering (Magenta:cluster 1, Cyan:cluster 2) (p_{sig}=%.2e)',clustTypeName,p_crit),'Interpreter','none','FontWeight','bold','FontSize',16)
fprintf('\nSaving STRFParams_kmeansClust_%s.png\n',clustTypeName)
saveas(h_k,sprintf('~/STRF/Figures/STRFParams_kmeansClust_%s_%s.png',clustTypeName,savetag),'png')



% Plot distribution of STRF parameters for core and belt
h_dist = figure;
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),900,900]);

numBin = 20;
for l = 1:size(strf_featArr,3)
	figure(h_dist)
	strf_featCol = reshape(strf_featArr(:,:,l),[numel(strf_featArr(:,:,l)),1]);
	strf_featCol_core = reshape(strf_featArr(idx_core,:,l),[numel(strf_featArr(idx_core,:,l)),1]);
	strf_featCol_belt = reshape(strf_featArr(idx_belt,:,l),[numel(strf_featArr(idx_belt,:,l)),1]);
	% Take out all nans
	strf_feat_nnan = strf_featCol(~isnan(strf_featCol));
	strf_feat_nnan_core = strf_featCol_core(~isnan(strf_featCol_core));
	strf_feat_nnan_belt = strf_featCol_belt(~isnan(strf_featCol_belt));
	% Take out all infs
	strf_feat_ninf = strf_feat_nnan(~isinf(strf_feat_nnan));
	strf_feat_ninf_core = strf_feat_nnan_core(~isinf(strf_feat_nnan_core));
	strf_feat_ninf_belt = strf_feat_nnan_belt(~isinf(strf_feat_nnan_belt));
	% Significant test
	[stat_dist(l).p, stat_dist(l).h, stat_dist(l).stat] = ranksum(strf_feat_ninf_core,strf_feat_ninf_belt);
	% Figure out bins
	%binMin = min(strf_featCol(~isinf(strf_featCol)),[],'omitnan');
	binMin = min(strf_feat_ninf);
	%binMax = max(strf_featCol(~isinf(strf_featCol)),[],'omitnan');
	binMax = max(strf_feat_ninf);
	binWidth = (binMax-binMin)/numBin;
	binEdges = [binMin:binWidth:binMax];
	figure(h_dist)
	subplot(3,3,l)
	%h1 = histogram(strf_featArr(idx_core,:,l),binEdges,'Normalization','probability','FaceColor','r');
	h1 = histogram(strf_feat_ninf_core,binEdges,'Normalization','probability','FaceColor','r');
	hold on
	%h2 = histogram(strf_featArr(idx_belt,:,l),binEdges,'Normalization','probability','FaceColor','g');
	h2 = histogram(strf_feat_ninf_belt,binEdges,'Normalization','probability','FaceColor','g');
	xlim([binMin binMax])
	%ylim([0 1])
	xlabel(sprintf('%s',feat_unit{l}))
	ylabel(sprintf('Proportion'))
	legend('core','belt')
	%legend([p2 p1],'core','belt')
	legend('location','northeast')
	if stat_dist(l).h == 1
		title(sprintf('%s (p**=%.2e)',feat_name{l},stat_dist(l).p))
	else
		title(sprintf('%s (p=%.2e)',feat_name{l},stat_dist(l).p))
	end
	drawnow
end
sgtitle(sprintf('%s (p_{sig}=%.2e)',clustTypeName,p_crit),'Interpreter','none','FontWeight','bold','FontSize',16)
fprintf('\nSaving STRFParams_Distribution_%s.png\n',clustTypeName)
saveas(h_dist,sprintf('~/STRF/Figures/STRFParams_Distribution_%s_%s.png',clustTypeName,savetag),'png')


%% Plot bar graph of STRF parameters
h_bar = figure;
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),900,900]);

for l = 1:size(strf_featArr,3)
	figure(h_bar)
	% For core
	subplot(3,3,l)
	data2plot_core = nonzeros(strf_featArr(idx_core,:,l));
	idx2plot_core = find(isfinite(data2plot_core));
	h_b(1) = bar(1,mean(data2plot_core(idx2plot_core)),'r');
	hold on
	errorbar(1,mean(data2plot_core(idx2plot_core)),std(data2plot_core(idx2plot_core))/sqrt(length(data2plot_core(idx2plot_core))),'k','LineWidth',2)
	% For belt
	data2plot_belt= nonzeros(strf_featArr(idx_belt,:,l));
	idx2plot_belt = find(isfinite(data2plot_belt));
	h_b(2) = bar(2,mean(data2plot_belt(idx2plot_belt)),'g');
	errorbar(2,mean(data2plot_belt(idx2plot_belt)),std(data2plot_belt(idx2plot_belt))/sqrt(length(data2plot_belt(idx2plot_belt))),'k','LineWidth',2)
	% Significance test
	[stat_mean(l).h, stat_mean(l).p, stat_mean(l).ci, stat_mean(l).stat] = ttest2(data2plot_core(idx2plot_core),data2plot_belt(idx2plot_belt));
	set(gca,'xticklabels',[])
	ylabel(sprintf('%s',feat_unit{l}))
	legend(h_b,'core','belt')
	legend('location','southeast')
	if stat_mean(l).h == 1
		title(sprintf('%s (p**=%.2e)',feat_name{l},stat_mean(l).p))
	else
		title(sprintf('%s (p=%.2e)',feat_name{l},stat_mean(l).p))
	end
	drawnow
end
sgtitle(sprintf('%s (p_{sig}=%.2e)',clustTypeName,p_crit),'Interpreter','none','FontWeight','bold','FontSize',16)
fprintf('\nSaving STRFParams_Bars_%s.png\n\n',clustTypeName)
saveas(h_bar,sprintf('~/STRF/Figures/STRFParams_Bars_%s_%s.png',clustTypeName,savetag),'png')



return

% The following plots are disabled

%% Plot STRF parameters vs channels
h_ch = figure;
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1000,900]);
for l = 1:size(strf_featArr,3)
	figure(h_ch)
	subplot(3,3,l)
	for j = 1:size(chanArr,2)
		scatter(strf_featArr(idx_core,j,l),chanArr(idx_core,j),20,'r','filled')
		hold on
		scatter(strf_featArr(idx_belt,j,l),chanArr(idx_belt,j),20,'g','filled')
		%hold on
		%drawnow
	end
	set(gca,'YDir','reverse')
	title(feat_label{l})
	ylabel('Channels')
	legend('core','belt')
	%legend('box','off')
	legend('location','northeast')
	drawnow
end
sgtitle(sprintf('%s',monkeyName,clustTypeName),'Interpreter','none','FontWeight','bold','FontSize',16)


%% Plot STRF parameters vs channels
Channels = 1:24;
h_ch_bar = figure;
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1000,900]);
for l = 1:size(strf_featArr,3)
    figure(h_ch_bar)
    subplot(3,3,l)
	featArr_core = strf_featArr(idx_core,:);
	chanArr_core = chanArr(idx_core,:);
	featArr_belt = strf_featArr(idx_belt,:);
	chanArr_belt = chanArr(idx_belt,:);
	for chan = Channels
		feat_chan_core = featArr_core(find(chanArr_core==chan));
		feat_chan_belt = featArr_belt(find(chanArr_belt==chan));
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
    title(feat_label{l})
	ylabel('Channels')
	drawnow
end
sgtitle(sprintf('%s, %s',monkeyName,clustTypeName),'Interpreter','none','FontWeight','bold','FontSize',16)


