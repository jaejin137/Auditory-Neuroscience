% Plot decoding accuracies as a function of the number of channels per band.
decodeWhat = input('Decode what?(r/t) ','s');
sessNum = input('Session number? ');
figNum = input('Figure number? ');

if ~exist('figNum') || isempty(figNum)
	figure
	figNum = get(gcf,'Number');
end

figure(figNum)
currPos = get(figNum,'Position'); set(figNum,'Position',[currPos(1),currPos(2),800,320]);

% Broad band
subplot(1,2,1)
if decodeWhat == 'r'
	h311 = errorbar(GrandDecodingAccuracy(sessNum).class_acc_select_broad_mean(:,end),GrandDecodingAccuracy(sessNum).class_acc_select_broad_std(:,end),'o-','LineWidth',2);
hold on
	h312 = errorbar(GrandDecodingAccuracy(sessNum).class_acc_select_broad_rnd_mean(:,end),GrandDecodingAccuracy(sessNum).class_acc_select_broad_std(:,end),'o--','LineWidth',2);
else
	h311 = errorbar(GrandDecodingAccuracyTNR(sessNum).class_acc_select_broad_mean(:,end),GrandDecodingAccuracyTNR(sessNum).class_acc_select_broad_std(:,end),'o-','LineWidth',2);
hold on
	h312 = errorbar(GrandDecodingAccuracyTNR(sessNum).class_acc_select_broad_rnd_mean(:,end),GrandDecodingAccuracyTNR(sessNum).class_acc_select_broad_std(:,end),'o--','LineWidth',2);
end
set(gca,'FontSize',16)
xlim([0 7])
ylim([0 100])
xlabel('No. of Channels'); ylabel('Class. Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('Broad Band')
drawnow
% Gamma band
figure(figNum)
subplot(1,2,2)
if decodeWhat == 'r'
	h351 = errorbar(GrandDecodingAccuracy(sessNum).class_acc_select_gamma_mean(:,end),GrandDecodingAccuracy(sessNum).class_acc_select_gamma_std(:,end),'o-','LineWidth',2);
	hold on
	h352 = errorbar(GrandDecodingAccuracy(sessNum).class_acc_select_gamma_rnd_mean(:,end),GrandDecodingAccuracy(sessNum).class_acc_select_gamma_std(:,end),'o--','LineWidth',2);
else
	h351 = errorbar(GrandDecodingAccuracyTNR(sessNum).class_acc_select_gamma_mean(:,end),GrandDecodingAccuracyTNR(sessNum).class_acc_select_gamma_std(:,end),'o-','LineWidth',2);
	hold on
	h352 = errorbar(GrandDecodingAccuracyTNR(sessNum).class_acc_select_gamma_rnd_mean(:,end),GrandDecodingAccuracyTNR(sessNum).class_acc_select_gamma_std(:,end),'o--','LineWidth',2);
end
set(gca,'FontSize',16)
xlim([0 7])
ylim([0 100])
xlabel('No. of Channels'); ylabel('Class. Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('Gamma band')
drawnow
