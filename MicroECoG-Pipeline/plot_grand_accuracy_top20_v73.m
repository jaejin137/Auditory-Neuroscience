% Plot decoding accuracies as a function of the number of channels per band.
decodeWhat = input('Decode what?(r/t) ','s');
sessNum = input('Session number? ');
figNum = input('Figure number? ');

if ~exist('figNum') || isempty(figNum)
	figure
	figNum = get(gcf,'Number');
end

figure(figNum)
currPos = get(figNum,'Position'); set(figNum,'Position',[currPos(1),currPos(2),800,1200]);

% Broad band
subplot(3,2,1)
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
xlabel('No. of Channels'); ylabel('CV Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('Broad Band')
drawnow
% Theta band
figure(figNum)
subplot(3,2,2)
if decodeWhat == 'r'
	h321 = errorbar(GrandDecodingAccuracy(sessNum).class_acc_select_theta_mean(:,end),GrandDecodingAccuracy(sessNum).class_acc_select_theta_std(:,end),'o-','LineWidth',2);
	hold on
	h322 = errorbar(GrandDecodingAccuracy(sessNum).class_acc_select_theta_rnd_mean(:,end),GrandDecodingAccuracy(sessNum).class_acc_select_theta_std(:,end),'o--','LineWidth',2);
else
	h321 = errorbar(GrandDecodingAccuracyTNR(sessNum).class_acc_select_theta_mean(:,end),GrandDecodingAccuracyTNR(sessNum).class_acc_select_theta_std(:,end),'o-','LineWidth',2);
	hold on
	h322 = errorbar(GrandDecodingAccuracyTNR(sessNum).class_acc_select_theta_rnd_mean(:,end),GrandDecodingAccuracyTNR(sessNum).class_acc_select_theta_std(:,end),'o--','LineWidth',2);
end
set(gca,'FontSize',16)
xlim([0 7])
ylim([0 100])
xlabel('No. of Channels'); ylabel('CV Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('Theta Band')
drawnow
% Alpha band
figure(figNum)
subplot(3,2,3)
if decodeWhat == 'r'
	h331 = errorbar(GrandDecodingAccuracy(sessNum).class_acc_select_alpha_mean(:,end),GrandDecodingAccuracy(sessNum).class_acc_select_alpha_std(:,end),'o-','LineWidth',2);
	hold on
	h332 = errorbar(GrandDecodingAccuracy(sessNum).class_acc_select_alpha_rnd_mean(:,end),GrandDecodingAccuracy(sessNum).class_acc_select_alpha_std(:,end),'o--','LineWidth',2);
else
	h331 = errorbar(GrandDecodingAccuracyTNR(sessNum).class_acc_select_alpha_mean(:,end),GrandDecodingAccuracyTNR(sessNum).class_acc_select_alpha_std(:,end),'o-','LineWidth',2);
	hold on
	h332 = errorbar(GrandDecodingAccuracyTNR(sessNum).class_acc_select_alpha_rnd_mean(:,end),GrandDecodingAccuracyTNR(sessNum).class_acc_select_alpha_std(:,end),'o--','LineWidth',2);
end
set(gca,'FontSize',16)
xlim([0 7])
ylim([0 100])
xlabel('No. of Channels'); ylabel('CV Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('Alpha Band')
drawnow
% Beta band
figure(figNum)
subplot(3,2,4)
if decodeWhat == 'r'
	h341 = errorbar(GrandDecodingAccuracy(sessNum).class_acc_select_beta_mean(:,end),GrandDecodingAccuracy(sessNum).class_acc_select_beta_std(:,end),'o-','LineWidth',2);
	hold on
	h342 = errorbar(GrandDecodingAccuracy(sessNum).class_acc_select_beta_rnd_mean(:,end),GrandDecodingAccuracy(sessNum).class_acc_select_beta_std(:,end),'o--','LineWidth',2);
else
	h341 = errorbar(GrandDecodingAccuracyTNR(sessNum).class_acc_select_beta_mean(:,end),GrandDecodingAccuracyTNR(sessNum).class_acc_select_beta_std(:,end),'o-','LineWidth',2);
	hold on
	h342 = errorbar(GrandDecodingAccuracyTNR(sessNum).class_acc_select_beta_rnd_mean(:,end),GrandDecodingAccuracyTNR(sessNum).class_acc_select_beta_std(:,end),'o--','LineWidth',2);end
set(gca,'FontSize',16)
xlim([0 7])
ylim([0 100])
xlabel('No. of Channels'); ylabel('CV Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('Beta band')
drawnow
% Gamma band
figure(figNum)
subplot(3,2,5)
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
xlabel('No. of Channels'); ylabel('CV Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('Gamma band')
drawnow
% HFO band
figure(figNum)
subplot(3,2,6)
if decodeWhat == 'r'
	h361 = errorbar(GrandDecodingAccuracy(sessNum).class_acc_select_hfo_mean(:,end),GrandDecodingAccuracy(sessNum).class_acc_select_hfo_std(:,end),'o-','LineWidth',2);
	hold on
	h362 = errorbar(GrandDecodingAccuracy(sessNum).class_acc_select_hfo_rnd_mean(:,end),GrandDecodingAccuracy(sessNum).class_acc_select_hfo_std(:,end),'o--','LineWidth',2);
else
	h361 = errorbar(GrandDecodingAccuracyTNR(sessNum).class_acc_select_hfo_mean(:,end),GrandDecodingAccuracyTNR(sessNum).class_acc_select_hfo_std(:,end),'o-','LineWidth',2);
	hold on
	h362 = errorbar(GrandDecodingAccuracyTNR(sessNum).class_acc_select_hfo_rnd_mean(:,end),GrandDecodingAccuracyTNR(sessNum).class_acc_select_hfo_std(:,end),'o--','LineWidth',2);
end
set(gca,'FontSize',16)
xlim([0 7])
ylim([0 100])
xlabel('No. of Channels'); ylabel('CV Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('HFO Band')
drawnow
