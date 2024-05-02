% Plot decoding accuracies as a function of the number of channels.
if decodeWhat == 'r';
	condIndUB = condIndUB_resp;
else
	condIndUB = condIndUB_tnr;
end

figure(3)
if sessNum == 1
	currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),1200,1200]);
end
%suptitle('Decoding Accuracies for Top 20 Percentile')
% Broad band
subplot(3,2,1)
h311 = errorbar(class_acc_select_broad_mean(:,length(condIndUB)),class_acc_select_broad_std(:,length(condIndUB)),'o-','LineWidth',2);
hold on
h312 = errorbar(class_acc_select_broad_rnd_mean(:,length(condIndUB)),class_acc_select_broad_std(:,length(condIndUB)),'o--','LineWidth',2);
ylim([0 100])
xlabel('No. of Channels'); ylabel('CV Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('Broad Band')
drawnow
% Theta band
figure(3)
subplot(3,2,2)
h321 = errorbar(class_acc_select_theta_mean(:,length(condIndUB)),class_acc_select_theta_std(:,length(condIndUB)),'o-','LineWidth',2);
hold on
h322 = errorbar(class_acc_select_theta_rnd_mean(:,length(condIndUB)),class_acc_select_theta_std(:,length(condIndUB)),'o--','LineWidth',2);
ylim([0 100])
xlabel('No. of Channels'); ylabel('CV Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('Theta Band')
drawnow
% Alpha band
figure(3)
subplot(3,2,3)
h331 = errorbar(class_acc_select_alpha_mean(:,length(condIndUB)),class_acc_select_alpha_std(:,length(condIndUB)),'o-','LineWidth',2);
hold on
h332 = errorbar(class_acc_select_alpha_rnd_mean(:,length(condIndUB)),class_acc_select_alpha_std(:,length(condIndUB)),'o--','LineWidth',2);
ylim([0 100])
xlabel('No. of Channels'); ylabel('CV Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('Alpha Band')
drawnow
% Beta band
figure(3)
subplot(3,2,4)
h341 = errorbar(class_acc_select_beta_mean(:,length(condIndUB)),class_acc_select_beta_std(:,length(condIndUB)),'o-','LineWidth',2);
hold on
h342 = errorbar(class_acc_select_beta_rnd_mean(:,length(condIndUB)),class_acc_select_beta_std(:,length(condIndUB)),'o--','LineWidth',2);
ylim([0 100])
xlabel('No. of Channels'); ylabel('CV Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('Beta band')
drawnow
% Gamma band
figure(3)
subplot(3,2,5)
ylim([0 100])
h351 = errorbar(class_acc_select_gamma_mean(:,length(condIndUB)),class_acc_select_gamma_std(:,length(condIndUB)),'o-','LineWidth',2);
hold on
h352 = errorbar(class_acc_select_gamma_rnd_mean(:,length(condIndUB)),class_acc_select_gamma_std(:,length(condIndUB)),'o--','LineWidth',2);
xlabel('No. of Channels'); ylabel('CV Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('Gamma band')
drawnow
% HFO band
figure(3)
subplot(3,2,6)
h361 = errorbar(class_acc_select_hfo_mean(:,length(condIndUB)),class_acc_select_hfo_std(:,length(condIndUB)),'o-','LineWidth',2);
hold on
h362 = errorbar(class_acc_select_hfo_rnd_mean(:,length(condIndUB)),class_acc_select_hfo_std(:,length(condIndUB)),'o--','LineWidth',2);
ylim([0 100])
xlabel('No. of Channels'); ylabel('CV Accuracy')
xticks([1 2 3 4 5])
xticklabels(numchpow)
title('HFO Band')
drawnow
