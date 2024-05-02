decodeWhat = input('Which do you want to decode between (r)esponse and (t)nr? ','s');

for idx_density=1:3
	chproc_given = decAccDensityCh{idx_density};
	decode_givench_v73_new
	decAccDensityResp(1,idx_density) = DecodingAccuracyResp.class_acc_given_broad_mean(end);
	decAccDensityResp(2,idx_density) = DecodingAccuracyResp.class_acc_given_theta_mean(end);
	decAccDensityResp(3,idx_density) = DecodingAccuracyResp.class_acc_given_alpha_mean(end);
	decAccDensityResp(4,idx_density) = DecodingAccuracyResp.class_acc_given_beta_mean(end);
	decAccDensityResp(5,idx_density) = DecodingAccuracyResp.class_acc_given_gamma_mean(end);
	decAccDensityResp(6,idx_density) = DecodingAccuracyResp.class_acc_given_hfo_mean(end);
end

% Plot the result
figure
plot(decAccDensityResp','o-','LineWidth',4)
xlim([.9 3.1]); ylim([50 100])
xticks([1 2 3])
xticklabels({'Low','Mid','High'})
xlabel('Channel Density'); ylabel('C.V. Accuracy'); set(gca,'FontSize',20)
legend on
legend({'Broad','Theta','Alpha','Beta','Gamma','HFO'},'FontSize',16)
legend('boxoff')
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),560,500]);

