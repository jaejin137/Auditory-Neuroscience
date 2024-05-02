m = input('Enter m: ');

% First split core, belt, supra-, and infra-layers.
split_htbl

for k = 1:14
%k = input('Enter k: ');
	tmparr_core_supra = htbl_core_supra(:,:,k);
	tmparr_core_infra = htbl_core_infra(:,:,k);
	tmparr_belt_supra = htbl_belt_supra(:,:,k);
	tmparr_belt_infra = htbl_belt_infra(:,:,k);
	core_supra = tmparr_core_supra(find(~isnan(htbl_core_supra(:,:,k))));
	core_infra = tmparr_core_infra(find(~isnan(htbl_core_infra(:,:,k))));
	belt_supra = tmparr_belt_supra(find(~isnan(htbl_belt_supra(:,:,k))));
	belt_infra = tmparr_belt_infra(find(~isnan(htbl_belt_infra(:,:,k))));
	cc_core = corrcoef(core_supra,core_infra);
	CompMetric(k).cc_core = cc_core(1,2);
	cc_belt = corrcoef(belt_supra,belt_infra);
	CompMetric(k).cc_belt = cc_belt(1,2);
	
	cv_core_supra = std(core_supra)/mean(core_supra);
	cv_core_infra = std(core_infra)/mean(core_infra);
	CompMetric(k).alpha_core = cv_core_supra/cv_core_infra;
	cv_belt_supra = std(belt_supra)/mean(belt_supra);
	cv_belt_infra = std(belt_infra)/mean(belt_infra);
	CompMetric(k).alpha_belt = cv_belt_supra/cv_belt_infra;
	
	core_supra_zscore = zscore(core_supra);
	core_infra_zscore = zscore(core_infra);
	bins = 10;
	[N_core_supra,~] = histcounts(core_supra_zscore,bins);
	[N_core_infra,~] = histcounts(core_infra_zscore,bins);
	minOfHists_core = min([N_core_infra; N_core_supra], [], 1);
	overlappedHist_core = sum(minOfHists_core);
	histogram_match_core = overlappedHist_core/sum(N_core_infra);
	%spaef_core(k,6) = 1- sqrt( (CompMetric(k).cc_core-1)^2 + (CompMetric(k).alpha_core-1)^2 + (histogram_match_core-1)^2 );
	CompMetric(k).spaef_core = 1- sqrt( (CompMetric(k).cc_core-1)^2 + (CompMetric(k).alpha_core-1)^2 );

	belt_supra_zscore = zscore(belt_supra);
	belt_infra_zscore = zscore(belt_infra);
	bins = 10;
	[N_belt_supra,~] = histcounts(belt_supra_zscore,bins);
	[N_belt_infra,~] = histcounts(belt_infra_zscore,bins);
	minOfHists_belt = min([N_belt_infra; N_belt_supra], [], 1);
	overlappedHist_belt = sum(minOfHists_belt);
	histogram_match_belt = overlappedHist_core/sum(N_belt_infra);
	%spaef_belt(k,6) = 1- sqrt( (CompMetric(k).cc_belt-1)^2 + (CompMetric(k).alpha_belt-1)^2 + (histogram_match_belt-1)^2 );
	CompMetric(k).spaef_belt = 1- sqrt( (CompMetric(k).cc_belt-1)^2 + (CompMetric(k).alpha_belt-1)^2 );
end
if m == 1
	CompMetric_Cassius = CompMetric;
else
	CompMetric_Miyagi = CompMetric;
end

cc_comb = [[CompMetric_Cassius.cc_core]' [CompMetric_Cassius.cc_belt]']
alpha_comb = [[CompMetric_Cassius.alpha_core]' [CompMetric_Cassius.alpha_belt]']
spaef_comb = [[CompMetric_Cassius.spaef_core]' [CompMetric_Cassius.spaef_belt]']

XTickLabels = {'Ratio','p_med','BF','BFHz','Delay','PD','t10','t50','Dur','PLI','BW','BWHz','bSMF','FmUpperCutoff'}
figure
subplot(3,1,1)
bar(cc_comb)
set(gca,'XtickLabel',XTickLabels,'FontSize',20)
xtickangle(-45)
title('Correlation Coefficient')

subplot(3,1,2)
bar(alpha_comb)
set(gca,'XtickLabel',XTickLabels,'FontSize',20)
xtickangle(-45)
title('Coefficient of Variance')

subplot(3,1,3)
bar(spaef_comb)
set(gca,'XtickLabel',XTickLabels,'FontSize',20)
xtickangle(-45)
title('Space Efficiency')
currPos = get(gcf,'Position'); set(gcf,'Position',[currPos(1),currPos(2),800,600])
