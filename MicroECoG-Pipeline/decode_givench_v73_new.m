%%============================== Decoding neural signal =======================
tic
%% For classification using SVM model

%decodeWhat = input('Which do you want to decode between (r)esponse and (t)nr? ','s');

% Get SVM models using increasing number of channels.
optLibSVM = 0;
itermax = 1;
numch = length(chproc_given);
%numchpow = [2.^PN(2.^PN<numch) numch];
numchbstr = numch;


% Define arrays for accuracy
class_acc_given_broad = nan(itermax,numel(numchbstr));
class_acc_given_broad_rnd = nan(itermax,numel(numchbstr));
%%if bandProc == 'm'
%	class_acc_given_theta = nan(itermax,numel(numchbstr));
%	class_acc_given_theta_rnd = nan(itermax,numel(numchbstr));
%	class_acc_given_alpha = nan(itermax,numel(numchbstr));
%	class_acc_given_alpha_rnd = nan(itermax,numel(numchbstr));
%	class_acc_given_beta = nan(itermax,numel(numchbstr));
%	class_acc_given_beta_rnd = nan(itermax,numel(numchbstr));
%	class_acc_given_gamma = nan(itermax,numel(numchbstr));
%	class_acc_given_gamma_rnd = nan(itermax,numel(numchbstr));
%	class_acc_given_hfo = nan(itermax,numel(numchbstr));
%	class_acc_given_hfo_rnd = nan(itermax,numel(numchbstr));
%%end

% Set SVM optimization options.
% (Use default options for linear kernel.)
%libsvmopt = '-s 0 -t 0 -c 10 -v 5 -q';
% RBF kernel with 5-fold crossvalidation.
%libsvmopt = '-t 2 -c 1e4 -g 1e4 -v 5 -q';


if decodeWhat == 'r'
	condIndUB = condIndUB_resp;
	classLabel = 9;
else
	condIndUB = condIndUB_tnr;
	classLabel = 10;
end

%for j = 1:numel(numchbstr)
for j = numel(numchbstr)
	fprintf('\nNo. of Channels = %i\n',numchbstr(j))
	for iternum = 1:itermax
		%chpick = sort(chproc_given(randi(numch,[numchbstr(j),1])));
		chpick = sort(chproc_given);
		% generate a label vector for the given condition
		labels_2d = repmat(cell2mat(behavSummary(condIndUB{end},classLabel)),numel(chpick),1);
		% construct 2D feature array
		features_broad_2d = [];
		%if bandProcPar{iternum} == 'm'
			features_theta_2d = [];
			features_alpha_2d = [];
			features_beta_2d = [];
			features_gamma_2d = [];
			features_hfo_2d = [];
		%end

		% Construct feature matrices and normalized it.
		for i=1:numel(chpick)
			features_broad_2d = [features_broad_2d; chstat_broad(condIndUB{end},1:2,chpick(i))];
			features_broad_2d_scaled = (features_broad_2d - repmat(min(features_broad_2d,[],1),size(features_broad_2d,1),1))*spdiags(1./(max(features_broad_2d,[],1)-min(features_broad_2d,[],1))',0,size(features_broad_2d,2),size(features_broad_2d,2));
			%if bandProcPar{iternum} == 'm'
				features_theta_2d = [features_theta_2d; chstat_theta(condIndUB{end},1:2,chpick(i))];
				features_theta_2d_scaled = (features_theta_2d - repmat(min(features_theta_2d,[],1),size(features_theta_2d,1),1))*spdiags(1./(max(features_theta_2d,[],1)-min(features_theta_2d,[],1))',0,size(features_theta_2d,2),size(features_theta_2d,2));
				features_alpha_2d = [features_alpha_2d; chstat_alpha(condIndUB{end},1:2,chpick(i))];
				features_alpha_2d_scaled = (features_alpha_2d - repmat(min(features_alpha_2d,[],1),size(features_alpha_2d,1),1))*spdiags(1./(max(features_alpha_2d,[],1)-min(features_alpha_2d,[],1))',0,size(features_alpha_2d,2),size(features_alpha_2d,2));
				features_beta_2d = [features_beta_2d; chstat_beta(condIndUB{end},1:2,chpick(i))];
				features_beta_2d_scaled = (features_beta_2d - repmat(min(features_beta_2d,[],1),size(features_beta_2d,1),1))*spdiags(1./(max(features_beta_2d,[],1)-min(features_beta_2d,[],1))',0,size(features_beta_2d,2),size(features_beta_2d,2));
				features_gamma_2d = [features_gamma_2d; chstat_gamma(condIndUB{end},1:2,chpick(i))];
				features_gamma_2d_scaled = (features_gamma_2d - repmat(min(features_gamma_2d,[],1),size(features_gamma_2d,1),1))*spdiags(1./(max(features_gamma_2d,[],1)-min(features_gamma_2d,[],1))',0,size(features_gamma_2d,2),size(features_gamma_2d,2));
				features_hfo_2d = [features_hfo_2d; chstat_hfo(condIndUB{end},1:2,chpick(i))];
				features_hfo_2d_scaled = (features_hfo_2d - repmat(min(features_hfo_2d,[],1),size(features_hfo_2d,1),1))*spdiags(1./(max(features_hfo_2d,[],1)-min(features_hfo_2d,[],1))',0,size(features_hfo_2d,2),size(features_hfo_2d,2));
			%end
		end
		
		%if iternum == 1
			% Calculate classification accuracy with cross-validation
			accCV_max = zeros(10,1);
			accCV_max_par = cell(10,1);
			parfor r=2:10
				for q=1:5
					for p=1:5
						accCV_curr = libsvmtrain(labels_2d,features_broad_2d_scaled,sprintf('-t 2 -c %d -g %d -v %i -q',10^(p),10^(q),r));
						if accCV_curr > accCV_max(r)
							accCV_max(r) = accCV_curr;
							accCV_max_par{r} = [p q r];
						end
					end
				end
			end
			class_acc_given_broad(iternum,j,k) = max(accCV_max(:));
			idx_opt = accCV_max_par(find(accCV_max==max(accCV_max(:))));
			libsvmopt_par(iternum) = {sprintf('-t 2 -c %d -g %d -v %i -q',10^(idx_opt{1}(1)),10^(idx_opt{1}(2)),idx_opt{1}(3))};
			display('Optimization done!')
		%end

%		fprintf('\nNo. of Channels = %i\n',numchbstr(j))
%		fprintf('\nIteration #%i\n',iternum)
		class_acc_given_broad(iternum,j,k) = libsvmtrain(labels_2d,features_broad_2d_scaled,libsvmopt_par{iternum});
%		class_acc_given_broad_rnd(iternum,j,k) = libsvmtrain(labels_rnd_2d,features_broad_2d_scaled,libsvmopt_par{iternum});
		%if bandProcPar{iternum} == 'm'
			class_acc_given_theta(iternum,j,k) = libsvmtrain(labels_2d,features_theta_2d_scaled,libsvmopt_par{iternum});
%			class_acc_given_theta_rnd(iternum,j,k) = libsvmtrain(labels_rnd_2d,features_theta_2d_scaled,libsvmopt_par{iternum});
			class_acc_given_alpha(iternum,j,k) = libsvmtrain(labels_2d,features_alpha_2d_scaled,libsvmopt_par{iternum});
%			class_acc_given_alpha_rnd(iternum,j,k) = libsvmtrain(labels_rnd_2d,features_alpha_2d_scaled,libsvmopt_par{iternum});
			class_acc_given_beta(iternum,j,k) = libsvmtrain(labels_2d,features_beta_2d_scaled,libsvmopt_par{iternum});
%			class_acc_given_beta_rnd(iternum,j,k) = libsvmtrain(labels_rnd_2d,features_beta_2d_scaled,libsvmopt_par{iternum});
			class_acc_given_gamma(iternum,j,k) = libsvmtrain(labels_2d,features_gamma_2d_scaled,libsvmopt_par{iternum});
%			class_acc_given_gamma_rnd(iternum,j,k) = libsvmtrain(labels_rnd_2d,features_gamma_2d_scaled,libsvmopt_par{iternum});
			class_acc_given_hfo(iternum,j,k) = libsvmtrain(labels_2d,features_hfo_2d_scaled,libsvmopt_par{iternum});
%			class_acc_given_hfo_rnd(iternum,j,k) = libsvmtrain(labels_rnd_2d,features_hfo_2d_scaled,libsvmopt_par{iternum});
		%end
	end
	class_acc_given_broad_mean(j,k) = nanmean(class_acc_given_broad(:,j,k),1);
	class_acc_given_broad_std(j,k) = nanstd(class_acc_given_broad(:,j,k),1);
	class_acc_given_broad_sem(j,k) = nanstd(class_acc_given_broad(:,j,k),1)/sqrt(itermax);
%	class_acc_given_broad_rnd_mean(j,k) = nanmean(class_acc_given_broad_rnd(:,j,k),1);
%	class_acc_given_broad_rnd_std(j,k) = nanstd(class_acc_given_broad_rnd(:,j,k),1);
%	class_acc_given_broad_rnd_sem(j,k) = nanstd(class_acc_given_broad_rnd(:,j,k),1)/sqrt(itermax);
	%if bandProc == 'm'
		class_acc_given_theta_mean(j,k) = nanmean(class_acc_given_theta(:,j,k),1);
		class_acc_given_theta_std(j,k) = nanstd(class_acc_given_theta(:,j,k),1);
		class_acc_given_theta_sem(j,k) = nanstd(class_acc_given_theta(:,j,k),1)/sqrt(itermax);
%		class_acc_given_theta_rnd_mean(j,k) = nanmean(class_acc_given_theta_rnd(:,j,k),1);
%		class_acc_given_theta_rnd_std(j,k) = nanstd(class_acc_given_theta_rnd(:,j,k),1);
%		class_acc_given_theta_rnd_sem(j,k) = nanstd(class_acc_given_theta_rnd(:,j,k),1)/sqrt(itermax);
		class_acc_given_alpha_mean(j,k) = nanmean(class_acc_given_alpha(:,j,k),1);
		class_acc_given_alpha_std(j,k) = nanstd(class_acc_given_alpha(:,j,k),1);
		class_acc_given_alpha_sem(j,k) = nanstd(class_acc_given_alpha(:,j,k),1)/sqrt(itermax);
%		class_acc_given_alpha_rnd_mean(j,k) = nanmean(class_acc_given_alpha_rnd(:,j,k),1);
%		class_acc_given_alpha_rnd_std(j,k) = nanstd(class_acc_given_alpha_rnd(:,j,k),1);
%		class_acc_given_alpha_rnd_sem(j,k) = nanstd(class_acc_given_alpha_rnd(:,j,k),1)/sqrt(itermax);
		class_acc_given_beta_mean(j,k) = nanmean(class_acc_given_beta(:,j,k),1);
		class_acc_given_beta_std(j,k) = nanstd(class_acc_given_beta(:,j,k),1);
		class_acc_given_beta_sem(j,k) = nanstd(class_acc_given_beta(:,j,k),1)/sqrt(itermax);
%		class_acc_given_beta_rnd_mean(j,k) = nanmean(class_acc_given_beta_rnd(:,j,k),1);
%		class_acc_given_beta_rnd_std(j,k) = nanstd(class_acc_given_beta_rnd(:,j,k),1);
%		class_acc_given_beta_rnd_sem(j,k) = nanstd(class_acc_given_beta_rnd(:,j,k),1)/sqrt(itermax);
		class_acc_given_gamma_mean(j,k) = nanmean(class_acc_given_gamma(:,j,k),1);
		class_acc_given_gamma_std(j,k) = nanstd(class_acc_given_gamma(:,j,k),1);
		class_acc_given_gamma_sem(j,k) = nanstd(class_acc_given_gamma(:,j,k),1)/sqrt(itermax);
%		class_acc_given_gamma_rnd_mean(j,k) = nanmean(class_acc_given_gamma_rnd(:,j,k),1);
%		class_acc_given_gamma_rnd_std(j,k) = nanstd(class_acc_given_gamma_rnd(:,j,k),1);
%		class_acc_given_gamma_rnd_sem(j,k) = nanstd(class_acc_given_gamma_rnd(:,j,k),1)/sqrt(itermax);
		class_acc_given_hfo_mean(j,k) = nanmean(class_acc_given_hfo(:,j,k),1);
		class_acc_given_hfo_std(j,k) = nanstd(class_acc_given_hfo(:,j,k),1);
		class_acc_given_hfo_sem(j,k) = nanstd(class_acc_given_hfo(:,j,k),1)/sqrt(itermax);
%		class_acc_given_hfo_rnd_mean(j,k) = nanmean(class_acc_given_hfo_rnd(:,j,k),1);
%		class_acc_given_hfo_rnd_std(j,k) = nanstd(class_acc_given_hfo_rnd(:,j,k),1);
%		class_acc_given_hfo_rnd_sem(j,k) = nanstd(class_acc_given_hfo_rnd(:,j,k),1)/sqrt(itermax);
	%end
end

% Save the decoding accuracy in a structure.
decAccVar = who('class_acc_given_*');
for i=1:length(decAccVar)
	if decodeWhat == 'r'
		DecodingAccuracyResp.(decAccVar{i}) = eval(decAccVar{i});
	else
		DecodingAccuracyTNR.(decAccVar{i}) = eval(decAccVar{i});
	end
end
clearvars class_acc_select*
%clearvars bandProcPar libsvmopt_par class_acc_select*


display('Decoding done.')

toc
