%%============================== Decoding neural signal =======================
decodeWhat = input('Which do you want to decode between (r)esponse and (t)nr? ','s');
%bandProc = input('Broad band only or multiband?(b/m) ','s');
bandProc = 'm';
tic

%% For classification using SVM model

% Get SVM models using increasing number of channels.
numch = length(chproc);
numchpow = [2.^PN(2.^PN<numch) numch];

% Build a table for all classes.
classInd = cell(length(tnrLabel),length(respLabel));
className = cell(length(tnrLabel),length(respLabel));
classNum = 0;
for j=1:length(respLabel)
    for i=1:length(tnrLabel)
        classInd{i,j} = find(strcmp(behavSummary(:,4),num2str(tnrLabel{i})).*strcmp(behavSummary(:,5),respLabel(j)));
        className{i,j} = {strcat(num2str(tnrLabel{i}),respLabel{j})};
        behavSummary(classInd{i,j},7) = {strcat(num2str(tnrLabel{i}),respLabel{j})};
        classNum = classNum+1;
        behavSummary(classInd{i,j},8) = {classNum};
    end
end

% number of classes
numClass_tnr = numel(unique(behavSummary(:,4)));
numClass_resp = numel(unique(behavSummary(:,5)));

% Convert label to integers
for i=1:numClass_resp
    behavSummary(trind{i},9) = {i};
end
for i=1:numClass_tnr
    behavSummary(tnrind{i},10) = {i};
end

%------------------------------------------------------------------------------
% Find trial indices unbiased for all the conditions that has at least two nonempty classes.
condIndUB_resp = cell(length(tnrLabel)+1,1);	% +1 for across all TNR
for i=1:numel(tnrLabel)
    if sum(~cellfun(@isempty,classInd(i,:))) >= 2
        minoccur = min(cellfun(@length,classInd(i,~cellfun(@isempty,classInd(i,:)))));
        if minoccur >= 2
            condIndUB_resp{i} = reshape(cell2mat(cellfun(@(x) x(1:minoccur),classInd(i,~cellfun(@isempty,classInd(i,:))),'UniformOutput',false)),[],1);
        end
    end
end

condIndUB_tnr = cell(length(respLabel)+1,1);	% +1 for all responses
j=1;		% only pick hit trials.
if sum(~cellfun(@isempty,classInd(:,j))) >= 2
    minoccur = min(cellfun(@length,classInd(~cellfun(@isempty,classInd(:,j)))));
    if minoccur >= 2
        condIndUB_tnr{j} = reshape(cell2mat(cellfun(@(x) x(1:minoccur),classInd(~cellfun(@isempty,classInd(:,j))),'UniformOutput',false)),[],1);
    end
end

% Collapse class indices across all TNR and create condition indices.
for j=1:length(respLabel)
    classIndAllTNR{1,j} = vertcat(classInd{:,j});
end
if sum(~cellfun(@isempty,classIndAllTNR)) >= 2	% check if at least two different responses to classify
	minoccur = min(cellfun(@length,classIndAllTNR(~cellfun(@isempty,classIndAllTNR))));
    condIndUB_resp{end} = reshape(cell2mat(cellfun(@(x) x(1:minoccur),classIndAllTNR(~cellfun(@isempty,classIndAllTNR)),'UniformOutput',false)),[],1);
end
% Collapse class indices across all responses and create condition indices.
for i=1:length(tnrLabel)
    classIndAllResp{i,1} = vertcat(classInd{i,:});
end
if sum(~cellfun(@isempty,classIndAllResp)) >= 2	% check if at least two different responses to classify
	minoccur = min(cellfun(@length,classIndAllResp(~cellfun(@isempty,classIndAllResp))));
    condIndUB_tnr{end} = reshape(cell2mat(cellfun(@(x) x(1:minoccur),classIndAllResp(~cellfun(@isempty,classIndAllResp)),'UniformOutput',false)),[],1);
end


% Remove null cell elements  (* Commented out to retain the shape)
%condIndUB_resp = condIndUB_resp(find(~cellfun(@isempty,condIndUB_resp)));
%condIndUB_tnr = condIndUB_tnr(find(~cellfun(@isempty,condIndUB_tnr)));

% Find all classes in each condition.
numClassUB_resp = cell(length(condIndUB_resp),1);
numClassUB_tnr = cell(length(condIndUB_tnr),1);

% Assign class label to use depending on what to decode.
if decodeWhat == 'r';	% decode responses
    classLabel = 9;
else					% decode TNRs
    classLabel = 10;
end
%------------------------------------------------------------------------------


% Define accuracy arrays

if decodeWhat == 'r'	% Decode responses.
	class_acc_single_resp_broad = nan(256,1);
	if bandProc == 'm'
		class_acc_single_resp_theta = nan(256,1);
		class_acc_single_resp_alpha = nan(256,1);
		class_acc_single_resp_beta = nan(256,1);
		class_acc_single_resp_gamma = nan(256,1);
		class_acc_single_resp_hfo = nan(256,1);
	end
	
	% Set SVM optimization options.
	%libsvmopt = '-t 2 -c 1000 -g 1000 -v 10 -q';
	
	%% Calculate decoding accuracy
	display('Starting classification...')
	% Change k depending on what you want to decode
	k = length(condIndUB_resp);
	
	parfor j=1:numch
		% generate a label vector for the given condition
		labels_2d = cell2mat(behavSummary(condIndUB_resp{k},classLabel));
		% construct 2D feature array
		features_broad_2d = [];
		if bandProc == 'm'
			features_theta_2d = [];
			features_alpha_2d = [];
			features_beta_2d = [];
			features_gamma_2d = [];
			features_hfo_2d = [];
		end
	
		features_broad_2d = chstat_broad(condIndUB_resp{k},1:2,chproc(j));
		features_broad_2d_scaled = (features_broad_2d - repmat(min(features_broad_2d,[],1),size(features_broad_2d,1),1))*spdiags(1./(max(features_broad_2d,[],1)-min(features_broad_2d,[],1))',0,size(features_broad_2d,2),size(features_broad_2d,2));
		if bandProc == 'm'
			features_theta_2d = chstat_theta(condIndUB_resp{k},1:2,chproc(j));
			features_theta_2d_scaled = (features_theta_2d - repmat(min(features_theta_2d,[],1),size(features_theta_2d,1),1))*spdiags(1./(max(features_theta_2d,[],1)-min(features_theta_2d,[],1))',0,size(features_theta_2d,2),size(features_theta_2d,2));
			features_alpha_2d = chstat_alpha(condIndUB_resp{k},1:2,chproc(j));
			features_alpha_2d_scaled = (features_alpha_2d - repmat(min(features_alpha_2d,[],1),size(features_alpha_2d,1),1))*spdiags(1./(max(features_alpha_2d,[],1)-min(features_alpha_2d,[],1))',0,size(features_alpha_2d,2),size(features_alpha_2d,2));
			features_beta_2d = chstat_beta(condIndUB_resp{k},1:2,chproc(j));
			features_beta_2d_scaled = (features_beta_2d - repmat(min(features_beta_2d,[],1),size(features_beta_2d,1),1))*spdiags(1./(max(features_beta_2d,[],1)-min(features_beta_2d,[],1))',0,size(features_beta_2d,2),size(features_beta_2d,2));
			features_gamma_2d = chstat_gamma(condIndUB_resp{k},1:2,chproc(j));
			features_gamma_2d_scaled = (features_gamma_2d - repmat(min(features_gamma_2d,[],1),size(features_gamma_2d,1),1))*spdiags(1./(max(features_gamma_2d,[],1)-min(features_gamma_2d,[],1))',0,size(features_gamma_2d,2),size(features_gamma_2d,2));
			features_hfo_2d = chstat_hfo(condIndUB_resp{k},1:2,chproc(j));
			features_hfo_2d_scaled = (features_hfo_2d - repmat(min(features_hfo_2d,[],1),size(features_hfo_2d,1),1))*spdiags(1./(max(features_hfo_2d,[],1)-min(features_hfo_2d,[],1))',0,size(features_hfo_2d,2),size(features_hfo_2d,2));
		end
	
	
		% Calculate classification accuracy with cross-validation
		class_acc_single_resp_max = zeros(10,1);
	    class_acc_single_resp_max_param = cell(10,1);
	    for r=2:10
	        for q=1:5
	            for p=1:5
	                class_acc_single_resp_curr = libsvmtrain(labels_2d,features_broad_2d_scaled,sprintf('-t 2 -c %d -g %d -v %i -q',10^(p),10^(q),r));
	                if class_acc_single_resp_curr > class_acc_single_resp_max(r)
	                    class_acc_single_resp_max(r) = class_acc_single_resp_curr;
	                    class_acc_single_resp_max_param{r} = [p q r];
	                end
	            end
	        end
	    end
		class_acc_single_resp_broad_par(j,1) = max(class_acc_single_resp_max(:));
	    idx_opt = class_acc_single_resp_max_param(find(class_acc_single_resp_max==max(class_acc_single_resp_max(:))));
	    libsvmopt = sprintf('-t 2 -c %d -g %d -v %i -q',10^(idx_opt{1}(1)),10^(idx_opt{1}(2)),idx_opt{1}(3))
	    display('Optimization done!')
		if bandProc == 'm'
			class_acc_single_resp_theta_par(j,1) = libsvmtrain(labels_2d,features_theta_2d_scaled,libsvmopt);
			class_acc_single_resp_alpha_par(j,1) = libsvmtrain(labels_2d,features_alpha_2d_scaled,libsvmopt);
			class_acc_single_resp_beta_par(j,1) = libsvmtrain(labels_2d,features_beta_2d_scaled,libsvmopt);
			class_acc_single_resp_gamma_par(j,1) = libsvmtrain(labels_2d,features_gamma_2d_scaled,libsvmopt);
			class_acc_single_resp_hfo_par(j,1) = libsvmtrain(labels_2d,features_hfo_2d_scaled,libsvmopt);
		end
	end
	
	for j=1:numch
		class_acc_single_resp_broad(chproc(j)) = class_acc_single_resp_broad_par(j);
		if bandProc == 'm'
			class_acc_single_resp_theta(chproc(j)) = class_acc_single_resp_theta_par(j);
			class_acc_single_resp_alpha(chproc(j)) = class_acc_single_resp_alpha_par(j);
			class_acc_single_resp_beta(chproc(j)) = class_acc_single_resp_beta_par(j);
			class_acc_single_resp_gamma(chproc(j)) = class_acc_single_resp_gamma_par(j);
			class_acc_single_resp_hfo(chproc(j)) = class_acc_single_resp_hfo_par(j);
		end
	end
	
	clearvars class_acc_single_resp_*_par

else	% Decode TNRs
	class_acc_single_tnr_broad = nan(256,1);
	if bandProc == 'm'
		class_acc_single_tnr_theta = nan(256,1);
		class_acc_single_tnr_alpha = nan(256,1);
		class_acc_single_tnr_beta = nan(256,1);
		class_acc_single_tnr_gamma = nan(256,1);
		class_acc_single_tnr_hfo = nan(256,1);
	end
	
	% Set SVM optimization options.
	%libsvmopt = '-t 2 -c 1000 -g 1000 -v 10 -q';
	
	%% Calculate decoding accuracy
	display('Starting classification...')
	% Change k depending on what you want to decode
	k = length(condIndUB_tnr);
	
	parfor j=1:numch
		% generate a label vector for the given condition
		labels_2d = cell2mat(behavSummary(condIndUB_tnr{k},classLabel));
		% construct 2D feature array
		features_broad_2d = [];
		if bandProc == 'm'
			features_theta_2d = [];
			features_alpha_2d = [];
			features_beta_2d = [];
			features_gamma_2d = [];
			features_hfo_2d = [];
		end
	
		features_broad_2d = chstat_broad(condIndUB_tnr{k},1:2,chproc(j));
		features_broad_2d_scaled = (features_broad_2d - repmat(min(features_broad_2d,[],1),size(features_broad_2d,1),1))*spdiags(1./(max(features_broad_2d,[],1)-min(features_broad_2d,[],1))',0,size(features_broad_2d,2),size(features_broad_2d,2));
		if bandProc == 'm'
			features_theta_2d = chstat_theta(condIndUB_tnr{k},1:2,chproc(j));
			features_theta_2d_scaled = (features_theta_2d - repmat(min(features_theta_2d,[],1),size(features_theta_2d,1),1))*spdiags(1./(max(features_theta_2d,[],1)-min(features_theta_2d,[],1))',0,size(features_theta_2d,2),size(features_theta_2d,2));
			features_alpha_2d = chstat_alpha(condIndUB_tnr{k},1:2,chproc(j));
			features_alpha_2d_scaled = (features_alpha_2d - repmat(min(features_alpha_2d,[],1),size(features_alpha_2d,1),1))*spdiags(1./(max(features_alpha_2d,[],1)-min(features_alpha_2d,[],1))',0,size(features_alpha_2d,2),size(features_alpha_2d,2));
			features_beta_2d = chstat_beta(condIndUB_tnr{k},1:2,chproc(j));
			features_beta_2d_scaled = (features_beta_2d - repmat(min(features_beta_2d,[],1),size(features_beta_2d,1),1))*spdiags(1./(max(features_beta_2d,[],1)-min(features_beta_2d,[],1))',0,size(features_beta_2d,2),size(features_beta_2d,2));
			features_gamma_2d = chstat_gamma(condIndUB_tnr{k},1:2,chproc(j));
			features_gamma_2d_scaled = (features_gamma_2d - repmat(min(features_gamma_2d,[],1),size(features_gamma_2d,1),1))*spdiags(1./(max(features_gamma_2d,[],1)-min(features_gamma_2d,[],1))',0,size(features_gamma_2d,2),size(features_gamma_2d,2));
			features_hfo_2d = chstat_hfo(condIndUB_tnr{k},1:2,chproc(j));
			features_hfo_2d_scaled = (features_hfo_2d - repmat(min(features_hfo_2d,[],1),size(features_hfo_2d,1),1))*spdiags(1./(max(features_hfo_2d,[],1)-min(features_hfo_2d,[],1))',0,size(features_hfo_2d,2),size(features_hfo_2d,2));
		end
	
	
		% Calculate classification accuracy with cross-validation
		class_acc_single_tnr_max = zeros(10,1);
	    class_acc_single_tnr_max_param = cell(10,1);
	    for r=2:10
	        for q=1:5
	            for p=1:5
	                class_acc_single_tnr_curr = libsvmtrain(labels_2d,features_broad_2d_scaled,sprintf('-t 2 -c %d -g %d -v %i -q',10^(p),10^(q),r));
	                if class_acc_single_tnr_curr > class_acc_single_tnr_max(r)
	                    class_acc_single_tnr_max(r) = class_acc_single_tnr_curr;
	                    class_acc_single_tnr_max_param{r} = [p q r];
	                end
	            end
	        end
	    end
		class_acc_single_tnr_broad_par(j,1) = max(class_acc_single_tnr_max(:));
	    idx_opt = class_acc_single_tnr_max_param(find(class_acc_single_tnr_max==max(class_acc_single_tnr_max(:))));
	    libsvmopt = sprintf('-t 2 -c %d -g %d -v %i -q',10^(idx_opt{1}(1)),10^(idx_opt{1}(2)),idx_opt{1}(3))
	    display('Optimization done!')
		if bandProc == 'm'
			class_acc_single_tnr_theta_par(j,1) = libsvmtrain(labels_2d,features_theta_2d_scaled,libsvmopt);
			class_acc_single_tnr_alpha_par(j,1) = libsvmtrain(labels_2d,features_alpha_2d_scaled,libsvmopt);
			class_acc_single_tnr_beta_par(j,1) = libsvmtrain(labels_2d,features_beta_2d_scaled,libsvmopt);
			class_acc_single_tnr_gamma_par(j,1) = libsvmtrain(labels_2d,features_gamma_2d_scaled,libsvmopt);
			class_acc_single_tnr_hfo_par(j,1) = libsvmtrain(labels_2d,features_hfo_2d_scaled,libsvmopt);
		end
	end
	
	for j=1:numch
		class_acc_single_tnr_broad(chproc(j)) = class_acc_single_tnr_broad_par(j);
		if bandProc == 'm'
			class_acc_single_tnr_theta(chproc(j)) = class_acc_single_tnr_theta_par(j);
			class_acc_single_tnr_alpha(chproc(j)) = class_acc_single_tnr_alpha_par(j);
			class_acc_single_tnr_beta(chproc(j)) = class_acc_single_tnr_beta_par(j);
			class_acc_single_tnr_gamma(chproc(j)) = class_acc_single_tnr_gamma_par(j);
			class_acc_single_tnr_hfo(chproc(j)) = class_acc_single_tnr_hfo_par(j);
		end
	end
	
	clearvars class_acc_single_tnr_*_par

end

display('Single channel decoding done.')

toc
