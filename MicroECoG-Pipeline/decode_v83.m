%%============================== Decoding neural signal =======================
% Use fake data?
useFake = input('Use fake data?(y/n) ','s');
% Equalize sizes of classes?
if useFake == 'y'
	equalFake = input('Equalize class size of fake data?(y/n) ','s');
	varyingNoise = input('Vary noise level?(y/n) ','s');
	if varyingNoise == 'n'
		nc = input('Enter noise coefficient: ');
	elseif varyingNoise == 'y'
		contJumpNoise = input('Continuous or jumping noise?(c/j) ','s');
	end
else
	% What band to process?
	bandProc = input('Broad band only or multiband?(b/m) ','s');
	useFeature = input('Use means or prominent peaks? (m/p) ','s');
end
% Sampling channels
sampleChannel = input('Random or linear sampling of channels?(r/l) ','s');
% Use ADASYN?
useADASYN = input('Use SMOTE?(y/n) ','s');
runParallel = input('Run in parallel mode?(y/n) ','s');

tic
%% Classification using SVM model with SMOTE (ADASYN implementation)
%%------------------------------- Bootstrapping -------------------------------




%% Preprocess data to conform libsvm format.
numClass = numel(unique(behavSummary(:,5)));
% Convert label to integers
for i=1:numClass
	behavSummary(trind{i},9) = {i};
end

%%=============== Setup fake data by manipulating indices. ====================
if useFake == 'y'
	display('Constructing a fake data ...')
	behavSummary_fake = behavSummary;
	trind_fake = trind;
	if equalFake == 'y'
		display('Equalizing class sizes ...')
		% Equalize the data sizes.
		numEqualInstances = floor(numtr/numClass);
		for i=1:numClass
			if i~=numClass
				trind_fake(i) = {[numEqualInstances*(i-1)+1:numEqualInstances*i]'};
			else
				trind_fake(i) = {[numEqualInstances*(i-1)+1:numtr]'};
			end
		end		
		for i=1:numClass
			behavSummary_fake(trind_fake{i},9) = {i};
		end
	end
		
	% Construct fake data
	chstat_fake = chstat_broad;
	for k=1:256
		% noise coeffient
		if varyingNoise == 'y'
			if contJumpNoise == 'c'
				nc = .01*k;
			elseif contJumpNoise == 'j'
				if k<=32
					nc = 1;
				else
					nc = 2;
				end
			end
		end
		% fake data
		for i=1:numClass
			chstat_fake(find(cell2mat(behavSummary_fake(:,9))==i),1,k) = (i-1)+nc*rand(length(trind_fake{i}),1);
			chstat_fake(find(cell2mat(behavSummary_fake(:,9))==i),2,k) = (i-1)+nc*rand(length(trind_fake{i}),1);
		end
	end
	% plot fake data
	display('Plotting fake data ...')
	figure(90)
	clf(90)
	for k=1:4
		for i=1:numClass
			subplot(2,2,k)
			plot(chstat_fake(find(cell2mat(behavSummary_fake(:,9))==i),1,64*k),chstat_fake(find(cell2mat(behavSummary_fake(:,9))==i),2,64*k),'.')
			xlabel('Mean')
			ylabel('Variance')
			title(sprintf('Ch=%i',64*k))
			hold on
		end
		legend show
	end
	drawnow
	
	% feed fake data
	display('Feeding a fake data ...')
	behavSummary_feed = behavSummary_fake;
	chstat_feed = cell(1,1);
	chstat_feed{1,1} = chstat_fake;
elseif useFake == 'n'
	% feed real data
	display('Feeding the real data ...')
	behavSummary_feed = behavSummary;
	if bandProc == 'b'
		if useFeature == 'p'
			% Backup original chstat arrays
			chstat_broad_orig = chstat_broad;
			% Replace channel mean with number of peaks
    		for j=1:numch
    		    for i=1:numtr
    		        [pk_broad pkloc_broad] = findpeaks(chdata_broad_dns{i,chproc(j)},'minPeakProminence',5);
    		        chstat_broad_numpeaks(i,chproc(j)) = numel(pk_broad);
    		    end
    		end
			chstat_broad(:,1,:) = chstat_broad_numpeaks;
		end
		chstat_feed = cell(1,1);
		chstat_feed{1,1} = chstat_broad;
	elseif bandProc == 'm'
		if useFeature == 'p'
			display('Detecting peaks ...')
			% Backup original chstat arrays
			chstat_broad_orig = chstat_broad;
			chstat_theta_orig = chstat_theta;
			chstat_alpha_orig = chstat_alpha;
			chstat_beta_orig = chstat_beta;
			chstat_gamma_orig = chstat_gamma;
			chstat_hfo_orig = chstat_hfo;
			% Replace channel mean with number of peaks
    		for j=1:numch
    		    for i=1:numtr
    		        [pk_broad pkloc_broad] = findpeaks(chdata_broad_dns{i,chproc(j)},'minPeakProminence',5);
    		        chstat_broad_numpeaks(i,chproc(j)) = numel(pk_broad);
    		        [pk_theta pkloc_theta] = findpeaks(chdata_theta_dns{i,chproc(j)},'minPeakProminence',5);
    		        chstat_theta_numpeaks(i,chproc(j)) = numel(pk_theta);
    		        [pk_alpha pkloc_alpha] = findpeaks(chdata_alpha_dns{i,chproc(j)},'minPeakProminence',5);
    		        chstat_alpha_numpeaks(i,chproc(j)) = numel(pk_alpha);
    		        [pk_beta pkloc_beta] = findpeaks(chdata_beta_dns{i,chproc(j)},'minPeakProminence',5);
    		        chstat_beta_numpeaks(i,chproc(j)) = numel(pk_beta);
    		        [pk_gamma pkloc_gamma] = findpeaks(chdata_gamma_dns{i,chproc(j)},'minPeakProminence',5);
    		        chstat_gamma_numpeaks(i,chproc(j)) = numel(pk_gamma);
    		        [pk_hfo pkloc_hfo] = findpeaks(chdata_hfo_dns{i,chproc(j)},'minPeakProminence',5);
    		        chstat_hfo_numpeaks(i,chproc(j)) = numel(pk_hfo);
    		    end
    		end
			chstat_broad(:,1,:) = chstat_broad_numpeaks;
			chstat_theta(:,1,:) = chstat_theta_numpeaks;
			chstat_alpha(:,1,:) = chstat_alpha_numpeaks;
			chstat_beta(:,1,:) = chstat_beta_numpeaks;
			chstat_gamma(:,1,:) = chstat_gamma_numpeaks;
			chstat_hfo(:,1,:) = chstat_hfo_numpeaks;
		end
		chstat_feed = cell(6,1);
		chstat_feed{1,1} = chstat_broad;
		chstat_feed{2,1} = chstat_theta;
		chstat_feed{3,1} = chstat_alpha;
		chstat_feed{4,1} = chstat_beta;
		chstat_feed{5,1} = chstat_gamma;
		chstat_feed{6,1} = chstat_hfo;
	else
		display('!! Cannot determine what band to process !!')
	end
		
else
	display('!! Cannot determine what data to use !!')
	return
end
%%=============================================================================

%% Setup bootstrap and SVM.

% Redefine channels to bootstrap over.
numch = length(chproc);
numchpow = [2.^PN(2.^PN<numch) numch];
numchbstr = numchpow;

% Get SVM models using increasing number of channels.
% Use optimization?
optLibSVM = 0;
% Set SVM optimization options.
% (Use default options for linear kernel.)
libsvmopt = '-s 0 -t 0 -c 1 -v 5 -q';
% Maximum number of iteration.
itermax = 100;

% Setup ADASYN parameters
if useADASYN == 'y'
	ADASYN_beta                     = 1;   %let ADASYN choose default(= [])
	ADASYN_kDensity                 = [];   %let ADASYN choose default
	ADASYN_kSMOTE                   = [];   %let ADASYN choose default
	ADASYN_featuresAreNormalized    = false;    %false lets ADASYN handle normalization
end

sigBand = {'broad','theta','alpha','beta','gamma','hfo'};

classAcc = cell(numel(chstat_feed),8);

for m = 1:numel(chstat_feed)
	fprintf('\n---------- Decoding %s band ----------\n',sigBand{m})
	% Reshape 3d channel stat array to 2d array by concatenating each 2d channel block.
	% Use fake data
	chstat_perm = permute(chstat_feed{m},[2 1 3]);
	chstat_2d = reshape(chstat_perm,size(chstat_perm,1),size(chstat_perm,2)*size(chstat_perm,3))';
	
	% Construct classification object.
	classObjPre = [repmat(cell2mat(behavSummary_feed(:,9)),[256,1]) chstat_2d];
	% Get indices for each class and plot them.
	display('Plotting all data points for classification ...')
	figure(91)
	clf(91)
	for i=1:numClass
		indClass(i) = {find(classObjPre(:,1)==i)};
		plot(classObjPre(indClass{i},2),classObjPre(indClass{i},3),'.')
		xlabel('Mean')
		ylabel('Variance')
		hold on
	end
	legend show
	suptitle('All channels')
	drawnow

	% Define arrays to store accuracy.
	class_acc = nan(itermax,numel(numchbstr));
	class_acc_rnd = nan(itermax,numel(numchbstr));
	class_acc_mean = nan(numel(numchbstr),1);
	class_acc_std = nan(numel(numchbstr),1);
	class_acc_sem = nan(numel(numchbstr),1);
	class_acc_rnd_mean = nan(numel(numchbstr),1);
	class_acc_rnd_std = nan(numel(numchbstr),1);
	class_acc_rnd_sem = nan(numel(numchbstr),1);
	
	%% Classify and bootstrap
	display('Starting classification and bootstrap...')
	k=1;
	chpick = [];
	
	for j=1:numel(numchpow)
		fprintf('\nNo. of Channels = %i\n',numchbstr(j))
		if runParallel == 'y'
			parfor iternum = 1:itermax
				classObj = classObjPre;
				if sampleChannel == 'r'
					chpick = sort(chproc(randperm(numch,numchbstr(j))))'
				elseif sampleChannel == 'l'
					chpick = chproc(1:numchbstr(j))'
				end
				indInstances = [];
				for i=1:numel(chpick)
					indInstances = [indInstances; [numtr*(chpick(i)-1)+1:numtr*chpick(i)]'];
				end
				if useADASYN == 'n'
					labels = classObj(indInstances,1);
					% Scale the data to normalize.
					classObj(indInstances,2:3) = (classObj(indInstances,2:3) - repmat(min(classObj(indInstances,2:3),[],1),size(classObj(indInstances,2:3),1),1))*spdiags(1./(max(classObj(indInstances,2:3),[],1)-min(classObj(indInstances,2:3),[],1))',0,size(classObj(indInstances,2:3),2),size(classObj(indInstances,2:3),2));
					features = classObj(indInstances,2:3);
				elseif useADASYN == 'y'
					% Define major and minor classes.
					% Assume class 1 is the majority.
					majClass = 1;
					majClassSize = nnz(classObj(indInstances,1)==1);
					indInstanceMaj = {indInstances(1)-1+find(classObj(indInstances,1)==majClass)};
					for i=2:numClass
						if nnz(classObj(indInstances,1)==i) > majClassSize
							majClass = i;
							majClassSize = nnz(classObj(indInstances,1)==i);
						end
					end
					% Redefine instance indices for true majority.
					indInstanceMaj = {indInstances(1)-1+find(classObj(indInstances,1)==majClass)};
	
					% Synthesize SV one minority by one minority against majority.
					classObjSyn = classObj(indInstances,:);
					for i=1:numClass
						if i ~= majClass
							indInstanceMin= {indInstances(1)-1+find(classObj(indInstances,1)==i)};
							[featuresSyn, labelsSyn] = ADASYN([classObj(indInstanceMaj{1},2:3);classObj(indInstanceMin{1},2:3)],[ones(size(indInstanceMaj{1}));zeros(size(indInstanceMin{1}))],ADASYN_beta, ADASYN_kDensity, ADASYN_kSMOTE, ADASYN_featuresAreNormalized);
							labelsSyn = i*ones(size(labelsSyn));
							classObjSyn = [classObjSyn; [labelsSyn, featuresSyn]];
						end
					end
					labels = classObjSyn(:,1);
					classObjSyn(:,2:3) = (classObjSyn(:,2:3) - repmat(min(classObjSyn(:,2:3),[],1),size(classObjSyn(:,2:3),1),1))*spdiags(1./(max(classObjSyn(:,2:3),[],1)-min(classObjSyn(:,2:3),[],1))',0,size(classObjSyn(:,2:3),2),size(classObjSyn(:,2:3),2));
					features = classObjSyn(:,2:3);
				else
					display('!! Cannot determine whether or not to use SMOTE !!')
				end
				% Classify
				class_acc(iternum,j,k) = libsvmtrain(labels,features,libsvmopt);
				labels_rnd = labels(randperm(length(labels)));
				class_acc_rnd(iternum,j,k) = libsvmtrain(labels_rnd,features,libsvmopt);
			end
		else
			figure(92)
			for iternum = 1:itermax
				classObj = classObjPre;
				if sampleChannel == 'r'
					chpick = sort(chproc(randperm(numch,numchbstr(j))))'
				elseif sampleChannel == 'l'
					chpick = chproc(1:numchbstr(j))'
				end
				indInstances = [];
				for i=1:numel(chpick)
					indInstances = [indInstances; [numtr*(chpick(i)-1)+1:numtr*chpick(i)]'];
				end
				if useADASYN == 'n'
					labels = classObj(indInstances,1);
					% Scale the data to normalize.
					classObj(indInstances,2:3) = (classObj(indInstances,2:3) - repmat(min(classObj(indInstances,2:3),[],1),size(classObj(indInstances,2:3),1),1))*spdiags(1./(max(classObj(indInstances,2:3),[],1)-min(classObj(indInstances,2:3),[],1))',0,size(classObj(indInstances,2:3),2),size(classObj(indInstances,2:3),2));
					features = classObj(indInstances,2:3);
				elseif useADASYN == 'y'
					% Define major and minor classes.
					% Assume class 1 is the majority.
					majClass = 1;
					majClassSize = nnz(classObj(indInstances,1)==1);
					indInstanceMaj = {indInstances(1)-1+find(classObj(indInstances,1)==majClass)};
					for i=2:numClass
						if nnz(classObj(indInstances,1)==i) > majClassSize
							majClass = i;
							majClassSize = nnz(classObj(indInstances,1)==i);
						end
					end
					% Redefine instance indices for true majority.
					indInstanceMaj = {indInstances(1)-1+find(classObj(indInstances,1)==majClass)};
	
					% Synthesize SV one minority by one minority against majority.
					classObjSyn = classObj(indInstances,:);
					for i=1:numClass
						if i ~= majClass
							indInstanceMin= {indInstances(1)-1+find(classObj(indInstances,1)==i)};
							[featuresSyn, labelsSyn] = ADASYN([classObj(indInstanceMaj{1},2:3);classObj(indInstanceMin{1},2:3)],[ones(size(indInstanceMaj{1}));zeros(size(indInstanceMin{1}))],ADASYN_beta, ADASYN_kDensity, ADASYN_kSMOTE, ADASYN_featuresAreNormalized);
							labelsSyn = i*ones(size(labelsSyn));
							classObjSyn = [classObjSyn; [labelsSyn, featuresSyn]];
						end
					end
					labels = classObjSyn(:,1);
					classObjSyn(:,2:3) = (classObjSyn(:,2:3) - repmat(min(classObjSyn(:,2:3),[],1),size(classObjSyn(:,2:3),1),1))*spdiags(1./(max(classObjSyn(:,2:3),[],1)-min(classObjSyn(:,2:3),[],1))',0,size(classObjSyn(:,2:3),2),size(classObjSyn(:,2:3),2));
					features = classObjSyn(:,2:3);
				else
					display('!! Cannot determine whether or not to use SMOTE !!')
				end
	
				% Plot data being classified at each iteration. (Only works in serial mode.)
				clf(92)
				figure(92)
				%subplot(numel(numchbstr),itermax,iternum*(j-1)+iternum)
				labels_uniq = sort(unique(labels));
				for l=1:length(labels_uniq)
					plot(features(find(labels==labels_uniq(l)),1),features(find(labels==labels_uniq(l)),2),'.')
					hold on
				end
				legend show
				title(sprintf('% band, numCh=%i, iterNum=%i',sigBand{m},numchbstr(j),iternum))
				xlabel('Mean')
				ylabel('Var')
				drawnow
				% Classify
				class_acc(iternum,j,k) = libsvmtrain(labels,features,libsvmopt);
				labels_rnd = labels(randperm(length(labels)));
				class_acc_rnd(iternum,j,k) = libsvmtrain(labels_rnd,features,libsvmopt);
			end
		end
		%% Calculate decoding accuracy for different number of channels.
		class_acc_mean(j,k) = nanmean(class_acc(:,j,k),1);
		class_acc_std(j,k) = nanstd(class_acc(:,j,k),1);
		class_acc_sem(j,k) = nanstd(class_acc(:,j,k),1)/sqrt(itermax);
		class_acc_rnd_mean(j,k) = nanmean(class_acc_rnd(:,j,k),1);
		class_acc_rnd_std(j,k) = nanstd(class_acc_rnd(:,j,k),1);
		class_acc_rnd_sem(j,k) = nanstd(class_acc_rnd(:,j,k),1)/sqrt(itermax);
	end
	classAcc{m,1} = class_acc;
	classAcc{m,2} = class_acc_mean;
	classAcc{m,3} = class_acc_std;
	classAcc{m,4} = class_acc_sem;
	classAcc{m,5} = class_acc_rnd;
	classAcc{m,6} = class_acc_rnd_mean;
	classAcc{m,7} = class_acc_rnd_std;
	classAcc{m,8} = class_acc_rnd_sem;
end

display('Decoding done.')


%%% Channel remapping.
%if numch == 256
%	% Redefine the data consistent with spatical channel distribution.
%	LUT = [78 250 249 229 230 186 177 165 119 111 100 56 55 39 38 82;
%	       164 245 241 234 238 180 185 171 117 121 107 48 53 47 42 120;
%	       174 247 243 232 236 178 175 169 122 113 104 108 59 51 40 112;
%	       242 244 239 235 176 181 168 172 114 115 106 110 52 49 43 44;
%	       240 251 231 228 184 187 173 170 116 118 101 102 58 57 41 46;
%	       70 246 237 233 182 179 183 167 123 103 109 105 50 54 45 204;
%	       33 158 66 124 31 36 25 37 199 205 197 202 196 206 200 203;
%	       35 1 7 63 3 61 5 29 207 254 217 252 215 209 211 201;
%	       62 2 60 6 94 88 98 0 213 140 194 210 198 208 192 212;
%	       96 90 92 91 95 93 97 193 20 143 139 138 142 132 144 136;
%	       99 89 65 67 125 71 127 248 27 147 148 151 149 141 145 133;
%	       64 69 126 160 155 163 157 255 21 75 68 76 72 150 146 154;
%	       153 161 159 189 135 191 137 162 26 83 87 81 77 79 73 74;
%	       131 128 129 188 130 190 134 156 16 32 86 30 84 24 85 28;
%	       216 224 220 226 214 222 218 166 22 11 10 14 8 80 12 34;
%	       223 219 227 221 225 195 253 152 18 23 19 9 17 13 15 4];
%	LUT = LUT + 1; %change to 1 base instead of 0 base
%	% For broad data
%	%for k=1:numch
%	%	chstat_broad_rearr(:,:,k) = chstat_broad(:,:,LUT(k));
%	%end
%	chstat_broad_newarr = permute(chstat_broad(:,:,LUT),[3 1 2]);
%	clear chstat_broad_rearr
%	for j=1:numtr
%		chstat_broad_mean_grid(:,:,j) = reshape(chstat_broad_newarr(:,j,1),16,16);
%		chstat_broad_var_grid(:,:,j) = reshape(chstat_broad_newarr(:,j,2),16,16);
%	end
%	clear chstat_broad_newarr
%	if ifFilter == 'y'
%		% For theta band
%		%for k=1:numch
%		%	chstat_theta_rearr(:,:,k) = chstat_theta(:,:,LUT(k));
%		%end
%		chstat_theta_newarr = permute(chstat_theta(:,:,LUT),[3 1 2]);
%		clear chstat_theta_rearr
%		for j=1:numtr
%			chstat_theta_mean_grid(:,:,j) = reshape(chstat_theta_newarr(:,j,1),16,16);
%			chstat_theta_var_grid(:,:,j) = reshape(chstat_theta_newarr(:,j,2),16,16);
%		end
%		clear chstat_theta_newarr
%		% For alpha band
%		%for k=1:numch
%		%	chstat_alpha_rearr(:,:,k) = chstat_alpha(:,:,LUT(k));
%		%end
%		chstat_alpha_newarr = permute(chstat_alpha(:,:,LUT),[3 1 2]);
%		clear chstat_alpha_rearr
%		for j=1:numtr
%			chstat_alpha_mean_grid(:,:,j) = reshape(chstat_alpha_newarr(:,j,1),16,16);
%			chstat_alpha_var_grid(:,:,j) = reshape(chstat_alpha_newarr(:,j,2),16,16);
%		end
%		clear chstat_alpha_newarr
%		% For beta band
%		%for k=1:numch
%		%	chstat_beta_rearr(:,:,k) = chstat_beta(:,:,LUT(k));
%		%end
%		chstat_beta_newarr = permute(chstat_beta(:,:,LUT),[3 1 2]);
%		clear chstat_beta_rearr
%		for j=1:numtr
%			chstat_beta_mean_grid(:,:,j) = reshape(chstat_beta_newarr(:,j,1),16,16);
%			chstat_beta_var_grid(:,:,j) = reshape(chstat_beta_newarr(:,j,2),16,16);
%		end
%		clear chstat_beta_newarr
%		% For gamma band
%		%for k=1:numch
%		%	chstat_gamma_rearr(:,:,k) = chstat_gamma(:,:,LUT(k));
%		%end
%		chstat_gamma_newarr = permute(chstat_gamma(:,:,LUT),[3 1 2]);
%		clear chstat_gamma_rearr
%		for j=1:numtr
%			chstat_gamma_mean_grid(:,:,j) = reshape(chstat_gamma_newarr(:,j,1),16,16);
%			chstat_gamma_var_grid(:,:,j) = reshape(chstat_gamma_newarr(:,j,2),16,16);
%		end
%		clear chstat_gamma_newarr
%		% For hfo band
%		%for k=1:numch
%		%	chstat_hfo_rearr(:,:,k) = chstat_hfo(:,:,LUT(k));
%		%end
%		chstat_hfo_newarr = permute(chstat_hfo(:,:,LUT),[3 1 2]);
%		clear chstat_hfo_rearr
%		for j=1:numtr
%			chstat_hfo_mean_grid(:,:,j) = reshape(chstat_hfo_newarr(:,j,1),16,16);
%			chstat_hfo_var_grid(:,:,j) = reshape(chstat_hfo_newarr(:,j,2),16,16);
%		end
%		clear chstat_hfo_newarr
%	end
%end

%%============================== End of bootstraplotpping =========================

toc
