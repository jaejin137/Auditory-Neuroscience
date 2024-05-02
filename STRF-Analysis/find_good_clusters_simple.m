%% This script simply access the STRF reliability indices and find good clusters.
function [goodClust] = find_good_clusters(monkeyName,sessionDate,p_crit)

%if nargin <6
%	NB = 100;
%elseif nargin < 5
	NBoot = 8;
	NB = 100;
%end
if nargin < 3
	p_crit = 1.e-8;
end


% Load RI data
target_RI = dir(sprintf('~/Research/Auditory/STRF/RI/RI_%s-%s_*_8-100.mat',monkeyName,sessionDate));
% Sanity check
if length(target_RI) < 1
	fprintf('\n\n!!!Error: No target RI data found. Skipped!!!\n')
elseif length(target_RI) > 1
	fprintf('\n\n!!!Error: Multiple target RI data found. First one used!!!\n')
end

load(target_RI(1).name);

% Initialize good cluster numbers.
goodClust = {};
% Find good clusters
i = 1;
for k = 1:length(STRFSig)
	if STRFSig(k).p < p_crit
		goodClust{i,1} = STRFSig(k).clustNum;
		goodClust{i,2} = STRFSig(k).clustType;
		i = i+1;
	end
end

