%% This script simply access the STRF reliability indices and find good clusters.
function [goodClust] = find_good_clusters(monkeyName,sessionDate,driveID,p_crit,NBoot,NB)

%if nargin <6
%	NB = 100;
%elseif nargin < 5
	NBoot = 8;
	NB = 100;
%end


% Load RI data
load(sprintf('~/Research/Auditory/STRF/RI/RI_%s-%s_%s_%i-%i.mat',monkeyName,sessionDate,driveID,NBoot,NB))
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

