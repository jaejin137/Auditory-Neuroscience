function spikeStruct = loadKSdir(ksDir,varargin)

if ~isempty(varargin)
    params = varargin{1};
else
    params = [];
end

if ~isfield(params, 'excludeNoise')
    params.excludeNoise = true;
end
if ~isfield(params, 'loadPCs')
    params.loadPCs = false;
end

% Find indices for ripple stimulus
%load('Z:\Parooa\Synapse\h\CC_sortingInfoTable.mat');

%idx_ripOn = sortingInfo.offSetStartSample(find((sortingInfo.blockName==blockName).*(sortingInfo.chamber=='AC').*strcmp(sortingInfo.recordGroupName,areaName)));
%len_ripple = sortingInfo.nSample(find((sortingInfo.blockName==blockName).*(sortingInfo.chamber=='AC').*strcmp(sortingInfo.recordGroupName,areaName)));
%idx_ripOff = idx_ripOn + len_ripple;

% load spike data

spikeStruct = loadParamsPy(fullfile(ksDir, 'params.py'));

% Load entire spike time and cluster spike time
ss = double(readNPY(fullfile(ksDir, 'spike_times.npy')));
%ss_clu = struct2array(load(fullfile(ksDir,'sorting_metrics','cluster_spike_output.mat')));

% Cut out only ripple part
%ss_ripple = ss(find((ss>=idx_ripOn).*(ss<=idx_ripOff)))-idx_ripOn;
%ss_clu_ripple = ss_clu(find((ss>=idx_ripOn).*(ss<=idx_ripOff)),:);
%ss_clu_ripple(:,2) = ss_clu_ripple(:,2) - idx_ripOn;
%clu_no = sort(unique(ss_clu_ripple(:,1)));
% Assign spike time
%spikeTimeRip = ss_ripple;
% Assign cluster spike time
%for i=1:length(clu_no)
%spikeTimeRipClus{i,1} = clu_no(i);
%spikeTimeRipClus{i,2} = ss_clu_ripple(find(ss_clu_ripple(:,1)==clu_no(i)),2);
%end
% Save ripple spike times
%save(fullfile(ksDir,'spike_times_ripple.mat'),'spikeTimeRip')
%save(fullfile(ksDir,'spike_times_ripple_clust.mat'),'spikeTimeRipClus')
%clearvars spikeTimeRip spikeTimeRipClus

%st = double(ss_ripple)/spikeStruct.sample_rate;
st = double(ss)/spikeStruct.sample_rate;
spikeTemplates = readNPY(fullfile(ksDir, 'spike_templates.npy')); % note: zero-indexed
%spikeTemplates = spikeTemplates(find((ss>=idx_ripOn).*(ss<=idx_ripOff)));

if exist(fullfile(ksDir, 'spike_clusters.npy'))
    clu = readNPY(fullfile(ksDir, 'spike_clusters.npy'));
%    clu = clu(find((ss>=idx_ripOn).*(ss<=idx_ripOff)));
else
    clu = spikeTemplates;
end

tempScalingAmps = readNPY(fullfile(ksDir, 'amplitudes.npy'));
%tempScalingAmps = tempScalingAmps(find((ss>=idx_ripOn).*(ss<=idx_ripOff)));

if params.loadPCs
    pcFeat = readNPY(fullfile(ksDir,'pc_features.npy')); % nSpikes x nFeatures x nLocalChannels
%    pcFeat = pcFeat(find((ss>=idx_ripOn).*(ss<=idx_ripOff)));
    pcFeatInd = readNPY(fullfile(ksDir,'pc_feature_ind.npy')); % nTemplates x nLocalChannels
else
    pcFeat = [];
    pcFeatInd = [];
end

cgsFile = '';
if exist(fullfile(ksDir, 'cluster_groups.csv')) 
    cgsFile = fullfile(ksDir, 'cluster_groups.csv');
end
if exist(fullfile(ksDir, 'cluster_group.tsv')) 
   cgsFile = fullfile(ksDir, 'cluster_group.tsv');
end 
if ~isempty(cgsFile)
    [cids, cgs] = readClusterGroupsCSV(cgsFile);

    if params.excludeNoise
        noiseClusters = cids(cgs==0);

        st = st(~ismember(clu, noiseClusters));
        spikeTemplates = spikeTemplates(~ismember(clu, noiseClusters));
        tempScalingAmps = tempScalingAmps(~ismember(clu, noiseClusters));        
        
        if params.loadPCs
            pcFeat = pcFeat(~ismember(clu, noiseClusters), :,:);
            %pcFeatInd = pcFeatInd(~ismember(cids, noiseClusters),:);
        end
        
        clu = clu(~ismember(clu, noiseClusters));
        cgs = cgs(~ismember(cids, noiseClusters));
        cids = cids(~ismember(cids, noiseClusters));
        
        
    end
    
else
    clu = spikeTemplates;
    
    cids = unique(spikeTemplates);
    cgs = 3*ones(size(cids));
end
    

coords = readNPY(fullfile(ksDir, 'channel_positions.npy'));
ycoords = coords(:,2); xcoords = coords(:,1);
temps = readNPY(fullfile(ksDir, 'templates.npy'));

winv = readNPY(fullfile(ksDir, 'whitening_mat_inv.npy'));

spikeStruct.st = st;
spikeStruct.spikeTemplates = spikeTemplates;
spikeStruct.clu = clu;
spikeStruct.tempScalingAmps = tempScalingAmps;
spikeStruct.cgs = cgs;
spikeStruct.cids = cids;
spikeStruct.xcoords = xcoords;
spikeStruct.ycoords = ycoords;
spikeStruct.temps = temps;
spikeStruct.winv = winv;
spikeStruct.pcFeat = pcFeat;
spikeStruct.pcFeatInd = pcFeatInd;