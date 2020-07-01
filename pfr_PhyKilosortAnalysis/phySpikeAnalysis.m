
%{ 
First Step Loading Data Sets
Obtain Spike Time and cluster ID from npy file from kilosort. 
You need npy-matlab in the path.
%}
sampFreq = 30000; % sampling rate in Hz
sampRate = 1/sampFreq;
spikeTime = readNPY('spike_times.npy'); % spike times from phy
spkIds = readNPY('spike_clusters.npy'); % cluster id from phy
clusterInfo = tdfread('cluster_info.tsv');
clustersId = clusterInfo.id;
clusterFR = clusterInfo.fr;
clInfo = [clusterInfo.id, clusterInfo.depth,clusterInfo.fr];
%% Second Step Filtering the 'noise' clusters
noiseLabelID= clusterInfo.id(regexp(clusterInfo.group(:,1)', 'n'),1); % identify noise cluster
noiseIdx = find(ismember(spkIds, noiseLabelID)); %
spkTimeClust = [spikeTime, spkIds]; % spike time and cluster ID
spkTimeClust(noiseIdx,:) = []; % eliminating noise from clusters and spikes
% spkTimeClustSort = double(sortrows(spkTimeClust,2)); % convert to double
% spkTimeClustSort(:,1) = spkTimeClustSort(:,1)./sampFreq; % sample to sec
% spkGoodMUAIdx = find(spikeLabels.group(:,1) == 'g'| ... 
%     spikeLabels.group(:,1) == 'm'); % good and mua clusters indices
spkNoiseIdx = find(clusterInfo.group(:,1) == 'n'); % indices noise cluester id
clustersId(spkNoiseIdx,:) = []; % remove noise clusters
clInfo(spkNoiseIdx,:) = [];

if isfieldRecursive(clusterInfo, 'neurontype') %Matthew Arthington (2020). Recursively check 
    ...fields of a structure exist  MATLAB Central File Exchange. Retrieved May 11, 2020.
    clusterInfo.neurontype(spkNoiseIdx,:) = [];
end
if isfieldRecursive(clusterInfo, 'region')
    clusterInfo.region(spkNoiseIdx,:) = [];
end
clusterInfo.group(spkNoiseIdx,:) = [];
% % 
%% Create a Struct with every Good or MUA cluster
idx = zeros(size(spkTimeClust,1),1);
for k = 1 : size(clustersId,1)
    clstIdx = find(spkTimeClust(:,2) == clustersId(k));
    spkClust(k).clusterId = spkTimeClust(clstIdx(1),2);
    spkClust(k).spkSamples = spkTimeClust(clstIdx,1); % Spike time samples
    spkClust(k).spkTime = double(spkClust(k).spkSamples) .* sampRate; % spike time
    spkClust(k).tetNum = clInfo(k,2);
    spkClust(k).clustFR = clInfo(k,3);
    spkClust(k).neurontype = clusterInfo.neurontype(k,:);
    spkClust(k).region = clusterInfo.region(k,:);
    spkClust(k).singleUnit = clusterInfo.group(k,:);
%     spkClust(k).channel = clInfo(k,2);
end

%%
save('units.mat', 'spkClust')
clear