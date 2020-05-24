%{ 
Function to import units/MUAs to MATLAB 
You need npy-matlab in matlab path. https://github.com/kwikteam/npy-matlab

Inputs:
    1- path to python files created from kilosort sorting procedure
    2- savePath - path to save the units/MUAS
    3- sampFreq(Default 30000) - sampling frequency. for open-ephys is 30KHz 
%}
% npy_dir = dir(npy_path);
function pfr_npy2spikes(npy_path, savePath, sampFreq)

if nargin < 3
    sampFreq = 30000;
end
sampRate = 1/sampFreq;
file_path = fullfile(npy_path,{'spike_times.npy'; 'spike_clusters.npy';...
    'cluster_info.tsv'}); 
spikeTime = readNPY(file_path{1}); % spike times from phy
spkIds = readNPY(file_path{2}); % cluster id from phy
clusterInfo = tdfread(file_path{3});
clustersId = clusterInfo.id;
clInfo = [clusterInfo.id, clusterInfo.depth];
%% Second Step Filtering the 'noise' clusters
noiseLabelID= clusterInfo.id(regexp(clusterInfo.group(:,1)', ...
    'n'),1); % identify noise cluster
noiseIdx = find(ismember(spkIds, noiseLabelID)); %
spkTimeClust = [spikeTime, spkIds]; % spike time and cluster ID
spkTimeClust(noiseIdx,:) = []; % eliminating noise from clusters and spikes
% spkGoodMUAIdx = find(spikeLabels.group(:,1) == 'g'| ... 
%     spikeLabels.group(:,1) == 'm'); % good and mua clusters indices
spkNoiseIdx = find(clusterInfo.group(:,1) == 'n'); % indices noise cluester id
clustersId(spkNoiseIdx,:) = []; % remove noise clusters
clInfo(spkNoiseIdx,:) = [];
%% Create a Struct with every Good or MUA cluster

for k = 1 : size(clustersId,1)
    clstIdx = find(spkTimeClust(:,2) == clustersId(k));
    spkClust(k).clusterId = spkTimeClust(clstIdx(1),2);
    spkClust(k).spkSamples = spkTimeClust(clstIdx,1); % Spike time samples
    spkClust(k).spkTime = double(spkClust(k).spkSamples) .* sampRate; % spike time
    spkClust(k).tetNum = clInfo(k,2);
%     spkClust(k).channel = clInfo(k,2);
end
%%
save([savePath filesep 'units.mat'], 'spkClust')