%% Binning spikes for firing rate calculation
Fs = 1000; % This value will depend
Tet = time; % just to calculate the duration of the recording
bins = (0 : 1/Fs : length(Tet) * (1/Fs) - (1/Fs))';
if exist('spkClust','var') == 0
    warning('spiking data spkClust is not in your workspace');
    return
end

%% CA1 
TetCa1 = [15, 2,10]; % choose the tetrode of interest
CA1Tetfind = cell(size(spkClust,2),1);
for kk = 1 : size(spkClust,2)
    CA1Tetfind{kk} = find(spkClust(kk).tetNum == TetCa1);
end
CA1TetIdx = find(~cellfun('isempty', CA1Tetfind));
CA1spiketrain = zeros(size(CA1TetIdx,1),length(bins));
for kk = 1 : size(CA1TetIdx,1)
    CA1spiketime = (spkClust(CA1TetIdx(kk)).spkTime)';
    CA1spiketrain(kk,:) = hist(CA1spiketime,bins);
end
%% Spike density estimation
CA1totalSpk = sum(CA1spiketrain,1);
sigma = 0.02;
edges = -3*sigma: 0.001: 3*sigma;
kernel = (normpdf(edges,0,sigma)).*(1/10);
% kernel = kernel*(1/Fs); % kernel muitiplied by bin size
CA1spkDensity = conv(CA1totalSpk, kernel);
center = ceil(length(edges)/2);
CA1spkDensity = CA1spkDensity(center:length(CA1totalSpk)+center-1);

%% CA3
TetCa3 = [7, 8, 9, 11]; % choose the tetrode of interest
CA3Tetfind = cell(size(spkClust,2),1);
for kk = 1 : size(spkClust,2)
    CA3Tetfind{kk} = find(spkClust(kk).tetNum == TetCa3);
end
CA3TetIdx = find(~cellfun('isempty', CA3Tetfind));
CA3spiketrain = zeros(size(CA3TetIdx,1),length(bins));
for kk = 1 : size(CA3TetIdx,1)
    CA3spiketime = (spkClust(CA3TetIdx(kk)).spkTime)';
    CA3spiketrain(kk,:) = hist(CA3spiketime,bins);
end
%% Spike density estimation
CA3totalSpk = sum(CA3spiketrain,1);
sigma = 0.02;
edges = -3*sigma: 0.001: 3*sigma;
kernel = (normpdf(edges,0,sigma)).*(1/10);
% kernel = kernel*(1/Fs); % kernel muitiplied by bin size
CA3spkDensity = conv(CA3totalSpk, kernel);
center = ceil(length(edges)/2);
CA3spkDensity = CA3spkDensity(center:length(CA3totalSpk)+center-1);


%% DG
TetDG = [5, 6]; % choose the tetrode of interest
DGTetfind = cell(size(spkClust,2),1);
for kk = 1 : size(spkClust,2)
    DGTetfind{kk} = find(spkClust(kk).tetNum == TetDG);
end
DGTetIdx = find(~cellfun('isempty', DGTetfind));
DGspiketrain = zeros(size(DGTetIdx,1),length(bins));
for kk = 1 : size(DGTetIdx,1)
    DGspiketime = (spkClust(DGTetIdx(kk)).spkTime)';
    DGspiketrain(kk,:) = hist(DGspiketime,bins);
end
%% Spike density estimation
DGtotalSpk = sum(DGspiketrain,1);
sigma = 0.02;
edges = -3*sigma: 0.001: 3*sigma;
kernel = (normpdf(edges,0,sigma)).*(1/10);
% kernel = kernel*(1/Fs); % kernel muitiplied by bin size
DGspkDensity = conv(DGtotalSpk, kernel);
center = ceil(length(edges)/2);
DGspkDensity = DGspkDensity(center:length(DGtotalSpk)+center-1);

%% Cortex
TetCtx = 13; % choose the tetrode of interest
CtxTetfind = cell(size(spkClust,2),1);
for kk = 1 : size(spkClust,2)
    CtxTetfind{kk} = find(spkClust(kk).tetNum == TetCtx);
end
CtxTetIdx = find(~cellfun('isempty', CtxTetfind));
Ctxspiketrain = zeros(size(CtxTetIdx,1),length(bins));
for kk = 1 : size(CtxTetIdx,1)
    Ctxspiketime = (spkClust(CtxTetIdx(kk)).spkTime)';
    Ctxspiketrain(kk,:) = hist(Ctxspiketime,bins);
end
%% Spike density estimation
CtxtotalSpk = sum(Ctxspiketrain,1);
sigma = 0.02;
edges = -3*sigma: 0.001: 3*sigma;
kernel = (normpdf(edges,0,sigma)).*(1/10);
% kernel = kernel*(1/Fs); % kernel muitiplied by bin size
CtxspkDensity = conv(CtxtotalSpk, kernel);
center = ceil(length(edges)/2);
CtxspkDensity = CtxspkDensity(center:length(CtxtotalSpk)+center-1);
