% Modified Pedro Feliciano from Jon Newman
% 2018.10.03
% Clusterless oKDE

%% Create training and test data sets
% load('../data/kde_decoder_args_amp_filter.mat');

% Tetrodes to decode
tetrodes = [2:8 10]';
source_filter = ismember(1:30,tetrodes);

% Bandwidth parameters
pos_bandwidth_m = 0.05;
amp_bandwidth_uV = 25;

%% Divided data into training and test sets
len = 3.59;                         % circumference of the maze in m
dx = 0.01;                          % position granularity in m
pos_grid = 0:dx:len;

n = size(RUN,1);
n_enc_segs = ceil(n/2);
n_dec_segs = n - n_enc_segs;
encode_epochs = RUN(1:n_enc_segs,:);
decode_epochs = RUN(n_enc_segs+1:n,:);

%% Examine the RUN segments to make sure they are sane
close all
figure
hold on
plot(T,X)
for i = 1:size(encode_epochs,1)
   b = T >= encode_epochs(i,1) & T <= encode_epochs(i,2);
   plot(T(b),X(b),'g')
end

for i = 1:size(decode_epochs,1)
   b = T >= decode_epochs(i,1) & T <= decode_epochs(i,2);
   plot(T(b),X(b),'r')
end

%% Decoding time bins 
dt = 0.250; % ms
decode_run_time_total = 0;
decode_time_total = decode_epochs(end,2) - decode_epochs(1,1);
decode_run_bins = [];
decode_run_bounds = [];

% Construct bins
for i = 1:size(decode_epochs,1)
    b0 = decode_epochs(i,1):dt:decode_epochs(i,2);
    b1 = b0 + dt;
    l0 = size(decode_run_bins,1)+1;
    decode_run_bins = [decode_run_bins; [b0' b1']];
    l1 = size(decode_run_bins,1);
    decode_run_bounds = [decode_run_bounds; [l0 l1]];
    
    decode_run_time_total = decode_run_time_total +  diff(decode_epochs(i,:));
end
n_dec_bins = size(decode_run_bins,1);

all_decode_bins = min(decode_run_bins(:,1)):dt:max(decode_run_bins(:,2))-dt;
all_decode_bins = [all_decode_bins' all_decode_bins' + dt];
n_all_decode_bins = size(all_decode_bins,1);

% Store time-series bin info
bins.dt = dt;
bins.all_decode_bins = all_decode_bins;
bins.decode_epochs =  decode_epochs;
bins.decode_bounds =  decode_run_bounds;
bins.decode_bins =  decode_run_bins;
bins.decode_run_time_total = decode_run_time_total;
bins.decode_time_total = decode_time_total;

%% Create kde_decoder
for i = 1:length(S)
   if (isempty(S{i}))
       S{i} = 1;
       XX{i} = 1;
       A{i} = [1 1 1 1];
   end
end


z = kde_decoder(T,X,S,XX,A, ...
    'encoding_segments', encode_epochs, ...
    'stimulus_variable_type', 'linear', ...
    'stimulus_grid', {pos_grid}, ...              
    'stimulus_kernel', 'gaussian', ...
    'stimulus_bandwidth', pos_bandwidth_m, ...
    'response_variable_type', 'linear', ...
    'response_kernel', 'gaussian', ...
    'response_bandwidth', amp_bandwidth_uV, ...
    'source_selection', source_filter, ...
    'parallel', 0);

%% Compute posterior
%P = posterior
%E = decoding performance

tic
[pos_decode, E_test, test_stim] = z.compute(all_decode_bins);
toc

%% Figure: all bin time series
close all
figure('units','centimeters','position', [1 1 30 10]) ; 
xlimit = [0 5200];
ylimit = [0 4];
%colormap(flipud(bone.^2))

% Full time series
mbins = mean(all_decode_bins,2);
binned_pos = zeros(n_all_decode_bins ,1);

% RUN portions of time series
mbins_run = [];
% pos_decode_run = [];
binned_pos_run = zeros(size(mbins));
error_run = [];

subplot(211)
hold on

imagesc(mbins,pos_grid,pos_decode);
set(gca,'YDir','normal','XTick',[],'Xlim',xlimit,'ylim',ylimit)

% All bins, even during non-RUN
b = T >= mbins(1) & T <= mbins(end);
plot(T(b),X(b),'b')
set(gca,'tickdir','out','ticklength',[0.005 0.01],'Box','off','XTick',[],'Xlim',xlimit,'ylim',ylimit)
ylabel('Distance (m)')

% RUN-only bins
t_vs_p = cell(1,size(decode_epochs,1));
for i = 1:size(decode_epochs,1)
   b = T >= decode_epochs(i,1) & T <= decode_epochs(i,2);
   t_vs_p{i} = [T(b), X(b)];
   plot(T(b),X(b),'r')
end

for j = 1:n_all_decode_bins ;
    
    [~, min_idx] = min(abs(T-mbins(j)));
    binned_pos(j) = X(min_idx);
end

for i = 1:size(decode_epochs,1)
    
    curr_idx = mbins >= decode_epochs(i,1) & mbins < decode_epochs(i,2);
    mbins_run = [mbins_run; mbins(curr_idx)];
    curr_P = pos_decode(:,curr_idx);
    
    [~, max_idx] = max(curr_P);
    mode_decode = pos_grid(max_idx)';
    err = abs(mode_decode - binned_pos(curr_idx));
    error_run = [error_run; err];

end

subplot(212)
hold on

for j = 1:numel(t_vs_p)
   
    t0 = t_vs_p{j}(1,1);
    t1 = t_vs_p{j}(end,1);
    
    plot([t0 t1],[3.5 3.5],'k','linewidth',4);
    
end

[~, max_idx] = max(pos_decode);
mode_decode = pos_grid(max_idx)';
err = abs(mode_decode - binned_pos);
plot(mbins, err, 'k')
set(gca,'tickdir','out','ticklength',[0.005 0.01],'Box','off','Xlim',xlimit,'ylim',ylimit)
ylabel('Distance (m)')
xlabel('Time (sec)')

%export_fig('decode_ts_01','-transparent','-pdf',gcf)

error.error_full = err;
error.error_run = error_run;