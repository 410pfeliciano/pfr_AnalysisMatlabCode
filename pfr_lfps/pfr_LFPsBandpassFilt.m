% Use the save path to save your data to a specific directory

savePath = 'G:\vGAT_CreCh2_4\ForcedD3FreeD1\LFPs'; % specify the saving path


%% Bandpass Filtering Delta 1-4Hz Filter

lfpFilter = 'delta'; % you can choose between 'ripple', 'lowGamma'
... 'fastGamma', 'theta', 'delta', 'spindle', 'deltaSpindle'

LFPtofilt = []; % Remove [] for the lfp of interest

if isempty(LFPtofilt)
    warning('You have not chosen the LFP to filter, LFPtofilt is empty');
    return
end

filtLFP = bandFilter (LFPtofilt, lfpFilter); % default Fs = 1000
prompt = 'Assign a new variable name to your data ? '; %e.g. for a lfp from
... ch1 put ch1LFP for ch60 put ch60LFP or (ch# tet# LFP) 
xx = input(prompt, 's');
assignin('base',xx,filtLFP) % assignin the new variable
file_name1 = sprintf('%s.mat',xx);
if exist('savePath','var') == 0
    warning('You have not chosen the path to save the data, Run the first part of this script');
    return
end
fprintf(1, 'Saving file\n'); 
save ([savePath filesep file_name1], xx,'time','-v7.3') % saving spike data
clear('filtLFP','LFPtofilt', 'lfpFilter', 'xx')
fprintf(1, 'Done!\n');


%% Bandpass Filtering Spindle 

lfpFilter = 'spindle'; % you can choose between 'ripple', 'lowGamma'
... 'fastGamma', 'theta', 'delta', 'spindle', 'deltaSpindle'

LFPtofilt = []; % Remove [] for the lfp of interest

if isempty(LFPtofilt)
    warning('You have not chosen the LFP to filter, LFPtofilt is empty');
    return
end

filtLFP = bandFilter (LFPtofilt, lfpFilter); % default Fs = 1000
prompt = 'Assign a new variable name to your data ? '; %e.g. for a lfp from
... ch1 put ch1LFP for ch60 put ch60LFP or (ch# tet# LFP) 
xx = input(prompt, 's');
assignin('base',xx,filtLFP) % assignin the new variable
file_name1 = sprintf('%s.mat',xx);
if exist('savePath','var') == 0
    warning('You have not chosen the path to save the data, Run the first part of this script');
    return
end
fprintf(1, 'Saving file\n'); 
save ([savePath filesep file_name1], xx,'time','-v7.3') % saving spike data
clear('filtLFP','LFPtofilt', 'lfpFilter', 'xx')
fprintf(1, 'Done!\n'); 


%% Bandpass Filtering DeltaSpindle

lfpFilter = 'deltaSpindle'; % you can choose between 'ripple', 'lowGamma'
... 'fastGamma', 'theta', 'delta', 'spindle', 'deltaSpindle'
% Choose the LFP data to Filter

LFPtofilt = ch42Tet13LFP - mean(ch42Tet13LFP); % Remove [] for the lfp of interest

if isempty(LFPtofilt)
    warning('You have not chosen the LFP to filter, LFPtofilt is empty');
    return
end

filtLFP = bandFilter (LFPtofilt, lfpFilter); % default Fs = 1000
prompt = 'Assign a new variable name to your data ? '; %e.g. for a lfp from
... ch1 put ch1LFP for ch60 put ch60LFP or (ch# tet# LFP) 
xx = input(prompt, 's');
assignin('base',xx,filtLFP) % assignin the new variable
file_name1 = sprintf('%s.mat',xx);
if exist('savePath','var') == 0
    warning('You have not chosen the path to save the data, Run the first part of this script');
    return
end
fprintf(1, 'Saving file\n'); 
save ([savePath filesep file_name1], xx,'time','-v7.3') % saving spike data
clear('filtLFP','LFPtofilt', 'lfpFilter', 'xx')
fprintf(1, 'Done!\n'); 



%% Bandpass Filtering theta
lfpFilter = 'theta'; % you can choose between 'ripple', 'lowGamma'
... 'fastGamma', 'theta', 'delta', 'spindle', 'deltaSpindle'
% Choose the LFP data to Filter

LFPtofilt = []; % Remove [] for the lfp of interest

if isempty(LFPtofilt)
    warning('You have not chosen the LFP to filter, LFPtofilt is empty');
    return
end

filtLFP = bandFilter (LFPtofilt, lfpFilter); % default Fs = 1000
prompt = 'Assign a new variable name to your data ? '; %e.g. for a lfp from
... ch1 put ch1LFP for ch60 put ch60LFP or (ch# tet# LFP) 
xx = input(prompt, 's');
assignin('base',xx,filtLFP) % assignin the new variable
file_name1 = sprintf('%s.mat',xx);
if exist('savePath','var') == 0
    warning('You have not chosen the path to save the data, Run the first part of this script');
    return
end
fprintf(1, 'Saving file\n'); 
save ([savePath filesep file_name1], xx,'time','-v7.3') % saving spike data
clear('filtLFP','LFPtofilt', 'lfpFilter', 'xx')
fprintf(1, 'Done!\n'); 


%% slow Gamma filtering
lfpFilter = 'lowGamma'; % you can choose between 'ripple', 'lowGamma'
... 'fastGamma', 'theta', 'delta', 'spindle', 'deltaSpindle'
% Choose the LFP data to Filter

LFPtofilt = []; % Remove [] for the lfp of interest

if isempty(LFPtofilt)
    warning('You have not chosen the LFP to filter, LFPtofilt is empty');
    return
end

filtLFP = bandFilter (LFPtofilt, lfpFilter); % default Fs = 1000
prompt = 'Assign a new variable name to your data ? '; %e.g. for a lfp from
... ch1 put ch1LFP for ch60 put ch60LFP or (ch# tet# LFP) 
xx = input(prompt, 's');
assignin('base',xx,filtLFP) % assignin the new variable
file_name1 = sprintf('%s.mat',xx);
if exist('savePath','var') == 0
    warning('You have not chosen the path to save the data, Run the first part of this script');
    return
end
fprintf(1, 'Saving file\n'); 
save ([savePath filesep file_name1], xx,'time','-v7.3') % saving spike data
clear('filtLFP','LFPtofilt', 'lfpFilter', 'xx')
fprintf(1, 'Done!\n'); 


%% fast Gamma Filtering
lfpFilter = 'fastGamma'; % you can choose between 'ripple', 'lowGamma'
... 'fastGamma', 'theta', 'spindle', 'delta', 'deltaSpindle'
% Choose the LFP data to Filter

LFPtofilt = ch42Tet13LFP; % Remove [] for the lfp of interest

if isempty(LFPtofilt)
    warning('You have not chosen the LFP to filter, LFPtofilt is empty');
    return
end

filtLFP = bandFilter (LFPtofilt, lfpFilter); % default Fs = 1000
prompt = 'Assign a new variable name to your data ? '; %e.g. for a lfp from
... ch1 put ch1LFP for ch60 put ch60LFP or (ch# tet# LFP) 
xx = input(prompt, 's');
assignin('base',xx,filtLFP) % assignin the new variable
file_name1 = sprintf('%s.mat',xx);
if exist('savePath','var') == 0
    warning('You have not chosen the path to save the data, Run the first part of this script');
    return
end
fprintf(1, 'Saving file\n'); 
save ([savePath filesep file_name1], xx,'time','-v7.3') % saving spike data
clear('filtLFP','LFPtofilt', 'lfpFilter', 'xx')
fprintf(1, 'Done!\n'); 


%% Ripple Filtering
lfpFilter = 'ripple'; % you can choose between 'ripple', 'lowGamma'
... 'fastGamma', 'theta', 'spindle', 'delta', 'deltaSpindle'
% Choose the LFP data to Filter

LFPtofilt = []; % Remove [] for the lfp of interest

if isempty(LFPtofilt)
    warning('You have not chosen the LFP to filter, LFPtofilt is empty');
    return
end

filtLFP = bandFilter (LFPtofilt, lfpFilter); % default Fs = 1000
prompt = 'Assign a new variable name to your data ? '; %e.g. for a lfp from
... ch1 put ch1LFP for ch60 put ch60LFP or (ch# tet# LFP) 
xx = input(prompt, 's');
assignin('base',xx,filtLFP) % assignin the new variable
file_name1 = sprintf('%s.mat',xx);
if exist('savePath','var') == 0
    warning('You have not chosen the path to save the data, Run the first part of this script');
    return
end
fprintf(1, 'Saving file\n'); 
save ([savePath filesep file_name1], xx,'time','-v7.3') % saving spike data
clear('filtLFP','LFPtofilt', 'lfpFilter', 'xx')
fprintf(1, 'Done!\n'); 




%% highFrequency filter for delta detection
lfpFilter = 'highfreq'; % you can choose between 'ripple', 'lowGamma'
... 'fastGamma', 'theta', 'spindle', 'delta', 'deltaSpindle'
% Choose the LFP data to Filter

LFPtofilt = ch42Tet13LFP; % Remove [] for the lfp of interest

if isempty(LFPtofilt)
    warning('You have not chosen the LFP to filter, LFPtofilt is empty');
    return
end

filtLFP = bandFilter (LFPtofilt, lfpFilter); % default Fs = 1000
prompt = 'Assign a new variable name to your data ? '; %e.g. for a lfp from
... ch1 put ch1LFP for ch60 put ch60LFP or (ch# tet# LFP) 
xx = input(prompt, 's');
assignin('base',xx,filtLFP) % assignin the new variable
file_name1 = sprintf('%s.mat',xx);
if exist('savePath','var') == 0
    warning('You have not chosen the path to save the data, Run the first part of this script');
    return
end
fprintf(1, 'Saving file\n'); 
save ([savePath filesep file_name1], xx,'time','-v7.3') % saving spike data
clear('filtLFP','LFPtofilt', 'lfpFilter', 'xx')
fprintf(1, 'Done!\n'); 