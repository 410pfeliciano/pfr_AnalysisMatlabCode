% Obtained the whole LFP. Add the all directory paths of the experiments
ephysPaths ={'G:\vGAT_CreCh2_4\ForcedD3FreeD1\LFPs\pre',...
    'G:\vGAT_CreCh2_4\ForcedD3FreeD1\LFPs\run1', ...
    'G:\vGAT_CreCh2_4\ForcedD3FreeD1\LFPs\run2', ...
    'G:\vGAT_CreCh2_4\ForcedD3FreeD1\LFPs\run3', ...
    'G:\vGAT_CreCh2_4\ForcedD3FreeD1\LFPs\run4', ...
    'G:\vGAT_CreCh2_4\ForcedD3FreeD1\LFPs\post1', ...
    'G:\vGAT_CreCh2_4\ForcedD3FreeD1\LFPs\post2'};
savePath = 'G:\vGAT_CreCh2_4\ForcedD3FreeD1\LFPs'; % specify the saving path


%% whole LFP function
% rawFs = 30000; % this value is default in the pfr_LFP function
% reduction = 30; % this value is default in the pfr_LFP function
tetChannel = 62;
[LFP, time] = pfr_LFP (tetChannel, ephysPaths);

prompt = 'Assign a new variable name to your data ? '; %e.g. for a lfp from
... ch1 put ch1LFP for ch60 put ch60LFP or (ch# tet# LFP)
    
xx = input(prompt, 's');
assignin('base',xx,LFP) % assignin the new variable
file_name1 = sprintf('LFPch#%d.mat',tetChannel);
fprintf(1, 'Saving file\n'); 
save ([savePath filesep file_name1], xx, 'time','-v7.3') % saving spike data
clear('LFP')
fprintf(1, 'Done!\n'); 



%% Bandpass Filtering Ripple
lfpFilter = 'ripple'; % you can choose between 'ripple', 'lowGamma'
... 'fastGamma', 'theta'
LFPtofilt = ch33Tet15LFP;

filtLFP = bandFilter (LFPtofilt, lfpFilter); % default Fs = 1000
prompt = 'Assign a new variable name to your data ? '; %e.g. for a lfp from
... ch1 put ch1LFP for ch60 put ch60LFP or (ch# tet# LFP) 
xx = input(prompt, 's');
assignin('base',xx,filtLFP) % assignin the new variable
file_name1 = sprintf('%s.mat',xx);
fprintf(1, 'Saving file\n'); 
save ([savePath filesep file_name1], xx,'time','-v7.3') % saving spike data
clear('filtLFP','LFPtofilt', 'lfpFilter', 'xx')
fprintf(1, 'Done!\n'); 


%% Bandpass Filtering Spindle
lfpFilter = 'spindle'; % you can choose between 'ripple', 'lowGamma'
... 'fastGamma', 'theta'
LFPtofilt = ch33Tet15LFP;

filtLFP = bandFilter (LFPtofilt, lfpFilter); % default Fs = 1000
prompt = 'Assign a new variable name to your data ? '; %e.g. for a lfp from
... ch1 put ch1LFP for ch60 put ch60LFP or (ch# tet# LFP) 
xx = input(prompt, 's');
assignin('base',xx,filtLFP) % assignin the new variable
file_name1 = sprintf('%s.mat',xx);
fprintf(1, 'Saving file\n'); 
save ([savePath filesep file_name1], xx,'time','-v7.3') % saving spike data
clear('filtLFP','LFPtofilt', 'lfpFilter', 'xx')
fprintf(1, 'Done!\n'); 





%% Delta Filtering
Fs = 1000;
LFPtofilt = ch42Tet13LFP;

filtLFP = pfr_DeltaFilter(LFPtofilt, Fs);

prompt = 'Assign a new variable name to your data ? '; %e.g. for a lfp from
... ch1 put ch1LFP for ch60 put ch60LFP or (ch# tet# LFP) 
xx = input(prompt, 's');
assignin('base',xx,filtLFP) % assignin the new variable
file_name1 = sprintf('%s.mat',xx);
fprintf(1, 'Saving file\n'); 
save ([savePath filesep file_name1], xx,'time','-v7.3') % saving spike data
clear('filtLFP','LFPtofilt', 'xx')
fprintf(1, 'Done!\n'); 

%%
