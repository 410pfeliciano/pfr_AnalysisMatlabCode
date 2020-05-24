% Obtained the whole LFP. Add the all directory paths of the experiments
ephysPaths ={'D:\Grik4_Ai40D_2\10-31-19\Grik4_Ai40D_2_Sleep_2019-10-31_17-38-34', ...
    'D:\Grik4_Ai40D_2\10-31-19\Grik4_Ai40D_2_Sleep_2019-10-31_17-44-41'};

savePath = 'D:\Grik4_Ai40D_2\10-31-19\lfps'; % specify the saving path
% Use the save path to save your data to a specific directory

%% whole LFP function
% rawFs = 30000; % this value is default in the pfr_LFP function
% reduction = 30; % this value is default in the pfr_LFP function
tetChannel = 32;
[LFP, time] = pfr_LFP(tetChannel, ephysPaths);

prompt = 'Assign a new variable name to your data ? '; %e.g. for a lfp from
... ch1 put ch1LFP for ch60 put ch60LFP or (ch# tet# LFP)
    
xx = input(prompt, 's');
assignin('base',xx,LFP) % assignin the new variable
file_name1 = sprintf('LFPch#%d.mat',tetChannel);
fprintf(1, 'Saving file\n'); 
save ([savePath filesep file_name1], xx, 'time','-v7.3') % saving spike data
clear('LFP')
fprintf(1, 'Done!\n'); 