
% Made with matlab R2018a by PFR 1/20/2020
% This function will merge all LFP recordings into a single file.
% Input:
% 1- tetChannel- raw LFP from a specific channel
% 2- ephyspaths- paths to raw LFPs. Cell Array containing all paths
% 3- reduction(Default 30), integer for downsampling the data. 
%   e.g. a reduction of 30 for a 30KHz recording will reduce the file into
%   a 1KHz or 1000Hz sampling rate file.
% 4- rawFs(Default 30000 or 30KHz) - raw sampling rate. Integer. 
% You need this to calculate the new sampling rate.

function [LFP, time] = pfr_LFP (tetChannel, ephysPaths,reduction, rawFs)

if size(ephysPaths,1) > size(ephysPaths,2) % ephysPaths can be a column/row cell array
 oe_dir = cell(size(ephysPaths,1));
    for k = 1 : max(size(ephysPaths))
    oe_dir{k} = dir([ephysPaths{k}]);
    end
else
   oe_dir = cell(size(ephysPaths,2),1);
    for k = 1 : max(size(ephysPaths))
    oe_dir{k} = dir([ephysPaths{k}]);
    end  
end

%%
oe_ch = cell(size(oe_dir,1),1);
oe_ch_filenames = cell(size(oe_dir,1),1);

for k = 1 : max(size(ephysPaths))    
    oe_ch{k} = cellfun(@any, regexp({oe_dir{k}.name},'\d.*?.continuous'));
    oe_ch_filenames{k} = sort(cellfun(@(x) [ephysPaths{k} filesep x], ...
    {oe_dir{k}(oe_ch{k}).name},'uni',false));
end
%%
formatSpec = 'CH%d.continuous';
st = sprintf(formatSpec,tetChannel);
% string2 = '\d.*?';
% st = strcat(string2,st1);
Idx = cell(1,size(ephysPaths,2));
for kk = 1 : max(size(ephysPaths))
    Idx{kk} = find(~cellfun(@isempty,regexp(oe_ch_filenames{kk}(1,:), st)));
end
%%
if isempty(Idx{1}) == 1
     warning('Make sure the LFP channel is located in one of the folders')
     return
else
    
end
%%
fprintf('Reading files!\n')
rawLFP = cell(size(oe_dir,1),1);
% rawTime = cell(size(oe_dir,1),1);
 for kk = 1 :max(size(ephysPaths))  
    [lfp] = load_open_ephys_data_faster(oe_ch_filenames{kk}{1,Idx{kk}});
    rawLFP{kk} = lfp;
%     rawTime{kk} = Time;
 end
 %%
rawLFP2 = cat(1,rawLFP{:});
 %%
 if nargin < 3
     rawFs = 30000;
     reduction = 30;
 end
LFP = downsample(rawLFP2,reduction);
newSR = 1/(rawFs/reduction);
time = 0: newSR : length(LFP) * newSR - newSR;
fprintf(1, 'Done!\n');

%         prompt = 'Assign a new variable name to your data ? ';
%         xx = input(prompt, 's');
%         assignin('base',xx,LFP) % assignin the new variable
%         fprintf(1, 'Saving Files!\n');
%         file_name1 = sprintf('LFPch#%d.mat',tetChannel);
%         save ([savePath filesep file_name1]) % saving spike data
%         % save ([savePath filesep file_name1], LFP, 'time','-v7.3') % saving spike data
%         fprintf(1, 'Done!\n');
end