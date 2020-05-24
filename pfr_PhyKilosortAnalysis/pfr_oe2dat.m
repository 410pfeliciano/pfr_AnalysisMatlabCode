function pfr_oe2dat(oe_path,save_path, tetConfig)
%{ 
Created by PFR2019 with code from https://github.com/cortex-lab/spikes
Converts OE format recorded from open-ephys into flat binary
Instructions:
1- oe_path - path from all the folders with open-ephys data to convert. 
    Important! oe_path needs to be a cell array
    e.g oe_path = {'C/PATH1', 'C/PATH2', etc....};
2- save_path - path to save converted data
3- Tet_config should be a vector with the tet/channel configuration from
   open-ephys mapping. 
%}

%% Set paths and memory
if size(oe_path,1) > size(oe_path,2) % oe_path can be a column/row cell array
 oe_dir = cell(size(oe_path,1));
    for k = 1 : max(size(oe_path))
    oe_dir{k} = dir([oe_path{k}]);
    end
else
   oe_dir = cell(size(oe_path,2));
    for k = 1 : max(size(oe_path))
    oe_dir{k} = dir([oe_path{k}]);
    end  
end
%
oe_ch = cell(size(oe_dir,1),1);
oe_ch_filenames = cell(size(oe_dir,1),1);

for k = 1 : max(size(oe_path))
    
    oe_ch{k} = cellfun(@any, regexp({oe_dir{k}.name},'CH\d.*?.continuous'));
    oe_ch_filenames{k} = sort(cellfun(@(x) [oe_path{k} filesep x], ...
    {oe_dir{k}(oe_ch{k}).name},'uni',false));

end

% Load in first channel to set memory restrictions
fprintf(1, 'reading first channels for initialization\n');
data = cell(size(oe_dir,1),1);
totalSamps = zeros(size(oe_ch,1),1);
% Samps = cell(size(oe_dir,1),1);
for k = 1 : max(size(oe_path))

    data{k} = load_open_ephys_data_faster(oe_ch_filenames{k}{1});
    totalSamps(k) = length(data{k});
    
end
  
dmem = memory; 
memToLeaveFree = 4 * 2^30; % num of GB to keep free
memToAllocate = dmem.MemAvailableAllArrays - memToLeaveFree;
memToAllocate = max(0, memToAllocate);
nint16s = memToAllocate/2;

%% Save electrophysiology (CH) channels
outFile = [save_path filesep 'ephys.dat']; % filesep = \
fid = fopen(outFile, 'w'); % opening for writing only
% DataArrange = cell(size(oe_dir,1),1);
for k = 1 : max(size(oe_path))

nCH = length(oe_ch_filenames{k});
chunkSizeSamps = nint16s/nCH;
nChunks = ceil(totalSamps(k)/chunkSizeSamps);

    for chunkInd = 1:nChunks
        fprintf(1, 'chunck %d/%d\n', k ,max(size(oe_path)) );
        chunkStart = (chunkInd-1)*chunkSizeSamps+1;
        chunkEnd = min(chunkInd*chunkSizeSamps, totalSamps(k));
        allData = zeros(nCH, chunkEnd-chunkStart+1, 'int16');

        for n = 1:nCH
            
            fprintf(1, 'reading CH %d/%d\n', n, nCH);
            data = load_open_ephys_data_faster(oe_ch_filenames{k}{n});        
            allData(n,:) = int16(data);
            dataArrange = allData(tetConfig,:); % sort channels by tetrodes

        end
        fprintf(1, 'Writing dat file\n');
        fwrite(fid, dataArrange, 'int16');
        fprintf(1, 'done with folder %d\n', k)
    end
    
end
%%
% fprintf(1, 'Concatenating Matrices\n');
% allData = cat(2,DataArrange{:});
% fprintf(1, 'Writing dat file\n');
% fwrite(fid, allData, 'int16');
fclose(fid);
fprintf(1, 'done\n')
fprintf(1, 'dat file finished\n');

end