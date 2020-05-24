function [OEspkWave, spkTimeStamps] = jclustSpikes (tetrodeNumber, SpkPaths, savePath)
if size(SpkPaths,1) > size(SpkPaths,2) % SpkPaths can be a column/row cell array
 oe_dir = cell(size(SpkPaths,1));
    for k = 1 : max(size(SpkPaths))
    oe_dir{k} = dir([SpkPaths{k}]);
    end
else
   oe_dir = cell(size(SpkPaths,2),1);
    for k = 1 : max(size(SpkPaths))
    oe_dir{k} = dir([SpkPaths{k}]);
    end  
end

%%
oe_ch = cell(size(oe_dir,1),1);
OEspikes_ch = cell(size(oe_dir,1),1);
oe_ch_filenames = cell(size(oe_dir,1),1);
OEspikes_ch_filenames = cell(size(oe_dir,1),1);
for k = 1 : max(size(SpkPaths))
    
    oe_ch{k} = cellfun(@any, regexp({oe_dir{k}.name},'\d.*?.continuous'));
    oe_ch_filenames{k} = sort(cellfun(@(x) [SpkPaths{k} filesep x], ...
    {oe_dir{k}(oe_ch{k}).name},'uni',false));
    OEspikes_ch{k} = cellfun(@any, regexp({oe_dir{k}.name},'0n\d.*?.spikes'));
    OEspikes_ch_filenames{k} = sort(cellfun(@(x) [SpkPaths{k} filesep x], ...
    {oe_dir{k}(OEspikes_ch{k}).name},'uni',false));

end
%% Calculating the recording time in each experiment
fprintf(1, 'reading continous channels for time calculation\n');
timeRecorded = cell(size(oe_dir,1),1);
totalTime = zeros(size(oe_ch,1),1);
% Samps = cell(size(oe_dir,1),1);
for k = 1 : max(size(SpkPaths))

    [~, timeRecorded{k}] = load_open_ephys_data_faster(oe_ch_filenames{k}{1});
    totalTime(k) = timeRecorded{k}(end);
    
end

%%
if tetrodeNumber == 1
        st = '\d.*?0n0.spikes';
    elseif tetrodeNumber == 2
        st = '\d.*?0n1.spikes';
    elseif tetrodeNumber == 3
        st = '\d.*?0n2.spikes';
    elseif tetrodeNumber == 4
        st = '\d.*?0n3.spikes';
    elseif tetrodeNumber == 5
        st = '\d.*?0n4.spikes';
    elseif tetrodeNumber == 6
        st = '\d.*?0n5.spikes';
    elseif tetrodeNumber == 7
        st = '\d.*?0n6.spikes';
    elseif tetrodeNumber == 8
        st = '\d.*?0n7.spikes';
    elseif tetrodeNumber == 9
        st = '\d.*?0n8.spikes';
    elseif tetrodeNumber == 10
        st = '\d.*?0n9.spikes';
    elseif tetrodeNumber == 11
        st = '\d.*?0n10.spikes';
    elseif tetrodeNumber == 12
        st = '\d.*?0n11.spikes';
    elseif tetrodeNumber == 13
        st = '\d.*?0n12.spikes';
    elseif tetrodeNumber == 14
        st = '\d.*?0n13.spikes';
    elseif tetrodeNumber == 15
        st = '\d.*?0n14.spikes';
    elseif tetrodeNumber == 16
        st = '\d.*?0n15.spikes';                
end
OEspikes = cell(size(SpkPaths,2),1);
timeStamps = cell(size(SpkPaths,2),1);
Idx = find(~cellfun(@isempty,regexp(OEspikes_ch_filenames{1,1}, st)));
%%
 for kk = 1 :max(size(SpkPaths))  
    [spikes, time_stamp] = load_open_ephys_data_faster(OEspikes_ch_filenames{kk}{Idx});
    OEspikes{kk} = permute(spikes, [3 2 1]);
    timeStamps{kk} = time_stamp;
end
%%
    totalTimeCumSum = cumsum(totalTime); % cumulative time vector
    spkTimeStp = cell(size(timeStamps));
    spkTimeStp{1} = timeStamps{1};
    for kk =2 : max(size(SpkPaths))
        spkTimeStp{kk} = timeStamps{kk} + totalTimeCumSum(kk-1);
    end
%% Merging the spiking and time stamp vector to a single file
OEspkWave = cat(3,OEspikes{:,1});
spkTimeStamps = cat(1,spkTimeStp{:});
fprintf(1, 'Finished!\n');
%%
fprintf(1, 'Saving Files!\n');
file_name1 = sprintf('spkWave&StampsTet#%d.mat',tetrodeNumber);
save ([savePath filesep file_name1], 'OEspkWave', 'spkTimeStamps','-v7.3') % saving spike data
fprintf(1, 'Done!\n');
end