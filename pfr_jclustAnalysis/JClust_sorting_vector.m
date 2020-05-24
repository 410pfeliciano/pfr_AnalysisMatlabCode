% Creating Vectors for Spike Sorting from Open Ephys Files
% This function depends on a particular Channel map configuration from the
% Open Ephys GUI.
% Channel map configuration:
% Tet1= 25, 27, 29, 31; Tet2= 17, 19, 21, 23; Tet3= 26, 28, 30, 32;
% Tet4= 18, 20, 22, 24; Tet5= 10, 12, 14, 16; Tet6= 2, 4, 6, 8;
% Tet7= 9, 11, 13, 15; Tet8= 1, 3, 5, 7; Tet9= 57, 59, 61, 63;
% Tet10= 49, 51, 53, 55; Tet11= 58, 60, 62, 64, Tet12= 34, 36, 38, 40;
% Tet13= 42, 44, 46, 48; Tet14= 41, 43, 45, 47;
% Tet15= 33, 35, 37, 39; Tet16 = 50, 52, 54, 56;
% 
% tetrode_num = tetrode number
% damage = vector with damaged channel(s) for the tetrode of interest
% Damaged = 0 is default
% If you have channels that were damaged, put the # of the channel(s):
%   e.g. If channel(s) 1, 3, 4 are damaged the vector should be [1, 3, 4]
%   e.g. If channel(s) 1 was damaged, just enter 1.

function [tet, spike, time_stamp] = JClust_sorting_vector(tetrode_num, damage)
if nargin < 2
      damage = 0;
end
if tetrode_num == 1
    if isfile('102_CH25.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH25.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH27.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH27.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH29.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH29.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH31.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH31.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channels')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n0.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   

    elseif tetrode_num == 2
    if isfile('102_CH17.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH17.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH19.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH19.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH21.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH21.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH23.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH23.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n1.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   

    elseif tetrode_num == 3
    if isfile('102_CH26.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH26.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH28.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH28.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH30.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH30.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH32.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH32.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
  
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n2.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   



    elseif tetrode_num == 4
    if isfile('102_CH18.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH18.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH20.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH20.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH22.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH22.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH24.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH24.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
    
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n3.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   elseif tetrode_num == 5
    if isfile('102_CH10.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH10.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH12.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH12.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH14.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH14.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH16.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH16.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
    
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n4.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   elseif tetrode_num == 6
    if isfile('102_CH2.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH2.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH4.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH4.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH6.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH6.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH8.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH8.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
    
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n5.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   elseif tetrode_num == 7
    if isfile('102_CH9.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH9.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH11.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH11.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH13.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH13.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH15.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH15.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
    
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n6.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   elseif tetrode_num == 8
    if isfile('102_CH1.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH1.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH3.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH3.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH5.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH5.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH7.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH7.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght,1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
        
    tet = [Ch1, Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n7.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   elseif tetrode_num == 9
    if isfile('102_CH9.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH9.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH11.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH11.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH13.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH13.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH15.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH15.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
    
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n8.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   elseif tetrode_num == 10
    if isfile('102_CH49.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH49.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH51.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH51.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH53.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH53.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH55.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH55.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
    
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n9.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   elseif tetrode_num == 11
    if isfile('102_CH58.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH58.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH60.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH60.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH62.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH62.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH64.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH64.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
    
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n10.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   elseif tetrode_num == 12
    if isfile('102_CH34.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH34.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH36.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH36.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH38.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH38.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH40.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH40.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
    
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n11.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   elseif tetrode_num == 13
    if isfile('102_CH42.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH42.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH44.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH44.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH46.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH46.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH48.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH48.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
    
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n12.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   elseif tetrode_num == 14
    if isfile('102_CH41.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH41.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH43.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH43.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH45.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH45.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH47.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH47.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
    
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n13.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   elseif tetrode_num == 15
    if isfile('102_CH33.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH33.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH35.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH35.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH37.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH37.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH39.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH39.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
    
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n14.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';

   elseif tetrode_num == 16
    if isfile('102_CH50.continuous')
    Ch1 = load_open_ephys_data_faster('102_CH50.continuous');
    tetlenght = length(Ch1);
    else
     % File does not exist.
    end
    if isfile('102_CH52.continuous')
    Ch2 = load_open_ephys_data_faster('102_CH52.continuous');
    tetlenght = length(Ch2);
    else
     % File does not exist.
    end
    if isfile('102_CH54.continuous')
    Ch3 = load_open_ephys_data_faster('102_CH54.continuous');
    tetlenght = length(Ch3);
    else
     % File does not exist.
    end
    if isfile('102_CH56.continuous')
    Ch4 = load_open_ephys_data_faster('102_CH56.continuous');
    tetlenght = length(Ch4);
    else
     % File does not exist.
    end
    if exist('Ch1','var')
        % If exist do nothing
    else
        Ch1 = zeros(tetlenght, 1);
    end
    if exist('Ch2','var')
        % If exist do nothing
    else
        Ch2 = zeros(tetlenght, 1);
    end
    if exist('Ch3','var')
        % If exist do nothing
    else
        Ch3 = zeros(tetlenght, 1);
    end
    if exist('Ch4','var')
        % If exist do nothing
    else
        Ch4 = zeros(tetlenght, 1);
    end
   
    tet = [Ch1 Ch2, Ch3, Ch4]'; % vector continuous data
    if isrow(damage) && damage > 0
     for c = damage % converting damaged channel(s) to zeros
          tet(c,:) = zeros(1,length(tet));
     end
     disp('Damaged channel(s) converted to 0.')
    else
     disp('Great! Tetrode with no damaged channel(s)')
    end
    [spike, time_stamp] = load_open_ephys_data_faster('TTp111.0n15.spikes');
    spike = permute(spike, [3 2 1]);
    time_stamp = time_stamp';
end

file_name1 = sprintf('spikes_wave_tet#%d.mat',tetrode_num);
file_name2 = sprintf('raw_signal_tet#%d.mat',tetrode_num);
save (file_name1, 'spike', 'time_stamp') % saving spike data
save(file_name2, 'tet','-v7.3') % Raw wideband data

end
