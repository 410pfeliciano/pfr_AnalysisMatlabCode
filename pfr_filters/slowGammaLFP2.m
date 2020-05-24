%{
tetNum = tetrode number from 1 to 16
ChNum = channel number of interest for that tetrode
reduction 1 if you want to downsampling the LFP signal
RedParm = Decrease the sample rate of the sequence by this factor
%}
function [slowGammaLFP] = slowGammaLFP2(tetNum, ChNum, reduction, RedParm)
        flag = 'scale';  % Sampling Flag
        Fs = 30000/RedParm;  % Sampling Frequency
        N    = 500;      % Order
        win = blackman(N+1);
        RFc1 = 25;      % First Cutoff Frequency
        RFc2 = 55;      % Second Cutoff Frequency
if tetNum == 1
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH25.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH27.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH29.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH31.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
 
elseif tetNum == 2
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH17.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH19.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH21.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH23.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
 
elseif tetNum == 3
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH26.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH28.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH30.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH32.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
  
elseif tetNum == 4
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH18.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH20.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH22.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH24.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
 
elseif tetNum == 5
     if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH10.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH12.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH14.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH16.continuous');
     end
     LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end

elseif tetNum == 6
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH2.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH4.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH6.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH8.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
 
elseif tetNum == 7
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH9.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH11.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH13.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH15.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
    
elseif tetNum == 8
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH1.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH3.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH5.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH7.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
 
elseif tetNum == 9
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH9.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH11.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH13.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH15.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
 
elseif tetNum == 10
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH49.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH51.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH53.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH55.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
  
elseif tetNum == 11
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH58.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH60.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH62.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH64.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
 
elseif tetNum == 12
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH34.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH36.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH38.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH40.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
   
elseif tetNum == 13
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH42.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH44.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH46.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH48.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
   
elseif tetNum == 14
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH41.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH43.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH45.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH47.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
   
elseif tetNum == 15
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH33.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH35.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH37.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH39.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
    
elseif tetNum == 16
    if ChNum == 1
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH50.continuous');
      elseif ChNum == 2
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH52.continuous');
      elseif ChNum == 3
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH54.continuous');
      elseif ChNum == 4
        [LFP, LFPt] = load_open_ephys_data_faster('102_CH56.continuous');
    end
    LFP = detrend(LFP,'linear',5); % detrend the original LFP
    if ~exist('reduction','var')
     % third parameter does not exist, so default it to something
      reduction = 42;
      RedParm = 42;
    elseif reduction == 1
        LFP = downsample(LFP,RedParm);
        LFPt = downsample(LFPt,RedParm);
        Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
        slowGammaHd = dfilt.dffir(Rb); % Ripple Filter
        slowGammaLFP = filtfilt(slowGammaHd.Numerator,1,LFP);
    end
 
end