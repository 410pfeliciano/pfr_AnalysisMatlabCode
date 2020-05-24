% Function to detect Ripples
function [SWRrate, meanSWRrateperTrials, HistoIRI, ripplesidx] = rippleDetect(LFP, Fs)
%% Filtering the data
Fs = 2000;  % Sampling Frequency
N    = 500;      % Order
TFc1  = 4;        % First Cutoff Frequency
TFc2  = 12;       % Second Cutoff Frequency
flag = 'scale';  % Sampling Flag
win = blackman(N+1);
% Calculate the coefficients using the FIR1 function.
Tb  = fir1(N, [TFc1 TFc2]/(Fs/2), 'bandpass', win, flag);
ThetaHd = dfilt.dffir(Tb); %Theta filter

RFc1  = 120;      % First Cutoff Frequency
RFc2  = 250;      % Second Cutoff Frequency
Rb  = fir1(N, [RFc1 RFc2]/(Fs/2), 'bandpass', win, flag);
RippleHd = dfilt.dffir(Rb); % Ripple Filter

DFc   = 4;
Db  = fir1(N, DFc/(Fs/2), 'low', win, flag);
DeltaHd = dfilt.dffir(Db);

sint = 1/Fs;
LFP = detrend(LFP,'linear',5); % detrend the original LFP


filtLFP = filtfilt(RippleHd.Numerator,1,LFP);
% ThetaLFP = filtfilt(ThetaHd.Numerator,1,LFP);
SharpwaveLFP = filtfilt(DeltaHd.Numerator,1,LFP);
time = 0:sint:length(LFP)*sint-sint;

%%

filtLFP = filtLFP - mean(filtLFP,1); % filtered LFP detrend
% SharpwaveLFP = SharpwaveLFP - mean(SharpwaveLFP,1); % filtered LFP detrend

    
%% Second Step, Calculation the mean power and standar deviation of the recording

hibfiltLFP = abs(hilbert(filtLFP)); % amplitude envelope of ripple-band filtered LFP
sqfiltLFP= hibfiltLFP.^2; % calculating the power of the filtered LFP
smsqfiltLFP = smooth(sqfiltLFP,10); % smoothing the power data
% m_sqfiltLFP = median(sqfiltLFP,1); % calculating the average power of smoothed data
% sd_sqfiltLFP = std(smsqfiltLFP,1); % calculating the SD of power of smoothed data
threshold = 2.5 * mean(median(smsqfiltLFP ./ .6745), 2);

%%  FIRST CRITERION: AT LEAST 3*STD OF FILTERED, SQUARED LFP


splog = smsqfiltLFP > threshold;
pos = find(splog == 1); % This vectors contain the indices of power larger the 3std


%% SECOND CRITERION: at least a difference of 15ms between the end of
...a ripple and the beginning of the next to be considered separate events

    difpos = find(diff(pos) >= (0.02/sint)); % events with inte4rvals > 20ms
    % end of all possible ripple events
    z = zeros(length(difpos)+1,1); % Creating a zero row vector
    z(1:end-1) = difpos;
    z(end) = length(pos); % End position of all potential ripples
    % Calculate the start
    st = difpos+1; % Start positions of all possible ripple events
    zst = zeros(length(difpos)+1,1);
    zst(2:end) = st; %
    zst(1)=1; % zst start position of all
    stepripples = [pos(zst), pos(z)]; % Start and End position(indices) 
        ... of Ripples Candidates


%% THIRD CRITERION: Eliminate 'Ripples' candidates with < sampling interval

duration = stepripples(:,2)- stepripples(:,1);
stepripples1 = stepripples((~(duration<=1)),:);

%% FOURTH CRITERION:  at least 20ms long

duration2 = (stepripples1(:,2) - stepripples1(:,1))*sint; %duration in sec.
ripplesidx = stepripples1((find(duration2<=0.15 & ...
    duration2 >= 0.015)),:); %Events < 150ms and > 15ms (Indices)
ripplessec = ripplesidx*sint; % Time(sec) End and Start of selected ripples
% startripplessec(kk,1) = zeros(length(ripplesidx{kk,1}) +1 ,1);
% startripplessec(kk,1)(2:end) = ripplessec(kk,1); % Onlt the Start of Ripple Candidates

%% Second Step
below = ripplesidx(:,1)-(0.01/sint);
abov = ripplesidx(:,2)+(0.01/sint);
below(below < 0) = 1; % Substitute start negative value with 1
abov(abov > length(LFP)) = length(LFP); % Making ripple candidates finishing at the end of the LFP
belowidx = find(linVel(below) < 0.03);
below = below(belowidx);
abov = abov(belowidx);
%% Visual Inspection
figure
subplot(211)
plot(linTime, linPos,'k', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('Position [m]')
axis([515 530 0 inf])
subplot(212)
plot(LFPt, LFP)
hold on
for kk = 1: length(below)
    plot(LFPt(below(kk):abov(kk)),LFP(below(kk):abov(kk)), 'r', 'LineWidth',2) 
end
plot(LFPt, filtLFP,'k','LineWidth',1)
plot(LFPt, smsqfiltLFP./80,'m', 'LineWidth',2)
hold off
axis([515 530 -500 600])


%% Visual Inspection all the ripple candidates LFPs

for jj = 1: size(below,1)
    figure
    tm = 0:sint:length(ripplesCandLFPs{kk,1}{jj,1})*sint-sint;
    plot(tm, ripplesCandLFPs{kk,1}{jj,1})
    axis tight
    hold on
    plot(tm, ripplesCandfiltLFPs{kk,1}{jj,1})
    plot(tm, ripplesCandSharpwaveLFPs{kk,1}{jj,1})
    plot(tm, ripplesCandRippleHibLFPs{kk,1}{jj,1}, 'LineWidth',2)
    
    xlabel('Time [s]');			%... and label axes.
    ylabel('Voltage [uV]');
    hold off
end

%% Gaussian fitting of ripple emvelop



%% Power Spectra Calulation in the ripple of Ripple Candidates
% FreqPowRipple = [];
% TimePowRipple = [];
% PowerRippleCand = [];
for kk = 1:length(LFP2)
    for jj = 1:length(ripplesCandfiltLFPs{kk,1}(:,1))
        interval = fix(length(ripplesCandfiltLFPs{kk,1}{jj,1})./ 4); %Specify the interval size.
        overlap = fix(interval * 0.98);	%Specify the overlap of intervals.
        nfft = round((length(ripplesCandfiltLFPs{kk,1}{jj,1}))/2); %Specify the FFT length.
        [S,F,T,P] = spectrogram(ripplesCandfiltLFPs{kk,1}{jj,1}...
            - mean(ripplesCandfiltLFPs{kk,1}{jj,1}), interval,overlap, nfft, Fs, 'psd');
        
        FreqPowRipple{kk,1}{jj,1} = F;
        TimePowRipple{kk,1}{jj,1}  = T;
        PowerRippleCand{kk,1}{jj,1} = P;
    end
end

%% Visual Inspection of Ripple Candidates Spectrograms

for kk = 1:length(LFP)
    for jj = 1:length(ripplesCandfiltLFPs{kk}(:,1))
        figure
        imagesc(TimePowRipple{kk,1}{jj,1},FreqPowRipple{kk,1}{jj,1},...
            10*log10(PowerRippleCand{kk,1}{jj,1}));	%... and plot it,
        colormap(copper(5));		%... with a colorbar,
        colormap;
        caxis([0 50]);
        axis xy	;					%... and origin in lower left, 
        ylim([0 500]);				%... set the frequency range,
        xlabel('Time [s]');			%... and label axes.
        % view(-45,65)
        ylabel('Frequency [Hz]');
    end
end

%% FIFTH CRITERION: A Frequency Peak Power larger than 125Hz

for kk = 1:length(LFP)
    for jj = 1:length(PowerRippleCand{kk,1}(:,1))
            maxPowerFreqDist{kk,1}{jj,1} =  max(PowerRippleCand{kk,1}{jj,1},[],2); % Find max value in rows
%         [val, idx] = max(maxPowerFreqDist);
    end
end

% visual inspection of Frequency Peak

for kk = 1:length(LFP)
    for jj = 1:length(maxPowerFreqDist{kk,1}(:,1))
        figure
        plot(FreqPowRipple{kk,1}{jj,1},maxPowerFreqDist{kk,1}{jj,1})
        xlim([0 350]);
        xlabel('Frequency [Hz]');
        ylabel('PSD [db/Hz]');
    end
end

% Extracting Max Frequency Peak for all Ripple Candidates

for kk = 1:length(LFP)
    for jj = 1:length(maxPowerFreqDist{kk,1}(:,1))
            [MaxVal{kk}{jj,1}, MaxIdx{kk}{jj,1}] =  max(maxPowerFreqDist{kk}{jj,1},[],1); % Find max value in column
            maxFreqCandidate{kk}{jj,1} = FreqPowRipple{kk,1}{jj,1}(MaxIdx{kk}{jj,1});
    end
end

% Selection of Ripples with Frequency Peak larger or less than 125-250Hz

for kk = 1:length(LFP)
    for jj = 1:length(maxPowerFreqDist{kk,1}(:,1))
        ripplesidx2{kk}{jj,1} = ripplesidx{kk}((maxFreqCandidate{kk}{jj,1} >= 140),:);  % selection of ripples with more the 15ms duration
    end
end

for kk = 1:length(LFP)
ripplesidx2{kk}(cellfun('isempty', ripplesidx2{kk})) = []; 
end
%% Inter Ripple Interval Histogram calculation

for kk = 1:length(LFP) 
    interRippleInterval{kk} = diff(startripplessec{kk}(:,1));
end
    % [nrows,ncols] = cellfun(@size, interRippleInterval); % getting the elemnts in the cell array
    % cellElements = sum(nrows); % total number of elements in the cell array
    % IRI = zeros(cellElements,1); % zero vector
    IRI = cell2mat(interRippleInterval'); % all Inter Ripple Interval of the cell into a single column vector
    bins = 0:0.1:max(IRI);
    HistoIRI = histogram(IRI, bins);
    figure(1)
    histogram(IRI, bins);
    set(gca,'fontsize', 15)
    title('Inter Ripple Interval')
    xlabel('Time sec')
    ylabel('Ripple Counts')
%     xlim([0 5])

%% Calculating the Rate of Ripple Candidates

    [nrows,ncols] = cellfun(@size, ripplessec); % getting the elemnts in the cell array
    RippleNumber = sum(nrows);
    
for kk = 1:length(LFP) 
    TimeTrials(kk,1) = time{kk}(end);
end

    TotalTime = sum(TimeTrials);
    SWRrate = RippleNumber./TotalTime;
    
for kk = 1:length(LFP) 
   SWRrateperTrials(kk,1) = length(ripplessec{kk})./time{kk}(end);
end
    meanSWRrateperTrials = mean(SWRrateperTrials,1);
    sdtSWRrateperTrials = std(SWRrateperTrials,1);
    

