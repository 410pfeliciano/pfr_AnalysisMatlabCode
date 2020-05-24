%% Peak Stimulation Data
clear
%%
LEDoutPK = load_open_ephys_data_faster('102_ADC6.continuous'); % Load voltage LED output
downsamplefactor = 30;
LEDPK = downsample(LEDoutPK,downsamplefactor); % Reduce sampling rate
Optodur = 0.01; % Duration of Stimulation
SampRate = 30000/downsamplefactor; % 30KHz is the old sampling rate
SampInt = 1/SampRate;
diffLEDoutPK = zeros(size(LEDPK,1),1); %Adding a sample to the derivative vector 
diffLEDoutPK(2:end, 1) = diff(LEDPK);
[LFPPK,LFPtPK,thetaLFPPK] = ThetaLFP(1, 1,1,downsamplefactor); 
phaseThetaPK =  wrapToPi(angle(hilbert(thetaLFPPK)));
pksStartPK = find(diffLEDoutPK > 0.4); % first stimulus sample indice
pksEndPK = pksStartPK + Optodur*SampRate;
%% Defining Stimulation Trials
diffPeakStart = diff(pksStartPK);
binaryTrials = zeros(size(pksStartPK,1),1);
binaryTrials(2:end,1) = diffPeakStart;
binaryTrials(binaryTrials<200)=0; 
binaryTrials(binaryTrials>200)=1;
StartTrials = find(binaryTrials);
EndTrials = zeros(size(StartTrials,1),1);
EndTrials(1,1) = StartTrials(2,1) - 1;
EndTrials(2:end-1,1) = StartTrials(3:end) - 1; 
EndTrials(end,1) = size(pksStartPK,1);
PeakTrials = cell(size(StartTrials,1),1);
rowElimination = find((EndTrials-StartTrials) < 5); %Eliminate Trials with less than 5 stimuli
StartTrials(rowElimination) = [];
EndTrials(rowElimination) = [];
for kk = 1: size(StartTrials,1)
    PeakTrials{kk} = pksStartPK(StartTrials(kk) : EndTrials(kk),1);
end
% PeakTrials(cellfun('length', PeakTrials)< 5) = []; %Eliminate Trials with less than 5 stimuli
PeakStartEnd = [pksStartPK(StartTrials,1), pksStartPK(EndTrials,1)];

LFPtrials = cell(size(StartTrials,1),1);
timetrials = cell(size(StartTrials,1),1);
for kk = 1: size(StartTrials,1)
    LFPtrials{kk,1} = detrend(LFPPK(PeakStartEnd(kk,1)-(0.1/SampInt): ...
    PeakStartEnd(kk,2)+(0.1/SampInt)));
    timetrials{kk,1} = LFPtPK(PeakStartEnd(kk,1)-(0.1/SampInt): ...
    PeakStartEnd(kk,2)+(0.1/SampInt));
end
%%
Thetatrials = cell(size(StartTrials,1),1);
for kk = 1: size(StartTrials,1)
    Thetatrials{kk,1} = thetaLFPPK(PeakStartEnd(kk,1)-(0.1/SampInt): ...
        PeakStartEnd(kk,2)+(0.1/SampInt));
end
%%
close all
for kk = 1 : size(StartTrials,1)
figure(kk)
plot(timetrials{kk}, LFPtrials{kk},'Color',[128/255 128/255 128/255])
hold on
plot(timetrials{kk}, Thetatrials{kk},'k','LineWidth',1)

for ii = 1 : size(PeakTrials{kk},1)
    x = [LFPtPK(PeakTrials{kk}(ii)) LFPtPK(PeakTrials{kk}(ii)+(0.01/SampInt))...
        LFPtPK(PeakTrials{kk}(ii)+(0.01/SampInt)) LFPtPK(PeakTrials{kk}(ii))];
    y = [-1000 -1000 1000 1000];
    patch(x,y,[210/255 105/255 30/255],'FaceAlpha',0.3,'EdgeColor','none')
end
   hold off 
   axis([-inf inf -2000 2000])
%    set(gca, 'visible', 'off');
   legend({'Raw LFP','Filt LFP(4-12Hz)','LED Driver output'},'Location',...
       'southwest','NumColumns',3)
   legend('boxoff')
end
%% Power Spectrum of LFP per Trials
% maxLength = max(cellfun('size',LFPtrials,1))*2;
maxLength = 10000;
for kk = 1 : size(LFPtrials,1)
    B{kk,1} = padarray(LFPtrials{kk},maxLength-length(LFPtrials{kk}),'pre');
%     B{kk,1} = padarray(B{kk,1},maxLength, 'pre');
end
newt = 0: SampInt : size(B{1},1)*SampInt-SampInt;

%% Power Spectrum with Taper
close all
TW = 5; % time-bandwidth product
LFPpspec = zeros(size(StartTrials,1), maxLength/2+1);
for kk = 1 : length(StartTrials)
    [LFPpspec(kk,:), freq] = pmtm(B{kk}, TW, maxLength, SampRate);
    
    hold on
    plot(freq, 10*log10(LFPpspec(kk,:)), 'Color',[192/255 192/255 192/255])
%     plot(faxis{kk}, Sxx{kk})
    axis([4 120 0 inf])
    xlabel('Frequency [Hz]')
    ylabel('Power [dB]')
end
% LFPpspecAvg = sum(LFPpspec,1)/size(StartTrials,1); % can also be done by using mean function
% LFPpspecSD = std(LFPpspec,1); % Std across trials
% LFPpspecSqt = LFPpspecSD/sqrt(size(StartTrials,1)); % std of the mean
% plot(freq, 10*log10(LFPpspecAvg),'Color',[210/255 105/255 30/255], 'LineWidth',2)
% plot(freq, 10*log10(LFPpspecAvg+2*LFPpspecSqt),':','Color',[210/255 105/255 30/255],'LineWidth',1)
% plot(freq, 10*log10(LFPpspecAvg-2*LFPpspecSqt),':','Color',[210/255 105/255 30/255],'LineWidth',1)
options.handle = figure(1);
options.color_area = [210/255 105/255 30/255]; 
options.error = 'c95';
options.color_line = [210/255 105/255 30/255];
options.line_width = 2;
options.alpha = 0.5; 
options.x_axis = freq;              
plot_areaerrorbar(10*log10(LFPpspec), options) % Function from matlab
hold off
title('LFP Power Spectrum Retrieval Tet7')
%%
clear % Clear the workspace


%% Trough Stimulation LFP
LEDoutTR = load_open_ephys_data_faster('102_ADC6.continuous'); % Load voltage LED output
downsamplefactor = 30;
LEDTR = downsample(LEDoutTR,downsamplefactor); % Reduce sampling rate
Optodur = 0.01; % Duration of Stimulation
SampRate = 30000/downsamplefactor; % 30KHz is the old sampling rate
SampInt = 1/SampRate;
diffLEDoutTR = zeros(size(LEDTR,1),1); %Adding a sample to the derivative vector 
diffLEDoutTR(2:end, 1) = diff(LEDTR);
[LFPTR,LFPtTR,thetaLFPTR] = ThetaLFP(1, 1,1,downsamplefactor); 
phaseThetaTR =  wrapToPi(angle(hilbert(thetaLFPTR)));
pksStartTR = find(diffLEDoutTR > 0.4); % first stimulus sample indice
pksEndTR = pksStartTR + Optodur*SampRate;
%% Defining Stimulation Trials
diffPeakStart = diff(pksStartTR);
binaryTrials = zeros(size(pksStartTR,1),1);
binaryTrials(2:end,1) = diffPeakStart;
binaryTrials(binaryTrials<200)=0; 
binaryTrials(binaryTrials>200)=1;
StartTrials = find(binaryTrials);
EndTrials = zeros(size(StartTrials,1),1);
EndTrials(1,1) = StartTrials(2,1) - 1;
EndTrials(2:end-1,1) = StartTrials(3:end) - 1; 
EndTrials(end,1) = size(pksStartTR,1);
rowElimination = find((EndTrials-StartTrials) < 5); %Eliminate Trials with less than 5 stimuli
StartTrials(rowElimination) = [];
EndTrials(rowElimination) = [];
PeakTrials = cell(size(StartTrials,1),1);
for kk = 1: size(StartTrials,1)
    PeakTrials{kk} = pksStartTR(StartTrials(kk) : EndTrials(kk),1);
end
% PeakTrials(cellfun('length', PeakTrials)< 5) = []; %Eliminate Trials with less than 5 stimuli
PeakStartEnd = [pksStartTR(StartTrials,1), pksStartTR(EndTrials,1)];
LFPtrials = cell(size(StartTrials,1),1);
timetrials = cell(size(StartTrials,1),1);
for kk = 1: size(StartTrials,1)
    LFPtrials{kk,1} = detrend(LFPTR(PeakStartEnd(kk,1)-(0.1/SampInt): ...
    PeakStartEnd(kk,2)+(0.1/SampInt)));
    timetrials{kk,1} = LFPtTR(PeakStartEnd(kk,1)-(0.1/SampInt): ...
    PeakStartEnd(kk,2)+(0.1/SampInt));
end
%%
Thetatrials = cell(size(StartTrials,1),1);
for kk = 1: size(StartTrials,1)
    Thetatrials{kk,1} = thetaLFPTR(PeakStartEnd(kk,1)-(0.1/SampInt): ...
        PeakStartEnd(kk,2)+(0.1/SampInt));
end
%%
close all
for kk = 1 : size(StartTrials,1)
figure(kk)
plot(timetrials{kk}, LFPtrials{kk},'Color',[128/255 128/255 128/255])
hold on
plot(timetrials{kk}, Thetatrials{kk},'k','LineWidth',1)

for ii = 1 : size(PeakTrials{kk},1)
    x = [LFPtTR(PeakTrials{kk}(ii)) LFPtTR(PeakTrials{kk}(ii)+(0.01/SampInt))...
        LFPtTR(PeakTrials{kk}(ii)+(0.01/SampInt)) LFPtTR(PeakTrials{kk}(ii))];
    y = [-1000 -1000 1000 1000];
    patch(x,y,[4/255 101/255 53/255],'FaceAlpha',0.3,'EdgeColor','none')
end
   hold off 
   axis([-inf inf -2000 2000])
%    set(gca, 'visible', 'off');
   legend({'Raw LFP','Filt LFP(4-12Hz)','LED Driver output'},'Location',...
       'southwest','NumColumns',3)
   legend('boxoff')
end
maxLength = 10000;
for kk = 1 : size(LFPtrials,1)
    B{kk,1} = padarray(LFPtrials{kk},maxLength-length(LFPtrials{kk}),'pre');
%     B{kk,1} = padarray(B{kk,1},maxLength, 'pre');
end
newt = 0: SampInt : size(B{1},1)*SampInt-SampInt;

%% Power Spectrum with Taper
close all
TW = 5; % time-bandwidth product
LFPpspec = zeros(size(StartTrials,1), maxLength/2+1);
for kk = 1 : length(StartTrials)
    [LFPpspec(kk,:), freq] = pmtm(B{kk}, TW, maxLength, SampRate);
    
    hold on
    plot(freq, 10*log10(LFPpspec(kk,:)), 'Color',[192/255 192/255 192/255])
%     plot(faxis{kk}, Sxx{kk})
    axis([4 100 0 inf])
    xlabel('Frequency [Hz]')
    ylabel('Power [dB]')
end
options.handle = figure(1);
options.color_area = [4/255 101/255 53/255]; 
options.error = 'c95';
options.color_line = [4/255 101/255 53/255];
options.line_width = 2;
options.alpha = 0.5; 
options.x_axis = freq;              
plot_areaerrorbar(10*log10(LFPpspec), options) % Function from matlab
hold off
title('LFP Power Spectrum Retrieval')