%{ 
    Code for analyzing optogenetic Closed-loop Theta-phase Stimulation
    We need the following data set:
        1- Position
        2- LFP(Raw and Filtered in theta)
        3- LED Driver voltage output
        4- Phase of Filtered LFP
%}

%% First Stimulus detection and timestamps
clear
LEDout = load_open_ephys_data_faster('102_ADC6.continuous'); % Load voltage LED output
% LFPin = load_open_ephys_data_faster('114_CH25.continuous'); % Load voltage LED output
downsamplefactor = 30;
LED = downsample(LEDout,downsamplefactor); % Reduce sampling rate
% LFPGUIin = downsample(LFPin,downsamplefactor); % Reduce sampling rate
Optodur = 0.01; % Duration of Stimulation
SampRate = 30000/downsamplefactor; % 30KHz is the old sampling rate
SampInt = 1/SampRate;
diffLEDout = zeros(size(LED,1),1); %Adding a sample to the derivative vector 
diffLEDout(2:end, 1) = diff(LED);
%% Theta, Low and High Gamma LFPs and Phase
[inputLFP,LFPt,thetaLFP] = ThetaLFP(1, 1,1,downsamplefactor); 
%% GUI Bandpass Filter for Comparison
... https://open-ephys.atlassian.net/wiki/spaces/OEW/pages/950415/Bandpass+Filter
lowcut = 6;
highcut = 10;
[b, a] = butter(2, [lowcut highcut]/(SampRate/2));
filtLFP = filter(b, a, inputLFP); % only filters in the forward direction

% slowGammaLFP = slowGammaLFP2(16, 4, 1, downsamplefactor);
% fastGammaLFP = fastGammaLFP2(16, 4, 1, downsamplefactor);
phaseTheta =  wrapToPi(angle(hilbert(thetaLFP)));
phaseTheta2 =  wrapToPi(angle(hilbert(filtLFP)));
% phaseLgamma =  wrapToPi(angle(hilbert(slowGammaLFP)));
% phaseFgamma =  wrapToPi(angle(hilbert(fastGammaLFP)));
%% Find the Events from the Open-ephys events file
[OEevents, OEtime, info] = load_open_ephys_data_faster('all_channels.events'); % load voltage LED output
events = find(OEevents < 1);
events1 = sub2ind(OEtime(events)*SampRate,1);
%% 
close all
figure
plot(LFPt, LED-5)
hold on
plot(LFPt, diffLEDout./2,'Linewidth',1)
% plot(LFPt, phaseTheta./25, 'm', 'Linewidth',1 )
% plot(LFPt, phaseTheta2./25, 'g', 'Linewidth',1 )
plot(OEtime(events), OEevents(events),'b.','MarkerSize',12)
% plot(LFPt, inputLFP/5000, 'k')
plot(LFPt, thetaLFP/2000, 'r', 'Linewidth',2 )
plot(LFPt, filtLFP/2000, 'b', 'Linewidth',2 )
% plot(LFPt, LFPGUIin/2000, 'g', 'Linewidth',1 )
hold off
axis([0 inf -1 0.4]);

%% 
close all
pksStart = find(diffLEDout > 0.4); % first stimulus sample indice
pksEnd = pksStart + Optodur*SampRate;
newt = (0:SampInt: size(diffLEDout,1)* SampInt-SampInt)';
figure
hold on
plot(LFPt, diffLEDout)
plot(LFPt(pksStart), pksStart./pksStart, 'bo')
plot(LFPt(pksEnd), pksEnd./pksEnd, 'ro')
hold off

%% Peak
StimThetaAnglePK = phaseTheta(pksStart); % Zero-phase Filter
degAnglePK = rad2deg(StimThetaAnglePK);
StimThetaAnglePK2 = phaseTheta2(pksStart); % GUI Bandpass filter
degAnglePK2 = rad2deg(StimThetaAnglePK2);

%% Trough
StimThetaAngleTR = phaseTheta(pksStart); % Zero-phase Filter
degAngleTR = rad2deg(StimThetaAngleTR);
StimThetaAngleTR2 = phaseTheta2(pksStart); % GUI Bandpass filter
degAngleTR2 = rad2deg(StimThetaAngleTR2);
% StimLowGammaAngle = phaseLgamma(pksEnd);
% StimLowGammaAngle2 = phaseLgamma(pksStart);
% StimHiGammaAngle = phaseFgamma(pksEnd);
% StimHiGammaAngle2 = phaseFgamma(pksStart);
%% Peak
close all
figure
subplot(221)
polarhistogram(StimThetaAnglePK,20,'Normalization','probability',...
    'FaceColor',[98/255 159/255 67/255],'FaceAlpha',1,'EdgeColor','none');
title('Zero-Phase Theta Filt Phase Angle Dist.')
subplot(222)
histogram(degAnglePK,20,'Normalization','probability',...
    'FaceColor',[98/255 159/255 67/255],'FaceAlpha',1,'EdgeColor','none');
hold on
x = rad2deg(-pi:0.01:pi);
amplred = 0.3;
plot(x,(amplred*cos(2*pi*.00135*x)).*(amplred*cos(2*pi*.00135*x))+.05,'k', 'LineWidth',2),
hold off
xlabel('Degrees')
ylabel('Probability')
axis([-180 180 0 0.2])
degreetick 'x'
subplot(223)
polarhistogram(StimThetaAnglePK2,20,'Normalization','probability',...
    'FaceColor',[98/255 159/255 67/255],'FaceAlpha',1,'EdgeColor','none');
title('GUI 6-10Hz Theta Filt Phase Angle Dist.')
subplot(224)
histogram(degAnglePK2,20,'Normalization','probability',...
    'FaceColor',[98/255 159/255 67/255],'FaceAlpha',1,'EdgeColor','none');
hold on
x = rad2deg(-pi:0.01:pi);
amplred = 0.35;
plot(x,(amplred*cos(2*pi*.00135*x)).*(amplred*cos(2*pi*.00135*x))+.05,'k', 'LineWidth',2),
hold off
sgtitle('4-12Hz Zero-Phase Filt vs 6-10Hz GUI filt  Tet1 vGATCreCh2-6')
xlabel('Degrees')
ylabel('Probability')
axis([-180 180 0 0.2])
degreetick 'x'

%% Trough
close all
figure
subplot(221)
polarhistogram(StimThetaAngleTR,20,'Normalization','probability',...
    'FaceColor',[98/255 159/255 67/255],'FaceAlpha',1,'EdgeColor','none');
title('Zero-Phase Theta Filt Phase Angle Dist.')
subplot(222)
histogram(degAngleTR,20,'Normalization','probability',...
    'FaceColor',[98/255 159/255 67/255],'FaceAlpha',1,'EdgeColor','none');
hold on
x = rad2deg(-pi:0.01:pi);
amplred = 0.3;
plot(x,(amplred*cos(2*pi*.00135*x)).*(amplred*cos(2*pi*.00135*x))+.05,'k', 'LineWidth',2),
hold off
xlabel('Degrees')
ylabel('Probability')
axis([-180 180 0 0.2])
degreetick 'x'
subplot(223)
polarhistogram(StimThetaAngleTR2,20,'Normalization','probability',...
    'FaceColor',[98/255 159/255 67/255],'FaceAlpha',1,'EdgeColor','none');
title('GUI 6-10Hz Theta Filt Phase Angle Dist.')
subplot(224)
histogram(degAngleTR2,20,'Normalization','probability',...
    'FaceColor',[98/255 159/255 67/255],'FaceAlpha',1,'EdgeColor','none');
hold on
x = rad2deg(-pi:0.01:pi);
amplred = 0.35;
plot(x,(amplred*cos(2*pi*.00135*x)).*(amplred*cos(2*pi*.00135*x))+.05,'k', 'LineWidth',2),
hold off
sgtitle('4-12Hz Zero-Phase Filt vs 6-10Hz GUI filt  Tet1 vGATCreCh2-6')
xlabel('Degrees')
ylabel('Probability')
axis([-180 180 0 0.2])
degreetick 'x'
%% 
save('peakDataTet.mat','StimThetaAngle','pksEnd','LFPt','inputLFP','thetaLFP','events')
clear
%% 
save('troughDataTet.mat','StimThetaAngle','pksEnd','LFPt','inputLFP','thetaLFP','events')
clear



