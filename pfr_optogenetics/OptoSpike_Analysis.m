%{ 
    Code for analyzing Spiking Activity after optogenetic Closed-loop Stimulation
    We need the following data set:
        1- MUA multiunit activity or
        2- SUA single unit activity
        3- LED Driver voltage output
%}
clear
%% First Stimulus detection and timestamps
[LEDout, LEDtime] = load_open_ephys_data_faster('102_ADC6.continuous'); % Load voltage LED output
Optodur = 0.01; % Duration of Stimulation
SampRate = 30000; % 30KHz is the raw sampling rate
SampInt = 1/SampRate; %Sampling interval
diffLEDout = zeros(size(LEDout,1),1); % Adding a sample to the derivative vector 
diffLEDout(2:end, 1) = diff(LEDout); % Adding a sample to the derivative vector
%%
plot(LEDtime, diffLEDout) % Verifying the LED output derivative
%%
pksStart = find(diffLEDout > 0.4); % First stimulus sample indice
%{
pksEnd can introduce variability because the last peak detection could be 
one sample larger or smaller than the real stimulus duration. How good
is the detection of the last peak will depends on the quality of the signal. 
Sometimes the quality of the recorded analog signal it's not so good. 
It's better to rely only in the first stimulus sample as a 
detection parameter of the Optogenetic stimulation and then extrapolate
for the stimulus duration that is already known. Usually 10ms or 250ms.
pksEnd = find(diffLEDout < -0.45 & diffLEDout > -2);
StimIdx = [pksStart, pksEnd]; 
StimDur = LEDtime(pksEnd) - LEDtime(pksStart);
%}
figure % Corroborate that the detection is correct
hold on
plot(LEDtime, diffLEDout)
plot(LEDtime(pksStart), pksStart./pksStart, 'ro')
hold off
% axis([0 100 -2 2])
%% Spiking Activity During Optogenetic Stimulation
bins = LEDtime;
SpikeTime = tet7_units.unit3.ts;
SpikeTrain = hist(SpikeTime,bins);
SpikeIdx = find(SpikeTrain);
trialVect = 1: 1: size(pksStart,1);
SpikeWinDur = 0.1; % 50ms window Spike analysis 
winTime = (-0.05 : SampInt : 0.05)*1000;
StimTrials = zeros(size(pksStart,1),SampRate*SpikeWinDur+1); 
for  kk = 1 : size(pksStart,1)
     StimTrials(kk,:) = SpikeTrain(pksStart(kk)- SampRate*0.05 : ...
         pksStart(kk)+ SampRate*0.05);
end
close all

figure 
subplot(211)
plot(winTime, StimTrials(1,:),'.','Color',[105/255 105/255 105/255])
hold on
for kk = 2 : size(pksStart,1)
   plot(winTime, StimTrials(kk,:)*kk,'.','Color',[105/255 105/255 105/255]) 
end
line([0,10],[size(pksStart,1)+5,size(pksStart,1)+5],...
    'Color', 'b', 'LineWidth', 2)
hold off
xlabel('Time(ms)')
ylabel('Stimulation Trial #')
title('Closed-Loop Simulation')
ymaxval = size(pksStart,1) + 10;
axis([-50 50 1 ymaxval])
subplot(212)
[spiketimes, spiketrials] = find(StimTrials');
SampFactor = 30; % 30 * Samppling interval would be 1ms
newSampling = SampInt * SampFactor; 
PSTH2 = hist(spiketimes, 1:SampFactor: size(StimTrials,2))...
    /size(pksStart,1)/SampInt*SampFactor;
bar((-0.05: newSampling : 0.05), PSTH2,...
    'FaceColor',[105/255 105/255 105/255],'EdgeColor','none','BarWidth',1)
xlabel('Time(sec)')
ylabel('Spike Rate(Spikes/sec)')

%% Analysinf LFP Trials during Optogenetic Stimulation
[inputLFP,LFPt,thetaLFP] = ThetaLFP(1,1,1,1);
%%
LFPTrials = zeros(size(pksStart,1),SampRate*SpikeWinDur+1); 
for  kk = 1 : size(pksStart,1)
     LFPTrials(kk,:) = inputLFP(pksStart(kk)- SampRate*0.05 : ...
         pksStart(kk)+ SampRate*0.05);
end
avgLFP = sum(LFPTrials)./size(pksStart,1);

thetaTrials = zeros(size(pksStart,1),SampRate*SpikeWinDur+1); 
for  kk = 1 : size(pksStart,1)
     thetaTrials(kk,:) = thetaLFP(pksStart(kk)- SampRate*0.05 : ...
         pksStart(kk)+ SampRate*0.05);
end
avgTheta = sum(thetaTrials)./size(pksStart,1);
ntime = (-0.05:SampInt:0.05)*1000;
%%
close all
figure
hold on
for kk = 1: size(pksStart,1)
   plot(ntime, LFPTrials(kk,:), 'Color',[192/255 192/255 192/255])
end
plot(ntime, avgLFP,'k')
hold off
ylabel('LFP voltage(uV)')
xlabel('Time(ms)')
figure
hold on
for kk = 1: size(pksStart,1)
   plot(ntime, thetaTrials(kk,:), 'Color',[192/255 192/255 192/255])
end
plot(ntime, avgTheta,'k')
hold off
ylabel('LFP voltage(uV)')
xlabel('Time(ms)')
%% Interpolation
VideoSR = 1/15;
postime = (0 :VideoSR : size(position,1) * VideoSR - VideoSR);
pos = interp1(postime, position, LFPt);
pos = fillmissing(pos,'linear');
%%
close all
plot(pos(:,1), pos(:,2),'Color',[128/255 128/255 128/255]);
hold on
p1 = plot(pos(pksStart,1), pos(pksStart,2), 'r.','MarkerSize',12 );
legend(p1, 'Retrieval Stimulation','Location',...
       'north')
legend('boxoff')
%%
close all
plot(pos(:,1), pos(:,2),'Color',[128/255 128/255 128/255]);
hold on
p1 = plot(pos(SpikeIdx,1), pos(SpikeIdx,2), 'r.','MarkerSize',12 );
legend(p1, 'Retrieval Stimulation','Location',...
       'north')
legend('boxoff')