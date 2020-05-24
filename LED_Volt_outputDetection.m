%{ 
    Code for analyzing optogenetic Closed-loop Theta-phase Stimulation
    We need the following data set:
        1- MUA
        2- SUA
        3- Position
        4- LFP(Raw and Filtered in theta)
        5- LED Driver voltage output
        6- Phase of Filtered LFP
%}

%% First Stimulus detection and timestamps
LEDout = load_open_ephys_data_faster('102_ADC6.continuous'); % Load voltage LED output
LED = downsample(LEDout,15); % Reduce sampling rate
[inputLFP,LFPt,thetaLFP] = ThetaLFP(5, 2,1,15); % This function downsample and filter the LFP signal
phaseTheta = angle(hilbert(thetaLFP));
Optodur = 0.01; % Duration of Stimulation
SampRate = 2000; 
SampInt = 1/SampRate;
diffLEDout = zeros(size(LED,1),1); %Adding a sample to the derivative vector 
diffLEDout(2:end, 1) = diff(LED); 
%%
[OEevents, OEtime, info] = load_open_ephys_data_faster('all_channels.events'); % load voltage LED output
events = find(OEevents < 1);
% events(2:2:end,:) = [];
plot(OEtime, OEevents*.0,'ro')
events1 = sub2ind(OEtime(events)*SampRate,1);
hold on
plot(OEtime(events), OEevents(events),'b*')
plot(OEtime, OEevents*.1,'k.')
axis([0 inf -20 20]);
hold off
%%
plot(LFPt, LED-5)
hold on
plot(LFPt, diffLEDout./5)
plot(LFPt, phaseTheta/25, 'm', 'Linewidth',2 )
plot(OEtime(events), OEevents(events),'b.','MarkerSize',12)
% plot(OEtime, OEevents*.0,'b.','MarkerSize',12)
plot(LFPt, inputLFP/5000, 'k')
plot(LFPt, thetaLFP/2000, 'r', 'Linewidth',2 )
hold off
axis([0 inf -1 0.4]);
%%
% save('LEDout.mat','LEDout','LED','diffLEDout')
%% pks = find(diffLEDout >=2|diffLEDout <=-2);
[pksStart,pksCol,pksVal] = find(diffLEDout >=0.5); % first stimulus sample indice
pksEnd = pksStart + Optodur * SampRate; % Last Stimulus sample indice
StimDur = [pksStart, pksEnd]; 
StimDurTime = StimDur * SampInt;
diffDur = diff(StimDurTime,1,2);

%% Interpolation
VideoSR = 1/15;
postime = (0 :VideoSR : size(pos,1) * VideoSR - VideoSR);
pos2 = interp1(postime, pos, LFPt);
pos2 = fillmissing(pos2,'linear');
%%
plot(pos2(:,1), pos2(:,2));
hold on
plot(pos2(pksEnd,1), pos2(pksEnd,2), 'r.','MarkerSize',12 );
% linVel = interp1(new_t, new_vel, LFPt);

%% Plotting Spikes and Linear Position
close all
spikeindex = pksStart;
figure
hold on
plot(linPos, LFPt, 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
plot(linPos(spikeindex),LFPt(spikeindex),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Retrieval Stimulation')
xlabel('Linear Position[m]')
ylabel('Time[sec]')
axis tight
%% Theta angle
close all
thetaAng = angle(hilbert(thetaLFP7));
H1 = plot(LFPt, thetaLFP7,'Color', [0.5,0.5,0.5])
hold on
H2 = plot(LFPt(spikeindex), thetaLFP7(spikeindex),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 5)
H3 = plot(OEtime(events), thetaLFP7(events),'o',...
    'MarkerEdgeColor','b','MarkerFaceColor','k','MarkerSize', 5)
hold off
title('Retrieval Stimulation')
xlabel('Time[sec]')
ylabel('Filtered Theta[uV]')
legend([H2 H3],{'LED Stimulation','OE Phase Detector'})
%%
close all
polarhistogram(thetaAng(spikeindex),20,'FaceColor','red','FaceAlpha',.5);
title('Theta Phase Stimulation Distribution')
histogram(thetaAng(spikeindex),20)
%%
save('FreeD19_PeakCL_analysis.mat')