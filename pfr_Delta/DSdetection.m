Fs = 1000;
sint = 1/Fs;
CA1LFP = [];
CA1LFPripple = [];
CA3LFP = [];
CA3LFPripple = [];
DGLFP = [];
DGLFPfgamma = [];
CTXLFP = [];
CTXLFPdeltaspin = [];
%% Plotting LFP and MUA
time1= 1605;
time2= 1609;
xt1= time1/sint;
xt2= time2/sint;
figure
% Cortex Tetrode 13
subplot(6,1,1)
plot(time(xt1:xt2),CTXLFP(xt1:xt2),'k') 
axis tight
legend('Ctx','Location', 'best')
legend boxoff
box off
set(gca,'xtickLabel',[])
ylabel('uV')

subplot(6,1,2)
plot(time(xt1:xt2),CTXLFPdeltaspin(xt1:xt2),'k','LineWidth',2) 
axis tight
legend('Delta 1-4Hz','Location', 'best')
legend boxoff
box off
set(gca,'xtickLabel',[])
ylabel('uV')

subplot(6,1,3)
plot(time(xt1:xt2),CA1LFP(xt1:xt2),'k') 
axis([-inf inf 0 200])
legend('100-450Hz')
legend boxoff
box off
set(gca,'xtickLabel',[])
ylabel('Ampl (uV)')

subplot(6,1,4)
plot(time(xt1:xt2),CA1LFPripple(xt1:xt2),'k') 
axis([-inf inf 0 200])
legend('100-450Hz')
legend boxoff
box off
set(gca,'xtickLabel',[])
ylabel('Ampl (uV)')

if exist('spkDensity','var') == 0
    warning('Spike density not in Workspace');
    return
end

subplot(6,1,5)
hold on
plot(spkClust(ctxTetIdx(1)).spkTime, ... 
    ones(length(spkClust(ctxTetIdx(1)).spkTime),1),'r.','MarkerSize',10)
for kk = 2: size(ctxTetIdx,1)
plot(spkClust(ctxTetIdx(kk)).spkTime, ...
    ones(length(spkClust(ctxTetIdx(kk)).spkTime),1) * kk,... 
    'r.','MarkerSize',10)
end
hold off
xlim([time1 time2])
ylim([0 inf])
set(gca,'xtickLabel',[])
xlabel('Time[sec]')
ylabel('MUA #')

subplot(6,1,6)
plot(time(xt1:xt2),spkDensity(xt1:xt2),'r','LineWidth',2) 
xlabel('Time[sec]')
ylabel('spks/s')
axis tight

%%
time1= 1600;
time2= 1609;
xt1= time1/sint;
xt2= time2/sint;
h= histogram(spkDensity(xt1:xt2),200,'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',1);
h.BinEdges = 0:30;
h.BinWidth = 0.01;
h.EdgeColor = 'k';
xlabel('MUA spike Density Hz')
ylabel('Probability Density')
legend('MUA spike Density during Behavior')
legend('boxoff')  
box off

%%
thr1 = 22.6; % 


%%  FIRST CRITERION: AT LEAST 3*STD OF FILTERED, SQUARED LFP
% splog = spkDensity' < thr1;

splog = highFreqAmpLFP <= thr1;
pos = find(splog == 1); 


%% SECOND CRITERION: at least a difference of 15ms between the end of
...a ripple and the beginning of the next to be considered separate events

    difpos = find(diff(pos) >= (0.05/sint)); % events with inte4rvals > 20ms
    z = zeros(length(difpos)+1,1); % Creating a zero row vector
    z(1:end-1) = difpos;
    z(end) = length(pos); % End position of all potential ripples
    % Calculate the start
    st = difpos+1; % Start positions of all possible ripple events
    zst = zeros(length(difpos)+1,1);
    zst(2:end) = st; %
    zst(1)=1; % zst start position of all
    stepDelta = [pos(zst), pos(z)]; % Start and End position(indices) 
        ... of Ripples Candidates


%% THIRD CRITERION: Eliminate 'Delta' candidates with < 1 sampling interval

duration = stepDelta(:,2)- stepDelta(:,1);
stepDelta1 = stepDelta((~(duration<=1)),:);

%% FOURTH CRITERION:  at least 20ms long

duration2 = (stepDelta1(:,2) - stepDelta1(:,1))*sint; %duration in sec.
deltaIdx = stepDelta1((find(duration2<=1 & ...
    duration2 >= 0.02)),:); %Events < 150ms and > 15ms (Indices)


%% Second Step
below = deltaIdx(:,1); % - (0.5/sint);
abov = deltaIdx(:,2); % + (0.5/sint);
below(below < 0) = 1; % Substitute start negative value with 1
abov(abov > length(LFP)) = length(LFP); % Making delta candidates finishing at the end of the LFP


%% Visual Inspection of LFP and Candidates
close all
time1= 1600;
time2= 1609;
xt1= time1/sint;
xt2= time2/sint;
figure
plot(time(xt1:xt2), LFP(xt1:xt2),'k')
hold on
for kk = 1: length(below)
    plot(time(below(kk):abov(kk)), LFP(below(kk):abov(kk)), 'r', 'LineWidth',1) 
end
hold off
axis([time1 time2 -500 600])


%%
deltaCanddur = abov - below;
h= histogram((deltaCanddur .* sint),'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',2);
h.BinEdges = 0:1;
h.BinWidth = 0.005;
h.EdgeColor = 'r';
xlabel('Down-state duration (s)')
ylabel('Probability Density')
legend('Down-State duriation')
legend('boxoff') 
xlim([0 0.5])
box off

%% Find Short < 100ms duration Down states

spindleLFP = ch42Tet13LFPspindle;
deltaSpinLFP = ch42Tet13LFPdeltaSpindle;
deltaSpindleAmpLFP = abs(hilbert(deltaSpinLFP));


%% Visual Inspection of Down-States candidates
close all
eventExt = 0.05 / sint;
for kk = 1:20 %: size(,1)
   figure
   subplot(211)
   plot(time(below(kk) - eventExt : abov(kk)+ eventExt), ...
       LFP(below(kk) - eventExt : abov(kk)+ eventExt))
   axis([-inf inf, -inf inf])
   
   subplot(212)
   plot(time(below(kk) - eventExt : abov(kk)+ eventExt),...
       deltaSpinLFP(below(kk) - eventExt : abov(kk)+ eventExt))
   axis([-inf inf, -inf inf])
   
end


%% correlation coefficient between the LFP candidates and the 1-16Hz filtered
%  Eliminate evemts that are lager than the LFP
if abov(end)+eventExt > length(LFP) 
    abov(end) = [];
    below(end) = [];
else
end
RR = zeros(size(below,1),1);
for kk = 1: size(below,1)-1
   R = corrcoef((LFP(below((kk)) - eventExt :... 
       abov((kk))+ eventExt)),... 
       (deltaSpinLFP(below((kk)) - eventExt :... 
       abov((kk))+ eventExt)));
   RR(kk,1) = R(2,1);
end
RRint = round(RR,4);


%% Distribution of the correlation
close all
h= histogram(RRint,'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',2);
h.BinEdges = -1:1;
h.BinWidth = 0.01;
h.EdgeColor = [213/300 29/300 91/300];
xlabel(' correlation coefficient')
ylabel('Probability Density')
legend('Short duration Down-States Candidates','Location','best')
legend('boxoff') 
xlim([0 1])
box off
%%
abovCorr = abov;
belowCorr = below;
noiseCandIdx = find(RR <= 0.6);
abovCorr((noiseCandIdx)) = [];
belowCorr((noiseCandIdx)) = [];


%% Down-state duration distribution visualization after correlation filtering

deltaCanddurRR = abovCorr - belowCorr;
h= histogram((deltaCanddurRR .* sint),'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',2);
h.BinEdges = 0:1;
h.BinWidth = 0.01;
h.EdgeColor = 'r';
xlabel('Down-state duration (s)')
ylabel('Probability Density')
legend('Down-State duriation')
legend('boxoff') 
xlim([0 0.6])
box off

%% Visual Inspection
close all
time1= 1535;
time2= 1540;
xt1= time1/sint;
xt2= time2/sint;
figure
plot(time, LFP,'k')
hold on
for kk = 1: length(belowCorr)
    plot(time(belowCorr(kk):abovCorr(kk)), LFP(belowCorr(kk):abovCorr(kk)), 'r', 'LineWidth',1) 
end
hold off
axis([time1 time2 -500 600])


%% Visual Inspection of Down-States candidates
close all
eventExt = 0.05 / sint;
for kk = 100:120 %: size(,1)
   figure
   subplot(211)
   plot(time(belowCorr((kk)) - eventExt :... 
       abovCorr((kk))+ eventExt),...
       LFP(belowCorr((kk)) - eventExt :... 
       abovCorr((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   ylabel('uV')
   legend('Raw LFP')
   legend('boxoff') 
   
   subplot(212)
   plot(time(belowCorr((kk)) - eventExt :... 
       abovCorr((kk))+ eventExt),...
       deltaSpinLFP(belowCorr((kk)) - eventExt :... 
       abovCorr((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   xlabel('Time (s)')
   ylabel('uV')
   legend('1-16Hz')
   legend('boxoff') 
end




%% Local Minima and local maxima from delta/spin lfp
maxt = 1:length(abovCorr);
maxMax = cell(size(maxt,1),1);
minMin = cell(size(maxt,1),1);
for kk = 1:size(maxt,2)
    ll = maxt(kk);
    candDelt = deltaSpinLFP(belowCorr((ll)) - eventExt : ... 
       abovCorr((ll))+ eventExt); % Local Min/Max Delta/Spin LFP
    maxDeltaAmpIdx = find(islocalmax(candDelt)); % find local max
    minDeltaAmpIdx = find(islocalmin(candDelt)); % find local min
   maxMax{kk,1} = max(candDelt(maxDeltaAmpIdx));
   minMin{kk,1} = min(candDelt(minDeltaAmpIdx));
end
%% Local Min/Max Amplitude of Delta/Spin LFP
deltaAmpRR = cellfun(@minus,maxMax,minMin,'Un',0);
empties = cellfun('isempty',deltaAmpRR);
deltaAmpRR(empties) = {NaN};
deltaAmpRR = cell2mat(deltaAmpRR); 
nanDeltaAmpRRIdx = find(isnan(deltaAmpRR));
abovCorr((nanDeltaAmpRRIdx)) = []; %Eliminate empty
belowCorr((nanDeltaAmpRRIdx)) = [];
%% Re-analaysis
maxt = 1:length(abovCorr);
maxMax = cell(size(maxt,1),1);
minMin = cell(size(maxt,1),1);
for kk = 1:size(maxt,2)
    ll = maxt(kk);
    candDelt = deltaSpinLFP(belowCorr((ll)) - eventExt : ... 
       abovCorr((ll))+ eventExt); % Local Min/Max Delta/Spin LFP
    maxDeltaAmpIdx = find(islocalmax(candDelt)); % find local max
    minDeltaAmpIdx = find(islocalmin(candDelt)); % find local min
   maxMax{kk,1} = max(candDelt(maxDeltaAmpIdx));
   minMin{kk,1} = min(candDelt(minDeltaAmpIdx));
end
deltaAmpRR = cellfun(@minus,maxMax,minMin,'Un',0);
empties = cellfun('isempty',deltaAmpRR);
deltaAmpRR(empties) = {NaN};
deltaAmpRR = cell2mat(deltaAmpRR); 
%% Peak Ampalitude
close all
h= histogram(deltaAmpRR(:,1),'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',2);
% h.BinEdges = -500:500;
% h.BinWidth = 1;
h.EdgeColor = [213/300 29/300 91/300];
xlabel(' 1-16Hz Peak Ampl. (uV)')
ylabel('Probability Density')
legend('Down-States Candidates','Location','best')
legend('boxoff') 
xlim([-inf inf])
box off

%%

abovCorrAmp = abovCorr;
belowCorrAmp = belowCorr;
noiseCandIdx = find(deltaAmpRR <= 50);
abovCorrAmp((noiseCandIdx)) = [];
belowCorrAmp((noiseCandIdx)) = [];


%% Visual Inspection
close all
time1= 1535;
time2= 1540;
xt1= time1/sint;
xt2= time2/sint;
figure
% subplot(211)
% plot(time, position,'k', 'LineWidth', 2)
% xlabel('Time [sec]')
% ylabel('Position [m]')
% axis([time1 time2 0 inf])
% subplot(212)
plot(time, LFP,'k')
hold on
for kk = 1: length(belowCorrAmp)
    plot(time(belowCorrAmp(kk):abovCorrAmp(kk)), ...
        LFP(belowCorrAmp(kk):abovCorrAmp(kk)),...
        'r', 'LineWidth',1) 
end
% plot(time, deltaLFP,'k','LineWidth',1)
% plot(time, (smsqfiltLFP./80)-400,'b', 'LineWidth',1)
hold off
axis([time1 time2 -500 600])

%% Visual Inspection of Down-States candidates
close all
eventExt = 0.05 / sint;
for kk = 1550:1560   %: size(,1)
   figure
   subplot(211)
   plot(time(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),...
       LFP(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   ylabel('uV')
   legend('Raw LFP')
   legend('boxoff') 
   
   subplot(212)
   plot(time(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),...
       deltaSpinLFP(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   xlabel('Time (s)')
   ylabel('Slope')
   legend('1-16Hz')
   legend('boxoff') 
end
%% Down-state duration distribution visualization after correlation filtering
diffdeltaCanddurRR = abovCorrAmp - belowCorrAmp;
h= histogram((diffdeltaCanddurRR .* sint),'Normalization','probability',... 
    'DisplayStyle','bar', 'LineWidth',2);
% h.BinEdges = 0:1;
h.BinWidth = 0.005;
h.EdgeColor = 'r';
xlabel('Down-state duration (s)')
ylabel('Probability Density')
legend('Down-State duriation')
legend('boxoff') 
xlim([0 0.6])
box off

%% Max Amplitude of diff delta/spin lfp 
diffdeltaSpinLFP = zeros(length(time),1);
diffdeltaSpinLFP(2:end) = diff(deltaSpinLFP);

deltaAmpDiff = zeros(size(abovCorrAmp,1),1);
for kk = 1 : size(abovCorrAmp,1)-1
   deltaAmpDiff(kk,1) = max(diffdeltaSpinLFP(belowCorrAmp(kk) - eventExt:... 
       abovCorrAmp(kk) + eventExt));
end

%% Differential Amplitude Threshold
abovCorrAmpDiff = abovCorrAmp;
belowCorrAmpDiff = belowCorrAmp;
noiseCandIdx = find(deltaAmpDiff <= 4);
abovCorrAmpDiff((noiseCandIdx)) = [];
belowCorrAmpDiff((noiseCandIdx)) = [];

%% Visual Inspection
close all
time1= 10000;
time2= 10010;
xt1= time1/sint;
xt2= time2/sint;
figure
plot(time, LFP,'k')
hold on
for kk = 1: length(belowCorrAmpDiff)
    plot(time(belowCorrAmpDiff(kk):abovCorrAmpDiff(kk)), ...
        LFP(belowCorrAmpDiff(kk):abovCorrAmpDiff(kk)),...
        'r', 'LineWidth',1) 
end
hold off
axis([time1 time2 -500 600])

%% Visual Inspection of Down-States candidates
close all
eventExt = 0.05 / sint;
for kk = 1550:1609    %: size(,1)
   figure
   subplot(211)
   plot(time(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),...
       LFP(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   ylabel('uV')
   legend('Raw LFP')
   legend('boxoff') 
   
   subplot(212)
   plot(time(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),...
       deltaSpinLFP(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   xlabel('Time (s)')
   ylabel('Slope')
   legend('1-16Hz')
   legend('boxoff') 
end
%% Down-state duration distribution visualization after correlation filtering

diffdeltaCanddurRR = abovCorrAmpDiff - belowCorrAmpDiff;
h= histogram((diffdeltaCanddurRR .* sint),'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',2);
% h.BinEdges = 0:1;
h.BinWidth = 0.007;
h.EdgeColor = 'r';
xlabel('Down-state duration (s)')
ylabel('Probability Density')
legend('Down-State duriation')
legend('boxoff') 
xlim([0 0.6])
box off
%% Hilbert-Huang transform
% 
% close all
% % hhTran = zeros(size(,1)-1,1);
% for kk = 19 %: size(,1)-1
%     emd(LFP(below((kk)) - eventExt :... 
%        abov((kk))+ eventExt),'Display',0)
%     [imf,residual,info] = emd(LFP(below((kk)) - eventExt :... 
%        abov((kk)) + eventExt),'Display',0,'MaxNumIMF',5);
%    infoA.energy{kk} = info.MeanEnvelopeEnergy;
%    infoA.maxNumIMF(kk) = max(info.NumIMF);
% %    hht(imf, Fs,'FrequencyLimits',[0 100]); 
% end


%% Plotting Down-states in the Maze
% Using Velocity to Detect delta wavess during inmobility
position = interp1(positiondata.time, positiondata.linearpos, time);
velocity = interp1(positiondata.time, positiondata.linVel, time);
pos2d = interp1(positiondata.time,positiondata.poscm, time);
%%
figure
% subplot(122)
hold on
h_y = plot(time,position, ...
    'Color', [0.5,0.5,0.5], 'LineWidth', 1);
h_Z = plot(time(belowCorrAmpDiff),position(belowCorrAmpDiff),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3);
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
title('"Down-State" Activity')
ylabel('Linear Position[cm]')
xlabel('Time[sec]')
axis tight
legend([h_Z, h_y],'"Down-States"', 'Linear Pos')
legend boxoff


% subplot(121)
figure
plot(pos2d(:,1),pos2d(:,2), 'Color', [0.5,0.5,0.5])
hold on
plot(pos2d(belowCorrAmpDiff,1),pos2d(belowCorrAmpDiff,2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
xlabel('x coordinate [cm]')
ylabel('y coordinate [cm]')
axis tight

%% Inter Event Histograms
deltaCanddur = abovCorrAmpDiff .*sint;
h= histogram(diff(deltaCanddur),'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',2);
% h.BinEdges = 0:1;
h.BinWidth = 0.01;
h.EdgeColor = 'r';
xlabel('Inter Down-state Interval (s)')
ylabel('Probability Density')
legend('Inter Down-state Interval')
legend('boxoff') 
xlim([0 2])
box off

%% Sleep Inter Down States Intervals
close all
preDownEvents = deltaCanddur(deltaCanddur>1000 & deltaCanddur<1609);
h= histogram(diff(preDownEvents),'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',2);
% h.BinEdges = 0:1;
h.BinWidth = 0.015;
h.EdgeColor = 'r';
xlabel('Inter Down-state Interval (s)')
ylabel('Probability Density')
legend('Pre "Sleeping Period"')
legend('boxoff') 
xlim([0 2])
ylim([0 .12])
box off
%%
close all
runDownEvents = deltaCanddur(deltaCanddur>1700 & deltaCanddur<7700);
h= histogram(diff(runDownEvents),'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',2);
% h.BinEdges = 0:1;
h.BinWidth = 0.015;
h.EdgeColor = 'k';
xlabel('Inter Down-state Interval (s)')
ylabel('Probability Density')
legend('Behavior')
legend('boxoff') 
xlim([0 2])
ylim([0 .12])
box off

%%
close all
postDownEvents = deltaCanddur(deltaCanddur>8000 & deltaCanddur<10000);
h= histogram(diff(postDownEvents),'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',2);
% h.BinEdges = 0:1;
h.BinWidth = 0.015;
h.EdgeColor = 'r';
xlabel('Inter Down-state Interval (s)')
ylabel('Probability Density')
legend('Post "Sleeping Period"')
legend('boxoff') 
xlim([0 2])
ylim([0 .12])
box off

%% cumulative histogram
close all
hold on
h1 = histogram(diff(preDownEvents),'Normalization','cdf',...
    'DisplayStyle','stairs','LineWidth',2,'BinWidth',0.001);
h2 = histogram(diff(runDownEvents),'Normalization','cdf',...
    'DisplayStyle','stairs','LineWidth',2,'BinWidth',0.001);
h3 = histogram(diff(postDownEvents),'Normalization','cdf',...
    'DisplayStyle','stairs','LineWidth',2,'BinWidth',0.001);
hold off
xlabel('Inter Down-state Interval (s)')
ylabel('Norm. Cumulative Probability')
legend([h1 h2 h3],'Pre "Sleeping Period"', 'Behavior', 'Post "Sleeping Period"')
legend('boxoff') 
xlim([0 .8])
box off

%% Local Maxima and minima
deltaCanddur = abovCorrAmpDiff .*sint;
interestEvents = find(deltaCanddur>8000 & deltaCanddur<10000);
eventExt = 0.05 / sint;
maxt = interestEvents;
maxMax = cell(size(maxt,1),1);
minMin = cell(size(maxt,1),1);
for kk = 1:size(maxt,1)
    ll = maxt(kk);
    candDelt = deltaSpinLFP(belowCorrAmpDiff((ll)) - eventExt :... 
       abovCorrAmpDiff((ll))+ eventExt);
    maxDeltaAmpIdx = find(islocalmax(candDelt)); % find local max
    minDeltaAmpIdx = find(islocalmin(candDelt)); % find local min
   maxMax{kk,1} = max(candDelt(maxDeltaAmpIdx));
   minMin{kk,1} = min(candDelt(minDeltaAmpIdx));
end
 %%  
deltaAmp = cellfun(@minus,maxMax,minMin,'Un',0);
close all
hold on
h1 = histogram(cell2mat(maxMax),'Normalization','probability',...
    'DisplayStyle','bar','LineWidth',1,'BinWidth',10,...
    'EdgeColor','c','FaceColor','c');
h2 = histogram(cell2mat(minMin),'Normalization','probability',...
    'DisplayStyle','bar','LineWidth',1,'BinWidth',10,...
    'EdgeColor','k','FaceColor','k');
h3 = histogram(cell2mat(deltaAmp),'Normalization','probability',...
    'DisplayStyle','bar','LineWidth',1,'BinWidth',10,...
    'EdgeColor','r','FaceColor','r');
hold off
xlabel('Delta/Spindle Amplitude (uV)')
ylabel('Probability')
legend([h1 h2 h3],'Local Maxima', 'Local Minima', 'Diff','Location','best')
legend('boxoff') 
% xlim([0 .8])
box off
%%
close all  
h1 = histogram(cell2mat(deltaAmp),'Normalization','probability',...
    'DisplayStyle','bar','LineWidth',1,'BinWidth',10,...
    'EdgeColor','r','FaceColor','r'); 
xlabel('Delta/Spindle Amplitude (uV)')
ylabel('Probability')
legend('Down-State Amplitude','Location','best') 
legend('boxoff') 
xlim([0 600])
box off

%%
deltaCanddur = abovCorrAmpDiff .*sint;
eventExt = 0.05 / sint;
maxt = 1:length(deltaCanddur);
maxMax = cell(size(maxt,1),1);
minMin = cell(size(maxt,1),1);
for kk = 1:size(maxt,2)-1
    ll = maxt(kk);
    candDelt = deltaSpinLFP(belowCorrAmpDiff((ll)) - eventExt :... 
       abovCorrAmpDiff((ll))+ eventExt);
    maxDeltaAmpIdx = find(islocalmax(candDelt)); % find local max
    minDeltaAmpIdx = find(islocalmin(candDelt)); % find local min
   maxMax{kk,1} = max(candDelt(maxDeltaAmpIdx));
   minMin{kk,1} = min(candDelt(minDeltaAmpIdx));
end
deltaAmp = cellfun(@minus,maxMax,minMin,'Un',0);   
empties = cellfun('isempty',deltaAmp);
deltaAmp(empties) = {NaN};
deltaAmp1 = cell2mat(deltaAmp);   

%%
close all
s = scatter(diffdeltaCanddurRR(1:length(deltaAmp)),deltaAmp1,... 
    'filled','SizeData',10,...
     'MarkerFaceColor','r',...
     'LineWidth',1);
ylabel('Down-State Amplitude (uV)')
xlabel('Down-State Duration (s)')
alpha(s,.1)
legend('Down-States','Location','best') 
legend('boxoff') 
xlim([0 600])
box off
   
%% Pre correlation
deltaCanddur = abovCorrAmpDiff .*sint;
interestEvents = find(deltaCanddur>10 & deltaCanddur<1609);
eventExt = 0.05 / sint;
maxt = interestEvents;
maxMax = cell(size(maxt,1),1);
minMin = cell(size(maxt,1),1);
for kk = 1:size(maxt,1)
    ll = maxt(kk);
    candDelt = deltaSpinLFP(belowCorrAmpDiff((ll)) - eventExt :... 
       abovCorrAmpDiff((ll))+ eventExt);
    maxDeltaAmpIdx = find(islocalmax(candDelt)); % find local max
    minDeltaAmpIdx = find(islocalmin(candDelt)); % find local min
   maxMax{kk,1} = max(candDelt(maxDeltaAmpIdx));
   minMin{kk,1} = min(candDelt(minDeltaAmpIdx));
end
deltaAmp = cellfun(@minus,maxMax,minMin,'Un',0);   
empties = cellfun('isempty',deltaAmp);
deltaAmp(empties) = {NaN};
deltaAmp1 = cell2mat(deltaAmp); 
close all
s = scatter(diffdeltaCanddurRR(interestEvents),deltaAmp1,... 
    'filled','SizeData',10,...
     'MarkerFaceColor','r',...
     'LineWidth',1);
ylabel('Down-State Amplitude (uV)')
xlabel('Down-State Duration (s)')
alpha(s,.3)
lgs = legend('Down-States','Location','best');
lgs.FontSize = 14;
lgs.TextColor ='k';
legend('boxoff') 
xlim([0 700])
ylim([0 700])
box off
figure
X = [diffdeltaCanddurRR(interestEvents),deltaAmp1];
hist3(X,'CDataMode','auto','FaceColor','interp'...
    ,'EdgeAlpha',0.5,'Edges',{0:10:700 0:10:700})
ylabel('Down-State Amplitude (uV)')
xlabel('Down-State Duration (s)')
colormap('hot') % Change color scheme 
view(2)
%% Running correlation amplitude and duration
deltaCanddur = abovCorrAmpDiff .*sint;
interestEvents = find(deltaCanddur>1700 & deltaCanddur<3700);
eventExt = 0.05 / sint;
maxt = interestEvents;
maxMax = cell(size(maxt,1),1);
minMin = cell(size(maxt,1),1);
for kk = 1:size(maxt,1)
    ll = maxt(kk);
    candDelt = deltaSpinLFP(belowCorrAmpDiff((ll)) - eventExt :... 
       abovCorrAmpDiff((ll))+ eventExt);
    maxDeltaAmpIdx = find(islocalmax(candDelt)); % find local max
    minDeltaAmpIdx = find(islocalmin(candDelt)); % find local min
   maxMax{kk,1} = max(candDelt(maxDeltaAmpIdx));
   minMin{kk,1} = min(candDelt(minDeltaAmpIdx));
end
deltaAmp = cellfun(@minus,maxMax,minMin,'Un',0);   
empties = cellfun('isempty',deltaAmp);
deltaAmp(empties) = {NaN};
deltaAmp1 = cell2mat(deltaAmp); 
close all
s = scatter(diffdeltaCanddurRR(interestEvents),deltaAmp1,... 
    'filled','SizeData',10,...
     'MarkerFaceColor','k',...
     'LineWidth',1);
ylabel('Down-State Amplitude (uV)')
xlabel('Down-State Duration (s)')
alpha(s,.3)
lgs = legend('Down-States','Location','best');
lgs.FontSize = 14;
lgs.TextColor ='k';
legend('boxoff') 
xlim([0 700])
ylim([0 700])
box off
figure
X = [diffdeltaCanddurRR(interestEvents),deltaAmp1];
hist3(X,'CDataMode','auto','FaceColor','interp'...
    ,'EdgeAlpha',0.5,'Edges',{0:10:700 0:10:700})
ylabel('Down-State Amplitude (uV)')
xlabel('Down-State Duration (s)')
colormap('hot') % Change color scheme 
view(2)
%% post sleep correlation amplitude and duration
deltaCanddur = abovCorrAmp .*sint;
interestEvents = find(deltaCanddur>8000 & deltaCanddur<10000);
eventExt = 0.05 / sint;
maxt = interestEvents;
maxMax = cell(size(maxt,1),1);
minMin = cell(size(maxt,1),1);
for kk = 1:size(maxt,1)
    ll = maxt(kk);
    candDelt = deltaSpinLFP(belowCorrAmp((ll)) - eventExt :... 
       abovCorrAmp((ll))+ eventExt);
    maxDeltaAmpIdx = find(islocalmax(candDelt)); % find local max
    minDeltaAmpIdx = find(islocalmin(candDelt)); % find local min
   maxMax{kk,1} = max(candDelt(maxDeltaAmpIdx));
   minMin{kk,1} = min(candDelt(minDeltaAmpIdx));
end
deltaAmp = cellfun(@minus,maxMax,minMin,'Un',0);   
empties = cellfun('isempty',deltaAmp);
deltaAmp(empties) = {NaN};
deltaAmp1 = cell2mat(deltaAmp); 
close all
s = scatter(diffdeltaCanddurRR(interestEvents),deltaAmp1,... 
    'filled','SizeData',10,...
     'MarkerFaceColor','r',...
     'LineWidth',1);
ylabel('Down-State Amplitude (uV)')
xlabel('Down-State Duration (s)')
alpha(s,.3)
lgs = legend('Down-States','Location','best');
lgs.FontSize = 14;
lgs.TextColor ='k';
legend('boxoff') 
xlim([0 700])
ylim([0 700])
box off
figure
X = [diffdeltaCanddurRR(interestEvents),deltaAmp1];
hist3(X,'CDataMode','auto','FaceColor','interp'...
    ,'EdgeAlpha',0.5,'Edges',{0:10:700 0:10:700})
ylabel('Down-State Amplitude (uV)')
xlabel('Down-State Duration (s)')
colormap('hot') % Change color scheme 
view(2)