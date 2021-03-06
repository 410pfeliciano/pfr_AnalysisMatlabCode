% Cortical LFP
ctxLFP = []; %Raw Cortical LFP
Fs = 1000; % samping rate
sint = 1/Fs; % sampling intervals
% bandfilter function needs to be in matlab path
ctxLFPdeltaspin = bandFilter (ctxLFP, 'deltaSpindle'); % Delta Spindle Filtered LFP
ctxDeltaLFP = bandFilter(ctxLFP, 'delta');
ctxSpindleLFP = bandFilter(ctxLFP, 'spindle');
ctxLFPmua = bandFilter (ctxLFP, 'highFreq'); % MUA envelop
% ctxLFPspindle = ch42Tet13LFPspindle;
% Smoothing the MUA envelop
ctxLFPmuaenvp = abs(hilbert(ctxLFPmua));
ctxLFPdeltaspinenvp = abs(hilbert(ctxLFPdeltaspin));
ctxLFPmuaenvpsmooth = smoothdata(ctxLFPmuaenvp,'gaussian',25); %smoothing

%% Estimating the Firing Rate from cortical MUA 
Tet = time; % just to calculate the duration of the recording
bins = (0 : 1/Fs : length(Tet) * (1/Fs) - (1/Fs))';

if exist('spkClust','var') == 0
    warning('spiking data spkClust is not in your workspace. You need to import the unit.mat file from kiloresult folder');
    return
end

TetCtx = []; % choose the tetrode of interest

CtxTetfind = cell(size(spkClust,2),1);
for kk = 1 : size(spkClust,2)
    CtxTetfind{kk} = find(spkClust(kk).tetNum == TetCtx);
end

CtxTetIdx = find(~cellfun('isempty', CtxTetfind));
Ctxspiketrain = zeros(size(CtxTetIdx,1),length(bins));
for kk = 1 : size(CtxTetIdx,1)
    Ctxspiketime = (spkClust(CtxTetIdx(kk)).spkTime)';
    Ctxspiketrain(kk,:) = hist(Ctxspiketime,bins);
end
% Spike density estimation
CtxtotalSpk = sum(Ctxspiketrain,1);
sigma = 0.02;
edges = -3*sigma: 0.001: 3*sigma;
kernel = (normpdf(edges,0,sigma)).*(1/10);
% kernel = kernel*(1/Fs); % kernel muitiplied by bin size
CtxspkDensity = conv(CtxtotalSpk, kernel);
center = ceil(length(edges)/2);
CtxspkDensity = CtxspkDensity(center:length(CtxtotalSpk)+center-1);
%% Plotting LFP and MUA envelop
time1= [];
time2= [];
xt1= time1/sint;
xt2= time2/sint;

subplot(5,1,1)
plot(time(xt1:xt2),ctxLFP(xt1:xt2),'k') 
axis tight
legend('Ctx','Location', 'best')
legend boxoff
box off
set(gca,'xtickLabel',[])
ylabel('uV')

subplot(5,1,2)
plot(time(xt1:xt2),ctxLFPmuaenvp(xt1:xt2),'k') 
axis tight
legend('MUA envp','Location', 'best')
legend boxoff
box off
set(gca,'xtickLabel',[])
ylabel('uV')

subplot(5,1,3)
plot(time(xt1:xt2),ctxLFPmuaenvpsmooth(xt1:xt2),'k') 
axis([-inf inf 0 inf])
legend('Smooth MUA','Location', 'best')
legend boxoff
box off
% set(gca,'xtickLabel',[])
ylabel('uV')

subplot(514)
hold on
plot(spkClust(CtxTetIdx(1)).spkTime, ... 
    ones(length(spkClust(CtxTetIdx(1)).spkTime),1),'r.','MarkerSize',5)
for kk = 2: size(CtxTetIdx,1)
plot(spkClust(CtxTetIdx(kk)).spkTime, ...
    ones(length(spkClust(CtxTetIdx(kk)).spkTime),1) * kk,... 
    'r.','MarkerSize',5)
end
hold off
xlim([time1 time2])
ylim([0 inf])
% set(gca,'xtickLabel',[])
xlabel('Time[sec]')
ylabel('MUA #')

subplot(515)
plot(time(xt1:xt2),CtxspkDensity(xt1:xt2),'r','LineWidth',2) 
xlabel('Time[sec]')
ylabel('spks/s')
axis tight

%% Distribution of MUA 
time1= [];
time2= [];
xt1= time1/sint;
xt2= time2/sint;
close all
h= histogram(ctxLFPmuaenvpsmooth(xt1:xt2),200,'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',1);
% h.BinEdges = 0:30;
% h.BinWidth = 0.01;
h.EdgeColor = 'k';
xlabel('MUA spike Density Hz')
ylabel('Probability Density')
legend('MUA spike Density')
legend('boxoff')  
box off
set(gca,'FontSize', 14)

%%  FIRST CRITERION: AT LEAST 3*STD OF FILTERED, SQUARED LFP
% splog = spkDensity' < thr1;
thr1 = []; % threshold
splog = ctxLFPmuaenvpsmooth <= thr1;
pos = find(splog == 1);

% SECOND CRITERION: at least a difference of 15ms between the end of
...a ripple and the beginning of the next to be considered separate events

difpos = find(diff(pos) >= (0.02/sint)); % events with inte4rvals > 20ms
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
            
% THIRD CRITERION: Eliminate 'Delta' candidates with < 1 sampling interval
duration = stepDelta(:,2)- stepDelta(:,1);
stepDelta1 = stepDelta((~(duration<=1)),:);

% FOURTH CRITERION:  at least 20ms long
duration2 = (stepDelta1(:,2) - stepDelta1(:,1))*sint; %duration in sec.
deltaIdx = stepDelta1((find(duration2<=1 & ...
    duration2 >= 0.025)),:); %Events < 150ms and > 15ms (Indices)
% Down-State intial and final idx from MUA threshold(First threshold)
below = deltaIdx(:,1); % - (0.5/sint);
abov = deltaIdx(:,2); % + (0.5/sint);
below(below < 0) = 1; % Substitute start negative value with 1
abov(abov > length(ctxLFP)) = length(ctxLFP); % Making delta candidates finishing at the end of the LFP


%% Visual Inspection of LFP and Down-State Candidates
close all
time1= [];
time2= [];
xt1= time1/sint;
xt2= time2/sint;
interBelowIdx = find(below.*sint >= time1-5 & below.*sint <= time2+5);
figure

subplot(611)
plot(time(xt1:xt2), ctxLFP(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        ctxLFP(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -500 600])

subplot(612)
plot(time(xt1:xt2), ctxLFPmuaenvpsmooth(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        ctxLFPmuaenvpsmooth(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 0 inf])
legend('Smooth MUA','Location', 'best')
legend boxoff
box off

subplot(613)
plot(time(xt1:xt2), ctxDeltaLFP(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        ctxDeltaLFP(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -inf inf])
legend('filt Delta','Location', 'best')
legend boxoff
box off

subplot(614)
plot(time(xt1:xt2), ctxSpindleLFP(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        ctxSpindleLFP(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -inf inf])
legend('filt Spindle','Location', 'best')
legend boxoff
box off

subplot(615)
hold on
plot(spkClust(CtxTetIdx(1)).spkTime, ... 
    ones(length(spkClust(CtxTetIdx(1)).spkTime),1),'r.','MarkerSize',5)
for kk = 2: size(CtxTetIdx,1)
plot(spkClust(CtxTetIdx(kk)).spkTime, ...
    ones(length(spkClust(CtxTetIdx(kk)).spkTime),1) * kk,... 
    'r.','MarkerSize',5)
end
hold off
xlim([time1 time2])
ylim([0 inf])
% set(gca,'xtickLabel',[])
xlabel('Time[sec]')
ylabel('MUA #')

subplot(616)
plot(time(xt1:xt2),CtxspkDensity(xt1:xt2),'r','LineWidth',2) 
xlabel('Time[sec]')
ylabel('spks/s')
axis tight



%% Visual Inspection of Down-States candidates
close all
eventExt = 0.05 / sint;
for kk = 2:250 %: size(,1)
   figure
   subplot(211)
   plot(time(below(kk) - eventExt : abov(kk)+ eventExt), ...
       ctxLFP(below(kk) - eventExt : abov(kk)+ eventExt))
   axis([-inf inf, -inf inf])
   
   subplot(212)
   plot(time(below(kk) - eventExt : abov(kk)+ eventExt),...
       ctxLFPdeltaspin(below(kk) - eventExt : abov(kk)+ eventExt))
   axis([-inf inf, -inf inf])
   
end



%% correlation coefficient between the LFP candidates and the 1-16Hz filtered
%  Eliminate evemts that are lager than the LFP
eventExt = 0.05 / sint;
if abov(end)+eventExt > length(ctxLFP) 
    abov(end) = [];
    below(end) = [];
elseif below(1)-eventExt < 1
     abov(1) = [];
    below(1) = [];
else
end
RR = zeros(size(below,1),1);
for kk = 1: size(below,1)-1
   R = corrcoef((ctxLFP(below((kk)) - eventExt :... 
       abov((kk))+ eventExt)),... 
       (ctxLFPdeltaspin(below((kk)) - eventExt :... 
       abov((kk))+ eventExt)));
   RR(kk,1) = R(2,1);
end
RRint = round(RR,4);


%% Distribution of the correlation
close all
h= histogram(RRint,'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',2);
% h.BinEdges = -1:1;
% h.BinWidth = 0.01;
h.EdgeColor = [213/300 29/300 91/300];
xlabel(' correlation coefficient')
ylabel('Probability Density')
legend('Short duration Down-States Candidates','Location','best')
legend('boxoff') 
xlim([0 1])
box off
%% Correlation Threhold
abovCorr = abov;
belowCorr = below;
noiseCandIdx = find(RR <= []); % choose the Th
abovCorr((noiseCandIdx)) = [];
belowCorr((noiseCandIdx)) = [];


%% Down-state duration distribution visualization after correlation filtering
close all
deltaCanddurRR = abovCorr - belowCorr;
h= histogram((deltaCanddurRR .* sint),'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',2);
% h.BinEdges = 0:1;
h.BinWidth = 0.01;
h.EdgeColor = 'r';
xlabel('Down-state duration (s)')
ylabel('Probability Density')
legend('Down-State duriation')
legend('boxoff') 
xlim([0 0.6])
box off
set(gca,'FontSize', 14)
%% Visual Inspection
close all
time1= [];
time2= [];
xt1= time1/sint;
xt2= time2/sint;
interBelowIdx = find(belowCorr.*sint >= time1-5 & belowCorr.*sint <= time2+5);
figure
subplot(211)
plot(time(xt1:xt2), ctxLFP(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(belowCorr(interBelowIdx(kk)):abovCorr(interBelowIdx(kk))),...
        ctxLFP(belowCorr(interBelowIdx(kk)):abovCorr(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -500 600])

subplot(212)
plot(time(xt1:xt2), ctxLFPdeltaspin(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(belowCorr(interBelowIdx(kk)):abovCorr(interBelowIdx(kk))),...
        ctxLFPdeltaspin(belowCorr(interBelowIdx(kk)):abovCorr(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -500 600])
%% Visual Inspection of Down-States candidates
close all
eventExt = 0.025 / sint;
for kk = 1:50 %: size(,1)
   figure
   subplot(311)
   plot(time(belowCorr((kk)) - eventExt :... 
       abovCorr((kk))+ eventExt),...
       ctxLFP(belowCorr((kk)) - eventExt :... 
       abovCorr((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   ylabel('uV')
   legend('Raw LFP')
   legend('boxoff') 
   
   subplot(312)
   plot(time(belowCorr((kk)) - eventExt :... 
       abovCorr((kk))+ eventExt),...
       ctxLFPdeltaspin(belowCorr((kk)) - eventExt :... 
       abovCorr((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   xlabel('Time (s)')
   ylabel('uV')
   legend('1-16Hz')
   legend('boxoff') 
   
   subplot(313)
   plot(time(belowCorr((kk)) - eventExt :... 
       abovCorr((kk))+ eventExt),...
       ctxLFPdeltaspinenvp(belowCorr((kk)) - eventExt :... 
       abovCorr((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   xlabel('Time (s)')
   ylabel('uV')
   legend('1-16Hz')
   legend('boxoff') 
end



%% Local Minima and local maxima from delta/spin lfp
% eventExt = 0.05; % extension of the candidate
maxt = 1:length(abovCorr);
maxMax = cell(size(maxt,1),1);
minMin = cell(size(maxt,1),1);
for kk = 1:size(maxt,2)
    ll = maxt(kk);
    candDelt = ctxLFPdeltaspin(belowCorr((ll)) - eventExt : ... 
       abovCorr((ll))+ eventExt); % Local Min/Max Delta/Spin LFP
    maxDeltaAmpIdx = find(islocalmax(candDelt)); % find local max
    minDeltaAmpIdx = find(islocalmin(candDelt)); % find local min
    maxMax{kk,1} = max(candDelt(maxDeltaAmpIdx));
    minMin{kk,1} = min(candDelt(minDeltaAmpIdx));
end
% Local Min/Max Amplitude of Delta/Spin LFP
deltaAmpRR = cellfun(@minus,maxMax,minMin,'Un',0);
empties = cellfun('isempty',deltaAmpRR);
deltaAmpRR(empties) = {NaN};
deltaAmpRR = cell2mat(deltaAmpRR); 
nanDeltaAmpRRIdx = find(isnan(deltaAmpRR));


%% Peak Ampalitude Distribution of Down-State Candidates
close all
hold on
h= histogram(deltaAmpRR(:,1),'Normalization','probability',... 
    'DisplayStyle','stairs','LineWidth',1,'BinWidth',10,...
    'EdgeColor','b');
h2 = histogram(cell2mat(minMin),'Normalization','probability',...
    'DisplayStyle','stairs','LineWidth',1,'BinWidth',10,...
    'EdgeColor','k');
h3 = histogram(cell2mat(maxMax),'Normalization','probability',...
    'DisplayStyle','stairs','LineWidth',1,'BinWidth',10,...
    'EdgeColor','r');
xlabel(' 1-16Hz Peak Ampl. (uV)')
ylabel('Probability Density')
legend('Down-States Candidates','Location','best')
legend('boxoff') 
xlim([-inf inf])
box off

%%
abovCorrAmp = abovCorr;
belowCorrAmp = belowCorr;
noiseCandIdx = find(deltaAmpRR <= []); % choose the Th
abovCorrAmp((noiseCandIdx)) = [];
belowCorrAmp((noiseCandIdx)) = [];


%% Visual Inspection
close all
time1= [];
time2= [];
xt1= time1/sint;
xt2= time2/sint;
interBelowIdx = find(belowCorrAmp.*sint >= time1-5 & belowCorrAmp.*sint <= time2+5);
figure
subplot(211)
plot(time(xt1:xt2), ctxLFP(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(belowCorrAmp(interBelowIdx(kk)):abovCorrAmp(interBelowIdx(kk))),...
        ctxLFP(belowCorrAmp(interBelowIdx(kk)):abovCorrAmp(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -500 600])
subplot(212)
plot(time(xt1:xt2), ctxLFPmuaenvpsmooth(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(belowCorrAmp(interBelowIdx(kk)):abovCorrAmp(interBelowIdx(kk))),...
        ctxLFPmuaenvpsmooth(belowCorrAmp(interBelowIdx(kk)):abovCorrAmp(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -inf inf])



%% Visual Inspection of Down-States candidates
close all
eventExt = 0.05 / sint;
for kk = 150:250   %: size(,1)
   figure
   subplot(311)
   plot(time(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),...
       ctxLFP(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   ylabel('uV')
   legend('Raw LFP')
   legend('boxoff') 
   
   subplot(312)
   plot(time(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),...
       ctxLFPdeltaspin(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   xlabel('Time (s)')
   ylabel('Slope')
   legend('1-16Hz')
   legend('boxoff') 
   
   subplot(313)
   plot(time(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),...
       ctxLFPdeltaspinenvp(belowCorrAmp((kk)) - eventExt :... 
       abovCorrAmp((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   xlabel('Time (s)')
   ylabel('Slope')
   legend('1-16Hz')
   legend('boxoff') 
end
%% Down-state duration distribution visualization after correlation filtering
close all
diffdeltaCanddurRR = abovCorrAmp - belowCorrAmp;
h= histogram((diffdeltaCanddurRR .* sint),'Normalization','probability',... 
    'DisplayStyle','stair', 'LineWidth',2);
% h.BinEdges = 0:1;
h.BinWidth = 0.01;
h.EdgeColor = 'r';
xlabel('Down-state duration (s)')
ylabel('Probability Density')
legend('Down-State duriation')
legend('boxoff') 
xlim([0 0.6])
box off
set(gca,'FontSize', 14)
%% Max Amplitude of envelop delta/spin lfp 
eventExt = 25;
maxt = 1:length(abovCorrAmp);
deltEnvp = zeros(length(abovCorrAmp),1);
for kk = 1:size(maxt,2)
    deltEnvpmax = max(ctxLFPdeltaspinenvp(belowCorrAmp((kk)) - eventExt : ... 
       abovCorrAmp((kk))+ eventExt)); % Local Min/Max Delta/Spin LFP
   deltEnvpmin = min(ctxLFPdeltaspinenvp(belowCorrAmp((kk)) - eventExt : ... 
       abovCorrAmp((kk))+ eventExt)); % Local Min/Max Delta/Spin LFP
   deltEnvp(kk,1) = deltEnvpmax - deltEnvpmin;
    
end

%% Peak Ampalitude Distribution of Down-State Candidates
close all
hold on
h= histogram(deltEnvp(:,1),'Normalization','probability',... 
    'DisplayStyle','stairs','LineWidth',1,'BinWidth',5,...
    'EdgeColor','b');
xlabel(' 1-16Hz Derivative')
ylabel('Probability Density')
legend('Down-States Candidates','Location','best')
legend('boxoff') 
xlim([-inf inf]) 
box off 

%% Differential Amplitude Threshold
abovCorrAmpEnvp = abovCorrAmp;
belowCorrAmpEnvp = belowCorrAmp;
noiseCandIdx = find(deltEnvp <= []); % choose the Th
abovCorrAmpEnvp((noiseCandIdx)) = [];
belowCorrAmpEnvp((noiseCandIdx)) = [];

%% Visual Inspection
close all
time1= [];
time2= [];
xt1= time1/sint;
xt2= time2/sint;
interBelowIdx = find(belowCorrAmpEnvp.*sint >= time1-5 & belowCorrAmpEnvp.*sint <= time2+5);
figure
subplot(211)
plot(time(xt1:xt2), ctxLFP(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(belowCorrAmpEnvp(interBelowIdx(kk)):abovCorrAmpEnvp(interBelowIdx(kk))),...
        ctxLFP(belowCorrAmpEnvp(interBelowIdx(kk)):abovCorrAmpEnvp(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -500 600])

subplot(212)
plot(time(xt1:xt2), ctxLFPmuaenvpsmooth(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(belowCorrAmpEnvp(interBelowIdx(kk)):abovCorrAmpEnvp(interBelowIdx(kk))),...
       ctxLFPmuaenvpsmooth(belowCorrAmpEnvp(interBelowIdx(kk)):abovCorrAmpEnvp(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -inf inf])

%% Visual Inspection of Down-States candidates
close all
eventExt = 0.05 / sint;
for kk = 1:50    %: size(,1)
   figure
   subplot(211)
   plot(time(belowCorrAmpEnvp((kk)) - eventExt :... 
       abovCorrAmpEnvp((kk))+ eventExt),...
       ctxLFP(belowCorrAmpEnvp((kk)) - eventExt :... 
       abovCorrAmpEnvp((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   ylabel('uV')
   legend('Raw LFP')
   legend('boxoff') 
   
   subplot(212)
   plot(time(belowCorrAmpEnvp((kk)) - eventExt :... 
       abovCorrAmpEnvp((kk))+ eventExt),...
       ctxLFPdeltaspin(belowCorrAmpEnvp((kk)) - eventExt :... 
       abovCorrAmpEnvp((kk))+ eventExt),'k')
   axis([-inf inf, -inf inf])
   xlabel('Time (s)')
   ylabel('Slope')
   legend('1-16Hz')
   legend('boxoff') 
end
%% Down-state duration distribution visualization after correlation filtering

diffdeltaCanddurRR = abovCorrAmpEnvp - belowCorrAmpEnvp;
h= histogram((diffdeltaCanddurRR .* sint),'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',2);
% h.BinEdges = 0:1;
h.BinWidth = 0.01;
h.EdgeColor = 'r';
xlabel('Down-state duration (s)')
ylabel('Probability Density')
legend('Down-State duriation')
legend('boxoff') 
xlim([0 0.6])
box off
set(gca,'FontSize', 18)
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
h_Z = plot(time(belowCorrAmpEnvp),position(belowCorrAmpEnvp),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3);
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
% title('"Down-State" Activity')
ylabel('Linear Position[pixels]')
xlabel('Time[sec]')
axis tight
% legend([h_Z, h_y],'"Down-States"', 'Linear Pos')
% legend boxoff
set(gca,'FontSize', 18)

% subplot(121)
figure
h_Z = plot(pos2d(:,1),pos2d(:,2), 'Color', [0.5,0.5,0.5]);
hold on
h_y = plot(pos2d(belowCorrAmpEnvp,1),pos2d(belowCorrAmpEnvp,2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3);
hold off
% legend([h_Z, h_y],'Position','"Down-States"')
% legend boxoff
xlabel('x coor. [pixels]')
ylabel('y coor. [pixels]')
axis tight
set(gca,'FontSize', 18)
figure
% subplot(122)
hold on
h_y = plot(time,velocity, ...
    'Color', [0.5,0.5,0.5], 'LineWidth', 1);
h_Z = plot(time(belowCorrAmpEnvp),velocity(belowCorrAmpEnvp),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3);
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
% title('"Down-State" Activity')
ylabel('Vel. pixels/sec')
xlabel('Time[sec]')
axis tight
% legend([h_Z, h_y],'"Down-States"', 'Velocity')
% legend boxoff
set(gca,'FontSize', 18)
%% Clasification ripples during sleep and run time
runtime1 = 1524; % Run start time 
runtime2 = 5797; % Run start time 2
sleepDownstateIdxbelow =belowCorrAmpEnvp;
sleepDownstateIdxabov = abovCorrAmpEnvp;
noise = find(belowCorrAmpEnvp.*sint > runtime1 & belowCorrAmpEnvp.*sint < runtime2);
sleepDownstateIdxbelow(noise)= [];
sleepDownstateIdxabov(noise) = [];
% Pre Sleep Ripple Candidates
PreDownstateIdxbelow = belowCorrAmpEnvp;
PreDownstateIdxabov = abovCorrAmpEnvp;
noise = find(belowCorrAmpEnvp.*sint > runtime1);
PreDownstateIdxbelow(noise)= [];
PreDownstateIdxabov(noise) = [];
% Post Sleep
PostDownstateIdxbelow = belowCorrAmpEnvp;
PostDownstateIdxabov = abovCorrAmpEnvp;
noise = find(belowCorrAmpEnvp.*sint < runtime2);
PostDownstateIdxbelow(noise)= [];
PostDownstateIdxabov(noise) = [];
% Run Ripple Candidates
RUNbelow = belowCorrAmpEnvp;
RUNabov = abovCorrAmpEnvp;
noise = find(belowCorrAmpEnvp.*sint < runtime1 | belowCorrAmpEnvp.*sint > runtime2);
RUNbelow(noise)= [];
RUNabov(noise) = [];
%% Cortex Data
%LFPs
ctxData.lfps.raw = ctxLFP;
ctxData.lfps.delta = ctxDeltaLFP;
ctxData.lfps.spindle = ctxSpindleLFP;
ctxData.lfps.muaEnvp.raw =  ctxLFPmua;
ctxData.lfps.muaEnvp.smooth = ctxLFPmuaenvpsmooth; 
% MUA
ctxData.spikes.Id = CtxTetIdx;
ctxData.spikes.density = CtxspkDensity;
% Down-State Candidates Idx
ctxData.dsIdx.ctxDownSlAbovIdx = sleepDownstateIdxabov;
ctxData.dsIdx.ctxDownSlBelowIdx = sleepDownstateIdxbelow;
ctxData.dsIdx.ctxDownPreSlAbovIdx = PreDownstateIdxabov;
ctxData.dsIdx.ctxDownPreSlBelowIdx = PreDownstateIdxbelow;
ctxData.dsIdx.ctxDownPostSlAbovIdx = PostDownstateIdxabov;
ctxData.dsIdx.ctxDownPostSlBelowIdx = PostDownstateIdxbelow;
ctxData.dsIdx.ctxDownRUNbelow = RUNbelow;
ctxData.dsIdx.ctxDownRUNRUNabov = RUNabov;
save('ctxData.mat', 'ctxData')
% save('ctxDownStateDetectionData.mat','ctxLFP','time','ctxLFPdeltaspin','ctxDeltaLFP',...
%     'ctxSpindleLFP','ctxLFPmuaenvpsmooth','ctxDownSlAbovIdx','ctxDownSlBelowIdx',...
%     'ctxDownPreSlAbovIdx','ctxDownPreSlBelowIdx',...
%     'ctxDownPostSlAbovIdx','ctxDownPostSlBelowIdx',...
%     'ctxDownRUNbelow','ctxDownRUNRUNabov')
% clear
%% Plotting Ripples in the Maze
linPos = interp1(positiondata.time, positiondata.linearpos, time);
poscm = interp1(positiondata.time, positiondata.poscm, time);
figure
subplot(421)
hold on
plot(time, linPos, 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
plot(time(sleepDownstateIdxbelow),linPos(sleepDownstateIdxbelow),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
title('Ripple Activity During Sleep')
xlabel('Linear Position[cm]')
ylabel('Time[sec]')
axis tight
subplot(422)
plot(poscm(:,1),poscm(:,2), 'Color', [0.5,0.5,0.5])
hold on
plot(poscm(sleepDownstateIdxbelow,1),poscm(sleepDownstateIdxbelow,2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
xlabel('x coordinate [cm]')
ylabel('y coordinate [cm]')

subplot(423)
hold on
plot(time, linPos, 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
plot(time(PreDownstateIdxbelow),linPos(PreDownstateIdxbelow),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Ripple Activity Pre-Sleep')
xlabel('Linear Position[cm]')
ylabel('Time[sec]')
axis tight
subplot(424)
plot(poscm(:,1),poscm(:,2), 'Color', [0.5,0.5,0.5])
hold on
plot(poscm(PreDownstateIdxbelow,1),poscm(PreDownstateIdxbelow,2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
xlabel('x coordinate [cm]')
ylabel('y coordinate [cm]')

subplot(425)
hold on
plot(time, linPos, 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
plot(time(PostDownstateIdxbelow),linPos(PostDownstateIdxbelow),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Ripple Activity Post-Sleep')
xlabel('Linear Position[cm]')
ylabel('Time[sec]')
axis tight
subplot(426)
plot(poscm(:,1),poscm(:,2), 'Color', [0.5,0.5,0.5])
hold on
plot(poscm(PostDownstateIdxbelow,1),poscm(PostDownstateIdxbelow,2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
xlabel('x coordinate [cm]')
ylabel('y coordinate [cm]')

subplot(427)
hold on
plot(time, linPos, 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
plot(time(RUNbelow),linPos(RUNbelow),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Ripple Activity RUN')
xlabel('Linear Position[cm]')
ylabel('Time[sec]')
axis tight
subplot(428)
plot(poscm(:,1),poscm(:,2), 'Color', [0.5,0.5,0.5])
hold on
plot(poscm(RUNbelow,1),poscm(RUNbelow,2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
xlabel('x coordinate [cm]')
ylabel('y coordinate [cm]')