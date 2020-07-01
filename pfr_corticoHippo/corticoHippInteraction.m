% Ripple visualization
Fs = 1000;
sint = 1/Fs;
%% append data to the detection CA1
% save('CA1RippDetectionData.mat', 'hippRippleenv','CA1Data.lfps.rippleEnvZscore', '-append')

%% Obtaining the Ripple peak
peakRipple = zeros(size(CA1Data.rippleIdx.CA1RipSlBelowIdx,1),1);
for kk = 1: size(CA1Data.rippleIdx.CA1RipSlBelowIdx,1)-1
   [R, peakRipple(kk)] = max(CA1Data.lfps.rippleEnvZscore(CA1Data.rippleIdx.CA1RipSlBelowIdx(kk): ...
       CA1Data.rippleIdx.CA1RipSlAbovIdx(kk)));
end

% this is to obtain the center indice based on the max ripple power peak
peakRippleIdx = zeros(size(peakRipple,1),1);
for kk = 1: size(CA1Data.rippleIdx.CA1RipSlBelowIdx,1)-1
   v = CA1Data.rippleIdx.CA1RipSlBelowIdx(kk): CA1Data.rippleIdx.CA1RipSlAbovIdx(kk);
   peakRippleIdx(kk) = v(peakRipple(kk));
end

%% Before and after distributions
% maxRipplepks = zeros(size(CA1Data.rippleIdx.CA1RipSlBelowIdx,1),1);
% maxHippZspkDen = zeros(size(CA1Data.rippleIdx.CA1RipSlBelowIdx,1),1); % New Addition
maxMUApks = zeros(size(CA1Data.rippleIdx.CA1RipSlBelowIdx,1),1);
maxMUActx = zeros(size(CA1Data.rippleIdx.CA1RipSlBelowIdx,1),1);
for kk = 1: size(CA1Data.rippleIdx.CA1RipSlBelowIdx,1)-1
%    maxRipplepks(kk) = CA1Data.lfps.rippleEnvZscore(peakRippleIdx(kk));
%    maxHippZspkDen(kk) = zHippSpkDensity(peakRippleIdx(kk)); 
   maxMUApks(kk) = CA1Data.spikes.density(peakRippleIdx(kk));
   maxMUActx(kk) = ctxData.lfps.muaEnvp.smooth(peakRippleIdx(kk));
end

%% Distribution of MUA 
close all
figure
% h= histogram(maxRipplepks,'Normalization','probability',... 
%     'DisplayStyle','stairs', 'LineWidth',1);
% h.BinEdges = 0:15;
% h.BinWidth = 1;
% h.EdgeColor = 'k';
% xlabel('Z-Score Ripple Power')
% ylabel('Probability Density')
% box off

hold on
h1= histogram(maxMUApks,'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',1);
h1.BinEdges = 0:500;
h1.BinWidth = 2;
h1.EdgeColor = 'r';
xlabel('Max MUA ripple activation (Spk/sec)')
ylabel('Probability Density')
box off

h2= histogram(maxMUActx,'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',1);
h2.BinEdges = 0:500;
h2.BinWidth = 1;
h2.EdgeColor = 'b';
xlabel('Max MUA ripple activation (Spk/sec)')
ylabel('Probability Density')

% h3= histogram(maxHippZspkDen ,'Normalization','probability',... 
%     'DisplayStyle','stairs', 'LineWidth',1);
% h3.BinEdges = 0:200;
% h3.BinWidth = 0.1;
% h3.EdgeColor = 'c';
% xlabel('Zscore MUA spiking Rate')
% ylabel('Probability Density')
legend([h1 h2],'CA1 MUA','Ctx MUA')
legend('boxoff')  
box off
hold off
%% Eliminating ripple candidates with a standar deaviation < 4 and MUA 
... activity < 20 spks/sec
medianMUAdensity = median(CA1Data.spikes.density);   
% stdMUAdensity = mad(CA1Data.spikes.density,1);
thMUA = medianMUAdensity; % median deviation flag = 1 mean deviation flag =0

pksRippleIdx = peakRippleIdx;   % ripple peak Idx                                                            
pksRippleIdx(maxMUApks < thMUA) = [];
ctxMUA = maxMUActx; % Ctx MUA matrix
ctxMUA(maxMUApks < thMUA) = [];
RippleBelowIdx = CA1Data.rippleIdx.CA1RipSlBelowIdx; % Ripple Start Idx
RippleBelowIdx(maxMUApks < thMUA) = [];
RippleAboveIdx = CA1Data.rippleIdx.CA1RipSlAbovIdx; % Ripple End Idx
RippleAboveIdx(maxMUApks < thMUA) = [];

%%   Use of Ripple Peak time to calculate the rate of SWR according to Down-States 
CA1pksRippleIdx = pksRippleIdx; % Ripple Peak Idx
% save('CA1pkssleeprippIdx.mat','CA1pksRippleIdx')
CA1ripple.indices.all.peak = pksRippleIdx;
CA1ripple.indices.all.start = RippleBelowIdx;
CA1ripple.indices.all.end = RippleAboveIdx;

%%
eventExt = 0.1 / sint;
tt = -0.1:sint:0.1;
trials = 1:length(ctxMUA);
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
ctxMUAmatrix = zeros(size(ctxMUA,1),Veclenght);
for kk = 1: size(ctxMUA,1)
  ctxMUAmatrix(kk,:) = ctxData.lfps.muaEnvp.smooth(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt);
end

% Cortex spiking activity
ctxSPKmatrix = zeros(size(ctxMUA,1),Veclenght);
for kk = 1: size(ctxMUA,1)
  ctxSPKmatrix(kk,:) = ctxData.spikes.density(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt);
end

% Hippocampal MUA activity
hippMUAmatrix = zeros(size(pksRippleIdx,1),Veclenght);
for kk = 1: size(pksRippleIdx,1)
  hippMUAmatrix(kk,:) = CA1Data.spikes.density(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt);
end
%% saving data for shuffle
ctxData.MUA.ctxMatrix = ctxMUAmatrix;
CA1Data.MUA.hippMatrix = hippMUAmatrix;
save('ctxData.mat', 'ctxData')
save('CA1Data.mat', 'CA1Data')

%% Cortex MUAenvp and Hipp MUA activity visualization
subplot(121)
imagesc(tt,trials,ctxMUAmatrix)
c = colorbar;
c.Label.String = 'MUA envelope (uV)';
c.FontSize = 11;
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
ylabel('Ripple #')
caxis([20 75])
colormap hot
set(gca,'FontSize',11)

subplot(122)
imagesc(tt,trials,hippMUAmatrix)
c = colorbar;
c.Label.String = 'CA1 Spike Rate (spk/sec)';
c.FontSize = 11;
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
caxis([0 25])
colormap hot
set(gca,'FontSize',11)
% colorbar('off')


%%
close all
% Adding the function jbfill to the path
% addpath(genpath('C:\Users\pdrfl\MATLAB Drive\pfr_AnalysisScripts\pfr_plots'));
ctxMUAavg = mean(ctxMUAmatrix,1);
ctxMUAerror = mad(ctxMUAmatrix,1);
hippMUAavg = mean(hippMUAmatrix,1);
hippMUAerror = mad(hippMUAmatrix,1);
subplot(211)
hold on
ctx = plot(tt, ctxMUAavg, 'r', 'LineWidth',2);
plot(tt, ctxMUAavg + ctxMUAerror, '--r')
plot(tt, ctxMUAavg - ctxMUAerror, '--r')
jbfill(tt,[ctxMUAavg + ctxMUAerror],[ctxMUAavg - ctxMUAerror],... 
    'r','r',0,0.1);
line([0 0],[0,80], 'Color','k','LineStyle','--')
legend(ctx, 'Cortex MUA')
legend boxoff
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA envp (uV)')
axis([-inf inf 0 75])
hold off
subplot(212)
hold on
hipp = plot(tt, hippMUAavg, 'b', 'LineWidth',2);
plot(tt, hippMUAavg + hippMUAerror, '--b')
plot(tt, hippMUAavg - hippMUAerror, '--b')
jbfill(tt,hippMUAavg + hippMUAerror,hippMUAavg - hippMUAerror,... 
    'b','b',0,0.1);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA (spk/sec)')
legend(hipp, 'Hipp MUA')
legend boxoff
axis([-inf inf 0 30])

figure
hold on
ctx = plot(tt, ctxMUAavg, 'r', 'LineWidth',2);
plot(tt, ctxMUAavg + ctxMUAerror, '--r')
plot(tt, ctxMUAavg - ctxMUAerror, '--r')
jbfill(tt,[ctxMUAavg + ctxMUAerror],[ctxMUAavg - ctxMUAerror],... 
    'r','r',0,0.1);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hipp = plot(tt, hippMUAavg, 'b', 'LineWidth',2);
plot(tt, hippMUAavg + hippMUAerror, '--b')
plot(tt, hippMUAavg - hippMUAerror, '--b')
jbfill(tt,hippMUAavg + hippMUAerror,hippMUAavg - hippMUAerror,... 
    'b','b',0,0.1);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'CA1 MUA')
legend boxoff
axis([-inf inf 0 50])
set(gca,'FontSize',11)

%% Calculate Cortical activity before, during, and after

ctxBelowmeddian = zeros(size(pksRippleIdx,1),1);
ctxCentermeddian = zeros(size(pksRippleIdx,1),1);
ctxAbovmeddian = zeros(size(pksRippleIdx,1),1);
for kk = 1:size(pksRippleIdx,1)
    ctxBelowmeddian(kk) = median(ctxMUAmatrix(kk,1:50));
    ctxCentermeddian(kk) = median(ctxMUAmatrix(kk,75:125));
    ctxAbovmeddian(kk) = median(ctxMUAmatrix(kk,150:201));
end

figure
hold on
h= histogram(ctxBelowmeddian,'Normalization','probability',... 
    'BinWidth', 2,'BinEdges',0:100,'DisplayStyle','stairs', 'LineWidth',2);
h1= histogram(ctxCentermeddian,'Normalization','probability',... 
    'BinWidth', 2,'BinEdges',0:100,'DisplayStyle','stairs', 'LineWidth',2);
h2= histogram(ctxAbovmeddian,'Normalization','probability',... 
    'BinWidth', 2,'BinEdges',0:100,'DisplayStyle','stairs', 'LineWidth',2);
hold off

h.BinWidth  = 2;
h1.BinWidth = 2;
h2.BinWidth = 2;
xlabel('Ctx MUA envp (uV)')
ylabel('Probability Density')
legend([h, h1, h2], 'Ctx MUA Before CA1 Ripple',...
    'Ctx MUA During CA1 Ripple','Ctx MUA After CA1 Ripple')
legend('boxoff')  
set(gca,'FontSize',11)
box off


%% Obtaining all the Ripple that occured during an up-state
close all
ripplePeakUpStateIdx = zeros(size(pksRippleIdx,1),1); % Ripple peak during up-state
rippleStartUpStateIdx = zeros(size(pksRippleIdx,1),1); % Ripple start during up-state
rippleEndUpStateIdx = zeros(size(pksRippleIdx,1),1); % Ripple end during up-state

for kk = 1: size(pksRippleIdx,1)
   if  ctxCentermeddian(kk) > 27 % if ripple occur during upstate
     if   ctxBelowmeddian(kk) < 27 &&  ctxAbovmeddian(kk) < 27
         ripplePeakUpStateIdx(kk) = pksRippleIdx(kk);
         rippleStartUpStateIdx(kk) = RippleBelowIdx(kk);
         rippleEndUpStateIdx(kk) = RippleAboveIdx(kk);
     elseif ctxBelowmeddian(kk) < 27 &&  ctxAbovmeddian(kk) > 27
         ripplePeakUpStateIdx(kk) = pksRippleIdx(kk);
         rippleStartUpStateIdx(kk) = RippleBelowIdx(kk);
         rippleEndUpStateIdx(kk) = RippleAboveIdx(kk);
     elseif ctxBelowmeddian(kk) > 27 &&  ctxAbovmeddian(kk) < 27
         ripplePeakUpStateIdx(kk) = pksRippleIdx(kk);
         rippleStartUpStateIdx(kk) = RippleBelowIdx(kk);
         rippleEndUpStateIdx(kk) = RippleAboveIdx(kk);
     elseif ctxBelowmeddian(kk) > 27 &&  ctxAbovmeddian(kk) > 27
         ripplePeakUpStateIdx(kk) = pksRippleIdx(kk); 
         rippleStartUpStateIdx(kk) = RippleBelowIdx(kk);
         rippleEndUpStateIdx(kk) = RippleAboveIdx(kk);
     else
     end
   end
end
ripplePeakUpStateIdx(ripplePeakUpStateIdx== 0) = [];
rippleStartUpStateIdx(rippleStartUpStateIdx== 0) = [];
rippleEndUpStateIdx(rippleEndUpStateIdx== 0) = [];

CA1ripple.indices.upstate.all.peak = ripplePeakUpStateIdx;
CA1ripple.indices.upstate.all.start = rippleStartUpStateIdx;
CA1ripple.indices.upstate.all.end = rippleEndUpStateIdx;
%% Classification of ripple occurence during "up-states configurations"
% Ripples that ocurred during "early up-state"
earlyUp = zeros(size(pksRippleIdx,1),1);
earlyUpRippleIdx = zeros(size(pksRippleIdx,1),2);
% Ripples that ocurred during "late up-state"
lateUp = zeros(size(pksRippleIdx,1),1);
lateUpRippleIdx = zeros(size(pksRippleIdx,1),2);
% Ripples that ocurred during "between up-state"
betweenUp = zeros(size(pksRippleIdx,1),1);
betweenUpRippleIdx = zeros(size(pksRippleIdx,1),2);
% Ripples that ocurred during "up up-state"
upUp = zeros(size(pksRippleIdx,1),1);
upUpRippleIdx = zeros(size(pksRippleIdx,1),2);



for kk = 1: size(pksRippleIdx,1)
   if  ctxCentermeddian(kk) > 27 % if ripple occur during upstate
     if   ctxBelowmeddian(kk) < 27 &&  ctxAbovmeddian(kk) < 27
         betweenUp(kk) = pksRippleIdx(kk);
         betweenUpRippleIdx(kk,1) = RippleBelowIdx(kk);
         betweenUpRippleIdx(kk,2) = RippleAboveIdx(kk);
     elseif ctxBelowmeddian(kk) < 27 &&  ctxAbovmeddian(kk) > 27
         earlyUp(kk) = pksRippleIdx(kk);
         earlyUpRippleIdx(kk,1) = RippleBelowIdx(kk);
         earlyUpRippleIdx(kk,2) = RippleAboveIdx(kk);
     elseif ctxBelowmeddian(kk) > 27 &&  ctxAbovmeddian(kk) < 27
         lateUp(kk) = pksRippleIdx(kk);
         lateUpRippleIdx(kk,1) = RippleBelowIdx(kk);
         lateUpRippleIdx(kk,2) = RippleAboveIdx(kk);
     elseif ctxBelowmeddian(kk) > 27 &&  ctxAbovmeddian(kk) > 27
         upUp(kk) = pksRippleIdx(kk); 
         upUpRippleIdx(kk,1) = RippleBelowIdx(kk);
         upUpRippleIdx(kk,2) = RippleAboveIdx(kk);
     end
   end
end
 
earlyUp(earlyUp == 0) = [];
earlyUpRippleIdx = (earlyUpRippleIdx(find(earlyUpRippleIdx(:,1)),:));
lateUp(lateUp == 0) = [];
lateUpRippleIdx = (lateUpRippleIdx(find(lateUpRippleIdx(:,1)),:));
betweenUp(betweenUp == 0) = [];
betweenUpRippleIdx = (betweenUpRippleIdx(find(betweenUpRippleIdx(:,1)),:));
upUp(upUp == 0) = [];
upUpRippleIdx = (upUpRippleIdx(find(upUpRippleIdx(:,1)),:));
% 
CA1ripple.indices.upstate.earlyUp.peak = earlyUp;
CA1ripple.indices.upstate.earlyUp.start = earlyUpRippleIdx(:,1);
CA1ripple.indices.upstate.earlyUp.end = earlyUpRippleIdx(:,2);

CA1ripple.indices.upstate.lateUp.peak = lateUp;
CA1ripple.indices.upstate.lateUp.start = lateUpRippleIdx(:,1);
CA1ripple.indices.upstate.lateUp.end = lateUpRippleIdx(:,2);

CA1ripple.indices.upstate.betweenUp.peak = betweenUp;
CA1ripple.indices.upstate.betweenUp.start = betweenUpRippleIdx(:,1);
CA1ripple.indices.upstate.betweenUp.end = betweenUpRippleIdx(:,2);

CA1ripple.indices.upstate.upUp.peak = upUp;
CA1ripple.indices.upstate.upUp.start = upUpRippleIdx(:,1);
CA1ripple.indices.upstate.upUp.end = upUpRippleIdx(:,2);

%% Matrix calculation for Ripples occuring during Up-states
eventExt = 0.1 / sint;
tt = -0.1:sint:0.1;
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
earlyUpMUAmat = zeros(size(earlyUp,1),Veclenght);
hippearlyUpMUAmat = zeros(size(earlyUp,1),Veclenght);
hippearlyUpRipplemat = zeros(size(earlyUp,1),Veclenght);
for kk = 1: size(earlyUp,1)
    earlyUpMUAmat(kk,:) = ctxData.lfps.muaEnvp.smooth(earlyUp(kk) - eventExt : earlyUp(kk)+ eventExt);
    hippearlyUpMUAmat(kk,:) = CA1Data.spikes.density(earlyUp(kk) - eventExt : earlyUp(kk)+ eventExt);
    hippearlyUpRipplemat(kk,:) = CA1Data.lfps.rippleEnvZscore(earlyUp(kk) - eventExt : earlyUp(kk)+ eventExt);
end

lateUpMUAmat= zeros(size(lateUp,1),Veclenght);
hipplateUpMUAmat = zeros(size(lateUp,1),Veclenght);
hipplateUpRipplemat = zeros(size(lateUp,1),Veclenght);
for kk = 1: size(lateUp,1)
    lateUpMUAmat(kk,:) = ctxData.lfps.muaEnvp.smooth(lateUp(kk) - eventExt : lateUp(kk)+ eventExt);
    hipplateUpMUAmat(kk,:) = CA1Data.spikes.density(lateUp(kk) - eventExt : lateUp(kk)+ eventExt);
    hipplateUpRipplemat(kk,:) = CA1Data.lfps.rippleEnvZscore(lateUp(kk) - eventExt : lateUp(kk)+ eventExt);
end

betweenUpMUAmat= zeros(size(betweenUp,1),Veclenght);
hippbetweenUpMUAmat = zeros(size(betweenUp,1),Veclenght);
hippbetweenUpRipplemat = zeros(size(betweenUp,1),Veclenght);
for kk = 1: size(betweenUp,1)
    betweenUpMUAmat(kk,:) = ctxData.lfps.muaEnvp.smooth(betweenUp(kk) - eventExt : betweenUp(kk)+ eventExt);
    hippbetweenUpMUAmat(kk,:) = CA1Data.spikes.density(betweenUp(kk) - eventExt : betweenUp(kk)+ eventExt);
    hippbetweenUpRipplemat(kk,:) = CA1Data.lfps.rippleEnvZscore(betweenUp(kk) - eventExt : betweenUp(kk)+ eventExt);
end

upUpMUAmat= zeros(size(upUp,1),Veclenght);
hippupUpMUAmat = zeros(size(upUp,1),Veclenght);
hippupUpRipplemat = zeros(size(upUp,1),Veclenght);
for kk = 1: size(upUp,1)
    upUpMUAmat(kk,:) = ctxData.lfps.muaEnvp.smooth(upUp(kk) - eventExt : upUp(kk)+ eventExt);
    hippupUpMUAmat(kk,:) = CA1Data.spikes.density(upUp(kk) - eventExt : upUp(kk)+ eventExt);
    hippupUpRipplemat(kk,:) = CA1Data.lfps.rippleEnvZscore(upUp(kk) - eventExt : upUp(kk)+ eventExt);
end

subplot(141)
trials = 1:length(earlyUp);
imagesc(tt,trials,earlyUpMUAmat)
c = colorbar;
axis xy
ylim([0 inf])
% xlabel('Time relative to SWR peak (sec)');
ylabel('Ripple event #')
caxis([20 75])
colormap hot
colorbar('off')
set(gca,'FontSize',11)

subplot(143)
trials = 1:length(lateUp);
imagesc(tt,trials,lateUpMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
% xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 75])
colormap hot
colorbar('off')
set(gca,'FontSize',11)

subplot(142)
trials = 1:length(betweenUp);
imagesc(tt,trials,betweenUpMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
% xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 75])
colormap hot
colorbar('off')
set(gca,'FontSize',11)

subplot(144)
trials = 1:length(upUp);
imagesc(tt,trials,upUpMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
c.FontSize = 11;
axis xy
ylim([0 inf])
% xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 75])
colormap hot
colorbar('off')
set(gca,'FontSize',11)
%%
corticoHippData.MUAmatrix.ctx.upstate.late = lateUpMUAmat;
corticoHippData.MUAmatrix.ctx.upstate.early = earlyUpMUAmat;
corticoHippData.MUAmatrix.ctx.upstate.upUp = upUpMUAmat;
corticoHippData.MUAmatrix.ctx.upstate.between = betweenUpMUAmat;
%%

lateUpMUAmatctxMUAmean = mean(lateUpMUAmat,1);
lateUpMUAmatctxMUAerror = mad(lateUpMUAmat,1);
hipplateUpMUAmatmean = mean(hipplateUpMUAmat,1);
hipplateUpMUAmaterror = mad(hipplateUpMUAmat,1);
hipplateUpRipplematmean = mean(hipplateUpRipplemat,1);
hipplateUpRipplematerror = mad(hipplateUpRipplemat,1);

figure
hold on
ctx = plot(tt, lateUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,lateUpMUAmatctxMUAmean + lateUpMUAmatctxMUAerror,... 
    lateUpMUAmatctxMUAmean - lateUpMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hipplateUpMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hipplateUpMUAmatmean + hipplateUpMUAmaterror,... 
    hipplateUpMUAmatmean - hipplateUpMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])
set(gca,'FontSize',11)

earlyUpMUAmatctxMUAmean = mean(earlyUpMUAmat,1);
earlyUpMUAmatctxMUAerror = mad(earlyUpMUAmat,1);
hippearlyUpMUAmatmean = mean(hippearlyUpMUAmat,1);
hippearlyUpMUAmaterror = mad(hippearlyUpMUAmat,1);
hippearlyUpRipplematmean = mean(hippearlyUpRipplemat,1);
hippearlyUpRipplematerror = mad(hippearlyUpRipplemat,1);
figure
hold on
ctx = plot(tt, earlyUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,earlyUpMUAmatctxMUAmean + earlyUpMUAmatctxMUAerror,... 
    earlyUpMUAmatctxMUAmean - earlyUpMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hippearlyUpMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippearlyUpMUAmatmean + hippearlyUpMUAmaterror,... 
    hippearlyUpMUAmatmean - hippearlyUpMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])
set(gca,'FontSize',11)

betweenUpMUAmatctxMUAmean = mean(betweenUpMUAmat,1);
betweenUpMUAmatctxMUAerror = mad(betweenUpMUAmat,1);
hippbetweenUpMUAmatmean = mean(hippbetweenUpMUAmat,1);
hippbetweenUpMUAmaterror = mad(hippbetweenUpMUAmat,1);
hippbetweenUpRipplematmean = mean(hippbetweenUpRipplemat,1);
hippbetweenUpRipplematerror = mad(hippbetweenUpRipplemat,1);
figure
hold on
ctx = plot(tt, betweenUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,betweenUpMUAmatctxMUAmean + betweenUpMUAmatctxMUAerror,... 
    betweenUpMUAmatctxMUAmean - betweenUpMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hippbetweenUpMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippbetweenUpMUAmatmean + hippbetweenUpMUAmaterror,... 
    hippbetweenUpMUAmatmean - hippbetweenUpMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])
set(gca,'FontSize',11)

upUpMUAmatctxMUAmean = mean(upUpMUAmat,1);
upUpMUAmatctxMUAerror = mad(upUpMUAmat,1);
hippupUpMUAmatmean = mean(hippupUpMUAmat,1);
hippupUpMUAmaterror = mad(hippupUpMUAmat,1);
hippupUpRipplematmean = mean(hippupUpRipplemat,1);
hippupUpRipplematerror = mad(hippupUpRipplemat,1);
figure
hold on
ctx = plot(tt, upUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,upUpMUAmatctxMUAmean + upUpMUAmatctxMUAerror,... 
    upUpMUAmatctxMUAmean - upUpMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hippupUpMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippupUpMUAmatmean + hippupUpMUAmaterror,... 
    hippupUpMUAmatmean - hippupUpMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])
set(gca,'FontSize',11)

%% Obtaining all the Ripple that occured during an Down-state
close all
ripplePeakDownStateIdx = zeros(size(pksRippleIdx,1),1); % Ripple peak during Down-state
rippleStartDownStateIdx = zeros(size(pksRippleIdx,1),1); % Ripple start during Down-state
rippleEndDownStateIdx = zeros(size(pksRippleIdx,1),1); % Ripple end during Down-state

for kk = 1: size(pksRippleIdx,1)
   if  ctxCentermeddian(kk) < 27 % if ripple occur during Downstate
     if   ctxBelowmeddian(kk) > 27 &&  ctxAbovmeddian(kk) > 27
         ripplePeakDownStateIdx(kk) = pksRippleIdx(kk);
         rippleStartDownStateIdx(kk) = RippleBelowIdx(kk);
         rippleEndDownStateIdx(kk) = RippleAboveIdx(kk);
     elseif ctxBelowmeddian(kk) > 27 &&  ctxAbovmeddian(kk) < 27
         ripplePeakDownStateIdx(kk) = pksRippleIdx(kk);
         rippleStartDownStateIdx(kk) = RippleBelowIdx(kk);
         rippleEndDownStateIdx(kk) = RippleAboveIdx(kk);
     elseif ctxBelowmeddian(kk) < 27 &&  ctxAbovmeddian(kk) > 27
         ripplePeakDownStateIdx(kk) = pksRippleIdx(kk);
         rippleStartDownStateIdx(kk) = RippleBelowIdx(kk);
         rippleEndDownStateIdx(kk) = RippleAboveIdx(kk);
     elseif ctxBelowmeddian(kk) < 27 &&  ctxAbovmeddian(kk) < 27
         ripplePeakDownStateIdx(kk) = pksRippleIdx(kk); 
         rippleStartDownStateIdx(kk) = RippleBelowIdx(kk);
         rippleEndDownStateIdx(kk) = RippleAboveIdx(kk);
     else
     end
   end
end
ripplePeakDownStateIdx(ripplePeakDownStateIdx== 0) = [];
rippleStartDownStateIdx(rippleStartDownStateIdx== 0) = [];
rippleEndDownStateIdx(rippleEndDownStateIdx== 0) = [];

CA1ripple.indices.downstate.all.peak = ripplePeakDownStateIdx;
CA1ripple.indices.downstate.all.start = rippleStartDownStateIdx;
CA1ripple.indices.downstate.all.end = rippleEndDownStateIdx;

% Classification of ripple occurence during "Down-states configurations"
% Ripples that ocurred during "early Down-state"
earlyDown = zeros(size(pksRippleIdx,1),1);
earlyDownRippleIdx = zeros(size(pksRippleIdx,1),2);
% Ripples that ocurred during "late Down-state"
lateDown = zeros(size(pksRippleIdx,1),1);
lateDownRippleIdx = zeros(size(pksRippleIdx,1),2);
% Ripples that ocurred during "between Down-state"
betweenDown = zeros(size(pksRippleIdx,1),1);
betweenDownRippleIdx = zeros(size(pksRippleIdx,1),2);
% Ripples that ocurred during "Down Down-state"
DownDown = zeros(size(pksRippleIdx,1),1);
DownDownRippleIdx = zeros(size(pksRippleIdx,1),2);



for kk = 1: size(pksRippleIdx,1)
   if  ctxCentermeddian(kk) < 27 % if ripple occur during Downstate
     if   ctxBelowmeddian(kk) > 27 &&  ctxAbovmeddian(kk) > 27
         betweenDown(kk) = pksRippleIdx(kk);
         betweenDownRippleIdx(kk,1) = RippleBelowIdx(kk);
         betweenDownRippleIdx(kk,2) = RippleAboveIdx(kk);
     elseif ctxBelowmeddian(kk) > 27 &&  ctxAbovmeddian(kk) < 27
         earlyDown(kk) = pksRippleIdx(kk);
         earlyDownRippleIdx(kk,1) = RippleBelowIdx(kk);
         earlyDownRippleIdx(kk,2) = RippleAboveIdx(kk);
     elseif ctxBelowmeddian(kk) < 27 &&  ctxAbovmeddian(kk) > 27
         lateDown(kk) = pksRippleIdx(kk);
         lateDownRippleIdx(kk,1) = RippleBelowIdx(kk);
         lateDownRippleIdx(kk,2) = RippleAboveIdx(kk);
     elseif ctxBelowmeddian(kk) < 27 &&  ctxAbovmeddian(kk) < 27
         DownDown(kk) = pksRippleIdx(kk); 
         DownDownRippleIdx(kk,1) = RippleBelowIdx(kk);
         DownDownRippleIdx(kk,2) = RippleAboveIdx(kk);
     else
     end
   end
end
 
earlyDown(earlyDown == 0) = [];
earlyDownRippleIdx = (earlyDownRippleIdx(find(earlyDownRippleIdx(:,1)),:));
lateDown(lateDown == 0) = [];
lateDownRippleIdx = (lateDownRippleIdx(find(lateDownRippleIdx(:,1)),:));
betweenDown(betweenDown == 0) = [];
betweenDownRippleIdx = (betweenDownRippleIdx(find(betweenDownRippleIdx(:,1)),:));
DownDown(DownDown == 0) = [];
DownDownRippleIdx = (DownDownRippleIdx(find(DownDownRippleIdx(:,1)),:));
% 
CA1ripple.indices.downstate.earlydown.peak = earlyDown;
CA1ripple.indices.downstate.earlydown.start = earlyDownRippleIdx(:,1);
CA1ripple.indices.downstate.earlydown.end = earlyDownRippleIdx(:,2);

CA1ripple.indices.downstate.latedown.peak = lateDown;
CA1ripple.indices.downstate.latedown.start = lateDownRippleIdx(:,1);
CA1ripple.indices.downstate.latedown.end = lateDownRippleIdx(:,2);

CA1ripple.indices.downstate.betweendown.peak = betweenDown;
CA1ripple.indices.downstate.betweendown.start = betweenDownRippleIdx(:,1);
CA1ripple.indices.downstate.betweendown.end = betweenDownRippleIdx(:,2);

CA1ripple.indices.downstate.downdown.peak = DownDown;
CA1ripple.indices.downstate.downdown.start = DownDownRippleIdx(:,1);
CA1ripple.indices.downstate.downdown.end = DownDownRippleIdx(:,2);

%% Matrix calculation for Ripples occuring during Down-states
eventExt = 0.1 / sint;
tt = -0.1:sint:0.1;
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
earlyDownMUAmat = zeros(size(earlyDown,1),Veclenght);
hippearlyDownMUAmat = zeros(size(earlyDown,1),Veclenght);
hippearlyDownRipplemat = zeros(size(earlyDown,1),Veclenght);
for kk = 1: size(earlyDown,1)
    earlyDownMUAmat(kk,:) = ctxData.lfps.muaEnvp.smooth(earlyDown(kk) - eventExt : earlyDown(kk)+ eventExt);
    hippearlyDownMUAmat(kk,:) = CA1Data.spikes.density(earlyDown(kk) - eventExt : earlyDown(kk)+ eventExt);
    hippearlyDownRipplemat(kk,:) = CA1Data.lfps.rippleEnvZscore(earlyDown(kk) - eventExt : earlyDown(kk)+ eventExt);
end

lateDownMUAmat= zeros(size(lateDown,1),Veclenght);
hipplateDownMUAmat = zeros(size(lateDown,1),Veclenght);
hipplateDownRipplemat = zeros(size(lateDown,1),Veclenght);
for kk = 1: size(lateDown,1)
    lateDownMUAmat(kk,:) = ctxData.lfps.muaEnvp.smooth(lateDown(kk) - eventExt : lateDown(kk)+ eventExt);
    hipplateDownMUAmat(kk,:) = CA1Data.spikes.density(lateDown(kk) - eventExt : lateDown(kk)+ eventExt);
    hipplateDownRipplemat(kk,:) = CA1Data.lfps.rippleEnvZscore(lateDown(kk) - eventExt : lateDown(kk)+ eventExt);
end

betweenDownMUAmat= zeros(size(betweenDown,1),Veclenght);
hippbetweenDownMUAmat = zeros(size(betweenDown,1),Veclenght);
hippbetweenDownRipplemat = zeros(size(betweenDown,1),Veclenght);
for kk = 1: size(betweenDown,1)
    betweenDownMUAmat(kk,:) = ctxData.lfps.muaEnvp.smooth(betweenDown(kk) - eventExt : betweenDown(kk)+ eventExt);
    hippbetweenDownMUAmat(kk,:) = CA1Data.spikes.density(betweenDown(kk) - eventExt : betweenDown(kk)+ eventExt);
    hippbetweenDownRipplemat(kk,:) = CA1Data.lfps.rippleEnvZscore(betweenDown(kk) - eventExt : betweenDown(kk)+ eventExt);
end

DownDownMUAmat= zeros(size(DownDown,1),Veclenght);
hippDownDownMUAmat = zeros(size(DownDown,1),Veclenght);
hippDownDownRipplemat = zeros(size(DownDown,1),Veclenght);
for kk = 1: size(DownDown,1)
    DownDownMUAmat(kk,:) = ctxData.lfps.muaEnvp.smooth(DownDown(kk) - eventExt : DownDown(kk)+ eventExt);
    hippDownDownMUAmat(kk,:) = CA1Data.spikes.density(DownDown(kk) - eventExt : DownDown(kk)+ eventExt);
    hippDownDownRipplemat(kk,:) = CA1Data.lfps.rippleEnvZscore(DownDown(kk) - eventExt : DownDown(kk)+ eventExt);
end
subplot(141)
trials = 1:length(earlyDown);
imagesc(tt,trials,earlyDownMUAmat)
c = colorbar;
axis xy
ylim([0 inf])
% xlabel('Time relative to SWR peak (sec)');
ylabel('Ripple event #')
caxis([20 75])
colormap hot
colorbar('off')
set(gca,'FontSize',11)

subplot(143)
trials = 1:length(lateDown);
imagesc(tt,trials,lateDownMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
% xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 75])
colormap hot
colorbar('off')
set(gca,'FontSize',11)

subplot(142)
trials = 1:length(betweenDown);
imagesc(tt,trials,betweenDownMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
ylim([0 inf])
% xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 75])
colormap hot
colorbar('off')
set(gca,'FontSize',11)

subplot(144)
trials = 1:length(DownDown);
imagesc(tt,trials,DownDownMUAmat)
c = colorbar;
c.Label.String = 'MUA envelope (uV)';
c.FontSize = 11;
axis xy
ylim([0 inf])
% xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 75])
ylim([1 inf])
colormap hot
colorbar('off')
set(gca,'FontSize',11)
%%
corticoHippData.MUAmatrix.ctx.downstate.late = lateDownMUAmat;
corticoHippData.MUAmatrix.ctx.downstate.early = earlyDownMUAmat;
corticoHippData.MUAmatrix.ctx.downstate.downDownUp = DownDownMUAmat;
corticoHippData.MUAmatrix.ctx.downstate.between = betweenDownMUAmat;

%%
lateDownMUAmatctxMUAmean = mean(lateDownMUAmat,1);
lateDownMUAmatctxMUAerror = mad(lateDownMUAmat,1);
hipplateDownMUAmatmean = mean(hipplateDownMUAmat,1);
hipplateDownMUAmaterror = mad(hipplateDownMUAmat,1);
hipplateDownRipplematmean = mean(hipplateDownRipplemat,1);
hipplateDownRipplematerror = mad(hipplateDownRipplemat,1);
close all
figure
hold on
ctx = plot(tt, lateDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,lateDownMUAmatctxMUAmean + lateDownMUAmatctxMUAerror,... 
    lateDownMUAmatctxMUAmean - lateDownMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hipplateDownMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hipplateDownMUAmatmean + hipplateDownMUAmaterror,... 
    hipplateDownMUAmatmean - hipplateDownMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])
set(gca,'FontSize',11)

%
earlyDownMUAmatctxMUAmean = mean(earlyDownMUAmat,1);
earlyDownMUAmatctxMUAerror = mad(earlyDownMUAmat,1);
hippearlyDownMUAmatmean = mean(hippearlyDownMUAmat,1);
hippearlyDownMUAmaterror = mad(hippearlyDownMUAmat,1);
hippearlyDownRipplematmean = mean(hippearlyDownRipplemat,1);
hippearlyDownRipplematerror = mad(hippearlyDownRipplemat,1);
figure
hold on
ctx = plot(tt, earlyDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,earlyDownMUAmatctxMUAmean + earlyDownMUAmatctxMUAerror,... 
    earlyDownMUAmatctxMUAmean - earlyDownMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hippearlyDownMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippearlyDownMUAmatmean + hippearlyDownMUAmaterror,... 
    hippearlyDownMUAmatmean - hippearlyDownMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])
set(gca,'FontSize',11)

%
betweenDownMUAmatctxMUAmean = mean(betweenDownMUAmat,1);
betweenDownMUAmatctxMUAerror = mad(betweenDownMUAmat,1);
hippbetweenDownMUAmatmean = mean(hippbetweenDownMUAmat,1);
hippbetweenDownMUAmaterror = mad(hippbetweenDownMUAmat,1);
hippbetweenDownRipplematmean = mean(hippbetweenDownRipplemat,1);
hippbetweenDownRipplematerror = mad(hippbetweenDownRipplemat,1);
figure
hold on
ctx = plot(tt, betweenDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,betweenDownMUAmatctxMUAmean + betweenDownMUAmatctxMUAerror,... 
    betweenDownMUAmatctxMUAmean - betweenDownMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hippbetweenDownMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippbetweenDownMUAmatmean + hippbetweenDownMUAmaterror,... 
    hippbetweenDownMUAmatmean - hippbetweenDownMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])
set(gca,'FontSize',11)

%
DownDownMUAmatctxMUAmean = mean(DownDownMUAmat,1);
DownDownMUAmatctxMUAerror = mad(DownDownMUAmat,1);
hippDownDownMUAmatmean = mean(hippDownDownMUAmat,1);
hippDownDownMUAmaterror = mad(hippDownDownMUAmat,1);
hippDownDownRipplematmean = mean(hippDownDownRipplemat,1);
hippDownDownRipplematerror = mad(hippDownDownRipplemat,1);
figure
hold on
ctx = plot(tt, DownDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,DownDownMUAmatctxMUAmean + DownDownMUAmatctxMUAerror,... 
    DownDownMUAmatctxMUAmean - DownDownMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')

hipp = plot(tt, hippDownDownMUAmatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippDownDownMUAmatmean + hippDownDownMUAmaterror,... 
    hippDownDownMUAmatmean - hippDownDownMUAmaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('MUA')
legend([ctx,hipp],'Ctx MUA', 'Hipp MUA')
legend boxoff
axis([-inf inf 0 75])
set(gca,'FontSize',11)

%% CA1 MUA during Up- and Down-States
hippMUAUpSatemat = [hipplateUpMUAmatmean; hippearlyUpMUAmatmean; ...
    hippbetweenUpMUAmatmean; hippupUpMUAmatmean];
hippMUAUpSatemean = mean(hippMUAUpSatemat,1);
hippMUAUpSateerror = std(hippMUAUpSatemat,1);

hippMUADownSatemat = [hipplateDownMUAmatmean; hippearlyDownMUAmatmean; ...
   hippbetweenDownMUAmatmean]; %hippDownDownMUAmatmean
hippMUADownSatemean = mean(hippMUADownSatemat,1);
hippMUADownSateerror = std(hippMUADownSatemat,1);
close all
figure
hold on
hipMUAUpState = plot(tt,hippMUAUpSatemean,... 
    'Color',[10/255 97/255 171/255], 'LineWidth',2);
jbfill(tt,hippMUAUpSatemean + hippMUAUpSateerror,... 
    hippMUAUpSatemean - hippMUAUpSateerror,... 
    [10/255 97/255 171/255],[10/255 97/255 171/255],0,0.2);

hipMUADownState = plot(tt, hippMUADownSatemean,...
    'Color',[213/255 29/255 91/255], 'LineWidth',2);
jbfill(tt,hippMUADownSatemean + hippMUADownSateerror,... 
    hippMUADownSatemean - hippMUADownSateerror,... 
    [213/255 29/255 91/255],[213/255 29/255 91/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('mean CA1 MUA (spk/sec) ')
legend([hipMUAUpState,hipMUADownState],'CA1 MUA during Up-States', ...
    'CA1 MUA during Down-States')
legend boxoff
axis([-inf inf 0 25])
set(gca,'FontSize',11)

%% Ripple power Up- and Down-States
hippRippleUpSatemat = [hipplateUpRipplematmean; hippearlyUpRipplematmean; ...
    hippbetweenUpRipplematmean; hippupUpRipplematmean];
hippRippleUpSatemean = mean(hippRippleUpSatemat,1);
hippRippleUpSateerror = std(hippRippleUpSatemat,1);

hippRippleDownSatemat = [hipplateDownRipplematmean; hippearlyDownRipplematmean; ...
   hippbetweenDownRipplematmean; hippDownDownRipplematmean];
hippRippleDownSatemean = mean(hippRippleDownSatemat,1);
hippRippleDownSateerror = std(hippRippleDownSatemat,1);

figure
hold on
hipRippleUpState = plot(tt,hippRippleUpSatemean,... 
    'Color',[10/255 97/255 171/255], 'LineWidth',2);
jbfill(tt,hippRippleUpSatemean + hippRippleUpSateerror,... 
    hippRippleUpSatemean - hippRippleUpSateerror,... 
    [10/255 97/255 171/255],[10/255 97/255 171/255],0,0.2);

hipRippleDownState = plot(tt, hippRippleDownSatemean,...
    'Color',[213/255 29/255 91/255], 'LineWidth',2);
jbfill(tt,hippRippleDownSatemean + hippRippleDownSateerror,... 
    hippRippleDownSatemean - hippRippleDownSateerror,... 
    [213/255 29/255 91/255],[213/255 29/255 91/255],0,0.2);
line([0 0],[0,8], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to SWR peak (sec)')
ylabel('mean Ripple(125-250Hz) Z-Score')
legend([hipRippleUpState,hipRippleDownState],'CA1 Ripple during Up-States', ...
    'CA1 Ripple during Down-States')
legend boxoff
axis([-inf inf -inf 20])
set(gca,'FontSize',11)

%%
totalRipdetected = size(lateUpMUAmat,1) + size(earlyUpMUAmat,1) ...
    + size(betweenUpMUAmat,1) + size(upUpMUAmat,1) + ...
    size(lateDownMUAmat,1) + size(earlyDownMUAmat,1) ...
    + size(betweenDownMUAmat,1) + size(DownDownMUAmat,1);
rippleUpStatePorcentage = (size(lateUpMUAmat,1) + size(earlyUpMUAmat,1) ...
    + size(betweenUpMUAmat,1) + size(upUpMUAmat,1))/totalRipdetected;
rippleDownStatePorcentage = (size(lateDownMUAmat,1) + size(earlyDownMUAmat,1) ...
    + size(betweenDownMUAmat,1) + size(DownDownMUAmat,1))/totalRipdetected;

corticoHippData.percentage.total.up = rippleUpStatePorcentage ;
corticoHippData.percentage.total.down = rippleDownStatePorcentage ;
%%
X = categorical({'"Up-State"','"Down-State"'});
X = reordercats(X,{'"Up-State"','"Down-State"'});
b = bar(X,[rippleUpStatePorcentage*100, rippleDownStatePorcentage*100]);
b.FaceColor = 'flat';
b.CData(2,:) = [213/255 29/255 91/255];
b.CData(1,:) = [10/255 1/255 171/255];
ylabel('CA1 Ripple Occurrence %')
ylim([0 100])
set(gca,'FontSize',11)
box off
%%
close all
subplot(121)
X = categorical({'Early Up-State','Late Up-State',...
    'Between Down-States', 'Up-States'});
X = reordercats(X,{'Early Up-State','Late Up-State',...
    'Between Down-States', 'Up-States'});
b = bar(X,[size(earlyUpMUAmat,1)/totalRipdetected*100,...
    size(lateUpMUAmat,1)/totalRipdetected*100,...
    size(betweenUpMUAmat,1)/totalRipdetected*100,...
    size(upUpMUAmat,1)/totalRipdetected*100],...
    'FaceColor', [10/255 1/255 171/255]);
corticoHippData.percentage.up.early = size(earlyUpMUAmat,1)/totalRipdetected*100;
corticoHippData.percentage.up.late = size(lateUpMUAmat,1)/totalRipdetected*100;
corticoHippData.percentage.up.between = size(betweenUpMUAmat,1)/totalRipdetected*100;
corticoHippData.percentage.up.up = size(upUpMUAmat,1)/totalRipdetected*100;
b.FaceColor = 'flat';
ylabel('Ripple Occurrence %')
ylim([0 65])
set(gca,'FontSize',11)
box off
subplot(122)
X = categorical({'Early Down-State','Late Down-State',...
    'Between Down-States', 'Down-States'});
X = reordercats(X,{'Early Down-State','Late Down-State',...
    'Between Down-States', 'Down-States'});
b = bar(X,[size(earlyDownMUAmat,1)/totalRipdetected*100,...
    size(lateDownMUAmat,1)/totalRipdetected*100,...
    size(betweenDownMUAmat,1)/totalRipdetected*100,...
    size(DownDownMUAmat,1)/totalRipdetected*100],...
    'FaceColor', [213/255 29/255 91/255]);
corticoHippData.percentage.down.early = size(earlyDownMUAmat,1)/totalRipdetected*100;
corticoHippData.percentage.down.late = size(lateDownMUAmat,1)/totalRipdetected*100;
corticoHippData.percentage.down.between = size(betweenDownMUAmat,1)/totalRipdetected*100;
corticoHippData.percentage.down.down = size(DownDownMUAmat,1)/totalRipdetected*100;
b.FaceColor = 'flat';
ylabel('Ripple Occurrence %')
set(gca,'FontSize',11)
ylim([0 65])
box off
%% verify the code

%%
save('CA1rippleClasification.mat', 'CA1ripple')
save('corticoHipp.mat', 'corticoHippData')
clear

%% Visual Inspection of Cortical activity based on ripple detection
close all
eventExt = 0.1 / sint;
for kk = 650:680 %: size(,1)
   figure
   subplot(611)
   plot(time(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt), ...
       CTXLFP(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),'k')
   axis([-inf inf, -inf inf])
   
   subplot(612)
   plot(time(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),...
       ctxData.lfps.muaEnvp.smooth(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),'k')
   axis([-inf inf, 0 inf])
   
   subplot(613)
   plot(time(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt), ...
       CA1LFP(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),'r')
   axis([-inf inf, -inf inf])
   
   subplot(614)
   plot(time(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt), ...
       CA1LFPripple(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),'r')
   axis([-inf inf, -inf inf])
   
   subplot(615)
   plot(time(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),...
       CA1Data.lfps.rippleEnvZscore(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),'r')
   axis([-inf inf, -inf inf])
   
   subplot(616)
    plot(time(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),... 
        CA1Data.spikes.density(pksRippleIdx(kk) - eventExt : pksRippleIdx(kk)+ eventExt),... 
        'r') 
    xlabel('Time[sec]')
    ylabel('spks/s')
    axis tight
   axis([-inf inf, 0 inf])
end










