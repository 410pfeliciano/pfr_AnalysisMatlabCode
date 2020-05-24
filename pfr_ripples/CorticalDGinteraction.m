% fastgamma visualization
Fs = 1000;
sint = 1/Fs;

%% Cortical activity
CTXLFP = ch42Tet13LFP;
ctxLFPmualfp = bandFilter (CTXLFP, 'highfreq'); % 100-450Hz Filter
ctxMUAenvp = abs(hilbert(ctxLFPmualfp));
ctxMUAenvpsm = smoothdata(ctxMUAenvp,'gaussian',20); 
%% Make sure the 
Tet = time; % just to calculate the duration of the recording
bins = (0 : 1/Fs : length(Tet) * (1/Fs) - (1/Fs))';
if exist('spkClust','var') == 0
    warning('spiking data spkClust is not in your workspace');
    return
end
TetCtx = 13; % choose the tetrode of interest
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
%% DG activity
smWin = 0.05; % smoothding window in seconds
hippfastgammaleenv = smoothdata(abs(hilbert(DGLFPfastgamma)),'gaussian',smWin*Fs); 
Zhippfastgammaleenv = zscore(hippfastgammaleenv);

%% Visual Inspection of Cortical activity based on fastgammale detection
% close all
% eventExt = 0.1 / sint;
% for kk = 1:10 %: size(,1)
%    figure
%    subplot(611)
%    plot(time(DGRipSlBelowIdx(kk) - eventExt : DGRipSlAbovIdx(kk)+ eventExt), ...
%        CTXLFP(DGRipSlBelowIdx(kk) - eventExt : DGRipSlAbovIdx(kk)+ eventExt),'k')
%    axis([-inf inf, -inf inf])
%    
%    subplot(612)
%    plot(time(DGRipSlBelowIdx(kk) - eventExt : DGRipSlAbovIdx(kk)+ eventExt),...
%        ctxMUAenvpsm(DGRipSlBelowIdx(kk) - eventExt : DGRipSlAbovIdx(kk)+ eventExt),'k')
%    axis([-inf inf, 0 inf])
%    
%    subplot(613)
%    plot(time(DGRipSlBelowIdx(kk) - eventExt : DGRipSlAbovIdx(kk)+ eventExt), ...
%        DGLFP(DGRipSlBelowIdx(kk) - eventExt : DGRipSlAbovIdx(kk)+ eventExt),'r')
%    axis([-inf inf, -inf inf])
%    
%    subplot(614)
%    plot(time(DGRipSlBelowIdx(kk) - eventExt : DGRipSlAbovIdx(kk)+ eventExt), ...
%        DGLFPfastgamma(DGRipSlBelowIdx(kk) - eventExt : DGRipSlAbovIdx(kk)+ eventExt),'r')
%    axis([-inf inf, -inf inf])
%    
%    subplot(615)
%    plot(time(DGRipSlBelowIdx(kk) - eventExt : DGRipSlAbovIdx(kk)+ eventExt),...
%        Zhippfastgammaleenv(DGRipSlBelowIdx(kk) - eventExt : DGRipSlAbovIdx(kk)+ eventExt),'r')
%    axis([-inf inf, -inf inf])
%    
%    subplot(616)
%     plot(time(DGRipSlBelowIdx(kk) - eventExt : DGRipSlAbovIdx(kk)+ eventExt),... 
%         DGspkDensity(DGRipSlBelowIdx(kk) - eventExt : DGRipSlAbovIdx(kk)+ eventExt),... 
%         'r','LineWidth',2) 
%     xlabel('Time[sec]')
%     ylabel('spks/s')
%     axis tight
%    axis([-inf inf, 0 inf])
% end


%% Obtaining the fastgammale peak
peakfastgamma = zeros(size(DGRipSlBelowIdx,1),1);
for kk = 1: size(DGRipSlBelowIdx,1)-1
   [R, peakfastgamma(kk)] = max(Zhippfastgammaleenv(DGRipSlBelowIdx(kk): DGRipSlAbovIdx(kk)));
end
% this is to obtain the center indice based on the max fastgamma power peak
peakfastgammaIdx = zeros(size(peakfastgamma,1),1);
for kk = 1: size(DGRipSlBelowIdx,1)-1
   v = DGRipSlBelowIdx(kk): DGRipSlAbovIdx(kk);
   peakfastgammaIdx(kk) = v(peakfastgamma(kk));
end
%% Before and after distributions
maxfastgammapks = zeros(size(DGRipSlBelowIdx,1),1);
maxMUApks = zeros(size(DGRipSlBelowIdx,1),1);
maxMUActx = zeros(size(DGRipSlBelowIdx,1),1);
for kk = 1: size(DGRipSlBelowIdx,1)-1
   maxfastgammapks(kk) = Zhippfastgammaleenv(peakfastgammaIdx(kk));
   maxMUApks(kk) = DGspkDensity(peakfastgammaIdx(kk));
   maxMUActx(kk) = ctxMUAenvpsm(peakfastgammaIdx(kk));
end
%% Distribution of MUA 
close all
figure
h= histogram(maxfastgammapks,'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',1);
h.BinEdges = 0:15;
h.BinWidth = 1;
h.EdgeColor = 'k';
xlabel('Z-Score fastgamma Power')
ylabel('Probability Density')
box off

hold on
h1= histogram(maxMUApks,'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',1);
h1.BinEdges = 0:100;
h1.BinWidth = 1;
h1.EdgeColor = 'r';
xlabel('Max MUA fastgamma activation (Spk/sec)')
ylabel('Probability Density')
box off

h2= histogram(maxMUActx,'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',1);
h2.BinEdges = 0:200;
h2.BinWidth = 1;
h2.EdgeColor = 'b';
xlabel('Max MUA fastgamma activation (Spk/sec)')
ylabel('Probability Density')
legend([h h1 h2],'fastgamma Peak Z-Score','DG MUA','Ctx MUA')
legend('boxoff')  
box off
hold off
%% Eliminating fastgamma candidates with a standar deaviation < 2

pksfastgammaIdx = peakfastgammaIdx;                                                               
pksfastgammaIdx(maxfastgammapks <= 2) = [];
ctxMUA = maxMUActx;
ctxMUA(maxfastgammapks <= 2) = [];
% to use to calculate the ripple rate according to down-states
DGpksGammaIdx = pksfastgammaIdx;
save('DGpkssleepGammaIdx.mat','DGpksGammaIdx')

%%
eventExt = 0.1 / sint;
tt = -0.1:sint:0.1;
trials = 1:length(ctxMUA);
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
ctxMUAmatrix = zeros(size(ctxMUA,1),Veclenght);
for kk = 1: size(ctxMUA,1)
  ctxMUAmatrix(kk,:) = ctxMUAenvpsm(pksfastgammaIdx(kk) - eventExt : pksfastgammaIdx(kk)+ eventExt);
end
% Cortex spiking activity
ctxSPKmatrix = zeros(size(ctxMUA,1),Veclenght);
for kk = 1: size(ctxMUA,1)
  ctxSPKmatrix(kk,:) = CtxspkDensity(pksfastgammaIdx(kk) - eventExt : pksfastgammaIdx(kk)+ eventExt);
end

% Hippocampal MUA activity
hippMUAmatrix = zeros(size(pksfastgammaIdx,1),Veclenght);
for kk = 1: size(pksfastgammaIdx,1)
  hippMUAmatrix(kk,:) = DGspkDensity(pksfastgammaIdx(kk) - eventExt : pksfastgammaIdx(kk)+ eventExt);
end

% Hippocampal MUA activity
hippFastGammamatrix = zeros(size(pksfastgammaIdx,1),Veclenght);
for kk = 1: size(pksfastgammaIdx,1)
  hippFastGammamatrix(kk,:) = Zhippfastgammaleenv(pksfastgammaIdx(kk) - eventExt : pksfastgammaIdx(kk)+ eventExt);
end


%%
% B = sort(hippMUAmatrix,1);
figure('DefaultAxesFontSize',10)
subplot(121)
imagesc(tt,trials,ctxMUAmatrix)
c = colorbar;
c.Label.String = 'MUA envelope(uV)';
c.FontSize = 11;
% colorbar('southoutside')
axis xy
ylim([0 inf])
xlabel('Time relative to FG peak (sec)');
ylabel('fastgamma #')
caxis([20 100])
colormap hot
set(gca,'FontSize',11)
% colorbar('off')
subplot(122)
imagesc(tt,trials,hippFastGammamatrix)
c = colorbar;
c.Label.String = 'DG Fast Gamma Z-Score';
c.FontSize = 11;
% colorbar('southoutside')
axis xy
ylim([0 inf])
xlabel('Time relative to FG peak (sec)');
% ylabel('fastgamma event #')
caxis([0 4])
colormap hot
set(gca,'FontSize',11)
% colorbar('off')


%%
close all
% Adding the function jbfill to the path
% addpath(genpath('C:\Users\pdrfl\MATLAB Drive\pfr_AnalysisScripts\pfr_plots'));
ctxMUAavg = mean(ctxMUAmatrix,1);
ctxMUAerror = mad(ctxMUAmatrix,1);
hippMUAavg = mean(hippFastGammamatrix,1);
hippMUAerror = mad(hippFastGammamatrix,1);
subplot(211)
hold on
ctx = plot(tt, ctxMUAavg, 'r', 'LineWidth',2);
plot(tt, ctxMUAavg + ctxMUAerror, '--r')
plot(tt, ctxMUAavg - ctxMUAerror, '--r')
jbfill(tt,[ctxMUAavg + ctxMUAerror],[ctxMUAavg - ctxMUAerror],... 
    'r','r',0,0.1);
line([0 0],[0,70], 'Color','k','LineStyle','--')
legend(ctx, 'Cortex MUA')
legend boxoff
xlabel('Time relative to FG peak (sec)')
ylabel('MUA envp (uV)')
axis([-inf inf 0 80])
set(gca,'FontSize',11)
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
xlabel('Time relative to FG peak (sec)')
ylabel('Z-Score')
legend(hipp, 'DG Fast-Gamma')
legend boxoff
axis([-inf inf 0 4])
set(gca,'FontSize',11)

% figure
% hold on
% ctx = plot(tt, ctxMUAavg, 'r', 'LineWidth',2);
% plot(tt, ctxMUAavg + ctxMUAerror, '--r')
% plot(tt, ctxMUAavg - ctxMUAerror, '--r')
% jbfill(tt,[ctxMUAavg + ctxMUAerror],[ctxMUAavg - ctxMUAerror],... 
%     'r','r',0,0.1);
% line([0 0],[0,100], 'Color','k','LineStyle','--')
% hipp = plot(tt, hippMUAavg, 'b', 'LineWidth',2);
% plot(tt, hippMUAavg + hippMUAerror, '--b')
% plot(tt, hippMUAavg - hippMUAerror, '--b')
% jbfill(tt,hippMUAavg + hippMUAerror,hippMUAavg - hippMUAerror,... 
%     'b','b',0,0.1);
% line([0 0],[0,100], 'Color','k','LineStyle','--')
% hold off
% xlabel('Time relative to FG peak (sec)')
% ylabel('')
% legend([ctx,hipp],'Ctx MUA', 'DG MUA')
% legend boxoff
% axis([-inf inf 0 80])
% set(gca,'FontSize',11)
%% Calculate Cortical activity before, during, and after
ctxBelowmeddian = zeros(size(pksfastgammaIdx,1),1);
ctxCentermeddian = zeros(size(pksfastgammaIdx,1),1);
ctxAbovmeddian = zeros(size(pksfastgammaIdx,1),1);
for kk = 1:size(pksfastgammaIdx,1)
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
% h2.BinEdges = 0:100;
h.BinWidth  = 2;
h1.BinWidth = 2;
h2.BinWidth = 2;
% h.EdgeColor = 'k';
xlabel('Ctx MUA envp (uV)')
ylabel('Probability Density')
ylim([0 .09])
legend([h, h1, h2], 'Ctx MUA Before DG fastgamma',...
    'Ctx MUA During DG fastgamma','Ctx MUA After DG fastgamma')
legend('boxoff')  
set(gca,'FontSize',11)
box off

%% Classification of fastgamma occurence during "up-states"
earlyUp = zeros(size(pksfastgammaIdx,1),1);
lateUp = zeros(size(pksfastgammaIdx,1),1);
betweenUp = zeros(size(pksfastgammaIdx,1),1);
upUp = zeros(size(pksfastgammaIdx,1),1);

for kk = 1: size(pksfastgammaIdx,1)
   if  ctxCentermeddian(kk) > 24 % if fastgamma occur during upstate
     if   ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) < 24
         betweenUp(kk) = pksfastgammaIdx(kk);
     elseif ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) > 24
         earlyUp(kk) = pksfastgammaIdx(kk);
     elseif ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) < 24
         lateUp(kk) = pksfastgammaIdx(kk);
     elseif ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) > 24
         upUp(kk) = pksfastgammaIdx(kk); 
     end
   end
end
 
earlyUp(earlyUp == 0) = [];
lateUp(lateUp == 0) = [];
betweenUp(betweenUp == 0) = [];
upUp(upUp == 0) = [];

%% Matrix calculation for fastgammas occuring during Up-states
eventExt = 0.1 / sint;
tt = -0.1:sint:0.1;
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
earlyUpMUAmat = zeros(size(earlyUp,1),Veclenght);
hippearlyUpMUAmat = zeros(size(earlyUp,1),Veclenght);
hippearlyUpfastgammamat = zeros(size(earlyUp,1),Veclenght);
for kk = 1: size(earlyUp,1)
    earlyUpMUAmat(kk,:) = ctxMUAenvpsm(earlyUp(kk) - eventExt : earlyUp(kk)+ eventExt);
    hippearlyUpMUAmat(kk,:) = DGspkDensity(earlyUp(kk) - eventExt : earlyUp(kk)+ eventExt);
    hippearlyUpfastgammamat(kk,:) = Zhippfastgammaleenv(earlyUp(kk) - eventExt : earlyUp(kk)+ eventExt);
end

lateUpMUAmat= zeros(size(lateUp,1),Veclenght);
hipplateUpMUAmat = zeros(size(lateUp,1),Veclenght);
hipplateUpfastgammamat = zeros(size(lateUp,1),Veclenght);
for kk = 1: size(lateUp,1)
    lateUpMUAmat(kk,:) = ctxMUAenvpsm(lateUp(kk) - eventExt : lateUp(kk)+ eventExt);
    hipplateUpMUAmat(kk,:) = DGspkDensity(lateUp(kk) - eventExt : lateUp(kk)+ eventExt);
    hipplateUpfastgammamat(kk,:) = Zhippfastgammaleenv(lateUp(kk) - eventExt : lateUp(kk)+ eventExt);
end

betweenUpMUAmat= zeros(size(betweenUp,1),Veclenght);
hippbetweenUpMUAmat = zeros(size(betweenUp,1),Veclenght);
hippbetweenUpfastgammamat = zeros(size(betweenUp,1),Veclenght);
for kk = 1: size(betweenUp,1)
    betweenUpMUAmat(kk,:) = ctxMUAenvpsm(betweenUp(kk) - eventExt : betweenUp(kk)+ eventExt);
    hippbetweenUpMUAmat(kk,:) = DGspkDensity(betweenUp(kk) - eventExt : betweenUp(kk)+ eventExt);
    hippbetweenUpfastgammamat(kk,:) = Zhippfastgammaleenv(betweenUp(kk) - eventExt : betweenUp(kk)+ eventExt);
end

upUpMUAmat= zeros(size(upUp,1),Veclenght);
hippupUpMUAmat = zeros(size(upUp,1),Veclenght);
hippupUpfastgammamat = zeros(size(upUp,1),Veclenght);
for kk = 1: size(upUp,1)
    upUpMUAmat(kk,:) = ctxMUAenvpsm(upUp(kk) - eventExt : upUp(kk)+ eventExt);
    hippupUpMUAmat(kk,:) = DGspkDensity(upUp(kk) - eventExt : upUp(kk)+ eventExt);
    hippupUpfastgammamat(kk,:) = Zhippfastgammaleenv(upUp(kk) - eventExt : upUp(kk)+ eventExt);
end

subplot(141)
trials = 1:length(earlyUp);
imagesc(tt,trials,earlyUpMUAmat)
c = colorbar;
axis xy
ylim([0 inf])
% xlabel('Time relative to FG peak (sec)');
ylabel('fastgamma event #')
caxis([20 100])
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
% xlabel('Time relative to FG peak (sec)');
% ylabel('fastgamma event #')
caxis([20 100])
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
% xlabel('Time relative to FG peak (sec)');
% ylabel('fastgamma event #')
caxis([20 100])
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
% xlabel('Time relative to FG peak (sec)');
% ylabel('fastgamma event #')
caxis([20 100])
colormap hot
colorbar('off')
set(gca,'FontSize',11)
%%
lateUpMUAmatctxMUAmean = mean(lateUpMUAmat,1);
lateUpMUAmatctxMUAerror = mad(lateUpMUAmat,1);
hipplateUpMUAmatmean = mean(hipplateUpMUAmat,1);
hipplateUpMUAmaterror = mad(hipplateUpMUAmat,1);
hipplateUpfastgammamatmean = mean(hipplateUpfastgammamat,1);
hipplateUpfastgammamaterror = mad(hipplateUpfastgammamat,1);

subplot(211)
hold on
ctx = plot(tt, lateUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,lateUpMUAmatctxMUAmean + lateUpMUAmatctxMUAerror,... 
    lateUpMUAmatctxMUAmean - lateUpMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
legend(ctx,'Ctx MUA')
ylabel('MUA envelope (uV)')
legend boxoff
set(gca,'FontSize',11)
hold off

subplot(212)
hold on
hipp = plot(tt, hipplateUpfastgammamatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hipplateUpfastgammamatmean + hipplateUpfastgammamaterror,... 
    hipplateUpfastgammamatmean - hipplateUpfastgammamaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,4], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to FG peak (sec)')
ylabel('z-score')
legend(hipp,'DG Fast-Gamma')
legend boxoff
axis([-inf inf 0 5])
set(gca,'FontSize',11)

%%
earlyUpMUAmatctxMUAmean = mean(earlyUpMUAmat,1);
earlyUpMUAmatctxMUAerror = mad(earlyUpMUAmat,1);
hippearlyUpMUAmatmean = mean(hippearlyUpMUAmat,1);
hippearlyUpMUAmaterror = mad(hippearlyUpMUAmat,1);
hippearlyUpfastgammamatmean = mean(hippearlyUpfastgammamat,1);
hippearlyUpfastgammamaterror = mad(hippearlyUpfastgammamat,1);

subplot(211)
hold on
ctx = plot(tt, earlyUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,earlyUpMUAmatctxMUAmean + earlyUpMUAmatctxMUAerror,... 
    earlyUpMUAmatctxMUAmean - earlyUpMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
legend(ctx,'Ctx MUA')
ylabel('MUA envelope (uV)')
legend boxoff
set(gca,'FontSize',11)
hold off

subplot(212)
hold on
hipp = plot(tt, hippearlyUpfastgammamatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippearlyUpfastgammamatmean + hippearlyUpfastgammamaterror,... 
    hippearlyUpfastgammamatmean - hippearlyUpfastgammamaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,4], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to FG peak (sec)')
ylabel('z-score')
legend(hipp,'DG Fast-Gamma')
legend boxoff
axis([-inf inf 0 5])
set(gca,'FontSize',11)

%%
betweenUpMUAmatctxMUAmean = mean(betweenUpMUAmat,1);
betweenUpMUAmatctxMUAerror = mad(betweenUpMUAmat,1);
hippbetweenUpMUAmatmean = mean(hippbetweenUpMUAmat,1);
hippbetweenUpMUAmaterror = mad(hippbetweenUpMUAmat,1);
hippbetweenUpfastgammamatmean = mean(hippbetweenUpfastgammamat,1);
hippbetweenUpfastgammamaterror = mad(hippbetweenUpfastgammamat,1);

subplot(211)
hold on
ctx = plot(tt, betweenUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,betweenUpMUAmatctxMUAmean + betweenUpMUAmatctxMUAerror,... 
    betweenUpMUAmatctxMUAmean - betweenUpMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
legend(ctx,'Ctx MUA')
ylabel('MUA envelope (uV)')
legend boxoff
set(gca,'FontSize',11)
hold off

subplot(212)
hold on
hipp = plot(tt, hippbetweenUpfastgammamatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippbetweenUpfastgammamatmean + hippbetweenUpfastgammamaterror,... 
    hippbetweenUpfastgammamatmean - hippbetweenUpfastgammamaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,4], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to FG peak (sec)')
ylabel('z-score')
legend(hipp,'DG Fast-Gamma')
legend boxoff
axis([-inf inf 0 5])
set(gca,'FontSize',11)

%%
upUpMUAmatctxMUAmean = mean(upUpMUAmat,1);
upUpMUAmatctxMUAerror = mad(upUpMUAmat,1);
hippupUpMUAmatmean = mean(hippupUpMUAmat,1);
hippupUpMUAmaterror = mad(hippupUpMUAmat,1);
hippupUpfastgammamatmean = mean(hippupUpfastgammamat,1);
hippupUpfastgammamaterror = mad(hippupUpfastgammamat,1);
subplot(211)
hold on
ctx = plot(tt,upUpMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,upUpMUAmatctxMUAmean +upUpMUAmatctxMUAerror,... 
   upUpMUAmatctxMUAmean -upUpMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
legend(ctx,'Ctx MUA')
ylabel('MUA envelope (uV)')
legend boxoff
set(gca,'FontSize',11)
hold off

subplot(212)
hold on
hipp = plot(tt, hippupUpfastgammamatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippupUpfastgammamatmean + hippupUpfastgammamaterror,... 
    hippupUpfastgammamatmean - hippupUpfastgammamaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,4], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to FG peak (sec)')
ylabel('z-score')
legend(hipp,'DG Fast-Gamma')
legend boxoff
axis([-inf inf 0 5])
set(gca,'FontSize',11)

%% Clasification of fastgamma occuridng during "Down-States"
earlyDown = zeros(size(pksfastgammaIdx,1),1);
lateDown = zeros(size(pksfastgammaIdx,1),1);
betweenDown = zeros(size(pksfastgammaIdx,1),1);
DownDown = zeros(size(pksfastgammaIdx,1),1);

for kk = 1: size(pksfastgammaIdx,1)
   if  ctxCentermeddian(kk) < 24 % if fastgamma occur during upstate
     if   ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) > 24
         betweenDown(kk) = pksfastgammaIdx(kk);
     elseif ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) < 24
         earlyDown(kk) = pksfastgammaIdx(kk);
     elseif ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) > 24
         lateDown(kk) = pksfastgammaIdx(kk);
     elseif ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) < 24
         DownDown(kk) = pksfastgammaIdx(kk); 
     end
   end
end
 
earlyDown(earlyDown == 0) = [];
lateDown(lateDown == 0) = [];
betweenDown(betweenDown == 0) = [];
DownDown(DownDown == 0) = [];

%% Matrix calculation for fastgammas occuring during Down-states
eventExt = 0.1 / sint;
tt = -0.1:sint:0.1;
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
earlyDownMUAmat = zeros(size(earlyDown,1),Veclenght);
hippearlyDownMUAmat = zeros(size(earlyDown,1),Veclenght);
hippearlyDownfastgammamat = zeros(size(earlyDown,1),Veclenght);
for kk = 1: size(earlyDown,1)
    earlyDownMUAmat(kk,:) = ctxMUAenvpsm(earlyDown(kk) - eventExt : earlyDown(kk)+ eventExt);
    hippearlyDownMUAmat(kk,:) = DGspkDensity(earlyDown(kk) - eventExt : earlyDown(kk)+ eventExt);
    hippearlyDownfastgammamat(kk,:) = Zhippfastgammaleenv(earlyDown(kk) - eventExt : earlyDown(kk)+ eventExt);
end

lateDownMUAmat= zeros(size(lateDown,1),Veclenght);
hipplateDownMUAmat = zeros(size(lateDown,1),Veclenght);
hipplateDownfastgammamat = zeros(size(lateDown,1),Veclenght);
for kk = 1: size(lateDown,1)
    lateDownMUAmat(kk,:) = ctxMUAenvpsm(lateDown(kk) - eventExt : lateDown(kk)+ eventExt);
    hipplateDownMUAmat(kk,:) = DGspkDensity(lateDown(kk) - eventExt : lateDown(kk)+ eventExt);
    hipplateDownfastgammamat(kk,:) = Zhippfastgammaleenv(lateDown(kk) - eventExt : lateDown(kk)+ eventExt);
end

betweenDownMUAmat= zeros(size(betweenDown,1),Veclenght);
hippbetweenDownMUAmat = zeros(size(betweenDown,1),Veclenght);
hippbetweenDownfastgammamat = zeros(size(betweenDown,1),Veclenght);
for kk = 1: size(betweenDown,1)
    betweenDownMUAmat(kk,:) = ctxMUAenvpsm(betweenDown(kk) - eventExt : betweenDown(kk)+ eventExt);
    hippbetweenDownMUAmat(kk,:) = DGspkDensity(betweenDown(kk) - eventExt : betweenDown(kk)+ eventExt);
    hippbetweenDownfastgammamat(kk,:) = Zhippfastgammaleenv(betweenDown(kk) - eventExt : betweenDown(kk)+ eventExt);
end

DownDownMUAmat= zeros(size(DownDown,1),Veclenght);
hippDownDownMUAmat = zeros(size(DownDown,1),Veclenght);
hippDownDownfastgammamat = zeros(size(DownDown,1),Veclenght);
for kk = 1: size(DownDown,1)
    DownDownMUAmat(kk,:) = ctxMUAenvpsm(DownDown(kk) - eventExt : DownDown(kk)+ eventExt);
    hippDownDownMUAmat(kk,:) = DGspkDensity(DownDown(kk) - eventExt : DownDown(kk)+ eventExt);
    hippDownDownfastgammamat(kk,:) = Zhippfastgammaleenv(DownDown(kk) - eventExt : DownDown(kk)+ eventExt);
end
subplot(141)
trials = 1:length(earlyDown);
imagesc(tt,trials,earlyDownMUAmat)
c = colorbar;
axis xy
ylim([0 inf])
% xlabel('Time relative to FG peak (sec)');
ylabel('fastgamma event #')
caxis([20 100])
colormap hot
colorbar('off')
set(gca,'FontSize',11)
ylim([1 inf])

subplot(143)
trials = 1:length(lateDown);
imagesc(tt,trials,lateDownMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
% xlabel('Time relative to FG peak (sec)');
% ylabel('fastgamma event #')
caxis([20 100])
colormap hot
colorbar('off')
set(gca,'FontSize',11)
ylim([1 inf])

subplot(142)
trials = 1:length(betweenDown);
imagesc(tt,trials,betweenDownMUAmat)
c = colorbar;
c.Label.String = 'MUA spk/sec';
axis xy
ylim([0 inf])
ylim([0 inf])
% xlabel('Time relative to FG peak (sec)');
% ylabel('fastgamma event #')
caxis([20 100])
colormap hot
colorbar('off')
set(gca,'FontSize',11)
ylim([1 inf])

subplot(144)
trials = 1:length(DownDown);
imagesc(tt,trials,DownDownMUAmat)
c = colorbar;
c.Label.String = 'MUA envelope (uV)';
c.FontSize = 11;
axis xy
ylim([1 inf])
% xlabel('Time relative to FG peak (sec)');
% ylabel('fastgamma event #')
caxis([20 100])
ylim([1 length(trials)])
colormap hot
colorbar('off')
set(gca,'FontSize',11)

%%
lateDownMUAmatctxMUAmean = mean(lateDownMUAmat,1);
lateDownMUAmatctxMUAerror = mad(lateDownMUAmat,1);
hipplateDownMUAmatmean = mean(hipplateDownMUAmat,1);
hipplateDownMUAmaterror = mad(hipplateDownMUAmat,1);
hipplateDownfastgammamatmean = mean(hipplateDownfastgammamat,1);
hipplateDownfastgammamaterror = mad(hipplateDownfastgammamat,1);

subplot(211)
hold on
ctx = plot(tt, lateDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,lateDownMUAmatctxMUAmean + lateDownMUAmatctxMUAerror,... 
    lateDownMUAmatctxMUAmean - lateDownMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
legend(ctx,'Ctx MUA')
ylabel('MUA envelope (uV)')
legend boxoff
set(gca,'FontSize',11)
hold off

subplot(212)
hold on
hipp = plot(tt, hipplateDownfastgammamatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hipplateDownfastgammamatmean + hipplateDownfastgammamaterror,... 
    hipplateDownfastgammamatmean - hipplateDownfastgammamaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,4], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to FG peak (sec)')
ylabel('z-score')
legend(hipp,'DG Fast-Gamma')
legend boxoff
axis([-inf inf 0 5])
set(gca,'FontSize',11)

%%
earlyDownMUAmatctxMUAmean = mean(earlyDownMUAmat,1);
earlyDownMUAmatctxMUAerror = mad(earlyDownMUAmat,1);
hippearlyDownMUAmatmean = mean(hippearlyDownMUAmat,1);
hippearlyDownMUAmaterror = mad(hippearlyDownMUAmat,1);
hippearlyDownfastgammamatmean = mean(hippearlyDownfastgammamat,1);
hippearlyDownfastgammamaterror = mad(hippearlyDownfastgammamat,1);

subplot(211)
hold on
ctx = plot(tt, earlyDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,earlyDownMUAmatctxMUAmean + earlyDownMUAmatctxMUAerror,... 
    earlyDownMUAmatctxMUAmean - earlyDownMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
legend(ctx,'Ctx MUA')
ylabel('MUA envelope (uV)')
legend boxoff
set(gca,'FontSize',11)
hold off

subplot(212)
hold on
hipp = plot(tt, hippearlyDownfastgammamatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippearlyDownfastgammamatmean + hippearlyDownfastgammamaterror,... 
    hippearlyDownfastgammamatmean - hippearlyDownfastgammamaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,4], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to FG peak (sec)')
ylabel('z-score')
legend(hipp,'DG Fast-Gamma')
legend boxoff
axis([-inf inf 0 5])
set(gca,'FontSize',11)
%%
betweenDownMUAmatctxMUAmean = mean(betweenDownMUAmat,1);
betweenDownMUAmatctxMUAerror = mad(betweenDownMUAmat,1);
hippbetweenDownMUAmatmean = mean(hippbetweenDownMUAmat,1);
hippbetweenDownMUAmaterror = mad(hippbetweenDownMUAmat,1);
hippbetweenDownfastgammamatmean = mean(hippbetweenDownfastgammamat,1);
hippbetweenDownfastgammamaterror = mad(hippbetweenDownfastgammamat,1);

subplot(211)
hold on
ctx = plot(tt, betweenDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,betweenDownMUAmatctxMUAmean + betweenDownMUAmatctxMUAerror,... 
    betweenDownMUAmatctxMUAmean - betweenDownMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
legend(ctx,'Ctx MUA')
ylabel('MUA envelope (uV)')
legend boxoff
set(gca,'FontSize',11)
hold off

subplot(212)
hold on
hipp = plot(tt, hippbetweenDownfastgammamatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippbetweenDownfastgammamatmean + hippbetweenDownfastgammamaterror,... 
    hippbetweenDownfastgammamatmean - hippbetweenDownfastgammamaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,4], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to FG peak (sec)')
ylabel('z-score')
legend(hipp,'DG Fast-Gamma')
legend boxoff
axis([-inf inf 0 5])
set(gca,'FontSize',11)

%%
DownDownMUAmatctxMUAmean = mean(DownDownMUAmat,1);
DownDownMUAmatctxMUAerror = mad(DownDownMUAmat,1);
hippDownDownMUAmatmean = mean(hippDownDownMUAmat,1);
hippDownDownMUAmaterror = mad(hippDownDownMUAmat,1);
hippDownDownfastgammamatmean = mean(hippDownDownfastgammamat,1);
hippDownDownfastgammamaterror = mad(hippDownDownfastgammamat,1);

subplot(211)
hold on
ctx = plot(tt, DownDownMUAmatctxMUAmean,'Color',[255/255 140/255 0/255], 'LineWidth',2);
jbfill(tt,DownDownMUAmatctxMUAmean + DownDownMUAmatctxMUAerror,... 
    DownDownMUAmatctxMUAmean - DownDownMUAmatctxMUAerror,... 
    [255/255 140/255 0/255],[255/255 140/255 0/255],0,0.2);
line([0 0],[0,100], 'Color','k','LineStyle','--')
legend(ctx,'Ctx MUA')
ylabel('MUA envelope (uV)')
legend boxoff
set(gca,'FontSize',11)
hold off

subplot(212)
hold on
hipp = plot(tt, hippDownDownfastgammamatmean, 'Color',[4/255 101/255 53/255], 'LineWidth',2);
jbfill(tt,hippDownDownfastgammamatmean + hippDownDownfastgammamaterror,... 
    hippDownDownfastgammamatmean - hippDownDownfastgammamaterror,... 
    [4/255 101/255 53/255],[4/255 101/255 53/255],0,0.2);
line([0 0],[0,4], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to FG peak (sec)')
ylabel('z-score')
legend(hipp,'DG Fast-Gamma')
legend boxoff
axis([-inf inf 0 5])
set(gca,'FontSize',11)

%% DG MUA during Up- and Down-States
% hippMUAUpSatemat = [hipplateUpMUAmatmean; hippearlyUpMUAmatmean; ...
%     hippbetweenUpMUAmatmean; hippupUpMUAmatmean];
% hippMUAUpSatemean = mean(hippMUAUpSatemat,1);
% hippMUAUpSateerror = mad(hippMUAUpSatemat,1);
% 
% hippMUADownSatemat = [hipplateDownMUAmatmean; hippearlyDownMUAmatmean; ...
%    hippbetweenDownMUAmatmean; hippDownDownMUAmatmean];
% hippMUADownSatemean = mean(hippMUADownSatemat,1);
% hippMUADownSateerror = mad(hippMUADownSatemat,1);
% 
% figure
% hold on
% hipMUAUpState = plot(tt,hippMUAUpSatemean,... 
%     'Color',[10/255 97/255 171/255], 'LineWidth',2);
% jbfill(tt,hippMUAUpSatemean + hippMUAUpSateerror,... 
%     hippMUAUpSatemean - hippMUAUpSateerror,... 
%     [10/255 97/255 171/255],[10/255 97/255 171/255],0,0.2);
% 
% hipMUADownState = plot(tt, hippMUADownSatemean,...
%     'Color',[213/255 29/255 91/255], 'LineWidth',2);
% jbfill(tt,hippMUADownSatemean + hippMUADownSateerror,... 
%     hippMUADownSatemean - hippMUADownSateerror,... 
%     [213/255 29/255 91/255],[213/255 29/255 91/255],0,0.2);
% line([0 0],[0,48], 'Color','k','LineStyle','--')
% hold off
% xlabel('Time relative to FG peak (sec)')
% ylabel('mean DG MUA (spk/sec) ')
% legend([hipMUAUpState,hipMUADownState],'DG MUA during Up-States', ...
%     'DG MUA during Down-States')
% legend boxoff
% axis([-inf inf 0 5])
% set(gca,'FontSize',11)

%% fastgamma power Up- and Down-States
hippfastgammaUpSatemat = [hipplateUpfastgammamatmean; hippearlyUpfastgammamatmean; ...
    hippbetweenUpfastgammamatmean; hippupUpfastgammamatmean];
hippfastgammaUpSatemean = mean(hippfastgammaUpSatemat,1);
hippfastgammaUpSateerror = mad(hippfastgammaUpSatemat,1);

hippfastgammaDownSatemat = [hipplateDownfastgammamatmean; hippearlyDownfastgammamatmean; ...
   hippbetweenDownfastgammamatmean; hippDownDownfastgammamatmean];
hippfastgammaDownSatemean = mean(hippfastgammaDownSatemat,1);
hippfastgammaDownSateerror = mad(hippfastgammaDownSatemat,1);

figure
hold on
hipfastgammaUpState = plot(tt,hippfastgammaUpSatemean,... 
    'Color',[10/255 97/255 171/255], 'LineWidth',2);
jbfill(tt,hippfastgammaUpSatemean + hippfastgammaUpSateerror,... 
    hippfastgammaUpSatemean - hippfastgammaUpSateerror,... 
    [10/255 97/255 171/255],[10/255 97/255 171/255],0,0.2);

hipfastgammaDownState = plot(tt, hippfastgammaDownSatemean,...
    'Color',[213/255 29/255 91/255], 'LineWidth',2);
jbfill(tt,hippfastgammaDownSatemean + hippfastgammaDownSateerror,... 
    hippfastgammaDownSatemean - hippfastgammaDownSateerror,... 
    [213/255 29/255 91/255],[213/255 29/255 91/255],0,0.2);
line([0 0],[0,8], 'Color','k','LineStyle','--')
hold off
xlabel('Time relative to FG peak (sec)')
ylabel('mean fastgamma(125-250Hz) Z-Score')
legend([hipfastgammaUpState,hipfastgammaDownState],'DG fastgamma during Up-States', ...
    'DG fastgamma during Down-States')
legend boxoff
axis([-inf inf -inf 5])
set(gca,'FontSize',11)

%%
totalRipdetected = size(lateUpMUAmat,1) + size(earlyUpMUAmat,1) ...
    + size(betweenUpMUAmat,1) + size(upUpMUAmat,1) + ...
    size(lateDownMUAmat,1) + size(earlyDownMUAmat,1) ...
    + size(betweenDownMUAmat,1) + size(DownDownMUAmat,1);
fastgammaUpStatePorcentage = (size(lateUpMUAmat,1) + size(earlyUpMUAmat,1) ...
    + size(betweenUpMUAmat,1) + size(upUpMUAmat,1))/totalRipdetected;
fastgammaDownStatePorcentage = (size(lateDownMUAmat,1) + size(earlyDownMUAmat,1) ...
    + size(betweenDownMUAmat,1) + size(DownDownMUAmat,1))/totalRipdetected;
%%
X = categorical({'"Up-State"','"Down-State"'});
X = reordercats(X,{'"Up-State"','"Down-State"'});
b = bar(X,[fastgammaUpStatePorcentage*100, fastgammaDownStatePorcentage*100]);
b.FaceColor = 'flat';
b.CData(2,:) = [213/255 29/255 91/255];
b.CData(1,:) = [10/255 1/255 171/255];
ylabel('DG fastgamma Occurrence %')
ylim([0 100])
set(gca,'FontSize',11)
box off
%%
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
b.FaceColor = 'flat';
ylabel('fastgamma Occurrence %')
ylim([0 70])
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
b.FaceColor = 'flat';
ylabel('fastgamma Occurrence %')
set(gca,'FontSize',11)
ylim([0 70])
box off











