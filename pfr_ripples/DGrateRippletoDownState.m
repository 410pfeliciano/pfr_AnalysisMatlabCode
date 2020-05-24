%% Matrix Presleep
Fs = 1000;
sint = 1/Fs;
tbins = time;
Gammats= DGpksGammaIdx*sint; % Sleep gamma idx converted to time stamps
%% Shuffle data to test random values
% startIdx = 1;
% endIdx= 1609/sint;
% presleepGammaIdx = sort(round((endIdx-startIdx).*rand(1000,1) + startIdx));
%% Pre sleep Before Down-State start
GammaTrain = hist(Gammats,tbins); % binning down-state time stamps
preBelowDSidx = ctxDownPreSlBelowIdx;
preAbovDSidx = ctxDownPreSlAbovIdx;
matTs = 0.2;% time to analyzed
preBelowtt = 0:-sint:-(matTs);
if preBelowDSidx(1)-(matTs/sint+1) < 1
    preBelowDSidx(1) = [];
end

%% Pre Below Down-States and Gamma activity
preBelowGammaMat = zeros(length(preBelowDSidx),length(preBelowtt));
for kk = 1:length(preBelowDSidx)
     preBelowGammaMat(kk,1) = GammaTrain(preBelowDSidx(kk,1));
    for ll = 2: size(preBelowGammaMat,2)
        preBelowGammaMat(kk,ll) = GammaTrain(1,preBelowDSidx(kk)-ll);
    end
end 
preBelowGamNum = size(find(preBelowGammaMat),1);
close all
figure
subplot(211)
plot(preBelowtt,preBelowGammaMat(1,:),'Color',[255/255  128/255 0/255])
hold on
for kk = 2:length(preBelowDSidx)
    plot(preBelowtt,preBelowGammaMat(kk,:)*kk,'.',...
        'Color',[255/255  128/255 0/255])
end
hold off
axis([-inf inf 1 inf])
ylabel('Down-State #')
set(gca,'FontSize',11)
box off
subplot(212)
preBelowGammaRate = sum(preBelowGammaMat)./length(preBelowDSidx)./sint;
% preBelowGammaRatehist = hist(preBelowGammaRate,100);
bar(preBelowtt,preBelowGammaRate,'FaceColor',[255/255  128/255 0/255],...
    'EdgeColor',[255/255  128/255 0/255])

set(gca,'FontSize',11)
hold on
startPoints = [0.5 -0.08 0.05 0.25];
% exclude = preBelowtt < -.06;
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
f = fit(preBelowtt',preBelowGammaRate',gaussEqn,'Start', startPoints)
p = plot(f,'r');
p.LineWidth = 2;
preBelowGammaRatesmooth = smoothdata(preBelowGammaRate, 'gaussian',50);
plot(preBelowtt,preBelowGammaRatesmooth, 'b', 'LineWidth',2)
ylabel('FG rate FGamma/sec')
xlabel('time rel. to Down-State(sec)')
ylim([0 2])
xlim([-inf 0])
box off
legend off
hold off


%% Pre Sleep Abov Down-State and Gamma Interaction 

preAbovtt = 0:sint:matTs;
preAbovGammaMat = zeros(length(preAbovDSidx),length(preAbovtt));
for kk = 1:length(preAbovDSidx)
     preAbovGammaMat(kk,1) = GammaTrain(preAbovDSidx(kk,1));
    for ll = 2: size(preAbovGammaMat,2)
        preAbovGammaMat(kk,ll) = GammaTrain(1,preAbovDSidx(kk)+ll);
    end
end 
preAbovGamNum = size(find(preAbovGammaMat),1);
close all
figure
subplot(211)
plot(preAbovtt,preAbovGammaMat(1,:),'Color',[255/255  128/255 0/255])
hold on
for kk = 2:length(preAbovDSidx)

    plot(preAbovtt,preAbovGammaMat(kk,:)*kk,'.', 'Color',[255/255  128/255 0/255])

end
hold off
axis([-inf inf 1 inf])

ylabel('Down-State #')
set(gca,'FontSize',11)
box off
subplot(212)
preAbovGammaRate  = sum(preAbovGammaMat)./length(preAbovDSidx)./sint;
bar(preAbovtt,preAbovGammaRate ,'FaceColor',[255/255  128/255 0/255],...
    'EdgeColor',[255/255  128/255 0/255])
set(gca,'FontSize',11)
hold on

startPoints = [0.5 .05 0.05 0.25];
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
f = fit(preAbovtt',preAbovGammaRate',gaussEqn,'Start', startPoints)
p = plot(f,'r');
p.LineWidth = 2;
preAbovGammaRatesmooth = smoothdata(preAbovGammaRate, 'gaussian',50);
plot(preAbovtt,preAbovGammaRatesmooth, 'b', 'LineWidth',2)
ylabel('FG rate fGamma/sec')
xlabel('time rel. to Down-State(sec)')
ylim([0 2])
xlim([0 inf])
box off
legend off
hold off















%% Post Above Cortico-Gamma activity
% Matrix Post sleep
matTs = 0.2;% time to analyzed
%% Post Sleep After Down-State
posttt = 0:sint:matTs;
posBelowDownStateIdx = ctxDownPostSlBelowIdx;
posAbovDownStateIdx = ctxDownPostSlAbovIdx;
if posAbovDownStateIdx(end)+200 > time(end)
    posAbovDownStateIdx(end) = [];
end
posAbovGammaMat= zeros(length(posAbovDownStateIdx),length(posttt));
for kk = 1:length(posAbovDownStateIdx)
     posAbovGammaMat(kk,1) = GammaTrain(posAbovDownStateIdx(kk,1));
    for ll = 2: size(posAbovGammaMat,2)
        posAbovGammaMat(kk,ll) = GammaTrain(1,posAbovDownStateIdx(kk)+ll);
    end
end 
posAbovGamNum = size(find(posAbovGammaMat),1);
close all
posAbovtt = 0:sint:matTs;
figure
subplot(211)
plot(posAbovtt,posAbovGammaMat(1,:),'Color',[255/255  128/255 0/255])
hold on
for kk = 2:length(posAbovDownStateIdx)

    plot(posAbovtt,posAbovGammaMat(kk,:)*kk,'.', 'Color',[255/255  128/255 0/255])

end
hold off
axis([-inf inf 1 inf])

ylabel('Down-State #')
set(gca,'FontSize',11)
box off
subplot(212)
posAbovGammaRate = sum(posAbovGammaMat)./length(posAbovDownStateIdx)./sint;
bar(posAbovtt,posAbovGammaRate,'FaceColor',[255/255  128/255 0/255],...
    'EdgeColor',[255/255  128/255 0/255])

set(gca,'FontSize',11)
hold on
startPoints = [0.5 .05 0.05 0.25];
% exclude = posAbovtt < -.06;
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
f = fit(posAbovtt',posAbovGammaRate',gaussEqn,'Start', startPoints)
p = plot(f,'r');
p.LineWidth = 2;
posAbovGammaRatesmooth = smoothdata(posAbovGammaRate, 'gaussian',50);
plot(preAbovtt,posAbovGammaRatesmooth, 'b', 'LineWidth',2)
ylabel('FG rate fGamma/sec')
xlabel('time rel. to Down-State(sec)')
ylim([0 2])
xlim([0 inf])
box off
legend off
hold off


%% Post Sleep Before(Below) Start Down State
posBelowGammaMat= zeros(length(posBelowDownStateIdx),length(posttt));
for kk = 1:length(posBelowDownStateIdx)
     posBelowGammaMat(kk,1) = GammaTrain(posBelowDownStateIdx(kk,1));
    for ll = 2: size(posBelowGammaMat,2)
        posBelowGammaMat(kk,ll) = GammaTrain(1,posBelowDownStateIdx(kk)-ll);
    end
end 
close all
posBelowGamNum = size(find(posBelowGammaMat),1);
posBelowtt = 0:-sint:-(matTs);
figure
subplot(211)
plot(posBelowtt,posBelowGammaMat(1,:),'Color',[255/255  128/255 0/255])
hold on
for kk = 2:length(posBelowDownStateIdx)

    plot(posBelowtt,posBelowGammaMat(kk,:)*kk,'.', 'Color',[255/255  128/255 0/255])

end
hold off
axis([-inf inf 1 inf])

ylabel('Down-State #')
set(gca,'FontSize',11)
box off
subplot(212)
posBelowGammaRate = sum(posBelowGammaMat)./length(posBelowDownStateIdx)./sint;
bar(posBelowtt,posBelowGammaRate,'FaceColor',[255/255  128/255 0/255],...
    'EdgeColor',[255/255  128/255 0/255])
set(gca,'FontSize',11)
hold on
startPoints = [0.5 -0.1 0.05 0.25];
% exclude = posBelowtt < -.06;
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
f = fit(posBelowtt',posBelowGammaRate',gaussEqn,'Start', startPoints)
p = plot(f,'r');
p.LineWidth = 2;
posBelowGammaRatesmooth = smoothdata(posBelowGammaRate, 'gaussian',50);
plot(posBelowtt,posBelowGammaRatesmooth, 'b', 'LineWidth',2)
ylabel('FG rate fGamma/sec')
xlabel('time rel. to Down-State(sec)')
ylim([0 2])
xlim([-inf 0])
box off
legend off
hold off

%% Putting all Down states and Gamma toguether
BelowGammaMat = [preBelowGammaMat; posBelowGammaMat];
ttmat = 0:-sint:-(matTs);
close all
figure
subplot(211)
plot(ttmat,BelowGammaMat(1,:),'Color',[255/255  128/255 0/255])
hold on
for kk = 2:length(BelowGammaMat)
    plot(ttmat,BelowGammaMat(kk,:)*kk,'.', 'Color',[255/255  128/255 0/255])
end
hold off
axis([-inf inf 10 inf])
ylabel('Down-State #')
set(gca,'FontSize',11)
box off
subplot(212)
BelowGammaRate = sum(BelowGammaMat)./length(BelowGammaMat)./sint;
bar(ttmat,BelowGammaRate,'FaceColor',[255/255  128/255 0/255],...
    'EdgeColor',[255/255  128/255 0/255])
ylim([2 inf])
set(gca,'FontSize',11)
hold on
startPoints = [0.25 -0.05 0.05 0.1];
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
f = fit(ttmat',BelowGammaRate',gaussEqn,'Start', startPoints)
p = plot(f,'r');
p.LineWidth = 2;
BelowGammaRatesmooth = smoothdata(BelowGammaRate, 'gaussian',50);
plot(ttmat,BelowGammaRatesmooth, 'b', 'LineWidth',2)
ylabel('FGamma rate Gamma/sec')
xlabel('time rel. to Down-State(sec)')
ylim([-inf 2.5])
xlim([-inf 0])
box off
legend off
hold off
%%
AbovGammaMat = [preAbovGammaMat; posAbovGammaMat];
ttmat = 0:sint:(matTs);
close all
figure
subplot(211)
plot(ttmat,AbovGammaMat(1,:),'Color',[255/255  128/255 0/255])
hold on
for kk = 2:length(AbovGammaMat)

    plot(ttmat,AbovGammaMat(kk,:)*kk,'.', 'Color',[255/255  128/255 0/255])

end
hold off
axis([-inf inf 10 inf])
ylabel('Down-State #')
set(gca,'FontSize',11)
box off
subplot(212)
abovGammaRate = sum(AbovGammaMat)./length(AbovGammaMat)./sint;
bar(ttmat,abovGammaRate,'FaceColor',[255/255  128/255 0/255],...
    'EdgeColor',[255/255  128/255 0/255])

set(gca,'FontSize',11)
hold on
startPoints = [0.5 0.1 0.1 0.1];
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
f = fit(ttmat',abovGammaRate',gaussEqn,'Start', startPoints)
p = plot(f,'r');
p.LineWidth = 2;
abovGammaRatesmooth = smoothdata(abovGammaRate, 'gaussian',50);
plot(ttmat,abovGammaRatesmooth, 'b', 'LineWidth',2)
ylabel('FGamma rate Gamma/sec')
xlabel('time rel. to Down-State(sec)')
ylim([-inf 2.5])
xlim([0 inf])
box off
legend off
hold off
