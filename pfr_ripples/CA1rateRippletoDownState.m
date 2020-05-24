% Matrix Presleep
Fs = 1000;
sint = 1/Fs;
ttbins = time;
sleepRipplets= CA1pksRippleIdx*sint;
rippletrain = hist(sleepRipplets,ttbins);
preBelowDSidx = ctxDownPreSlBelowIdx;
preAbovDSidx = ctxDownPreSlAbovIdx;
%% Pre Below Ripple and Cortico Interaction
matTs = 0.2;% time to analyzed
tt = -(matTs):sint:0;
if preBelowDSidx(1)-length(tt) < ttbins(1)
    preBelowDSidx(1) = [];
end
preBelowRippleMat = zeros(length(preBelowDSidx),length(tt));
for kk = 1:length(preBelowDSidx)
     preBelowRippleMat(kk,1) = rippletrain(preBelowDSidx(kk,1));
    for ll = 2: size(preBelowRippleMat,2)
        preBelowRippleMat(kk,ll) = rippletrain(1,preBelowDSidx(kk)-ll);
    end
end 
%%
ttmat = 0:-sint:-(matTs);
figure
subplot(211)
plot(ttmat,preBelowRippleMat(1,:),'Color',[128/255  128/255 128/255])
hold on
for kk = 2:length(preBelowDSidx)

    plot(ttmat,preBelowRippleMat(kk,:)*kk,'.', 'Color',[128/255  128/255 128/255])

end
preBelowRippNum = size(find(preBelowRippleMat),1);
hold off
axis([-inf inf 1 inf])
ylabel('Down-State #')
set(gca,'FontSize',11)
box off
subplot(212)
preBelowRippleRate = sum(preBelowRippleMat)./length(preBelowDSidx)./sint;
bar(ttmat,preBelowRippleRate,'FaceColor',[128/255  128/255 128/255],...
    'EdgeColor',[128/255  128/255 128/255])
set(gca,'FontSize',11)
hold on
startPoints = [1.2 -0.03 0.01 0.35];
% exclude = ttmat < -.06;
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
f = fit(ttmat',preBelowRippleRate',gaussEqn,'Start', startPoints)
p = plot(f,'r');
p.LineWidth = 2;
preBelowRippRatesmooth = smoothdata(preBelowRippleRate, 'gaussian',50);
plot(ttmat,preBelowRippRatesmooth, 'b', 'LineWidth',2)
ylabel('SWR rate ripple/sec')
xlabel('time rel. to Down-State(sec)')
ylim([0 3])
xlim([-inf 0])
box off
legend off
hold off
%% Pre Above Ripple Cortico Interaction
Prettmat = 0:sint:matTs;
preAbovRippleMat = zeros(length(preAbovDSidx),length(Prettmat));
for kk = 1:length(preAbovDSidx)
     preAbovRippleMat(kk,1) = rippletrain(preAbovDSidx(kk,1));
    for ll = 2: size(preAbovRippleMat,2)
        preAbovRippleMat(kk,ll) = rippletrain(1,preAbovDSidx(kk)+ll);
    end
end 
preAbovRippNum = size(find(preAbovRippleMat),1);
close all
figure
subplot(211)
plot(Prettmat,preAbovRippleMat(1,:),'Color',[128/255  128/255 128/255])
hold on
for kk = 2:length(preAbovDSidx)

    plot(Prettmat,preAbovRippleMat(kk,:)*kk,'.', 'Color',[128/255  128/255 128/255])

end
hold off
axis([-inf inf 1 inf])

ylabel('Down-State #')
set(gca,'FontSize',11)
box off
subplot(212)
preAbovRippleRate = sum(preAbovRippleMat)./length(preAbovDSidx)./sint;
bar(Prettmat,preAbovRippleRate,'FaceColor',[128/255  128/255 128/255],...
    'EdgeColor',[128/255  128/255 128/255])

set(gca,'FontSize',11)
hold on
startPoints = [0.5 .1 .25 0.025];
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
f = fit(Prettmat',preAbovRippleRate',gaussEqn,'Start', startPoints)
p = plot(f,'r');
p.LineWidth = 2;
preAbovRippRatesmooth = smoothdata(preAbovRippleRate, 'gaussian',50);
plot(Prettmat,preAbovRippRatesmooth, 'b', 'LineWidth',2)
ylabel('SWR rate ripple/sec')
xlabel('time rel. to Down-State(sec)')
ylim([0 3])
xlim([0 inf])
box off
legend off
hold off







%%
posBelowDownStateIdx = ctxDownPostSlBelowIdx;
posAbovDownStateIdx = ctxDownPostSlAbovIdx;
%%
matTs = 0.2;% time to analyzed
postt = -(matTs):sint:0;
if posBelowDownStateIdx(1)-length(postt) < ttbins(1)
    posBelowDownStateIdx(1) = [];
end
posBelowRippleMat = zeros(length(posBelowDownStateIdx),length(postt));
for kk = 1:length(posBelowDownStateIdx)
     posBelowRippleMat (kk,1) = rippletrain(posBelowDownStateIdx(kk,1));
    for ll = 2: size(posBelowRippleMat ,2)
        posBelowRippleMat (kk,ll) = rippletrain(1,posBelowDownStateIdx(kk)-ll);
    end
end 
posBelowRippNum = size(find(posBelowRippleMat),1);
close all
ttmat = 0:-sint:-(matTs);
figure
subplot(211)
plot(ttmat,posBelowRippleMat (1,:),'.','Color',[128/255  128/255 128/255])
hold on
for kk = 2:length(posBelowDownStateIdx)

    plot(ttmat,posBelowRippleMat (kk,:)*kk,'.','Color',[128/255  128/255 128/255])

end
hold off
axis([-inf inf 1 inf])
ylabel('Down-State #')
set(gca,'FontSize',11)
box off
subplot(212)
posBelowRippleRate = sum(posBelowRippleMat )./length(posBelowDownStateIdx)./sint;
bar(ttmat,posBelowRippleRate, 'FaceColor',[128/255  128/255 128/255],...
    'EdgeColor',[128/255  128/255 128/255])
hold on
startPoints = [1 -0.05 0.02 0.2];
% exclude = ttmat < -.06;
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
f = fit(ttmat',posBelowRippleRate',gaussEqn,'Start', startPoints)
p = plot(f,'r');
p.LineWidth = 2;
posBelowRippRatesmooth = smoothdata(posBelowRippleRate, 'gaussian',50);
plot(ttmat,posBelowRippRatesmooth, 'b', 'LineWidth',2)
ylabel('SWR rate ripple/sec')
xlabel('time rel. to Down-State(sec)')
ylim([0 3])
xlim([-inf 0])
set(gca,'FontSize',11)
box off
legend off
hold off


%%
posAbovRippleMat = zeros(length(posAbovDownStateIdx),length(postt));
for kk = 1:length(posAbovDownStateIdx)
     posAbovRippleMat(kk,1) = rippletrain(posAbovDownStateIdx(kk,1));
    for ll = 2: size(posAbovRippleMat,2)
        posAbovRippleMat(kk,ll) = rippletrain(1,posAbovDownStateIdx(kk)+ll);
    end
end 
posAbovRippNum = size(find(posAbovRippleMat),1);
close all
posttmat = 0:sint:matTs;
figure
subplot(211)
plot(posttmat,posAbovRippleMat(1,:),'Color',[128/255  128/255 128/255])
hold on
for kk = 2:length(posAbovDownStateIdx)

    plot(posttmat,posAbovRippleMat(kk,:)*kk,'.', 'Color',[128/255  128/255 128/255])

end
hold off
axis([-inf inf 1 inf])
ylabel('Down-State #')
set(gca,'FontSize',11)
box off
subplot(212)
posAbovRippleRate= sum(posAbovRippleMat)./length(posAbovDownStateIdx)./sint;
bar(posttmat,posAbovRippleRate,'FaceColor',[128/255  128/255 128/255],...
    'EdgeColor',[128/255  128/255 128/255])

set(gca,'FontSize',11)
hold on
% startPoints = [0.5 .1 0.1 0.25];
% % exclude = posttmat < -.06;
% gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
% f = fit(posttmat',posAbovRippleRate',gaussEqn,'Start', startPoints)
% p = plot(f,'r');
% p.LineWidth = 2;
posAbovRippRatesmooth = smoothdata(posAbovRippleRate, 'gaussian',50);
plot(posttmat,posAbovRippRatesmooth, 'r', 'LineWidth',2)
ylabel('SWR rate ripple/sec')
xlabel('time rel. to Down-State(sec)')
ylim([0 3])
xlim([0 0.2])
box off
legend off
hold off

%% Putting all Down states and Ripple toguther
BelowRippleMat = [preBelowRippleMat; posBelowRippleMat];
ttmat = 0:-sint:-(matTs);
figure
subplot(211)
plot(ttmat,BelowRippleMat(1,:),'Color',[128/255  128/255 128/255])
hold on
for kk = 2:length(BelowRippleMat)

    plot(ttmat,BelowRippleMat(kk,:)*kk,'.', 'Color',[128/255  128/255 128/255])

end
hold off
axis([-inf inf 10 inf])
ylabel('Down-State #')
set(gca,'FontSize',11)
box off
subplot(212)
BelowRippleRate = sum(BelowRippleMat)./length(BelowRippleMat)./sint;
bar(ttmat,BelowRippleRate,'FaceColor',[128/255  128/255 128/255],...
    'EdgeColor',[128/255  128/255 128/255])

set(gca,'FontSize',11)
hold on
% startPoints = [1.2 -0.02 0.1 0.01];
% % exclude = ttmat < -.06;
% gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
% f = fit(ttmat',BelowRippleRate',gaussEqn,'Start', startPoints)
% p = plot(f,'r');
% p.LineWidth = 2;
BelowRippRatesmooth = smoothdata(BelowRippleRate, 'gaussian',50);
plot(ttmat,BelowRippRatesmooth, 'r', 'LineWidth',2)
ylabel('SWR rate ripple/sec')
xlabel('time rel. to Down-State(sec)')
ylim([-inf 3])
xlim([-inf 0])
box off
legend off
hold off
%%
AbovRippleMat = [preAbovRippleMat; posAbovRippleMat];
ttmat = 0:sint:(matTs);
figure
subplot(211)
plot(ttmat,AbovRippleMat(1,:),'Color',[128/255  128/255 128/255])
hold on
for kk = 2:length(AbovRippleMat)

    plot(ttmat,AbovRippleMat(kk,:)*kk,'.', 'Color',[128/255  128/255 128/255])

end
hold off
axis([-inf inf 10 inf])
ylabel('Down-State #')
set(gca,'FontSize',11)
box off
subplot(212)
abovRippleRate = sum(AbovRippleMat)./length(AbovRippleMat)./sint;
bar(ttmat,abovRippleRate,'FaceColor',[128/255  128/255 128/255],...
    'EdgeColor',[128/255  128/255 128/255])
set(gca,'FontSize',11)
hold on
% startPoints = [0.5 0.05 0.2 0.1];
% % exclude = ttmat < -.06;
% gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
% f = fit(ttmat',abovRippleRate',gaussEqn,'Start', startPoints)
% p = plot(f,'r');
% p.LineWidth = 2;
abovRippRatesmooth = smoothdata(abovRippleRate, 'gaussian',50);
plot(ttmat,abovRippRatesmooth, 'r', 'LineWidth',2)
ylabel('SWR rate ripple/sec')
xlabel('time rel. to Down-State(sec)')
ylim([-inf 2.5])
xlim([0 inf])
box off
legend off
hold off
