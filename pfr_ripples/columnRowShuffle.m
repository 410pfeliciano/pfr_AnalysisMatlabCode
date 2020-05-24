Fs = 1000;
sint = 1/Fs;
%%
[ Y  ] = shuffle( ctxMUAreltoripp , 0 );

%% Cortical activity
CTXLFP = ch42Tet13LFP;
ctxLFPmualfp = bandFilter (CTXLFP, 'highfreq'); % 100-450Hz Filter
ctxMUAenvp = abs(hilbert(ctxLFPmualfp));
ctxMUAenvpsm = smoothdata(ctxMUAenvp,'gaussian',20); 
%%
eventExt = 0.1 / sint;
tt = -0.1:sint:0.1;
trials = 1:length(ctxMUAreltoripp);
figure('DefaultAxesFontSize',10)
subplot(121)
imagesc(tt,trials,ctxMUAreltoripp)
c = colorbar;
c.Label.String = 'MUA envelope (uV)';
c.FontSize = 11;
% colorbar('southoutside')
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
ylabel('Ripple #')
caxis([20 100])
colormap hot
set(gca,'FontSize',11)
% colorbar('off')
subplot(122)
imagesc(tt,trials,Y)
c = colorbar;
c.Label.String = 'Shuffle MUA envelope (uV)';
c.FontSize = 11;
% colorbar('southoutside')
axis xy
ylim([0 inf])
xlabel('Time relative to SWR peak (sec)');
% ylabel('Ripple event #')
caxis([20 100])
colormap hot
set(gca,'FontSize',11)
% colorbar('off')
%%
shuff = 1000;
shuffleIdxUpStatePorcentage = zeros(shuff,1); 
shuffleIdxDownStatePorcentage = zeros(shuff,1); 
earlyUpPorcentage = zeros(shuff,1); 
lateUpPorcentage = zeros(shuff,1); 
betweenUpPorcentage  = zeros(shuff,1); 
upUpPorcentage = zeros(shuff,1); 
earlyDownPorcentage = zeros(shuff,1); 
lateDownPorcentage = zeros(shuff,1); 
betweenDownPorcentage  = zeros(shuff,1); 
DownDownPorcentage = zeros(shuff,1); 
for xxx = 1:shuff
shuffleIdx = shuffle(ctxMUAreltoripp , 0 );
%
ctxBelowmeddian = zeros(size(shuffleIdx,1),1);
ctxCentermeddian = zeros(size(shuffleIdx,1),1);
ctxAbovmeddian = zeros(size(shuffleIdx,1),1);
for kk = 1:size(shuffleIdx,1)
    ctxBelowmeddian(kk) = median(shuffleIdx(kk,1:50));
    ctxCentermeddian(kk) = median(shuffleIdx(kk,75:125));
    ctxAbovmeddian(kk) = median(shuffleIdx(kk,150:201));
end
% separating up-states configurations based on corical MUA 
earlyUp = zeros(size(shuffleIdx,1),1);
lateUp = zeros(size(shuffleIdx,1),1);
betweenUp = zeros(size(shuffleIdx,1),1);
upUp = zeros(size(shuffleIdx,1),1);
for kk = 1: size(shuffleIdx,1)
   if  ctxCentermeddian(kk) > 24 % if shuffleIdx occur during upstate
     if   ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) < 24
         betweenUp(kk) = kk;
     elseif ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) > 24
         earlyUp(kk) = kk;
     elseif ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) < 24
         lateUp(kk) = kk;
     elseif ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) > 24
         upUp(kk) = kk; 
     end
   end
end
earlyUp(earlyUp == 0) = [];
lateUp(lateUp == 0) = [];
betweenUp(betweenUp == 0) = [];
upUp(upUp == 0) = [];

% Matrix creation for Up-States
eventExt = 0.1 / sint;
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
earlyUpMUAmat = zeros(size(earlyUp,1),Veclenght);
for kk = 1: size(earlyUp,1)
    earlyUpMUAmat(kk,:) = shuffleIdx(earlyUp(kk),:);
end

lateUpMUAmat= zeros(size(lateUp,1),Veclenght);
for kk = 1: size(lateUp,1)
    lateUpMUAmat(kk,:) = shuffleIdx(lateUp(kk),:);
end

betweenUpMUAmat= zeros(size(betweenUp,1),Veclenght);
for kk = 1: size(betweenUp,1)
    betweenUpMUAmat(kk,:) = shuffleIdx(betweenUp(kk),:);
end

upUpMUAmat= zeros(size(upUp,1),Veclenght);
for kk = 1: size(upUp,1)
    upUpMUAmat(kk,:) = shuffleIdx(upUp(kk),:);
end

% Down state clasification based on cortical mua
earlyDown = zeros(size(shuffleIdx,1),1);
lateDown = zeros(size(shuffleIdx,1),1);
betweenDown = zeros(size(shuffleIdx,1),1);
DownDown = zeros(size(shuffleIdx,1),1);
for kk = 1: size(shuffleIdx,1)
   if  ctxCentermeddian(kk) < 24 % if shuffleIdx occur during Downstate
     if   ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) > 24
         betweenDown(kk) = kk;
     elseif ctxBelowmeddian(kk) > 24 &&  ctxAbovmeddian(kk) < 24
         earlyDown(kk) = kk;
     elseif ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) > 24
         lateDown(kk) = kk;
     elseif ctxBelowmeddian(kk) < 24 &&  ctxAbovmeddian(kk) < 24
         DownDown(kk) = kk; 
     end
   end
end
earlyDown(earlyDown == 0) = [];
lateDown(lateDown == 0) = [];
betweenDown(betweenDown == 0) = [];
DownDown(DownDown == 0) = [];

% Matrix creation for Down-States
eventExt = 0.1 / sint;
Veclenght = length(1:eventExt) + length(1:eventExt)+1;
earlyDownMUAmat = zeros(size(earlyDown,1),Veclenght);
for kk = 1: size(earlyDown,1)
    earlyDownMUAmat(kk,:) = shuffleIdx(earlyDown(kk),:);
end

lateDownMUAmat= zeros(size(lateDown,1),Veclenght);
for kk = 1: size(lateDown,1)
    lateDownMUAmat(kk,:) = shuffleIdx(lateDown(kk),:);
end

betweenDownMUAmat= zeros(size(betweenDown,1),Veclenght);
for kk = 1: size(betweenDown,1)
    betweenDownMUAmat(kk,:) = shuffleIdx(betweenDown(kk),:);
end

DownDownMUAmat= zeros(size(DownDown,1),Veclenght);
for kk = 1: size(DownDown,1)
    DownDownMUAmat(kk,:) = shuffleIdx(DownDown(kk),:);
end

% Porcentage calculation
totalRipdetected = size(lateUpMUAmat,1) + size(earlyUpMUAmat,1) ...
    + size(betweenUpMUAmat,1) + size(upUpMUAmat,1) + ...
    size(lateDownMUAmat,1) + size(earlyDownMUAmat,1) ...
    + size(betweenDownMUAmat,1) + size(DownDownMUAmat,1);
%total up- and down-states
shuffleIdxUpStatePorcentage(xxx) = (size(lateUpMUAmat,1) + size(earlyUpMUAmat,1) ...
    + size(betweenUpMUAmat,1) + size(upUpMUAmat,1))/totalRipdetected*100;
shuffleIdxDownStatePorcentage(xxx) = (size(lateDownMUAmat,1) + size(earlyDownMUAmat,1) ...
    + size(betweenDownMUAmat,1) + size(DownDownMUAmat,1))/totalRipdetected*100;
%Up states
earlyUpPorcentage(xxx) = size(earlyUpMUAmat,1)/totalRipdetected*100;
lateUpPorcentage(xxx) = size(lateUpMUAmat,1)/totalRipdetected*100;
betweenUpPorcentage(xxx)  = size(betweenUpMUAmat,1)/totalRipdetected*100;
upUpPorcentage(xxx) = size(upUpMUAmat,1)/totalRipdetected*100;
%Down states
earlyDownPorcentage(xxx) = size(earlyDownMUAmat,1)/totalRipdetected*100;
lateDownPorcentage(xxx) = size(lateDownMUAmat,1)/totalRipdetected*100;
betweenDownPorcentage(xxx)  = size(betweenDownMUAmat,1)/totalRipdetected*100;
DownDownPorcentage(xxx) = size(DownDownMUAmat,1)/totalRipdetected*100;
end
%% figure
hold on
nbins = 10;
h= histogram(shuffleIdxUpStatePorcentage,nbins,...
    'Normalization','probability','DisplayStyle','bar',...
    'EdgeColor','none','FaceColor','b');
h1= histogram(shuffleIdxDownStatePorcentage,nbins,...
    'Normalization','probability','DisplayStyle','bar',...
    'EdgeColor','none','FaceColor','r');
hold off
xlabel('Shuffle Up- and Down-States %')
ylabel('Probability Density')
ylim([0 inf])
xlim([0 100])
legend([h, h1], 'Up-States Shuffle','Down-State Shuffle')
legend('boxoff')  
set(gca,'FontSize',11)
box off
% pd = fitdist(shuffleIdxUpStatePorcentage,'Normal')
% x_values = shuffleIdxUpStatePorcentage;
% y = pdf(pd,x_values);
% plot(x_values,y,'.')
% ci = paramci(pd)
%%
figure
hold on
nbins = 10;
h= histogram(earlyUpPorcentage,nbins,...
    'Normalization','probability','DisplayStyle','bar',...
    'EdgeColor','none','FaceColor','r');
h1= histogram(lateUpPorcentage,nbins,...
    'Normalization','probability','DisplayStyle','bar',...
    'EdgeColor','none','FaceColor','b');
h2= histogram(betweenUpPorcentage,nbins,...
    'Normalization','probability','DisplayStyle','bar',...
    'EdgeColor','none','FaceColor','k');
h3= histogram(upUpPorcentage,nbins,...
    'Normalization','probability','DisplayStyle','bar',...
    'EdgeColor','none','FaceColor','g');
hold off
xlabel('Shuffle Up-States %')
ylabel('Probability Density')
ylim([0 inf])
xlim([0 100])
legend([h, h1, h2 h3], 'earlyUp','lateUp','between','up-state')
legend('boxoff')  
set(gca,'FontSize',11)
box off
%% figure
hold on
nbins = 10;
h= histogram(earlyDownPorcentage,nbins,...
    'Normalization','probability','DisplayStyle','bar',...
    'EdgeColor','none','FaceColor','r');
h1= histogram(lateDownPorcentage,nbins,...
    'Normalization','probability','DisplayStyle','bar',...
    'EdgeColor','none','FaceColor','b');
h2= histogram(betweenDownPorcentage,nbins,...
    'Normalization','probability','DisplayStyle','bar',...
    'EdgeColor','none','FaceColor','k');
h3= histogram(DownDownPorcentage,nbins,...
    'Normalization','probability','DisplayStyle','bar',...
    'EdgeColor','none','FaceColor','g');
hold off
xlabel('Shuffle Down-States %')
ylabel('Probability Density')
ylim([0 inf])
xlim([0 100])
legend([h, h1, h2 h3], 'earlyDown','lateDown','between','Down-state')
legend('boxoff')  
set(gca,'FontSize',11)
box off
