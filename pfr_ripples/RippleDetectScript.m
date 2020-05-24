% LFP of Interest
hippLFP = ch9Tet7LFP;
rippleLFP = ch9Tet7LFPripple; % filtered LFP detrend
%% Hippocampal Spike Density

%%
bandFilter (hippLFP, 'sharpwave');

%%
sharpwaveLFP = ch9Tet7LFPsharpwave;
Fs = 1000;
sint = 1/Fs;
%% Second Step, Calculation the mean power and standar deviation of the recording
filtWin = 0.02; % filter window in in seconds
hibfiltLFP = abs(hilbert(rippleLFP)); % amplitude envelope of ripple-band filtered LFP
sharpLFPenv = abs(hilbert(sharpwaveLFP));
sqfiltLFP= hibfiltLFP.^2; % calculating the power of the filtered LFP
sqsharpwaveEnvp = sharpLFPenv.^2;
smsqfiltLFP = smoothdata(sqfiltLFP,'gaussian',filtWin/sint); % smoothing the power data
smoothsharpwaveenvp = smoothdata(sqsharpwaveEnvp,'gaussian',filtWin/sint);
%% Plotting LFP and MUA envelop
time1= 1595;
time2= 1602;
xt1= time1/sint;
xt2= time2/sint;

subplot(5,1,1)
plot(time(xt1:xt2),hippLFP(xt1:xt2),'k') 
axis([-inf inf -inf inf])
legend('Ctx','Location', 'best')
legend boxoff
box off
set(gca,'xtickLabel',[])
ylabel('uV')

subplot(5,1,2)
plot(time(xt1:xt2),rippleLFP(xt1:xt2),'r') 
axis tight
legend('filt. ripple','Location', 'best')
legend boxoff
box off
set(gca,'xtickLabel',[])
ylabel('uV')

subplot(5,1,3)
plot(time(xt1:xt2),sharpwaveLFP(xt1:xt2),'b') 
axis tight
legend('Sharpwave','Location', 'best')
legend boxoff
box off
% set(gca,'xtickLabel',[])
ylabel('uV')

subplot(5,1,4)
plot(time(xt1:xt2),smsqfiltLFP(xt1:xt2),'r') 
axis tight
legend('Smooth Ripple Envp.','Location', 'best')
legend boxoff
box off
% set(gca,'xtickLabel',[])
ylabel('uV')

subplot(5,1,5)
plot(time(xt1:xt2),smoothsharpwaveenvp(xt1:xt2),'b') 
axis tight
legend('Smooth Sharpwave Envp','Location', 'best')
legend boxoff
box off
% set(gca,'xtickLabel',[])
ylabel('uV')
%% Distribution of MUA 
time1= 1959;
time2= 6433;
xt1= time1/sint;
xt2= time2/sint;
h= histogram(smsqfiltLFP(xt1:xt2),200,'Normalization','probability',... 
    'DisplayStyle','stairs', 'LineWidth',1);
h.BinEdges = 0:6000;
h.BinWidth = 10;
h.EdgeColor = 'k';
xlabel('Ripple Power Density')
ylabel('Probability Density')
legend('Ripple Power')
legend('boxoff')  
box off

%%  FIRST CRITERION: AT LEAST 3*STD OF FILTERED, SQUARED LFP
time1= 4959;
time2= 6433;
% time1= 1000;
% time2= 1600;
xt1= time1/sint;
xt2= time2/sint;
thMad = mad(smsqfiltLFP(xt1:xt2),0); % median deviation flag = 1 mean deviation flag =0
th = median(smsqfiltLFP(xt1:xt2)) + (3*thMad);
% thMad2 = mad(smoothsharpwaveenvp(xt1:xt2),1); % median deviation flag = 1
% th2 = median(smoothsharpwaveenvp(xt1:xt2)) + (5*thMad2);
splog = smsqfiltLFP > th;
pos = find(splog == 1); % This vectors contain the indices of power larger the 3std


%% SECOND CRITERION: at least a difference of 15ms between the end of
...a ripple and the beginning of the next to be considered separate events

    difpos = find(diff(pos) >= (0.02/sint)); % events with inte4rvals > 20ms
    % end of all possible ripple evnts
    z = zeros(length(difpos)+1,1); % Creating a zero row vector
    z(1:end-1) = difpos;
    z(end) = length(pos); % End position of all potential ripples
    % Calculate the start
    st = difpos+1; % Start positions of all possible ripple events
    zst = zeros(length(difpos)+1,1);
    zst(2:end) = st; %
    zst(1)=1; % zst start position of all
    stepripples = [pos(zst), pos(z)]; % Start and End position(indices) 
        ... of Ripples Candidates


%% THIRD CRITERION: Eliminate 'Ripples' candidates with < sampling interval

duration = stepripples(:,2)- stepripples(:,1);
stepripples1 = stepripples((~(duration<=1)),:);

%% FOURTH CRITERION:  at least 20ms long

duration2 = (stepripples1(:,2) - stepripples1(:,1))*sint; %duration in sec.
ripplesidx = stepripples1((find(duration2<=0.15 & ...
    duration2 >= 0.025)),:); %Events < 150ms and > 20ms (Indices)
ripplessec = ripplesidx*sint; % Time(sec) End and Start of selected ripples
% startripplessec(kk,1) = zeros(length(ripplesidx{kk,1}) +1 ,1);
% startripplessec(kk,1)(2:end) = ripplessec(kk,1); % Onlt the Start of Ripple Candidates

%% Second Step
below = ripplesidx(:,1)-(0.01/sint);
abov = ripplesidx(:,2)+(0.01/sint);
below(below < 0) = 1; % Substitute start negative value with 1
abov(abov > length(time)) = length(time); % Making ripple candidates finis
linVel = abs(interp1(positiondata.time, positiondata.linVel, time));
VelIdx = find(linVel(below) < 3);
below = below(VelIdx);
abov = abov(VelIdx);
MUAIdx = find(spkDensity(abov) > 15);
below = below(MUAIdx);
abov = abov(MUAIdx);

%% Visual Inspection of LFP and Ripple Candidates
% close all
time1= 1110;
time2= 1122;
xt1= time1/sint;
xt2= time2/sint;
interBelowIdx = find(below.*sint >= time1-5 & below.*sint <= time2+5);
% figure
subplot(711)
plot(time(xt1:xt2), hippLFP(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        hippLFP(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -inf inf])

subplot(712)
plot(time(xt1:xt2), rippleLFP(xt1:xt2),'Color', [93/254, 109/254, 126/254])
hold on
for kk = 1: length(interBelowIdx)
    plot(time(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        rippleLFP(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -inf inf])

subplot(713)
plot(time(xt1:xt2), sharpwaveLFP(xt1:xt2),'Color', 'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        sharpwaveLFP(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
axis([time1 time2 -inf inf])
% line([time1 time2],[th th],'Color','b','LineStyle','--','LineWidth',2)
hold off

subplot(714)
plot(time(xt1:xt2), smsqfiltLFP(xt1:xt2),'Color', [93/254, 109/254, 126/254])
hold on
for kk = 1: length(interBelowIdx)
    plot(time(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        smsqfiltLFP(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
axis([time1 time2 -5000 inf])
line([time1 time2],[th th],'Color','b','LineStyle','--','LineWidth',2)
hold off

subplot(715)
plot(time(xt1:xt2), smoothsharpwaveenvp(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        smoothsharpwaveenvp(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
axis([time1 time2 -5000 inf])
% line([time1 time2],[th th],'Color','b','LineStyle','--','LineWidth',2)
hold off

TetIdx = CA3TetIdx; % make sure you have the tetrode of interest

subplot(716)
hold on
plot(spkClust(TetIdx(1)).spkTime, ... 
    ones(length(spkClust(TetIdx(1)).spkTime),1),'r.','MarkerSize',5)
for kk = 2: size(TetIdx,1)
plot(spkClust(TetIdx(kk)).spkTime, ...
    ones(length(spkClust(TetIdx(kk)).spkTime),1) * kk,... 
    'r.','MarkerSize',5)
end
hold off
xlim([time1 time2])
ylim([0 inf])
% set(gca,'xtickLabel',[])
xlabel('Time[sec]')
ylabel('MUA #')

subplot(717)

spkDensity = CA3spkDensity; % make sure you have the tetrode of interest

plot(time(xt1:xt2),spkDensity(xt1:xt2),'r','LineWidth',2) 
xlabel('Time[sec]')
ylabel('spks/s')
axis tight

% %%
% noisetime1 = 4959;
% noisetime2 = 6433;
%% ripples candidates before and after noisy recording
% below1 = below;
% abov1 = abov;
% noise = find(below.*sint > noisetime1 & below.*sint < noisetime2);
% below1(noise)= [];
% abov1(noise) = [];
% %% ripple candidates during noisy recordig period
% noisetime1 = 5000;
% noisetime2 = 6433;
% below2 = below;
% abov2 = abov;
% noise2 = find(below.*sint < noisetime1 | below.*sint > noisetime2);
% below2(noise2) = [];
% abov2(noise2) = [];
% 
% %%
% abovIdx = zeros(length(abov1)+length(abov2),1);
% abovIdx(1:length(abov1)) =abov1;
% abovIdx(length(abov1)+1:end) =abov2;
% abovIdx = sort(abovIdx);
% 
% belowIdx = zeros(length(below1)+length(below2),1);
% belowIdx(1:length(below1)) =below1;
% belowIdx(length(below1)+1:end) =below2;
% belowIdx = sort(belowIdx);

%% If no noise was detetcted
abovIdx = abov;
belowIdx = below;
%% Plotting Ripples in the Maze
linPos = interp1(positiondata.time, positiondata.linearpos, time);
poscm = interp1(positiondata.time, positiondata.poscm, time);
figure
subplot(121)
hold on
plot(time, linPos, 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
plot(time(belowIdx),linPos(belowIdx),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Ripple Activity Tet#2')
xlabel('Linear Position[cm]')
ylabel('Time[sec]')
axis tight
subplot(122)
plot(poscm(:,1),poscm(:,2), 'Color', [0.5,0.5,0.5])
hold on
plot(poscm(belowIdx,1),poscm(belowIdx,2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
xlabel('x coordinate [cm]')
ylabel('y coordinate [cm]')

%%













