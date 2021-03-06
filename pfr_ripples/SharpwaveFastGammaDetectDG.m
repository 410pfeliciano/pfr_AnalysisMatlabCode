% LFP of Interest
DGLFP = ch10Tet5LFP;
% DGLFPfastgamma = []; % filtered LFP detrend

%% Filter the LFP to have the SharpWave
DGLFPsharpwave = bandFilter (DGLFP, 'sharpwave');
DGLFPfastgamma = bandFilter (DGLFP, 'fastgamma');
%% Hippocampal Spike Density
Fs = 1000; % This value will depend
sint = 1/Fs;
Tet = time; % just to calculate the duration of the recording
bins = (0 : 1/Fs : length(Tet) * (1/Fs) - (1/Fs))';
if exist('spkClust','var') == 0
    warning('spiking data spkClust is not in your workspace');
    return
end

TetDG = [5, 6]; % choose the tetrode of interest
DGTetfind = cell(size(spkClust,2),1);
for kk = 1 : size(spkClust,2)
    DGTetfind{kk} = find(spkClust(kk).tetNum == TetDG);
end
DGTetIdx = find(~cellfun('isempty', DGTetfind));
DGspiketrain = zeros(size(DGTetIdx,1),length(bins));
for kk = 1 : size(DGTetIdx,1)
    DGspiketime = (spkClust(DGTetIdx(kk)).spkTime)';
    DGspiketrain(kk,:) = hist(DGspiketime,bins);
end
% Spike density estimation
DGtotalSpk = sum(DGspiketrain,1);
sigma = 0.02;
edges = -3*sigma: 0.001: 3*sigma;
kernel = (normpdf(edges,0,sigma)).*(1/10);
% kernel = kernel*(1/Fs); % kernel muitiplied by bin size
DGspkDensity = conv(DGtotalSpk, kernel);
center = ceil(length(edges)/2);
DGspkDensity = DGspkDensity(center:length(DGtotalSpk)+center-1);

%% Second Step, Calculation the mean power and standar deviation of the recording
filtWin = 0.02; % filter window in in seconds
hibfiltLFP = abs(hilbert(DGLFPfastgamma)); % amplitude envelope of ripple-band filtered LFP
sharpLFPenv = abs(hilbert(DGLFPsharpwave));
sqfiltLFP= hibfiltLFP.^2; % calculating the power of the filtered LFP
sqsharpwaveEnvp = sharpLFPenv.^2;
smsqfiltLFP = smoothdata(sqfiltLFP,'gaussian',filtWin/sint); % smoothing the power data
smoothsharpwaveenvp = smoothdata(sqsharpwaveEnvp,'gaussian',filtWin/sint);
%% Plotting LFP and MUA envelop
time1= 1255;
time2= 1260;
xt1= time1/sint;
xt2= time2/sint;

subplot(5,1,1)
plot(time(xt1:xt2),DGLFP(xt1:xt2),'k') 
axis([-inf inf -inf inf])
legend('Ctx','Location', 'best')
legend boxoff
box off
set(gca,'xtickLabel',[])
ylabel('uV')

subplot(5,1,2)
plot(time(xt1:xt2),DGLFPfastgamma(xt1:xt2),'r') 
axis tight
legend('Fast Gamma','Location', 'best')
legend boxoff
box off
set(gca,'xtickLabel',[])
ylabel('uV')

subplot(5,1,3)
plot(time(xt1:xt2),DGLFPsharpwave(xt1:xt2),'b') 
axis tight
legend('Sharpwave','Location', 'best')
legend boxoff
box off
% set(gca,'xtickLabel',[])
ylabel('uV')

subplot(5,1,4)
plot(time(xt1:xt2),smsqfiltLFP(xt1:xt2),'r') 
axis tight
legend('Smooth Fast-Gamma Envp.','Location', 'best')
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
h.BinEdges = 0:600000;
h.BinWidth = 10;
h.EdgeColor = 'k';
xlabel('Fast Gamma Power Density')
ylabel('Probability Density')
legend('Fast Gamma Power')
legend('boxoff')  
box off

%%  FIRST CRITERION: AT LEAST 3*STD OF FILTERED, SQUARED LFP
time1= 1;
time2= 1600;
% time1= 1000;
% time2= 1600;
xt1= time1/sint;
xt2= time2/sint;
thMad = mad(smsqfiltLFP(xt1:xt2),0); % median deviation flag = 1 mean deviation flag =0
th = median(smsqfiltLFP(xt1:xt2)) + (4*thMad);
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
%% Second Step
% load('PosVelSites.mat')
below = ripplesidx(:,1)-(0.01/sint);
abov = ripplesidx(:,2)+(0.01/sint);
below(below < 0) = 1; % Substitute start negative value with 1
abov(abov > length(time)) = length(time); % Making ripple candidates finis
linVel = abs(interp1(positiondata.time, positiondata.linVel, time));
VelIdx = find(linVel(below) < 3);
below = below(VelIdx);
abov = abov(VelIdx);
%% Optional Putting a threhold here
% MUAthMad = mad(DGspkDensity,0); % median deviation flag = 1 mean deviation flag =0
% MUAth = median(DGspkDensity) + (MUAthMad);
% MUAIdx = find(DGspkDensity(abov-10) > MUAth);
% below = below(MUAIdx);
% abov = abov(MUAIdx);

%% Visual Inspection of LFP and Ripple Candidates
% close all
time1= 1590; 
time2= 1610; 
xt1= time1/sint;
xt2= time2/sint;
interBelowIdx = find(below.*sint >= time1-5 & below.*sint <= time2+5);
% figure
subplot(711)
plot(time(xt1:xt2), DGLFP(xt1:xt2),'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        DGLFP(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -inf inf])

subplot(712)
plot(time(xt1:xt2), DGLFPfastgamma(xt1:xt2),'Color', [93/254, 109/254, 126/254])
hold on
for kk = 1: length(interBelowIdx)
    plot(time(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        DGLFPfastgamma(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        'r', 'LineWidth',1)
end
hold off
axis([time1 time2 -inf inf])

subplot(713)
plot(time(xt1:xt2), DGLFPsharpwave(xt1:xt2),'Color', 'k')
hold on
for kk = 1: length(interBelowIdx)
    plot(time(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
        DGLFPsharpwave(below(interBelowIdx(kk)):abov(interBelowIdx(kk))),...
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


subplot(716)
hold on
plot(spkClust(DGTetIdx(1)).spkTime, ... 
    ones(length(spkClust(DGTetIdx(1)).spkTime),1),'r.','MarkerSize',5)
for kk = 2: size(DGTetIdx,1)
plot(spkClust(DGTetIdx(kk)).spkTime, ...
    ones(length(spkClust(DGTetIdx(kk)).spkTime),1) * kk,... 
    'r.','MarkerSize',5)
end
hold off
xlim([time1 time2])
ylim([0 inf])
% set(gca,'xtickLabel',[])
xlabel('Time[sec]')
ylabel('MUA #')

subplot(717)
plot(time(xt1:xt2),DGspkDensity(xt1:xt2),'r','LineWidth',2) 
xlabel('Time[sec]')
ylabel('spks/s')
axis tight

%% Clasification ripples during sleep and run time
runtime1 = 1609.6; % Run start time 
runtime2 = 7700; % Run start time 2
Sleepbelow = below;
Sleepabov = abov;
noise = find(below.*sint > runtime1 & below.*sint < runtime2);
Sleepbelow(noise)= [];
Sleepabov(noise) = [];
% Pre Sleep Ripple Candidates
PreSleepbelow = below;
PreSleepabov = abov;
noise = find(below.*sint > runtime1);
PreSleepbelow(noise)= [];
PreSleepabov(noise) = [];
% Post Sleep
PostSleepbelow = below;
PostSleepabov = abov;
noise = find(below.*sint < runtime2);
PostSleepbelow(noise)= [];
PostSleepabov(noise) = [];
% Run Ripple Candidates
RUNbelow1 = below;
RUNabov1 = abov;
noise = find(below.*sint < runtime1 | below.*sint > runtime2);
RUNbelow1(noise)= [];
RUNabov1(noise) = [];
%% Plotting Ripples in the Maze
linPos = interp1(positiondata.time, positiondata.linearpos, time);
poscm = interp1(positiondata.time, positiondata.poscm, time);
figure
subplot(421)
hold on
plot(time, linPos, 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
plot(time(Sleepbelow),linPos(Sleepbelow),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 1)
hold off
title('Ripple Activity During Sleep')
xlabel('Linear Position[cm]')
ylabel('Time[sec]')
axis tight
subplot(422)
plot(poscm(:,1),poscm(:,2), 'Color', [0.5,0.5,0.5])
hold on
plot(poscm(Sleepbelow,1),poscm(Sleepbelow,2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 1)
hold off
xlabel('x coordinate [cm]')
ylabel('y coordinate [cm]')

subplot(423)
hold on
plot(time, linPos, 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
plot(time(PreSleepbelow),linPos(PreSleepbelow),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 1)
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Ripple Activity Pre-Sleep')
xlabel('Linear Position[cm]')
ylabel('Time[sec]')
axis tight
subplot(424)
plot(poscm(:,1),poscm(:,2), 'Color', [0.5,0.5,0.5])
hold on
plot(poscm(PreSleepbelow,1),poscm(PreSleepbelow,2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 1)
hold off
xlabel('x coordinate [cm]')
ylabel('y coordinate [cm]')

subplot(425)
hold on
plot(time, linPos, 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
plot(time(PostSleepbelow),linPos(PostSleepbelow),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 1)
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Ripple Activity Post-Sleep')
xlabel('Linear Position[cm]')
ylabel('Time[sec]')
axis tight
subplot(426)
plot(poscm(:,1),poscm(:,2), 'Color', [0.5,0.5,0.5])
hold on
plot(poscm(PostSleepbelow,1),poscm(PostSleepbelow,2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 1)
hold off
xlabel('x coordinate [cm]')
ylabel('y coordinate [cm]')

subplot(427)
hold on
plot(time, linPos, 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
plot(time(RUNbelow1),linPos(RUNbelow1),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 1)
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Ripple Activity RUN')
xlabel('Linear Position[cm]')
ylabel('Time[sec]')
axis tight
subplot(428)
plot(poscm(:,1),poscm(:,2), 'Color', [0.5,0.5,0.5])
hold on
plot(poscm(RUNbelow1,1),poscm(RUNbelow1,2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 1)
hold off
xlabel('x coordinate [cm]')
ylabel('y coordinate [cm]')

%% Sleep Ripple Candidates
DGRipSlAbovIdx = Sleepabov;
DGRipSlBelowIdx = Sleepbelow;
DGRipPreSlAbovIdx = PreSleepabov;
DGRipPreSlBelowIdx = PreSleepbelow;
DGRipPostSlAbovIdx = PostSleepabov;
DGRipPostSlBelowIdx = PostSleepbelow;
DGRunRipBelowIdx = RUNbelow1;
DGRunRipAbovIdx = RUNabov1;
% %% Second detection in case of noisy recordings segments
% % Eliminating ripples during noise time
% runtime1 = 4961;
% runtime2 = 6433;
% noise = find(RUNbelow1.*sint > runtime1 & RUNbelow1.*sint < runtime2);
% RUNbelow1(noise)= [];
% RUNabov1(noise) = [];
% figure
% hold on
% plot(time, linPos, 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
% plot(time(RUNbelow1),linPos(RUNbelow1),'o',...
%     'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 1)
% hold off
% % legend('Linear Position','Single Unit Spikes','location','southoutside')
% title('Ripple Activity Tet#2')
% xlabel('Linear Position[cm]')
% ylabel('Time[sec]')
% axis tight
% %% Second threshold for regions in the recording that are more noisy
% time1= 4961; % Start noisy region
% time2= 6433; % Finish noisy region
% xt1= time1/sint;
% xt2= time2/sint;
% thMad2 = mad(smsqfiltLFP(xt1:xt2),0); % median deviation flag = 1 mean deviation flag =0
% th2 = median(smsqfiltLFP(xt1:xt2)) + (7*thMad2);
% splog = smsqfiltLFP > th2;
% pos = find(splog == 1); % This vectors contain the indices of power larger the 3std
% difpos = find(diff(pos) >= (0.02/sint)); % events with inte4rvals > 20ms
% z = zeros(length(difpos)+1,1); % Creating a zero row vector
% z(1:end-1) = difpos;
% z(end) = length(pos); % End position of all potential ripples
% st = difpos+1; % Start positions of all possible ripple events
% zst = zeros(length(difpos)+1,1);
% zst(2:end) = st;
% zst(1)=1; % zst start position of all
% stepripples = [pos(zst), pos(z)]; % Start and End position(indices) 
% duration = stepripples(:,2)- stepripples(:,1);
% stepripples1 = stepripples((~(duration<=1)),:);
% duration2 = (stepripples1(:,2) - stepripples1(:,1))*sint; %duration in sec.
% ripplesidx = stepripples1((find(duration2<=0.15 & ...
%     duration2 >= 0.025)),:); %Events < 150ms and > 20ms (Indices)
% below = ripplesidx(:,1)-(0.01/sint);
% abov = ripplesidx(:,2)+(0.01/sint);
% below(below < 0) = 1; % Substitute start negative value with 1
% abov(abov > length(time)) = length(time); % Making ripple candidates finis
% linVel = abs(interp1(positiondata.time, positiondata.linVel, time));
% VelIdx = find(linVel(below) < 3);
% below = below(VelIdx);
% abov = abov(VelIdx);
% 
% %% Eliminating ripples during sleep time
% noisy1= 4961; % Start noisy region
% noisy2= 6433; % Finish noisy region
% % ripples candidates before and after noisy recording
% RUNbelow2 = below;
% RUNabov2 = abov;
% noise = find(below.*sint < noisy1 | below.*sint > noisy2);
% RUNbelow2(noise)= [];
% RUNabov2(noise) = [];
% %%
% figure
% hold on
% plot(time, linPos, 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
% plot(time(RUNbelow2),linPos(RUNbelow2),'o',...
%     'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 1)
% hold off
% % legend('Linear Position','Single Unit Spikes','location','southoutside')
% title('Ripple Activity Tet#2')
% xlabel('Linear Position[cm]')
% ylabel('Time[sec]')
% axis tight
% 
% %% Mergin all ripple candidates duribng running
% DGRunRipBelowIdx = [RUNbelow1;RUNbelow2];
% DGRunRipAbovIdx = [RUNabov1;RUNabov2];
% DGRunRipBelowIdx = sort(DGRunRipBelowIdx);
% DGRunRipAbovIdx = sort(DGRunRipAbovIdx);
% %%
% figure
% hold on
% plot(time, linPos, 'Color', [0.5,0.5,0.5], 'LineWidth', 1)
% plot(time(DGRunRipBelowIdx),linPos(DGRunRipBelowIdx),'o',...
%     'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 1)
% hold off
% % legend('Linear Position','Single Unit Spikes','location','southoutside')
% title('Ripple Activity Tet#2')
% xlabel('Linear Position[cm]')
% ylabel('Time[sec]')
% axis tight

%% Saving DG Parameters
save('DGgammaDetectionData.mat','DGLFP','DGLFPsharpwave','DGLFPfastgamma',...
    'DGTetIdx','DGspkDensity','DGRipSlAbovIdx','DGRipSlBelowIdx',...
    'DGRipPreSlAbovIdx','DGRipPreSlBelowIdx','DGRipPostSlAbovIdx',...
    'DGRipPostSlBelowIdx','DGRunRipBelowIdx','DGRunRipAbovIdx')



