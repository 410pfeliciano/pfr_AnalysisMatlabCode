%% Bin Spiking Data to Visualize Spiking in time and Position
close all
% Increase the sample rate of position and the input signal by a factor of r.
videoSR = 1/15;
int_factor = 30;
linear_pos = interp(positiondata.linearpos,int_factor);
new_SR = videoSR / int_factor;
new_t = (0 : new_SR : length(linear_pos) * new_SR - new_SR)';
velocity = interp(positiondata.linVel,int_factor);
figure
subplot(2,1,1)
plot(new_t, linear_pos, 'r', 'LineWidth',1)
legend('Raw Position')
xlabel('Time [sec]')
ylabel('Position [cm]')
axis tight
subplot(2,1,2)
plot(new_t, velocity, 'r')
xlabel('Time [sec]')
ylabel('Speed [cm * sec^{-1}]')
legend('Interpolated Velocity')
axis tight
%% Load Spiking Information
close all
x_pos = positiondata.poscm(:,1);
new_x_pos = interp(x_pos,int_factor);
y_pos = positiondata.poscm(:,2);
new_y_pos = interp(y_pos,int_factor);
bins = (0 : new_SR : length(linear_pos) * new_SR - new_SR);
for kk = 1  : size(spkClust,2)
spiketime = spkClust(kk).spkTime;
spiketrain = histcounts(spiketime,bins);
% spiketrain = h.Values;
% Spike train when animnalis moving
spiketrain2 = zeros(length(spiketrain),1);
for qq = 1:length(spiketrain2)
    if velocity(qq) > 4
         
spiketrain2(qq) = spiketrain(qq);
    else
        spiketrain2(qq) = 0;
    end
end
% Plotting Spikes and Linear Position
spikeindex = find(spiketrain);
figure
subplot(3,1,1)
hold on
plot(positiondata.linearpos, positiondata.time, 'Color', [0.7,0.7,0.7])
plot(linear_pos(spikeindex),new_t(spikeindex),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Single Unit Activity')
xlabel('Linear Position [cm]')
ylabel('Time [sec]')
axis([0 inf 1500 8000])

subplot(3,1,2)
spikeindex2 = find(spiketrain2);
hold on
plot(positiondata.linearpos, positiondata.time, 'Color', [0.7,0.7,0.7])
plot(linear_pos(spikeindex2),new_t(spikeindex2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Single Unit Activity During Movement')
xlabel('Linear Position [m]')
ylabel('Time [sec]')
axis([0 inf 1500 8000])

subplot(3,1,3)
% runStartTime = 1609.93;
% runEndTime = runStartTime + 6023.79;
plot(new_x_pos, new_y_pos,'k') 
axis tight
xlabel('x coordinate [cm]')
ylabel('y coordinate [cm]')
hold on;
%  Add dots where spikes occured. 
plot(new_x_pos(spikeindex2), new_y_pos(spikeindex2),'or', ...
    'MarkerSize', 2, 'MarkerFaceColor', 'r')
title('Position and Spikes')
hold off
end
% subplot(2,2,4)
% hold on
% plot( new_t, linear_pos, 'Color', [0.7,0.7,0.7],'LineWidth',2)
% plot(new_t(spikeindex2),linear_pos(spikeindex2),'o',...
%     'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 3)
% hold off
% % legend('Linear Position','Single Unit Spikes','location','southoutside')
% title('Single Unit Activity During Movement')
% xlabel('Time [sec]')
% ylabel('Linear Position [m]')
% axis tight
%% First plot is the rat's path + dots at actual spikes 
close all
x_pos = positiondata.poscm(:,1);
new_x_pos = interp(x_pos,int_factor);
y_pos = positiondata.poscm(:,2);
new_y_pos = interp(y_pos,int_factor);

figure
subplot(1,4,1)
plot(new_x_pos, new_y_pos,'k') 
axis tight
xlabel('x coordinate [cm]')
ylabel('y coordinate [cm]')
hold on;
%  Add dots where spikes occured. 
plot(new_x_pos(spikeindex), new_y_pos(spikeindex),'or', ...
    'MarkerSize', 2, 'MarkerFaceColor', 'r')
title('Position and Spikes')
hold off
%%
subplot(1,4,4)
Xlimit = (0:.01:.75);
Ylimit = (0:.008:.70);
pos_hist= histogram2(new_x_pos, new_y_pos, Xlimit,Ylimit, ...
    'DisplayStyle','tile','ShowEmptyBins','off');
colormap(bone)
title('Position Hist')
xlabel('x coordinate [m]')
ylabel('y coordinate [m]')
% total_time = sum(sum(pos_hist.BinCounts));

subplot(1,4,3)
spikeHist_pos = histogram2(new_x_pos(spikeindex),new_y_pos(spikeindex), ...
    Xlimit,Ylimit,'DisplayStyle','tile','ShowEmptyBins','on');
title('Spikes Hist')
xlabel('x coordinate [m]')
ylabel('y coordinate [m]')
pos_occupancy = pos_hist.BinCounts.*new_SR; % changing occupanncy to seconds
pos_spiketime = spikeHist_pos.BinCounts;
pos_FiringRate = pos_spiketime./pos_occupancy;
pos_FiringRate(isnan(pos_FiringRate))=0;

subplot(1,4,2)
histogram2('XBinEdges',Xlimit,'YBinEdges',Ylimit,'BinCounts', ...
    pos_FiringRate,'DisplayStyle','tile','ShowEmptyBins','on')
title('Ocupanccy Norm')
xlabel('x coordinate [m]')
ylabel('y coordinate [m]')
c = colorbar;
c.Label.String = 'Occupancy Norm. (Spikes/sec)';
% caxis auto
caxis([2 10])

figure
subplot(1,4,1)
plot(new_x_pos, new_y_pos,'k') 
axis([0 .75 0 .70])
xlabel('x coordinate [m]')
ylabel('y coordinate [m]')
hold on;
%  Add dots where spikes occured. 
plot(new_x_pos(spikeindex2), new_y_pos(spikeindex2),'or', ...
    'MarkerSize', 2, 'MarkerFaceColor', 'r')
title('Position and Spikes')
hold off

subplot(1,4,2)
Xlimit = (0:.01:.75);
Ylimit = (0:.008:.70);
pos_hist= histogram2(new_x_pos, new_y_pos, Xlimit,Ylimit, ...
    'DisplayStyle','tile','ShowEmptyBins','off');
colormap(bone)
title('Position Hist')
xlabel('x coordinate [m]')
ylabel('y coordinate [m]')
total_time = sum(sum(pos_hist.BinCounts));

subplot(1,4,3)
spikeHist_pos = histogram2(new_x_pos(spikeindex2),new_y_pos(spikeindex2), ...
    Xlimit,Ylimit,'DisplayStyle','tile','ShowEmptyBins','on');
title('Spikes Hist')
xlabel('x coordinate [m]')
ylabel('y coordinate [m]')
pos_occupancy = pos_hist.BinCounts.*new_SR;
pos_spiketime = spikeHist_pos.BinCounts;
pos_FiringRate = pos_spiketime./pos_occupancy;
pos_FiringRate(isnan(pos_FiringRate))=0;

subplot(1,4,4)
histogram2('XBinEdges',Xlimit,'YBinEdges',Ylimit,'BinCounts', ...
    pos_FiringRate,'DisplayStyle','tile','ShowEmptyBins','on')
title('Occup. During Locomotion')
xlabel('x coordinate [m]')
ylabel('y coordinate [m]')
c = colorbar;
c.Label.String = 'Occupancy Norm. (Spikes/sec)';
% caxis auto
caxis([2 20])

figure
subplot(1,2,1)
plot(new_x_pos, new_y_pos,'k') 
axis([0 .75 0 .70])
xlabel('x coordinate [m]')
ylabel('y coordinate [m]')
hold on;
%  Add dots where spikes occured. 
plot(new_x_pos(spikeindex2), new_y_pos(spikeindex2),'or', ...
    'MarkerSize', 2, 'MarkerFaceColor', 'r')
title('Position and Spikes')
hold off

subplot(1,2,2)
histogram2('XBinEdges',Xlimit,'YBinEdges',Ylimit,'BinCounts', ...
    pos_FiringRate,'DisplayStyle','tile','ShowEmptyBins','on')
title('Occup. During Locomotion')
xlabel('x coordinate [m]')
ylabel('y coordinate [m]')
c = colorbar;
c.Label.String = 'Occupancy Norm. (Spikes/sec)';
colormap(bone)
% caxis auto
caxis([0 30])

%% Gaussian Estimation
close all
lpos = linear_pos;
pos_bins = (0:.02:1.30)';
x_bins = (0:.025:.70);
y_bins = (0:.025:.70);
spikehist =(hist(lpos(spikeindex),pos_bins))';%Histogram positions @ spikes.
occupancy =(hist(lpos,pos_bins)* new_SR)';	%Convert occupancy to seconds.
spikehist2 =(hist(lpos(spikeindex2),pos_bins))';%Histogram positions @ spikes.
spikes_per_sec = spikehist2./occupancy;
spikes_per_sec(isnan(spikes_per_sec)) = 0;
start_arm = find(lpos <= .35);
center_arm = find(lpos > .35 & lpos < .90);
reward_arm = find(lpos >= .90);

% Gaussian fits Practice
f1 = fit(lpos ,spiketrain2,'gauss1');
figure % figure
subplot(2,1,1)
plot(f1,lpos,spiketrain2)
hold on
a1 = f1.a1/new_SR;
b1 = f1.b1;
c1 = f1.c1;
f1_x = a1*exp(-((lpos-b1)/c1).^2);
plot(lpos, f1_x)
axis([0 1.30 0 inf])
hold off
subplot(2,1,2)
plot(f1,lpos,spiketrain2,'Residuals')
axis([0 1.30 -inf inf])

[f2, f2gof, f2out] = fit(pos_bins ,spikes_per_sec,'gauss1');
figure % Figure
subplot(2,1,1)
bar(pos_bins, spikes_per_sec)
hold on
hplot1 = plot(f2,pos_bins,spikes_per_sec,'predfunc');
set(hplot1, 'LineWidth',2)
hold off
axis([0 1.30 0 inf])
subplot(2,1,2)
plot(f2,pos_bins,spikes_per_sec,'Residuals')
axis([0 1.30 -inf inf])

figure % Figure 
subplot(2,1,2)
bar(pos_bins, spikes_per_sec)
hold on
% hplot2 = plot(f2,'r','predfunc');
% set(hplot2, 'LineWidth',2)
hold off
xlabel('Position [m]')			%Label the axes.
ylabel('Occup. Norm(Spikes/sec)')
legend('off')
axis([0 1.30 0 inf])

subplot(2,1,1)
plot(lpos, new_t, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
hold on
plot(lpos(center_arm), new_t(center_arm),'.','Color', [0.55,0.55,0.55])
plot(lpos(start_arm), new_t(start_arm),'.','Color', [0.7,0.7,0.7])
plot(lpos(reward_arm), new_t(reward_arm),'.','Color', [0.7,0.7,0.7])
plot(lpos(spikeindex2),new_t(spikeindex2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 4)
title('Single Unit Activity')
xlabel('Linear Position [m]')
ylabel('Time [sec]')
axis([0 1.30 0 1009])
hold off

figure % Figure 
subplot(2,3,1)
hold on
plot(lpos, new_t, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
plot(lpos(spikeindex),new_t(spikeindex),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 4)
hold off
title('Single Unit Activity(Quite/Mov)')
xlabel('Linear Position [m]')
ylabel('Time [sec]')
axis([0 1.30 0 1009])

subplot(2,3,4)
bar(pos_bins,spikehist./occupancy); %Plot results as bars.
xlabel('Position [m]')			%Label the axes.
ylabel('Occup. Norm Counts (spikes/s)')
axis([0 1.30 0 inf])

subplot(2,3,2)
hold on
plot(lpos, new_t, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
plot(lpos(spikeindex2),new_t(spikeindex2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 4)
hold off
% legend('Linear Position','Single Unit Spikes','location','northoutside')
title('Single Unit Activity(Mov)')
xlabel('Linear Position [m]')
ylabel('Time [sec]')
axis([0 1.30 0 1009]) 

subplot(2,3,5)
bar(pos_bins,spikes_per_sec); %Plot results as bars.
xlabel('Position [m]')			%Label the axes.
ylabel('Occup. Norm Counts (spikes/s)')
axis([0 1.30 0 inf])

subplot(2,3,3)
hold on
plot(lpos, new_t, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
plot(lpos(spikeindex2),new_t(spikeindex2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 4)
hold off
% legend('Linear Position','Single Unit Spikes','location','northoutside')
title('Single Unit Activity(Mov)')
xlabel('Linear Position [m]')
ylabel('Time [sec]')
axis([0 1.30 0 1009])      
subplot(2,3,6)
bar(pos_bins,spikehist2); %Plot results as bars.
xlabel('Position [m]')			%Label the axes.
ylabel('Spike Counts')
axis([0 1.30 0 inf])

%% Fit Model 3 to the spike train data (omitting last input).
close all 
% 'Weights' vector
w_data = zeros(length(spiketrain2),1);
for qq = 1:length(w_data)
    if lpos(qq) > 0.7 && lpos(qq) < 1.1
        w_data(qq) = spiketrain2(qq) + 50;
    else
        w_data(qq) = 1;
    end
end
w_data1 = zeros(length(spiketrain2),1);
for qq = 1:length(w_data1)
    if lpos(qq) > 0 && lpos(qq) < .6
        w_data1(qq) = spiketrain2(qq) + 10;
    else
        w_data1(qq) = 1;
    end
end
% Exluding logical vector (1= observation to exclude; 0 = observation to
% include
excl = zeros(length(spiketrain2),1);
for qq = 1:length(excl)
    if lpos(qq) > 0 && lpos(qq) < 0.6
        excl(qq) = 1;
    else
        excl(qq) = 0;
    end
end
predictor = [lpos lpos.^2 lpos.^4 lpos.^6 lpos.^8 lpos.^10];
predictor2 = [lpos lpos.^2];
[beta3, dev3, stats3] = glmfit(predictor,spiketrain2,... 
    'poisson');
beta = fitglm(predictor2,spiketrain2,'linear','distr','poisson',... 
    'Exclude',logical(excl), 'Weights', w_data);
beta1 = table2array(beta.Coefficients(:,1));
lambda_estimates = exp([ones(size(spiketrain2,1),1) predictor] * beta3);

xs=((0:.875:130.5)./100)';	%Define interval of positions,
%Evaluate Model 4 in direction 0 (X decreases).
newpred = [xs xs.^2 xs.^4 xs.^6 xs.^8 xs.^10];
newpred2 = [xs xs.^2];
[lambda3_0, up0, low0]=glmval(beta3,newpred,'log',stats3);
[lambda3_1]=glmval(beta1,newpred2,'log');

figure, 
bar(pos_bins,spikes_per_sec); %Plot occupancy norm. hist.
hold on							%...freeze graphics,\
plot(xs,lambda3_0./new_SR,'r', 'LineWidth', 4);			%Plot Model 4, X decreasing,
plot(xs,lambda3_1./new_SR,'m', 'LineWidth', 4);			%Plot Model 4, X decreasing,
plot(xs,lambda3_0./new_SR+up0./new_SR,'b--')	%...add upper CI,
plot(xs,lambda3_0./new_SR-low0./new_SR,'b--')	%...and lower CI.
axis([0 1.3 0 inf])
hold off

figure, 
hold on
bar(pos_bins,spikes_per_sec); %Plot occupancy norm. hist.
plot(lpos, lambda_estimates./new_SR,'r'),
plot(xs,lambda3_0./new_SR,'b','LineWidth',4);			%Plot Model 4, X decreasing,
plot(xs,lambda3_0./new_SR+up0./new_SR,'b--')	%...add upper CI,
plot(xs,lambda3_0./new_SR-low0./new_SR,'b--')	%...and lower CI.
axis([0 1.30 0 20]), hold off

figure, subplot(3,1,1)
bar(pos_bins,spikes_per_sec); %Plot occupancy norm. hist.
hold on
plot(lpos, lambda_estimates./new_SR,'r')
plot(xs,lambda3_0./new_SR,'b','LineWidth',2);
hold off              
xlabel('Position [m]')		
ylabel('Occupancy Norm. (spikes/s)')
axis([0 1.30 0 20])
subplot(3,1,2)
hold on
plot(lpos, new_t, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
plot(lpos(spikeindex2),new_t(spikeindex2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 5)
hold off
title('Single Unit Activity')
xlabel('Linear Position [m]')
ylabel('Time [sec]')
axis([0 1.30 0 1009])
subplot(3,1,3)
plot(lpos, spiketrain2, 'o')
xlabel('Position [m]')		
ylabel('Spikes')
axis([0 1.30 0.25 2-.25])

%% GLM model # 4
close all
dir = [0; diff(lpos)>0];
w_data = zeros(length(spiketrain),1);
for qq = 1:length(w_data)
    if linear_pos(qq) > 0 && linear_pos(qq) < 0.6
        w_data(qq) = spiketrain2(qq) + 1;
    else
        w_data(qq) = 1;
    end
end
predictor = [lpos lpos.^2 lpos.^4 lpos.^6 ... 
    dir_RSRR dir_RRLR dir_LSLR dir_LLRL];
[b4,dev4,stats4]=glmfit(predictor,spiketrain2,'poisson','Weights',w_data);

% figure
% subplot(3,1,1)
% bar(pos_bins,spikes_per_sec); %Plot occupancy norm. hist.
% hold on;   
% yfit4 = glmval(b4,predictor,'log');
% yfit4 = yfit4./new_SR;
% plot(lpos, yfit4)
% hold off              
% xlabel('Position [m]')		
% ylabel('Occupancy Norm. (spikes/s)')
% axis([0 1.30 0 inf])
% subplot(3,1,2)
% hold on
% plot(lpos, new_t, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
% plot(lpos(spikeindex2),new_t(spikeindex2),'o',...
%     'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 5)
% hold off
% % legend('Linear Position','Single Unit Spikes','location','northoutside')
% title('Single Unit Activity')
% xlabel('Linear Position [m]')
% ylabel('Time [sec]')
% axis([0 1.30 0 1009])
% subplot(3,1,3)
% plot(lpos, spiketrain2, 'o')
% xlabel('Position [m]')		
% ylabel('Spikes')
% axis([0 1.30 0.25 2-.25])

xs=((0:.875:130.5)./100)';	%Define interval of positions
Ns = [size(xs,1),1];
newpredict = [xs xs.^2 xs.^4 xs.^6];


figure
%Evaluate Model 4 in direction 0 (X decreases).
[lambda4_0, up0, low0]=glmval(b4,[newpredict ... 
    zeros(Ns) ones(Ns) zeros(Ns) zeros(Ns)],'log',stats4);
%Evaluate Model 4 in direction 1 (X increases).
[lambda4_1, up1, low1]=glmval(b4,[newpredict ... 
    ones(Ns) zeros(Ns) zeros(Ns) zeros(Ns)],'log',stats4);
bar(posbins2, RSRRSpikesPerSec,'r')
hold on							%...freeze graphics,
% plot(lpos, yfit4)
bar(posbins2, RRLRSpikesPerSec, 'g')
bar(pos_bins,spikes_per_sec); %Plot occupancy norm. hist.
hhh1 = plot(xs,lambda4_0./new_SR,'g','LineWidth',4);			%Plot Model 4, X decreasing,
plot(xs,lambda4_0./new_SR+up0./new_SR,'g--')	%...add upper CI,
plot(xs,lambda4_0./new_SR-low0./new_SR,'g--')	%...and lower CI.
hhh2 = plot(xs,lambda4_1./new_SR,'r','LineWidth',4);			%Plot Model 4, X increasing,
plot(xs,lambda4_1./new_SR+up1./new_SR,'r--')	%...add upper CI,
plot(xs,lambda4_1./new_SR-low1./new_SR,'r--');	%...add lower CI.
legend([hhh1 hhh2],'Reward to Start','Start to Reward')
axis tight
hold off

figure
plot(new_t, yfit4)
hold on
plot(new_t,lpos/10, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
plot(new_t(spikeindex2), lpos(spikeindex2)/10,'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 5)
plot(new_t, lpos./0.1)
hold off

%%
w_data = zeros(length(spiketrain),1);
for qq = 1:length(w_data)
    if linear_pos(qq) > 0.1 && linear_pos(qq) < 0.4
        w_data(qq) = spiketrain2(qq) + 1000;
    else
        w_data(qq) = 1;
    end
end
% excluding vector
excl = (find(lpos>.4 & lpos<.8))';
predictor = [lpos lpos.^2 dir_RSRR1 dir_LSLR1 dir_RRLR dir_LLRL];
[b5,dev5,stats5]=glmfit(predictor,spiketrain2,'poisson');
beta5 = fitglm(predictor,spiketrain2,'linear','Distribution','poisson','Exclude',excl,'weights', w_data);
beta51 = table2array(beta5.Coefficients(:,1));


figure
bar(pos_bins,spikes_per_sec); %Plot occupancy norm. hist.
hold on;   
yfit5 = glmval(b5,predictor,'log'); 
yfit51 = glmval(beta51,predictor,'log'); 
yfit5 = yfit5./new_SR;  
yfitInd5 = find(yfit5>=0);  
plot(lpos, yfit5,'r')  
hold off              
xlabel('Position [m]')		
ylabel('Occupancy Norm. (spikes/s)')    
axis([0 1.30 0 inf]) 


figure
plot(new_t, yfit5)
hold on
plot(new_t,lpos/10, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
plot(new_t(spikeindex2), lpos(spikeindex2)/10,'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 5)
plot(new_t, lpos./0.1)
hold off


xs=((0:.875:130.5)./100)';					%Define interval of positions,
Ns = [size(xs,1),1];

%Evaluate Model 4 in direction 0 (X decreases).
% dir_RSRR,dir_LSLR,dir_RRRS,dir_RRLS,dir_LRLS,dir_LRRS
[lambda5_0, up0, low0]=glmval(b5,[xs xs.^2 ones(Ns) zeros(Ns) ...
    zeros(Ns) zeros(Ns)],'log',stats5);
[lambda5_1, up1, low1]=glmval(b5,[xs xs.^2 zeros(Ns) ones(Ns) ...
    zeros(Ns) zeros(Ns)],'log',stats5);
[lambda5_2, up2, low2]=glmval(b5,[xs xs.^2 zeros(Ns) zeros(Ns) ...
    ones(Ns) zeros(Ns)],'log',stats5);
[lambda5_3, up3, low3]=glmval(b5,[xs xs.^2 zeros(Ns) zeros(Ns) ...
    zeros(Ns) ones(Ns)],'log',stats5);


figure
subplot(1,2,1)
bar(pos_bins, spikes_per_sec)
hold on	
hhh1 = plot(xs,lambda5_0./new_SR,'b','LineWidth',2); %Plot Model 4, X decreasing
plot(xs,(lambda5_0./new_SR)+up0./new_SR,'b--')	%...add upper CI,
plot(xs,(lambda5_0./new_SR)-low0./new_SR,'b--')	%...and lower CI.
hhh2 = plot(xs,lambda5_1./new_SR,'r','LineWidth',2); %Plot Model 4
plot(xs,(lambda5_1./new_SR)+up1./new_SR,'r--')	%...add upper CI,
plot(xs,(lambda5_1./new_SR)-low1./new_SR,'r--');	%...add lower CI.
hold off
legend([hhh1 hhh2],'RS-RR','LS-LR')
title('Choice Routes')
axis([0 1.30 0 35])

subplot(1,2,2)
bar(posbins2, llrlSpikesPerSec)
hold on
hhh3 = plot(xs,lambda5_2./new_SR,'k','LineWidth',2); %Plot Model 4
plot(xs,lambda5_2./new_SR+up2./new_SR,'k--')	%...add upper CI,
plot(xs,lambda5_2./new_SR-low2./new_SR,'k--');	%...add lower CI.
hhh4 = plot(xs,lambda5_3./new_SR,'g','LineWidth',2); %Plot Model 4
plot(xs,lambda5_3./new_SR+up3./new_SR,'g--')	%...add upper CI,
plot(xs,lambda5_3./new_SR-low3./new_SR,'g--');	%...add lower CI.
legend([hhh3 hhh4],'RR-LR','LL-RL')
title('Forced Routes')
axis([0 1.30 0 35])
hold off
%% 
% close all
figure
subplot(2,3,1)
plot(new_t, dir_RSRR,'k', 'LineWidth',2)
hold on
plot(new_t, lpos)
hold off
legend('Right Start to Right Reward Vector','Location', 'northoutside')
axis tight

subplot(2,3,2)
plot(new_t, dir_LSLR,'k', 'LineWidth',2)
hold on
plot(new_t, lpos)
hold off
legend('Left Start to Left Reward Vector','Location', 'northoutside')
axis tight

subplot(2,3,3)
plot(new_t, dir_RRRS,'k', 'LineWidth',2)
hold on
plot(new_t, lpos)
hold off
legend('Right Reward to Right Start Vector','Location', 'northoutside')
axis tight

subplot(2,3,4)
plot(new_t, dir_RRLS,'k', 'LineWidth',2)
hold on
plot(new_t, lpos)
hold off
legend('Right Reward to Left Start Vector','Location', 'northoutside')
axis tight

subplot(2,3,5)
plot(new_t, dir_LRLS,'k', 'LineWidth',2)
hold on
plot(new_t, lpos)
hold off
legend('Left Reward to Left Start Vector','Location', 'northoutside')
axis tight

subplot(2,3,6)
plot(new_t, dir_LRRS,'k', 'LineWidth',2)
hold on
plot(new_t, lpos)
hold off
legend('Left Reward to Right Start Vector','Location', 'northoutside')
axis tight
%%
% close all
w_data = zeros(length(spiketrain),1);
for qq = 1:length(w_data)
    if linear_pos(qq) > 0.1 && linear_pos(qq) < 0.2
        w_data(qq) = spiketrain2(qq) + 1000;
    else
        w_data(qq) = 1;
    end
end
predictor = [lpos lpos.^2 dir_RSRR dir_LSLR dir_RRRS dir_RRLS ...
     dir_LRLS dir_LRRS];
[b6,dev6,stats6]=glmfit(predictor,spiketrain2,'poisson', 'weights', w_data);


figure
bar(pos_bins,spikes_per_sec); %Plot occupancy norm. hist.
hold on;   
yfit6 = glmval(b6,predictor,'log'); 
yfit6 = yfit6./new_SR;  
yfitInd6 = find(yfit6>=0);  
[N,edges] = histcounts(yfit6, pos_bins); 
plot(lpos, yfit6,'r')  
hold off              
xlabel('Position [m]')		
ylabel('Occupancy Norm. (spikes/s)')    
axis([0 1.30 0 inf]) 


figure
plot(new_t, yfit6)
hold on
plot(new_t,lpos/10, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
plot(new_t(spikeindex2), lpos(spikeindex2)/10,'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 5)
plot(new_t, lpos./0.1)
hold off


xs=((0:.875:130.5)./100)';					%Define interval of positions,
Ns = [size(xs,1),1];

%Evaluate Model 4 in direction 0 (X decreases).
% dir_RSRR,dir_LSLR,dir_RRRS,dir_RRLS,dir_LRLS,dir_LRRS
[lambda6_0, up0, low0]=glmval(b6,[xs xs.^2 ones(Ns) zeros(Ns) ...
    zeros(Ns) zeros(Ns) zeros(Ns) zeros(Ns)],'log',stats6);
[lambda6_1, up1, low1]=glmval(b6,[xs xs.^2 zeros(Ns) ones(Ns) ...
    zeros(Ns) zeros(Ns) zeros(Ns) zeros(Ns)],'log',stats6);
[lambda6_2, up2, low2]=glmval(b6,[xs xs.^2 zeros(Ns) zeros(Ns) ...
    ones(Ns) zeros(Ns) zeros(Ns) zeros(Ns)],'log',stats6);
[lambda6_3, up3, low3]=glmval(b6,[xs xs.^2 zeros(Ns) zeros(Ns) ...
    zeros(Ns) ones(Ns) zeros(Ns) zeros(Ns)],'log',stats6);
[lambda6_4, up4, low4]=glmval(b6,[xs xs.^2 zeros(Ns) zeros(Ns) ...
    zeros(Ns) zeros(Ns) ones(Ns) zeros(Ns)],'log',stats6);
[lambda6_5, up5, low5]=glmval(b6,[xs xs.^2 zeros(Ns) zeros(Ns) ...
    zeros(Ns) zeros(Ns) zeros(Ns) ones(Ns)],'log',stats6);

figure
subplot(1,2,1)
bar(pos_bins, spikes_per_sec)
hold on	
hhh1 = plot(xs,lambda6_0./new_SR,'b','LineWidth',2); %Plot Model 4, X decreasing
plot(xs,(lambda6_0./new_SR)+up0./new_SR,'b--')	%...add upper CI,
plot(xs,(lambda6_0./new_SR)-low0./new_SR,'b--')	%...and lower CI.
hhh2 = plot(xs,lambda6_1./new_SR,'r','LineWidth',2); %Plot Model 4
plot(xs,(lambda6_1./new_SR)+up1./new_SR,'r--')	%...add upper CI,
plot(xs,(lambda6_1./new_SR)-low1./new_SR,'r--');	%...add lower CI.
hold off
legend([hhh1 hhh2],'RS-RR','LS-LR')
title('Choice Routes')
axis([0 1.30 0 inf])
subplot(1,2,2)
bar(pos_bins, spikes_per_sec)
hold on
hhh3 = plot(xs,lambda6_2./new_SR,'k','LineWidth',2); %Plot Model 4
plot(xs,lambda6_2./new_SR+up2./new_SR,'k--')	%...add upper CI,
plot(xs,lambda6_2./new_SR-low2./new_SR,'k--');	%...add lower CI.
hhh4 = plot(xs,lambda6_3./new_SR,'g','LineWidth',2); %Plot Model 4
plot(xs,lambda6_3./new_SR+up3./new_SR,'g--')	%...add upper CI,
plot(xs,lambda6_3./new_SR-low3./new_SR,'g--');	%...add lower CI.
hhh5 = plot(xs,lambda6_4./new_SR,'c','LineWidth',2); %Plot Model 4
plot(xs,lambda6_4./new_SR+up4./new_SR,'c--')	%...add upper CI,
plot(xs,lambda6_4./new_SR-low4./new_SR,'c--');	%...add lower CI.
hhh6 = plot(xs,lambda6_5./new_SR,'m','LineWidth',2); %Plot Model 4
plot(xs,lambda6_5./new_SR+up5./new_SR,'m--')	%...add upper CI,
plot(xs,lambda6_5./new_SR-low5./new_SR,'m--');	%...add lower CI.
legend([hhh3 hhh4 hhh5 hhh6],'RR-RS','RR-LS','LR-LS','LR-RS')
title('Forced Routes')
axis([0 1.30 0 inf])
hold off

%% Spike Histogram based on


%% Video play to separate trials

implay('vGAT_CreCh2_4_ForcedD32018-04-13T18_09_04.avi')