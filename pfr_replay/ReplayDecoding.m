
%% MUA activity not for decoding just for detecting regions of high activity
NewSr = 1/2000;
MUAvect2 = sum(MUAvect, 2);
sigma = 0.01; % Standard deviation of the kernel 15ms
edges = -3*sigma: NewSr : 3*sigma; % Evaluate ranges form -3*sd to 3*sd dev
kernel = normpdf(edges,0,sigma); % Evaluate the Gaussian kernel
MUAdensity = conv(MUAvect2, kernel); % Convolved MUA with the kernel
MUAdensity = MUAdensity(1:length(MUAdensity)-length(edges)+1);
new_tt = (0 : NewSr : length(MUAdensity) * NewSr - NewSr);
plot(edges,kernel)
%% Visual inspection of MUA activity
figure
subplot(411)
plot(linTime, linPos,'k', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('Position [m]')
axis([0 20 0 1.4])
subplot(412)
plot(linTime, linPos*12,'Color',[.7 .7 .7], 'LineWidth', 2)
hold on
plot(linTime(MUAidx{1}), MUAvect(MUAidx{1},1), '.K');
for kk = 2: size(MUAvect,2)
    plot(linTime(MUAidx{kk}), MUAvect(MUAidx{kk},kk)+kk-1, '.K');
end
hold off
xlabel('Time [sec]')
ylabel('Unit #')
axis([0 20 0 size(MUAvect,2)+3])
subplot(413)
plot(new_tt, MUAdensity,'b', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('MUA [Hz]')
axis([0 20 0 inf])
subplot(414)
plot(linTime, linVel,'k', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('m/sec')
axis([0 20 0 inf])
%% Finding MUA "reactivation"
threshold = 3.5 * mean(median(MUAdensity ./ .6745), 2);
MUAmax = MUAdensity >= threshold;
MUAmax = find(MUAmax == 1);
difMUAmax = find(diff(MUAmax) >= (0.02/NewSr)); 
z = zeros(length(difMUAmax)+1,1); 
z(1:end-1) = difMUAmax;
z(end) = length(MUAmax); 
st = difMUAmax+1; 
zst = zeros(length(difMUAmax)+1,1);
zst(2:end) = st;
zst(1)=1; 
MUAsteps = [MUAmax(zst), MUAmax(z)];  
duration = MUAsteps(:,2)- MUAsteps(:,1);
MUAsteps1 = MUAsteps((~(duration<=1)),:);
duration2 = (MUAsteps1(:,2) - MUAsteps1(:,1))*NewSr; %duration in sec.
MUAsidx = MUAsteps1((find(duration2<=0.550 & ...
    duration2 >= 0.01)),:); %Events < 150ms and > 15ms (Indices)
MUAsec = MUAsidx*NewSr;
below = MUAsidx(:,1)-(0.25/NewSr);
abov = MUAsidx(:,2)+(0.25/NewSr);
below(below < 0) = 1; % Substitute start negative value with 1
abov(abov > length(LFPt)) = length(LFPt); % Making MUA candidates finishing at the end of the LFP
belowidx = find(linVel(below) < 0.025);
below = below(belowidx);
abov = abov(belowidx);
% MUAregions = [below abov];
% MUAregionsSec = (MUAregions(:,2) .*NewSr) - (MUAregions(:,1) .*NewSr);
%% Visual Inspection of MUA "reactivation" detection
subplot(211)
plot(linTime,linPos)
axis([0 100 0 inf])
subplot(212)
plot(new_tt, MUAdensity)
hold on
for kk = 1: length(below)
    plot(new_tt(below(kk):abov(kk)),MUAdensity(below(kk):abov(kk)), 'r', 'LineWidth',2) 
end
% plot(MUAdensityt, filtMUAdensity,'k','LineWidth',1)
% plot(MUAdensityt, smsqfiltMUAdensity./80,'m', 'LineWidth',2)
hold off
axis([0 100 0 inf])
%% Creating a matrix with spikes, poistion, velocity and LFP based on MUA
MUAreact = cell(1,size(below,1));
for jj = 1:size(below,1)
        ind1 = below(jj,1);
        ind2 = abov(jj,1);
        MUAreact{1,jj}(:,1) = linTime(ind1: ind2); % Time
        MUAreact{1,jj}(:,2) = linPos(ind1: ind2); % Position
        MUAreact{1,jj}(:,3) = linVel(ind1: ind2); % Velocity
        MUAreact{1,jj}(:,4) = LFP15(ind1: ind2); %LFP
        MUAreact{1,jj}(:,5) = swrLFP15(ind1: ind2); %Theta filtered LFP
        MUAreact{1,jj}(:,6) = LFP2(ind1: ind2); %LFP
        MUAreact{1,jj}(:,7) = swrLFP2(ind1: ind2); %Theta filtered LFP
        for kk = 1: size(MUAvect,2) % Spikes(0 or 1)
            MUAreact{1,jj}(:,7+kk) = MUAvect(ind1:ind2,kk);
        end
end
%% Finding Spikes in the MUA
MUAreactIdx = cell (1, length(below));
for kk = 1: length(below)
    for jj = 1: size(MUAvect,2)
    MUAreactIdx{1, kk}{1,jj} = find(MUAreact{1,kk}(:, 7+jj));
    end
end
%% Visual Inspection
close all
for kk = 1: length(below)
    figure(kk)
    subplot(411)
    yyaxis right
    plot(MUAreact{1,kk}(:,1), MUAreact{1,kk}(:,2),'Color',[.7 .7 .7], 'LineWidth', 2)
    hold on
    yyaxis left
        for jj = 1: size(MUAvect,2)
            plot(MUAreact{1,kk}(MUAreactIdx{1,kk}{1,jj},1), ...
                MUAreact{1,kk}(MUAreactIdx{1,kk}{1,jj},jj+7)+jj-1, '.K');
        end
    hold off
    yyaxis left
    xlabel('Time [sec]')
    ylabel('Unit #')
    axis([-inf inf 0 size(MUAvect,2)+3])
    yyaxis right
    
    ylabel('Position[m]')
    xlabel('Time [sec]')
    axis([-inf inf 0 1.3])
    subplot(412)
    plot(new_tt(below(kk):abov(kk)), MUAdensity(below(kk):abov(kk)),'b', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('MUA [Hz]')
    axis([-inf inf 0 600])
    subplot(413)
    plot(linTime(below(kk):abov(kk)), LFP15(below(kk):abov(kk)),'b')
    hold on
    plot(linTime(below(kk):abov(kk)), swrLFP15(below(kk):abov(kk))-300,'k')
    hold off
    xlabel('Time [sec]')
    ylabel('uV')
    axis([-inf inf -inf inf])
    subplot(414)
    plot(linTime(below(kk):abov(kk)), LFP2(below(kk):abov(kk)),'b')
    hold on
    plot(linTime(below(kk):abov(kk)), swrLFP2(below(kk):abov(kk))-300,'k')
    hold off
    xlabel('Time [sec]')
    ylabel('uV')
    axis([-inf inf -inf inf])
end
%% Creating Time bins
BinTime = 0.02;
Timebins = cell(1,size(MUAreact,2));
TimeTrialbins = cell(1,size(MUAreact,2));
for kk = 1: size(MUAreactIdx,2)
    % Timebins to bin Spikes over time
    Timebins{1,kk} = (0: BinTime: (MUAreact{1,kk}(end,1) - MUAreact{1,kk}(1,1)))'; 
    % Timebins for bin Spikes over time/Trials
    TimeTrialbins{1,kk} = (MUAreact{1,kk}(1,1): BinTime : MUAreact{1,kk}(end,1)); 
end

%% Histogram Spikes over time/trials
MUAreactTimeHist = cell(1, size(MUAreactIdx,2));
for jj = 1: size(MUAreactIdx,2)
    for kk = 1: size(MUAreactIdx{1,1},2)
     MUAreactTimeHist{1,jj}(:,kk) =(hist(MUAreact{1,jj}(MUAreactIdx{1,jj}{1,kk},1), ...
         TimeTrialbins{1,jj}))';
    end
end
%% Tuning Curves
SumFRExp = exp(-(BinTime).* sum(rrTuningCurveSm,2)); 
SumFRExp1 = exp(-(BinTime).* sum(llTuningCurveSm,2));
% SumFRExp2 = exp(-(BinTime).* sum(frrrTuningCurveSm,2));
% SumFRExp3 = exp(-(BinTime).* sum(frllTuningCurveSm,2));
% SumFRExp4 = exp(-(BinTime).* sum(frrlTuningCurveSm,2));
% SumFRExp5 = exp(-(BinTime).* sum(frlrTuningCurveSm,2));

%% Firing rate elevated by time spikes
% Replay decoding using Choice RR tuning Curve
MUAreact_eSpikesRR = cell(1, size(MUAreactIdx,2));
for jj = 1: size(MUAreactIdx,2)
    for tt = 1: length(MUAreactTimeHist{1,jj}(:,1))
        MUAreact_eSpikesRR{1, jj}(tt,:) = (prod(rrTuningCurveSm(:,:).^ ...
            MUAreactTimeHist{1,jj}(tt,:),2))';
    end
end
% Replay decoding using Choice LL tuning Curve
MUAreact_eSpikesLL = cell(1, size(MUAreactIdx,2));
for jj = 1: size(MUAreactIdx,2)
    for tt = 1: length(MUAreactTimeHist{1,jj}(:,1))
        MUAreact_eSpikesLL{1, jj}(tt,:) = (prod(llTuningCurveSm(:,:).^ ...
            MUAreactTimeHist{1,jj}(tt,:),2))';
    end
end
% % Replay decoding using Forced RR tuning Curve
% MUAreact_eSpikesFrRR = cell(1, size(MUAreactIdx,2));
% for jj = 1: size(MUAreactIdx,2)
%     for tt = 1: length(MUAreactTimeHist{1,jj}(:,1))
%         MUAreact_eSpikesFrRR{1, jj}(tt,:) = (prod(frrrTuningCurveSm(:,:).^ ...
%             MUAreactTimeHist{1,jj}(tt,:),2))';
%     end
% end
% % Replay decoding using Forced LL tuning Curve
% MUAreact_eSpikesFrLL = cell(1, size(MUAreactIdx,2));
% for jj = 1: size(MUAreactIdx,2)
%     for tt = 1: length(MUAreactTimeHist{1,jj}(:,1))
%         MUAreact_eSpikesFrLL{1, jj}(tt,:) = (prod(frllTuningCurveSm(:,:).^ ...
%             MUAreactTimeHist{1,jj}(tt,:),2))';
%     end
% end
% % Replay decoding using Forced RL tuning Curve
% MUAreact_eSpikesFrRL = cell(1, size(MUAreactIdx,2));
% for jj = 1: size(MUAreactIdx,2)
%     for tt = 1: length(MUAreactTimeHist{1,jj}(:,1))
%         MUAreact_eSpikesFrRL{1, jj}(tt,:) = (prod(frrlTuningCurveSm(:,:).^ ...
%             MUAreactTimeHist{1,jj}(tt,:),2))';
%     end
% end
% % Replay decoding using Forced LR tuning Curve
% MUAreact_eSpikesFrLR = cell(1, size(MUAreactIdx,2));
% for jj = 1: size(MUAreactIdx,2)
%     for tt = 1: length(MUAreactTimeHist{1,jj}(:,1))
%         MUAreact_eSpikesFrLR{1, jj}(tt,:) = (prod(frlrTuningCurveSm(:,:).^ ...
%             MUAreactTimeHist{1,jj}(tt,:),2))';
%     end
% end
%% Likelihood of Position over Time
MUAreactProbRR = cell(1, size(MUAreactIdx,2));
SumMUAreactProbRR = cell(1, size(MUAreactIdx,2));
MUAreactLikelihoodRR = cell(1, size(MUAreactIdx,2));
for kk = 1: size(MUAreactIdx,2)
    MUAreact_eSpikesRR{1,kk} = (MUAreact_eSpikesRR{1,kk})';
    MUAreactProbRR{1,kk} = SumFRExp.* MUAreact_eSpikesRR{1,kk};
    SumMUAreactProbRR{1,kk} = sum(MUAreactProbRR{1,kk},1);
    MUAreactLikelihoodRR{1,kk} = MUAreactProbRR{1,kk}./SumMUAreactProbRR{1,kk};
end
%% Visual Inspection
close all
for kk = 1: size(MUAreactIdx,2)
    figure(kk)
    hold on
    s = plot3(MUAreact{1,kk}(:,1), MUAreact{1,kk}(:,2), ...
        ones(size(MUAreact{1,kk}(:,2),1),1),'-.r','LineWidth', 0.5);
    alpha(s,0.1)
    s1 = surf(TimeTrialbins{1,kk}, tuningbins,MUAreactLikelihoodRR{1,kk});
    s1.EdgeColor = 'none';
    c = colorbar;
    c.Label.String = 'Probability';
%     colormap(flipud(hot))
    colormap(flipud(hot))
    % caxis auto
    axis([TimeTrialbins{1,kk}(1) TimeTrialbins{1,kk}(end) 0 1.3])
    caxis([0.0 .25])
    xlabel('Time [sec]')
    ylabel('Linear Position [m]')
    zlabel('Probability') 
    legend(s,'Actual Position','Location','northwest')
    box on
    hold off
end
%% Likelihood of Position Choice Left to Left
MUAreactProbLL = cell(1, size(MUAreactIdx,2));
SumMUAreactProbLL = cell(1, size(MUAreactIdx,2));
MUAreactLikelihoodLL = cell(1, size(MUAreactIdx,2));
for kk = 1: size(MUAreactIdx,2)
    MUAreact_eSpikesLL{1,kk} = (MUAreact_eSpikesLL{1,kk})';
    MUAreactProbLL{1,kk} = SumFRExp1.* MUAreact_eSpikesLL{1,kk};
    SumMUAreactProbLL{1,kk} = sum(MUAreactProbLL{1,kk},1);
    MUAreactLikelihoodLL{1,kk} = MUAreactProbLL{1,kk}./SumMUAreactProbLL{1,kk};
end

%% Visual Inspection
close all
for kk = 1: size(MUAreactIdx,2)
    figure(kk)
    hold on
    s = plot3(MUAreact{1,kk}(:,1), MUAreact{1,kk}(:,2), ...
        ones(size(MUAreact{1,kk}(:,2),1),1),'-.r','LineWidth', 0.5);
    alpha(s,0.1)
    s1 = surf(TimeTrialbins{1,kk}, tuningbins,MUAreactLikelihoodLL{1,kk});
    s1.EdgeColor = 'none';
    c = colorbar;
    c.Label.String = 'Probability';
    colormap(flipud(hot))
    % caxis auto
    axis([TimeTrialbins{1,kk}(1) TimeTrialbins{1,kk}(end) 0 1.3])
    caxis([0.0 .15])
    xlabel('Time [sec]')
    ylabel('Linear Position [m]')
    zlabel('Probability') 
    legend(s,'Actual Position','Location','northwest')
    box on
    hold off
end

%% Likelihood of Position Forced Right to Right
MUAreactProbFrRR = cell(1, size(MUAreactIdx,2));
SumMUAreactProbFrRR = cell(1, size(MUAreactIdx,2));
MUAreactLikelihoodFrRR = cell(1, size(MUAreactIdx,2));
for kk = 1: size(MUAreactIdx,2)
    MUAreact_eSpikesFrRR{1,kk} = (MUAreact_eSpikesFrRR{1,kk})';
    MUAreactProbFrRR{1,kk} = SumFRExp2.* MUAreact_eSpikesFrRR{1,kk};
    SumMUAreactProbFrRR{1,kk} = sum(MUAreactProbFrRR{1,kk},1);
    MUAreactLikelihoodFrRR{1,kk} = MUAreactProbFrRR{1,kk}./SumMUAreactProbFrRR{1,kk};
end

%% Visual Inspection
close all
for kk = 1: size(MUAreactIdx,2)
    figure(kk)
    hold on
    s = plot3(MUAreact{1,kk}(:,1), MUAreact{1,kk}(:,2), ...
        ones(size(MUAreact{1,kk}(:,2),1),1),'-.r','LineWidth', 0.5);
    alpha(s,0.1)
    s1 = surf(TimeTrialbins{1,kk}, tuningbins,MUAreactLikelihoodFrRR{1,kk});
    s1.EdgeColor = 'none';
    c = colorbar;
    c.Label.String = 'Probability';
    colormap(flipud(hot))
    % caxis auto
    axis([TimeTrialbins{1,kk}(1) TimeTrialbins{1,kk}(end) 0 1.3])
    caxis([0.0 .15])
    xlabel('Time [sec]')
    ylabel('Linear Position [m]')
    zlabel('Probability') 
    legend(s,'Actual Position','Location','northwest')
    box on
    hold off
end
%% Forced Left to Left
MUAreactProbFrLL = cell(1, size(MUAreactIdx,2));
SumMUAreactProbFrLL = cell(1, size(MUAreactIdx,2));
MUAreactLikelihoodFrLL = cell(1, size(MUAreactIdx,2));
for kk = 1: size(MUAreactIdx,2)
    MUAreact_eSpikesFrLL{1,kk} = (MUAreact_eSpikesFrLL{1,kk})';
    MUAreactProbFrLL{1,kk} = SumFRExp3.* MUAreact_eSpikesFrLL{1,kk};
    SumMUAreactProbFrLL{1,kk} = sum(MUAreactProbFrLL{1,kk},1);
    MUAreactLikelihoodFrLL{1,kk} = MUAreactProbFrLL{1,kk}./SumMUAreactProbFrLL{1,kk};
end

%% Visual Inspection
close all
for kk = 1: size(MUAreactIdx,2)
    figure(kk)
    hold on
    s = plot3(MUAreact{1,kk}(:,1), MUAreact{1,kk}(:,2), ...
        ones(size(MUAreact{1,kk}(:,2),1),1),'-.r','LineWidth', 0.5);
    alpha(s,0.1)
    s1 = surf(TimeTrialbins{1,kk}, tuningbins,MUAreactLikelihoodFrLL{1,kk});
    s1.EdgeColor = 'none';
    c = colorbar;
    c.Label.String = 'Probability';
    colormap(flipud(hot))
    % caxis auto
    axis([TimeTrialbins{1,kk}(1) TimeTrialbins{1,kk}(end) 0 1.3])
    caxis([0.005 .15])
    xlabel('Time [sec]')
    ylabel('Linear Position [m]')
    zlabel('Probability') 
    legend(s,'Actual Position','Location','northwest')
    box on
    hold off
end

%% Forced Right to Left
MUAreactProbFrRL = cell(1, size(MUAreactIdx,2));
SumMUAreactProbFrRL = cell(1, size(MUAreactIdx,2));
MUAreactLikelihoodFrRL = cell(1, size(MUAreactIdx,2));
for kk = 1: size(MUAreactIdx,2)
    MUAreact_eSpikesFrRL{1,kk} = (MUAreact_eSpikesFrRL{1,kk})';
    MUAreactProbFrRL{1,kk} = SumFRExp4.* MUAreact_eSpikesFrRL{1,kk};
    SumMUAreactProbFrRL{1,kk} = sum(MUAreactProbFrRL{1,kk},1);
    MUAreactLikelihoodFrRL{1,kk} = MUAreactProbFrRL{1,kk}./SumMUAreactProbFrRL{1,kk};
end

%% Visual Inspection
close all
for kk = 1: size(MUAreactIdx,2)
    figure(kk)
    hold on
    s = plot3(MUAreact{1,kk}(:,1), MUAreact{1,kk}(:,2), ...
        ones(size(MUAreact{1,kk}(:,2),1),1),'-.r','LineWidth', 0.5);
    alpha(s,0.1)
    s1 = surf(TimeTrialbins{1,kk}, tuningbins,MUAreactLikelihoodFrRL{1,kk});
    s1.EdgeColor = 'none';
    c = colorbar;
    c.Label.String = 'Probability';
    colormap(flipud(hot))
    % caxis auto
    axis([TimeTrialbins{1,kk}(1) TimeTrialbins{1,kk}(end) 0 1.3])
    caxis([0.0 .15])
    xlabel('Time [sec]')
    ylabel('Linear Position [m]')
    zlabel('Probability') 
    legend(s,'Actual Position','Location','northwest')
    box on
    hold off
end

%% Forced Left TO rIGHT
MUAreactProbFrLR = cell(1, size(MUAreactIdx,2));
SumMUAreactProbFrLR = cell(1, size(MUAreactIdx,2));
MUAreactLikelihoodFrLR = cell(1, size(MUAreactIdx,2));
for kk = 1: size(MUAreactIdx,2)
    MUAreact_eSpikesFrLR{1,kk} = (MUAreact_eSpikesFrLR{1,kk})';
    MUAreactProbFrLR{1,kk} = SumFRExp5.* MUAreact_eSpikesFrLR{1,kk};
    SumMUAreactProbFrLR{1,kk} = sum(MUAreactProbFrLR{1,kk},1);
    MUAreactLikelihoodFrLR{1,kk} = MUAreactProbFrLR{1,kk}./SumMUAreactProbFrLR{1,kk};
end

%% Visual Inspection
close all
for kk = 1: size(MUAreactIdx,2)
    figure(kk)
    hold on
    s = plot3(MUAreact{1,kk}(:,1), MUAreact{1,kk}(:,2), ...
        ones(size(MUAreact{1,kk}(:,2),1),1),'-.r','LineWidth', 0.5);
    alpha(s,0.1)
    s1 = surf(TimeTrialbins{1,kk}, tuningbins,MUAreactLikelihoodFrLR{1,kk});
    s1.EdgeColor = 'none';
    c = colorbar;
    c.Label.String = 'Probability';
    colormap(flipud(hot))
    % caxis auto
    axis([TimeTrialbins{1,kk}(1) TimeTrialbins{1,kk}(end) 0 1.3])
    caxis([0 .15])
    xlabel('Time [sec]')
    ylabel('Linear Position [m]')
    zlabel('Probability') 
    legend(s,'Actual Position','Location','northwest')
    box on
    hold off
end











