NewSr = 1/2000;
%%
semiMUAtrains = MUAvect;
MUAvectActive = zeros(size(semiMUAtrains));
for kk = 1 : size(semiMUAtrains,2)
    for qq = 1:size(semiMUAtrains,1)
        if linVel(qq) > 0.03
        MUAvectActive(qq,kk) = semiMUAtrains(qq,kk);
        else
        MUAvectActive(qq,kk) = 0;
        end
    end
end

%  %% Eliminating Posible Interneurons
% RecTime = LFPt(end);
% for kk = 1: size(semiMUAtrains,2)
%     if size(semiMUAtimes2{kk},1)/RecTime > 5
%        MUAvectActive(:,kk) = [];
%     elseif size(semiMUAtimes2{kk},1)/RecTime < 0.1
%         MUAvectActive(:,kk) = [];
%     else
% 
%     end
% end
%% 
rrFR = cell(1,size(ch_run,1)); 
for jj = 1:size(ch_run,1)
    if ch_run{jj,3} == 'R' && ch_run{jj,4} == 'R'
        ind1 = ch_run{jj,1} * (1/NewSr);
        ind2 = ch_run{jj,2} * (1/NewSr);
        rrFR{1,jj}(:,1) = linTime(ind1: ind2); % Time
        rrFR{1,jj}(:,2) = linPos(ind1: ind2); % Position
        rrFR{1,jj}(:,3) = linVel(ind1: ind2); % Velocity
        rrFR{1,jj}(:,4) = LFP7(ind1: ind2); % LFP
        rrFR{1,jj}(:,5) = thetaLFP7(ind1: ind2); % Theta filtered LFP
        for kk = 1: size(MUAvectActive,2) % Units and MUA spikes
            rrFR{1,jj}(:,5+kk) = MUAvectActive(ind1:ind2,kk);
        end
    end
end
% Removing [] cell arrays
rrFR = rrFR(~cellfun(@isempty,rrFR));

%% Finding Spikes index
rrFRidx = cell(1, size(rrFR,2));
for jj = 1 : size(rrFR,2)
    for kk = 1 : size(MUAvectActive,2)
        rrFRidx{1,jj}{:,kk} = find(rrFR{1,jj}(:, kk+5));%Column 5 and up are spikes
    end
end

%% Creating Encoding Model
EncTrials = [ 4 5 6 7 ]; % Training Data
binFac = 0.04;
tuningbins = (0.0:binFac:1.3)';
rrTuning = cell(1, size(MUAvectActive,2));
for kk = 1 : size(EncTrials,2)
   for jj = 1 : size(MUAvectActive,2)
       rrTuning{1,jj}(:,kk) = (hist(rrFR{1,kk}(rrFRidx{1,kk}{1,jj},2),tuningbins))';
   end
end
rrTuningSum = cell(1, size(rrTuning,2));
for kk = 1: size(rrTuning,2)
    rrTuningSum{1,kk} = sum(rrTuning{1,kk},2);
end
rrOccup = zeros(size(tuningbins,1), size(rrFR,2));
for kk = 1: size(rrFR,2)
   rrOccup(:,kk) = (hist(rrFR{1,kk}(:,2),tuningbins).*NewSr)';
end
rrOccup = sum(rrOccup,2);
rrTuningCurve = zeros(size(tuningbins,1), size(rrTuningSum,2));
for kk = 1: size(rrTuningSum,2)
    rrTuningCurve(:,kk) = rrTuningSum{1,kk}./rrOccup;
end
rrTuningCurve(isnan(rrTuningCurve)) = 0;
rrTuningCurveSm = zeros((size(rrTuningCurve,1)), size(rrTuningCurve,2));
for kk = 1: size(rrTuningCurve,2)
   rrTuningCurveSm(:,kk) = (smoothdata(rrTuningCurve(:,kk),'gaussian',6))+.01; % Adding 0.01Hz
end
rrnSR = 1.3/(length(rrTuningCurve(:,1))+10);
rrPosbins = (0 : rrnSR : 1.3)';
close all
for kk = 1: size(rrTuningCurve,2)
   s(kk) = subplot(size(MUAvectActive,2),1,kk);
   hold on
   area(tuningbins, rrTuningCurveSm(:,kk),'FaceColor','k')
%     line('XData', [.9 .9], 'YData', [0 55], 'LineStyle', '-.', ...
%     'LineWidth', 2, 'Color',[0.5 0.5 0.5])
%     line('XData', [.35 .35], 'YData', [0 55], 'LineStyle', '-.', ...
%     'LineWidth', 2, 'Color',[0.5 0.5 0.5])
    hold off
   axis([0.0 1.3 0 inf])
end
xlabel('Position [m]')
ylabel('FR[Hz]')
title(s(1),'Forced/Choice Right to Right','FontSize',12)
txt1 = 'Start ';
text(.2,115,txt1,'HorizontalAlignment','Center','FontSize', 8)
txt2 = 'Center ';
text(.6,115,txt2,'HorizontalAlignment','Center','FontSize', 8)
txt3 = 'Reward ';
text(1.1,115,txt3,'HorizontalAlignment','Center','FontSize', 8)
txt4 = '\rightarrow';
text(.1,240,txt4,'HorizontalAlignment','Center','FontSize', 40)
% axis(s(1),[0 1.3 0 25])
% axis(s(2),[0 1.3 0 30])
% axis(s(3),[0 1.3 0 53])
% axis(s(4),[0 1.3 0 20])
% axis(s(5),[0 1.3 0 10])
% axis(s(6),[0 1.3 0 25])
% axis(s(7),[0 1.3 0 20])
% axis(s(8),[0 1.3 0 10])
% axis(s(9),[0 1.3 0 50])
% axis(s(10),[0 1.3 0 45])
% axis(s(11),[0 1.3 0 20])
% axis(s(12),[0 1.3 0 30])
% axis(s(13),[0 1.3 0 45])
% axis(s(14),[0 1.3 0 10])
%% Position Time histogram
BinTime = 0.2;
Timebins = cell(1,size(rrFR,2));
TimeTrialbins = cell(1,size(rrFR,2));
for kk = 1: size(rrFRidx,2)
    % Timebins to bin Spikes over time
    Timebins{1,kk} = (0: BinTime: (rrFR{1,kk}(end,1) - rrFR{1,kk}(1,1)))'; 
    % Timebins for bin Spikes over time/Trials
    TimeTrialbins{1,kk} = (rrFR{1,kk}(1,1): BinTime : rrFR{1,kk}(end,1)); 
end
%% Testing Data/Spikes
rrFR2 = cell(1,size(ch_run,1)); 
for jj = 1:size(ch_run,1)
    if ch_run{jj,3} == 'R' && ch_run{jj,4} == 'R'
        ind1 = ch_run{jj,1} * (1/NewSr);
        ind2 = ch_run{jj,2} * (1/NewSr);
        rrFR2{1,jj}(:,1) = linTime(ind1: ind2); % Time
        rrFR2{1,jj}(:,2) = linPos(ind1: ind2); % Position
        rrFR2{1,jj}(:,3) = linVel(ind1: ind2); % Velocity
        rrFR2{1,jj}(:,4) = LFP7(ind1: ind2); % LFP
        rrFR2{1,jj}(:,5) = thetaLFP7(ind1: ind2); % Theta filtered LFP
        for kk = 1: size(MUAvectActive,2) % Units and MUA spikes
            rrFR2{1,jj}(:,5+kk) = semiMUAtrains(ind1:ind2,kk);
        end
    end
end
% Removing [] cell arrays
rrFR2 = rrFR2(~cellfun(@isempty,rrFR2));
rrFRidx2 = cell(1, size(rrFR2,2));
for jj = 1 : size(rrFR2,2)
    for kk = 1 : size(semiMUAtrains,2)
        rrFRidx2{1,jj}{:,kk} = find(rrFR2{1,jj}(:, kk+5));%Column 5 and up are spikes
    end
end
%% Histogram Spikes over time/trials
rrTimeHist = cell(1, size(rrFRidx2,2));
for jj = 1: size(rrFRidx2,2)
    for kk = 1: size(rrFRidx2{1,1},2)
     rrTimeHist{1,jj}(:,kk) =(hist(rrFR{1,jj}(rrFRidx2{1,jj}{1,kk},1), ...
         TimeTrialbins{1,jj}))';
    end
end

%% Summing the firing rates per position
rrSumFRate = (sum(rrTuningCurveSm,2));
rrSumFRExp = exp(-(BinTime).* sum(rrTuningCurveSm,2));-0.2.*(exp(rrSumFRate));
figure
% subplot(2,1,1)
plot(tuningbins, rrSumFRate, '-k','LineWidth', 2)
xlabel('Position [m]')
ylabel('\Sigma FR[Hz]')
axis([0.0 1.3 0 inf])
title('Population Firing Rate')
% subplot(2,1,2)
% plot(tuningbins, rrSumFRExp)
% xlabel('Position [m]')
% ylabel('e^{-\Deltat \Sigma FR[Hz]}')
% axis([0.0 1.2 0 inf])
%% Firing rate elevated by time spikes
rrFR_eSpikes = cell(1, size(rrFRidx,2));
for jj = 1: size(rrFRidx,2)
    for tt = 1: length(rrTimeHist{1,jj}(:,1))
        rrFR_eSpikes{1, jj}(tt,:) = (prod(rrTuningCurveSm(:,:).^ ...
            rrTimeHist{1,jj}(tt,:),2))';
    end
end

%%
rrProb = cell(1, size(rrFRidx,2));
SumrrProb = cell(1, size(rrFRidx,2));
rrLikelihood = cell(1, size(rrFRidx,2));
for kk = 1: size(rrFRidx,2)
    rrFR_eSpikes{1,kk} = (rrFR_eSpikes{1,kk})';
    rrProb{1,kk} = rrSumFRExp.* rrFR_eSpikes{1,kk};
    SumrrProb{1,kk} = sum(rrProb{1,kk},1);
    rrLikelihood{1,kk} = rrProb{1,kk}./SumrrProb{1,kk};
end

%% Visual Inspection
close all
for kk = 1: size(rrFRidx,2)
    figure(kk)
    hold on
    s = plot3(rrFR{1,kk}(:,1), rrFR{1,kk}(:,2), ...
        ones(size(rrFR{1,kk}(:,2),1),1),'-.r','LineWidth', 0.5);
    alpha(s,0.1)
    s1 = surf(TimeTrialbins{1,kk}, tuningbins,rrLikelihood{1,kk});
    s1.EdgeColor = 'none';
    c = colorbar;
    c.Label.String = 'Probability';
    colormap(flipud(gray))
    % caxis auto
    axis([TimeTrialbins{1,kk}(1) TimeTrialbins{1,kk}(end) 0.0 1.3])
    caxis([0 1])
    xlabel('Time [sec]')
    ylabel('Linear Position [m]')
    zlabel('Probability') 
    legend(s,'Actual Position','Location','northwest')
    box on
    hold off
end
%%
rrIndLikeli = cell(1, size(rrFRidx,2));
for kk = 1 : 3%size(rrFRidx,2)
    [MaxLikeli, rrIndLikeli{kk}] = max(rrLikelihood{kk});
end
rrMaxLikeliPos = cell(1, size(rrFRidx,2));
for kk = 1:3%size(rrFRidx,2)
    rrMaxLikeliPos{kk} = tuningbins(rrIndLikeli{kk},1);
end
close all
for kk = 1 :3% size(rrFRidx,2)
    figure
    plot(TimeTrialbins{kk}, rrMaxLikeliPos{kk}, 'o')
    hold on
    plot(rrFR{1,kk}(:,1),rrFR{1,kk}(:,2))
    hold off
end
%% Position histogram over time
rrPosTimeHist = cell(1,size(rrFRidx,2)); % Memory Preallocation
rrVelTimeHist = cell(1,size(rrFRidx,2));
for jj = 1 : 3%size(rrFRidx,2)
    [a ,b] = hist(rrFR{1,jj}(1:Timebins{jj}(1+1)*2000, 2),1); % Fist bin hist
    [a1 ,b1] = hist(rrFR{1,jj}(1:Timebins{jj}(1+1)*2000, 3),1); % Fist bin hist
    rrPosTimeHist{jj}(1,1) = b;
    rrVelTimeHist{jj}(1,1) = b1;
    for kk = 2 : length(Timebins{jj})-1
        [a ,b] = hist(rrFR{1,jj}(Timebins{jj}(kk)*2000:Timebins{jj}(kk+1)*2000,2),1);
        rrPosTimeHist{jj}(kk,1) = b;
        [a1 ,b1] = hist(rrFR{1,jj}(Timebins{jj}(kk)*2000:Timebins{jj}(kk+1)*2000,3),1);
        rrVelTimeHist{jj}(kk,1) = b1;
    end
    [a ,b] = hist(rrFR{1,jj}(Timebins{jj}(size(TimeTrialbins{jj},2))*2000:Timebins{jj}(end)*2000, 2),1); % Fist bin hist
    rrPosTimeHist{jj}(size(TimeTrialbins{jj},2),1) = b;
    [a1 ,b1] = hist(rrFR{1,jj}(Timebins{jj}(size(TimeTrialbins{jj},2))*2000:Timebins{jj}(end)*2000, 3),1); % Fist bin hist
    rrVelTimeHist{jj}(size(TimeTrialbins{jj},2),1) = b1;
end
close all
for kk = 1 : 3%size(rrFRidx,2)
    figure
    plot(TimeTrialbins{kk}, rrMaxLikeliPos{kk}, 'ob')
    hold on
    plot(rrFR{1,kk}(:,1),rrFR{1,kk}(:,2))
    plot(TimeTrialbins{kk}, rrPosTimeHist{kk}, '-r')
    plot(TimeTrialbins{kk}, rrVelTimeHist{kk}, '-m')
    hold off
end
%% Calculating Decoding Error
rrLikeliError = cell(1, size(rrFRidx,2));
for kk = 1 : 3%size(rrFRidx,2)
    rrLikeliError{kk} = abs(rrPosTimeHist{kk}- rrMaxLikeliPos{kk});
end

close all
for kk = 1 : 3%size(rrFRidx,2)
    figure
    plot(TimeTrialbins{kk}, rrLikeliError{kk}, '-k', 'LineWidth',2)
    xlabel('Time [sec]')
    ylabel('Error [m]')
    axis tight
end
%% Eliminating non running samples
for kk = 1 : 3%size(rrFRidx,2)
    for jj = 1: size(TimeTrialbins{kk},2)
        if rrVelTimeHist{kk}(jj,1) > 0.03
            rrPosTimeHist1{kk}(jj,1) = rrPosTimeHist{kk}(jj,1);
            rrMaxLikeliPos1{kk}(jj,1) = rrMaxLikeliPos{kk}(jj,1);
        else
             rrPosTimeHist1{kk}(jj,1) = NaN;
             rrMaxLikeliPos1{kk}(jj,1) = NaN;
        end
    end
end
%% Calculating Decoding Error from running
rrLikeliError1 = cell(1, 3);
for kk = 1 :3 %size(rrFRidx,2)
    rrLikeliError1{kk} = abs(rrPosTimeHist1{kk}- rrMaxLikeliPos1{kk});
end

close all
for kk = 1 :3 %size(rrFRidx,2)
    figure
    plot(TimeTrialbins{kk}, rrLikeliError1{kk}, '-b')
end
%%
close all
for kk = 1 :3 %size(rrFRidx,2)
    figure
    h = histogram(rrLikeliError1{kk},1000,'Normalization','cdf');
    ErrorCum = (h.Values)';
    cumbins = h.BinEdges;
    MaxError = max(ErrorCum);
    rrErrorCum(:,kk) = ErrorCum./MaxError;
end
%%
close all
rrAvgCumError = zeros(size(rrErrorCum,1)+1,1);
rrAvgCumError(2:end) = mean(rrErrorCum,2);
cumbins(1,1) = 0;
figure
plot(cumbins,rrAvgCumError,'-k','LineWidth',2)
xlabel('Linear Position [m]')
ylabel('Error Cum. Prob')
axis([0 1.2 0 1])
%%
close all
for kk = 1:3 %size(rrFRidx,2)
    labels = unique(rrPosTimeHist{kk});
    C{kk} = confusionmat(rrPosTimeHist{kk},rrMaxLikeliPos{kk});
%     figure
% %     subplot(211)
%     heatmap(C{kk});
%     axis([0 1.3 0 1.3])
%     subplot(212)
    figure
    plot(rrPosTimeHist{kk},rrMaxLikeliPos{kk},'sk')
    axis([0 1.3 0 1.3])
%     figure
%     plotConfMat(rrPosTimeHist{kk},rrMaxLikeliPos{kk})
end
%% Visualization
close all
for i=1:size(rrFR,2)
    figure(i)
    subplot(411)
    plot(rrFR{1,i}(:,1), rrFR{1,i}(:,2),'k', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('Position [m]')
    title('Forced/Choice Right to Right')
    axis tight
    
    subplot(412)
    plot(rrFR2{1,i}(:,1), rrFR2{1,i}(:,2)*(size(rrFR2{1,i},2)-20),'LineWidth',2,'Color', [0.7 0.7 0.7])
    for kk = 1: size(rrFRidx2{1,1},2)
        hold on
         plot(rrFR2{1,i}(rrFRidx2{1,i}{1,kk},1), rrFR2{1,i}(rrFRidx2{1,i}{1,kk},kk+5)+kk-1, '.K');
         hold off
    end
    axis([rrFR2{1,i}(1) rrFR2{1,i}(end,1) 0 kk+2])
    xlabel('Time [sec]')
    ylabel('MUA #')
    
    subplot(413)
    plot(rrFR{1,i}(:,1), rrFR{1,i}(:,3),'r', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('m/sec')
    axis tight
    
    subplot(414)
    hold on
%     plot(rrFR{1,i}(:,1), rrFR{1,i}(:,4),'b')
    plot(rrFR{1,i}(:,1), rrFR{1,i}(:,5),'k', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('LPF[uV]')
    axis tight
    
end
%%
