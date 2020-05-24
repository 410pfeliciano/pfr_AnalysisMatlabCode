NewSr = 1/2000;
%%
llFR = cell(1,size(ch_run,1)); 
for jj = 1:size(ch_run,1)
    if ch_run{jj,3} == 'L' && ch_run{jj,4} == 'L'
        ind1 = ch_run{jj,1} * (1/NewSr);
        ind2 = ch_run{jj,2} * (1/NewSr);
        llFR{1,jj}(:,1) = linTime(ind1: ind2); % Time
        llFR{1,jj}(:,2) = linPos(ind1: ind2); % Position
        llFR{1,jj}(:,3) = linVel(ind1: ind2); % Velocity
        llFR{1,jj}(:,4) = LFP7(ind1: ind2); % LFP
        llFR{1,jj}(:,5) = thetaLFP7(ind1: ind2); % Theta filtered LFP
        for kk = 1: size(semiMUAtrains,2) % Units and MUA spikes
            llFR{1,jj}(:,5+kk) = semiMUAtrains(ind1:ind2,kk);
        end
    end
end
% Removing [] cell arrays
llFR = llFR(~cellfun(@isempty,llFR));

%% Finding Spikes index
llFRidx = cell(1, size(llFR,2));
for jj = 1 : size(llFR,2)
    for kk = 1 : size(semiMUAtrains,2)
        llFRidx{1,jj}{:,kk} = find(llFR{1,jj}(:, kk+5));%Column 5 and up are spikes
    end
end

%% Creating Encoding Model
EncTrials = [4 5 6]; % Training Data
binFac = 0.04;
tuningbins = (0.0:binFac:1.3)';
llTuning = cell(1, size(semiMUAtrains,2));
for kk = 1 : size(EncTrials,2)
   for jj = 1 : size(semiMUAtrains,2)
       llTuning{1,jj}(:,kk) = (hist(llFR{1,kk}(llFRidx{1,kk}{1,jj},2),tuningbins))';
   end
end
llTuningSum = cell(1, size(llTuning,2));
for kk = 1: size(llTuning,2)
    llTuningSum{1,kk} = sum(llTuning{1,kk},2);
end
llOccup = zeros(size(tuningbins,1), size(llFR,2));
for kk = 1: size(llFR,2)
   llOccup(:,kk) = (hist(llFR{1,kk}(:,2),tuningbins).*NewSr)';
end
llOccup = sum(llOccup,2);
llTuningCurve = zeros(size(tuningbins,1), size(llTuningSum,2));
for kk = 1: size(llTuningSum,2)
    llTuningCurve(:,kk) = llTuningSum{1,kk}./llOccup;
end
llTuningCurve(isnan(llTuningCurve)) = 0;
llTuningCurveSm = zeros((size(llTuningCurve,1)), size(llTuningCurve,2));
for kk = 1: size(llTuningCurve,2)
   llTuningCurveSm(:,kk) = (smoothdata(llTuningCurve(:,kk),'gaussian',5))+.01; % Adding 0.01Hz
end
llnSR = 1.3/(length(llTuningCurve(:,1))+10);
llPosbins = (0 : llnSR : 1.3)';
close all
for kk = 1: size(llTuningCurve,2)
   s(kk) = subplot(size(semiMUAtrains,2),1,kk);
   hold on
   area(tuningbins, llTuningCurveSm(:,kk),'FaceColor','R')
%     line('XData', [.9 .9], 'YData', [0 55], 'LineStyle', '-.', ...
%     'LineWidth', 2, 'Color',[0.5 0.5 0.5])
%     line('XData', [.35 .35], 'YData', [0 55], 'LineStyle', '-.', ...
%     'LineWidth', 2, 'Color',[0.5 0.5 0.5])
    hold off
   axis([0.0 1.3 0 inf])
end
xlabel('Position [m]')
ylabel('FR[Hz]')
title(s(1),'Forced/Choice Left to Left','FontSize',12)
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
Timebins = cell(1,size(llFR,2));
TimeTrialbins = cell(1,size(llFR,2));
for kk = 1: size(llFRidx,2)
    % Timebins to bin Spikes over time
    Timebins{1,kk} = (0: BinTime: (llFR{1,kk}(end,1) - llFR{1,kk}(1,1)))'; 
    % Timebins for bin Spikes over time/Trials
    TimeTrialbins{1,kk} = (llFR{1,kk}(1,1): BinTime : llFR{1,kk}(end,1)); 
end
%% Testing Data/Spikes
llFR2 = cell(1,size(ch_run,1)); 
for jj = 1:size(ch_run,1)
    if ch_run{jj,3} == 'L' && ch_run{jj,4} == 'L'
        ind1 = ch_run{jj,1} * (1/NewSr);
        ind2 = ch_run{jj,2} * (1/NewSr);
        llFR2{1,jj}(:,1) = linTime(ind1: ind2); % Time
        llFR2{1,jj}(:,2) = linPos(ind1: ind2); % Position
        llFR2{1,jj}(:,3) = linVel(ind1: ind2); % Velocity
        llFR2{1,jj}(:,4) = LFP7(ind1: ind2); % LFP
        llFR2{1,jj}(:,5) = thetaLFP7(ind1: ind2); % Theta filtered LFP
        for kk = 1: size(semiMUAtrains,2) % Units and MUA spikes
            llFR2{1,jj}(:,5+kk) = semiMUAtrains(ind1:ind2,kk);
        end
    end
end
% Removing [] cell allays
llFR2 = llFR2(~cellfun(@isempty,llFR2));
llFRidx2 = cell(1, size(llFR2,2));
for jj = 1 : size(llFR2,2)
    for kk = 1 : size(semiMUAtrains,2)
        llFRidx2{1,jj}{:,kk} = find(llFR2{1,jj}(:, kk+5));%Column 5 and up are spikes
    end
end
%% Histogram Spikes over time/trials
llTimeHist = cell(1, size(llFRidx2,2));
for jj = 1: size(llFRidx2,2)
    for kk = 1: size(llFRidx2{1,1},2)
     llTimeHist{1,jj}(:,kk) =(hist(llFR{1,jj}(llFRidx2{1,jj}{1,kk},1), ...
         TimeTrialbins{1,jj}))';
    end
end

%% Summing the firing rates per position
llSumFRate = (sum(llTuningCurveSm,2));
llSumFRExp = exp(-(BinTime).* sum(llTuningCurveSm,2));-0.2.*(exp(llSumFRate));
figure
% subplot(2,1,1)
plot(tuningbins, llSumFRate,'-r', 'LineWidth',2)
xlabel('Position [m]')
ylabel('\Sigma FR[Hz]')
axis([0.0 1.3 0 inf])
% subplot(2,1,2)
% plot(tuningbins, llSumFRExp)
% xlabel('Position [m]')
% ylabel('e^{-\Deltat \Sigma FR[Hz]}')
% axis([0.0 1.3 0 inf])
%% Firing rate elevated by time spikes
llFR_eSpikes = cell(1, size(llFRidx,2));
for jj = 1: size(llFRidx,2)
    for tt = 1: length(llTimeHist{1,jj}(:,1))
        llFR_eSpikes{1, jj}(tt,:) = (prod(llTuningCurveSm(:,:).^ ...
            llTimeHist{1,jj}(tt,:),2))';
    end
end

%%
llProb = cell(1, size(llFRidx,2));
SumllProb = cell(1, size(llFRidx,2));
llLikelihood = cell(1, size(llFRidx,2));
for kk = 1: size(llFRidx,2)
    llFR_eSpikes{1,kk} = (llFR_eSpikes{1,kk})';
    llProb{1,kk} = llSumFRExp.* llFR_eSpikes{1,kk};
    SumllProb{1,kk} = sum(llProb{1,kk},1);
    llLikelihood{1,kk} = llProb{1,kk}./SumllProb{1,kk};
end

%% Visual Inspection
close all
for kk = 1: size(llFRidx,2)
    figure(kk)
    hold on
    s = plot3(llFR{1,kk}(:,1), llFR{1,kk}(:,2), ...
        ones(size(llFR{1,kk}(:,2),1),1),'-.r','LineWidth', 0.5);
    alpha(s,0.1)
    s1 = surf(TimeTrialbins{1,kk}, tuningbins,llLikelihood{1,kk});
    s1.EdgeColor = 'none';
    c = colorbar;
    c.Label.String = 'Probability';
    colormap(flipud(bone))
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
llIndLikeli = cell(1, size(llFRidx,2));
for kk = 1: 3%size(llFRidx,2)
    [MaxLikeli, llIndLikeli{kk}] = max(llLikelihood{kk});
end
llMaxLikeliPos = cell(1, size(llFRidx,2));
for kk = 1:3%size(llFRidx,2)
    llMaxLikeliPos{kk} = tuningbins(llIndLikeli{kk},1);
end
close all
for kk = 1:3% size(llFRidx,2)
    figure
    plot(TimeTrialbins{kk}, llMaxLikeliPos{kk}, 'o')
    hold on
    plot(llFR{1,kk}(:,1),llFR{1,kk}(:,2))
    hold off
end
%% Position histogram over time
llPosTimeHist = cell(1,size(llFRidx,2)); % Memory Preallocation
llVelTimeHist = cell(1,size(llFRidx,2));
for jj = 1: 3%size(llFRidx,2)
    [a ,b] = hist(llFR{1,jj}(1:Timebins{jj}(1+1)*2000, 2),1); % Fist bin hist
    [a1 ,b1] = hist(llFR{1,jj}(1:Timebins{jj}(1+1)*2000, 3),1); % Fist bin hist
    llPosTimeHist{jj}(1,1) = b;
    llVelTimeHist{jj}(1,1) = b1;
    for kk = 2: length(Timebins{jj})-1
        [a ,b] = hist(llFR{1,jj}(Timebins{jj}(kk)*2000:Timebins{jj}(kk+1)*2000,2),1);
        llPosTimeHist{jj}(kk,1) = b;
        [a1 ,b1] = hist(llFR{1,jj}(Timebins{jj}(kk)*2000:Timebins{jj}(kk+1)*2000,3),1);
        llVelTimeHist{jj}(kk,1) = b1;
    end
    [a ,b] = hist(llFR{1,jj}(Timebins{jj}(size(TimeTrialbins{jj},2))*2000:Timebins{jj}(end)*2000, 2),1); % Fist bin hist
    llPosTimeHist{jj}(size(TimeTrialbins{jj},2),1) = b;
    [a1 ,b1] = hist(llFR{1,jj}(Timebins{jj}(size(TimeTrialbins{jj},2))*2000:Timebins{jj}(end)*2000, 3),1); % Fist bin hist
    llVelTimeHist{jj}(size(TimeTrialbins{jj},2),1) = b1;
end
close all
for kk = 1: 3%size(llFRidx,2)
    figure
    plot(TimeTrialbins{kk}, llMaxLikeliPos{kk}, 'ob')
    hold on
    plot(llFR{1,kk}(:,1),llFR{1,kk}(:,2))
    plot(TimeTrialbins{kk}, llPosTimeHist{kk}, '-r')
    plot(TimeTrialbins{kk}, llVelTimeHist{kk}, '-m')
    hold off
end
%% Calculating Decoding Error
llLikeliError = cell(1, size(llFRidx,2));
for kk = 1: 3%size(llFRidx,2)
    llLikeliError{kk} = abs(llPosTimeHist{kk}- llMaxLikeliPos{kk});
end

close all
for kk = 1: 3%size(llFRidx,2)
    figure
    plot(TimeTrialbins{kk}, llLikeliError{kk}, '-r', 'LineWidth',2)
    xlabel('Time [sec]')
    ylabel('Error [m]')
    axis tight
end
%% Eliminating non running samples
for kk = 1: 3%size(llFRidx,2)
    for jj = 1: size(TimeTrialbins{kk},2)
        if llVelTimeHist{kk}(jj,1) > 0.03
            llPosTimeHist1{kk}(jj,1) = llPosTimeHist{kk}(jj,1);
            llMaxLikeliPos1{kk}(jj,1) = llMaxLikeliPos{kk}(jj,1);
        else
             llPosTimeHist1{kk}(jj,1) = NaN;
             llMaxLikeliPos1{kk}(jj,1) = NaN;
        end
    end
end
%% Calculating Decoding Error from running
llLikeliError1 = cell(1, 3);
for kk = 1:3 %size(llFRidx,2)
    llLikeliError1{kk} = abs(llPosTimeHist1{kk}- llMaxLikeliPos1{kk});
end

close all
for kk = 1:3 %size(llFRidx,2)
    figure
    plot(TimeTrialbins{kk}, llLikeliError1{kk}, '-b')
end
%%
close all
for kk = 1:3 %size(llFRidx,2)
    figure
    h = histogram(llLikeliError1{kk},1000,'Normalization','cdf');
    ErrorCum = (h.Values)';
    cumbins = h.BinEdges;
    MaxError = max(ErrorCum);
    llErrorCum(:,kk) = ErrorCum./MaxError;
end
%%
close all
llAvgCumError = zeros(size(llErrorCum,1)+1,1);
llAvgCumError(2:end) = mean(llErrorCum,2);
cumbins(1,1) = 0;
figure
plot(cumbins,llAvgCumError,'-r','LineWidth',2)
xlabel('Linear Position [m]')
ylabel('Error Cum. Prob')
axis([0 1.2 0 1])
%%
for kk = 1:3 %size(llFRidx,2)
    C{kk} = confusionmat(llPosTimeHist{kk},llMaxLikeliPos{kk});
    figure
    subplot(211)
    imagesc(llPosTimeHist{kk},llMaxLikeliPos{kk},C{kk});
    axis([0 1.3 0 1.3])
    subplot(212)
    plot(llPosTimeHist{kk},llMaxLikeliPos{kk},'or')
    axis([0 1.3 0 1.3])
%     figure
%     plotConfMat(llPosTimeHist{kk},llMaxLikeliPos{kk})
end
%% Visualization
close all
for i=1:size(llFR,2)
    figure(i)
    subplot(411)
    plot(llFR{1,i}(:,1), llFR{1,i}(:,2),'k', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('Position [m]')
    title('Forced/Choice Right to Right')
    axis([-inf inf 0 1.3])
    
    subplot(412)
    plot(llFR2{1,i}(:,1), llFR2{1,i}(:,2)*(size(llFR2{1,i},2)-20),'LineWidth',2,'Color', [0.7 0.7 0.7])
    for kk = 1: size(llFRidx2{1,1},2)
        hold on
         plot(llFR2{1,i}(llFRidx2{1,i}{1,kk},1), llFR2{1,i}(llFRidx2{1,i}{1,kk},kk+5)+kk-1, '.r');
         hold off
    end
    axis([llFR2{1,i}(1) llFR2{1,i}(end,1) 0 kk+2])
    xlabel('Time [sec]')
    ylabel('Units #')
    
    subplot(413)
    plot(llFR{1,i}(:,1), llFR{1,i}(:,3),'r', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('m/sec')
    axis tight
    
    subplot(414)
    hold on
%     plot(llFR{1,i}(:,1), llFR{1,i}(:,4),'b')
    plot(llFR{1,i}(:,1), llFR{1,i}(:,5),'k', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('LPF[uV]')
    axis tight
    
end