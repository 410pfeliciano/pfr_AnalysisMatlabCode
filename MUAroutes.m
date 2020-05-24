NewSr = 1/2000; % 2Khz sampling rate after downsampling the LFP signal
linPos = interp1(new_t, lpos, LFPt);
linVel = interp1(new_t, new_vel, LFPt);
linTime = 0:NewSr: length(linPos)* NewSr - NewSr;
%%
save('PosVel.mat', 'linPos', 'linVel', 'linTime')
%%
MUAvectActive = zeros(size(MUAvect));
for kk = 1 : size(MUAvect,2)
    for qq = 1:size(MUAvect,1)
        if linVel(qq) > 0.03
        MUAvectActive(qq,kk) = MUAvect(qq,kk);
        else
        MUAvectActive(qq,kk) = 0;
        end
    end
end

%%
rrFR = cell(1,size(ch_run,1)); 
llFR = cell(1,size(ch_run,1));
for jj = 1:size(ch_run,1)
    if ch_run{jj,3} == 'R' && ch_run{jj,4} == 'R'
        ind1 = ch_run{jj,1} * (1/NewSr);
        ind2 = ch_run{jj,2} * (1/NewSr);
        rrFR{1,jj}(:,1) = linTime(ind1: ind2); % Time
        rrFR{1,jj}(:,2) = linPos(ind1: ind2); % Position
        rrFR{1,jj}(:,3) = linVel(ind1: ind2); % Velocity
        rrFR{1,jj}(:,4) = LFP7(ind1: ind2); %LFP
        rrFR{1,jj}(:,5) = thetaLFP7(ind1: ind2); %Theta filtered LFP
        for kk = 1: size(MUAvectActive,2) % Spikes(0 or 1)
            rrFR{1,jj}(:,5+kk) = MUAvectActive(ind1:ind2,kk);
        end
    elseif ch_run{jj,3} == 'L' && ch_run{jj,4} == 'L'
        ind1 = ch_run{jj,1} * (1/NewSr);
        ind2 = ch_run{jj,2} * (1/NewSr);
        llFR{1,jj}(:,1) = linTime(ind1: ind2);
        llFR{1,jj}(:,2) = linPos(ind1: ind2);
        llFR{1,jj}(:,3) = linVel(ind1: ind2);
        llFR{1,jj}(:,4) = LFP7(ind1: ind2);
        llFR{1,jj}(:,5) = thetaLFP7(ind1: ind2);
        for kk = 1: size(MUAvectActive,2)
            llFR{1,jj}(:,5+kk) = MUAvectActive(ind1:ind2,kk);
        end
    end
end
% Removing [] cell arrays
rrFR = rrFR(~cellfun(@isempty,rrFR)); % Right to Right Choice
llFR = llFR(~cellfun(@isempty,llFR)); % Left to Left Choice

%% Finding Spikes
rrFRidx = cell(1, size(rrFR,2));
for jj = 1 : size(rrFR,2)
    for kk = 1 : size(MUAvectActive,2)
        rrFRidx{1,jj}{:,kk} = find(rrFR{1,jj}(:, kk+5));
    end
end
llFRidx = cell(1, size(llFR,2));
for jj = 1 : size(llFR,2)
    for kk = 1 : size(MUAvectActive,2)
        llFRidx{1,jj}{:,kk} = find(llFR{1,jj}(:, kk+5));
    end
end

%% Tuning Curve for Choice Right to Right
binFac = 0.025;
tuningbins = (0:binFac:1.3)';
rrTuning = cell(1, size(MUAvectActive,2));
for kk = 1 : size(rrFR,2)
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
   s(kk) = subplot(size(MUAvect,2),1,kk);
   hold on
   area(tuningbins, rrTuningCurveSm(:,kk),'FaceColor','k')
    line('XData', [.9 .9], 'YData', [0 55], 'LineStyle', '-.', ...
    'LineWidth', 2, 'Color',[0.5 0.5 0.5])
    line('XData', [.35 .35], 'YData', [0 55], 'LineStyle', '-.', ...
    'LineWidth', 2, 'Color',[0.5 0.5 0.5])
    hold off
   axis([0 1.3 0 inf])
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
axis(s(1),[0 inf 0 25])
axis(s(2),[0 inf 0 30])
axis(s(3),[0 inf 0 53])
axis(s(4),[0 inf 0 20])
axis(s(5),[0 inf 0 10])
axis(s(6),[0 inf 0 25])
axis(s(7),[0 inf 0 20])
axis(s(8),[0 inf 0 10])
axis(s(9),[0 inf 0 50])
axis(s(10),[0 inf 0 45])
axis(s(11),[0 inf 0 20])
axis(s(12),[0 inf 0 30])
axis(s(13),[0 inf 0 45])
axis(s(14),[0 inf 0 10])
%%
save('rrTuningCurves.mat','rrFR', 'rrFRidx', 'tuningbins', 'rrTuningCurve', ...
    'rrTuningCurveSm') 

%%
llTuning = cell(1, size(MUAvectActive,2));
for kk = 1 : size(llFR,2)
   for jj = 1 : size(MUAvectActive,2)
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
   llTuningCurveSm(:,kk) = (smoothdata(llTuningCurve(:,kk),'gaussian',6))+.01; % Adding 0.01Hz
end
llnSR = 1.3/(length(llTuningCurve(:,1))+5);
llPosbins = (0 : llnSR : 1.3)';
close all
for kk = 1: size(llTuningCurve,2)
   s(kk) = subplot(size(MUAvect,2),1,kk);
   hold on
   area(tuningbins, llTuningCurveSm(:,kk),'FaceColor','r')
%    plot(tuningbins, llTuningCurve(:,kk), 'LineWidth' ,1)
    line('XData', [.9 .9], 'YData', [0 55], 'LineStyle', '-.', ...
    'LineWidth', 2, 'Color',[0.5 0.5 0.5])
    line('XData', [.35 .35], 'YData', [0 55], 'LineStyle', '-.', ...
    'LineWidth', 2, 'Color',[0.5 0.5 0.5])
   hold off
   axis([0 1.3 0 inf])
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
axis(s(1),[0 inf 0 25])
axis(s(2),[0 inf 0 30])
axis(s(3),[0 inf 0 53])
axis(s(4),[0 inf 0 20])
axis(s(5),[0 inf 0 10])
axis(s(6),[0 inf 0 25])
axis(s(7),[0 inf 0 20])
axis(s(8),[0 inf 0 10])
axis(s(9),[0 inf 0 50])
axis(s(10),[0 inf 0 45])
axis(s(11),[0 inf 0 20])
axis(s(12),[0 inf 0 30])
axis(s(13),[0 inf 0 45])
axis(s(14),[0 inf 0 10])

%% Visual Inspection
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
    plot(rrFR{1,i}(:,1), rrFR{1,i}(:,2)*(size(rrFR{1,i},2)-8),'LineWidth',2,'Color', [0.7 0.7 0.7])
    for kk = 1: size(rrFRidx{1,1},2)
        hold on
         plot(rrFR{1,i}(rrFRidx{1,i}{1,kk},1), rrFR{1,i}(rrFRidx{1,i}{1,kk},kk+5)+kk-1, '.K');
         hold off
    end
    axis([rrFR{1,i}(1) rrFR{1,i}(end,1) 0 kk+2])
    xlabel('Time [sec]')
    ylabel('Units #')
    
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
save('llTuningCurves.mat','llFR', 'llFRidx', 'tuningbins', 'llTuningCurve', ...
    'llTuningCurveSm')
%%
close all
for i=1:size(llFR,2)
    figure(i)
    subplot(411)
    plot(llFR{1,i}(:,1), llFR{1,i}(:,2),'k', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('Position [m]')
    axis tight
    title('Forced/Choice Left to Left')
    subplot(412)
    plot(llFR{1,i}(:,1), llFR{1,i}(:,2)*(size(llFR{1,i},2)-8),'LineWidth',2,'Color', [0.7 0.7 0.7])
    for kk = 1: size(llFRidx{1,1},2)
        hold on
         plot(llFR{1,i}(llFRidx{1,i}{1,kk},1), llFR{1,i}(llFRidx{1,i}{1,kk},kk+5)+kk-1, '.K');
         hold off
    end
    axis([llFR{1,i}(1) llFR{1,i}(end,1) 0 kk+2])
    xlabel('Time [sec]')
    ylabel('Units #')
    subplot(413)
    plot(llFR{1,i}(:,1), llFR{1,i}(:,3),'r', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('m/sec')
    axis tight
    
    subplot(414)
    hold on
    plot(llFR{1,i}(:,1), llFR{1,i}(:,4),'b')
    hold on
    plot(llFR{1,i}(:,1), llFR{1,i}(:,5),'k', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('LPF[uV]')
    axis tight
end


%%
frrrFR = cell(1,size(fr_run,1)); 
frllFR = cell(1,size(fr_run,1));
frrlFR = cell(1,size(fr_run,1)); 
frlrFR = cell(1,size(fr_run,1));
for jj = 1:size(fr_run,1)
    if fr_run{jj,3} == 'R' && fr_run{jj,4} == 'R'
        ind1 = fr_run{jj,1} * (1/NewSr);
        ind2 = fr_run{jj,2} * (1/NewSr);
        frrrFR{1,jj}(:,1) = linTime(ind1: ind2);
        frrrFR{1,jj}(:,2) = linPos(ind1: ind2);
        frrrFR{1,jj}(:,3) = linVel(ind1: ind2);
        frrrFR{1,jj}(:,4) = LFP7(ind1: ind2);
        frrrFR{1,jj}(:,5) = thetaLFP7(ind1: ind2);
        for kk = 1: size(MUAvectActive,2)
            frrrFR{1,jj}(:,5+kk) = MUAvectActive(ind1:ind2,kk);
        end
    elseif fr_run{jj,3} == 'L' && fr_run{jj,4} == 'L'
        ind1 = fr_run{jj,1} * (1/NewSr);
        ind2 = fr_run{jj,2} * (1/NewSr);
        frllFR{1,jj}(:,1) = linTime(ind1: ind2);
        frllFR{1,jj}(:,2) = linPos(ind1: ind2);
        frllFR{1,jj}(:,3) = linVel(ind1: ind2);
        frllFR{1,jj}(:,4) = LFP7(ind1: ind2);
        frllFR{1,jj}(:,5) = thetaLFP7(ind1: ind2);
        for kk = 1: size(MUAvectActive,2)
            frllFR{1,jj}(:,5+kk) = MUAvectActive(ind1:ind2,kk);
        end
    elseif fr_run{jj,3} == 'R' && fr_run{jj,4} == 'L'
        ind1 = fr_run{jj,1} * (1/NewSr);
        ind2 = fr_run{jj,2} * (1/NewSr);
        frrlFR{1,jj}(:,1) = linTime(ind1: ind2);
        frrlFR{1,jj}(:,2) = linPos(ind1: ind2);
        frrlFR{1,jj}(:,3) = linVel(ind1: ind2);
        frrlFR{1,jj}(:,4) = LFP7(ind1: ind2);
        frrlFR{1,jj}(:,5) = thetaLFP7(ind1: ind2);
        for kk = 1: size(MUAvectActive,2)
            frrlFR{1,jj}(:,5+kk) = MUAvectActive(ind1:ind2,kk);
        end
    elseif fr_run{jj,3} == 'L' && fr_run{jj,4} == 'R'
        ind1 = fr_run{jj,1} * (1/NewSr);
        ind2 = fr_run{jj,2} * (1/NewSr);
        frlrFR{1,jj}(:,1) = linTime(ind1: ind2);
        frlrFR{1,jj}(:,2) = linPos(ind1: ind2);
        frlrFR{1,jj}(:,3) = linVel(ind1: ind2);
        frlrFR{1,jj}(:,4) = LFP7(ind1: ind2);
        frlrFR{1,jj}(:,5) = thetaLFP7(ind1: ind2);
        for kk = 1: size(MUAvectActive,2)
            frlrFR{1,jj}(:,5+kk) = MUAvectActive(ind1:ind2,kk);
        end
    end
end
% Removing [] cell arrays
frrrFR = frrrFR(~cellfun(@isempty,frrrFR));
frllFR = frllFR(~cellfun(@isempty,frllFR));
frrlFR = frrlFR(~cellfun(@isempty,frrlFR));
frlrFR = frlrFR(~cellfun(@isempty,frlrFR));

%% Finding Spikes
frrrFRidx = cell(1, size(frrrFR,2));
for jj = 1 : size(frrrFR,2)
    for kk = 1 : size(MUAvectActive,2)
        frrrFRidx{1,jj}{:,kk} = find(frrrFR{1,jj}(:, kk+5));
    end
end
frllFRidx = cell(1, size(frllFR,2));
for jj = 1 : size(frllFR,2)
    for kk = 1 : size(MUAvectActive,2)
        frllFRidx{1,jj}{:,kk} = find(frllFR{1,jj}(:, kk+5));
    end
end
frrlFRidx = cell(1, size(frrlFR,2));
for jj = 1 : size(frrlFR,2)
    for kk = 1 : size(MUAvectActive,2)
        frrlFRidx{1,jj}{:,kk} = find(frrlFR{1,jj}(:, kk+5));
    end
end
frlrFRidx = cell(1, size(frlrFR,2));
for jj = 1 : size(frlrFR,2)
    for kk = 1 : size(MUAvectActive,2)
        frlrFRidx{1,jj}{:,kk} = find(frlrFR{1,jj}(:, kk+5));
    end
end
llTuning = cell(1, size(MUAvectActive,2));
for kk = 1 : size(llFR,2)
   for jj = 1 : size(MUAvectActive,2)
       llTuning{1,jj}(:,kk) = (hist(llFR{1,kk}(llFRidx{1,kk}{1,jj},2),tuningbins))';
   end
end

%% Tuning Curve Forced Right Right
frrrTuning = cell(1, size(MUAvectActive,2));
for kk = 1 : size(frrrFR,2)
   for jj = 1 : size(MUAvectActive,2)
       frrrTuning{1,jj}(:,kk) = (hist(frrrFR{1,kk}(frrrFRidx{1,kk}{1,jj},2),tuningbins))';
   end
end
frrrTuningSum = cell(1, size(frrrTuning,2));
for kk = 1: size(frrrTuning,2)
    frrrTuningSum{1,kk} = sum(frrrTuning{1,kk},2);
end
frrrOccup = zeros(size(tuningbins,1), size(frrrFR,2));
for kk = 1: size(frrrFR,2)
   frrrOccup(:,kk) = (hist(frrrFR{1,kk}(:,2),tuningbins).*NewSr)';
end
frrrOccup = sum(frrrOccup,2);
frrrTuningCurve = zeros(size(tuningbins,1), size(frrrTuningSum,2));
for kk = 1: size(frrrTuningSum,2)
    frrrTuningCurve(:,kk) = frrrTuningSum{1,kk}./frrrOccup;
end
frrrTuningCurve(isnan(frrrTuningCurve)) = 0;
frrrTuningCurveSm = zeros((size(frrrTuningCurve,1)), size(frrrTuningCurve,2));
for kk = 1: size(frrrTuningCurve,2)
   frrrTuningCurveSm(:,kk) = (smoothdata(frrrTuningCurve(:,kk),'gaussian',6))+.01; % Adding 0.01Hz
end
close all
for kk = 1: size(frrrTuningCurve,2)
   s(kk) = subplot(size(MUAvect,2),1,kk);
   hold on
   area(tuningbins, frrrTuningCurveSm(:,kk),'FaceColor','b')
    line('XData', [.9 .9], 'YData', [0 55], 'LineStyle', '-.', ...
    'LineWidth', 2, 'Color',[0.5 0.5 0.5])
    line('XData', [.35 .35], 'YData', [0 55], 'LineStyle', '-.', ...
    'LineWidth', 2, 'Color',[0.5 0.5 0.5])
   hold off
   axis([0 1.3 0 inf])
end
xlabel('Position [m]')
ylabel('FR[Hz]')
title(s(1),'Forced Right to Right','FontSize',12)
txt1 = 'Start ';
text(.2,115,txt1,'HorizontalAlignment','Center','FontSize', 8)
txt2 = 'Center ';
text(.6,115,txt2,'HorizontalAlignment','Center','FontSize', 8)
txt3 = 'Reward ';
text(1.1,115,txt3,'HorizontalAlignment','Center','FontSize', 8)
txt4 = '\leftarrow';
text(1.23,120,txt4,'HorizontalAlignment','Center','FontSize', 40)
axis(s(1),[0 inf 0 25])
axis(s(2),[0 inf 0 30])
axis(s(3),[0 inf 0 53])
axis(s(4),[0 inf 0 40])
axis(s(5),[0 inf 0 10])
axis(s(6),[0 inf 0 25])
axis(s(7),[0 inf 0 20])
axis(s(8),[0 inf 0 10])
axis(s(9),[0 inf 0 50])
axis(s(10),[0 inf 0 45])
axis(s(11),[0 inf 0 20])
axis(s(12),[0 inf 0 30])
axis(s(13),[0 inf 0 45])
axis(s(14),[0 inf 0 5])

%%
frllTuning = cell(1, size(MUAvectActive,2));
for kk = 1 : size(frllFR,2)
   for jj = 1 : size(MUAvectActive,2)
       frllTuning{1,jj}(:,kk) = (hist(frllFR{1,kk}(frllFRidx{1,kk}{1,jj},2),tuningbins))';
   end
end
frllTuningSum = cell(1, size(frllTuning,2));
for kk = 1: size(frllTuning,2)
    frllTuningSum{1,kk} = sum(frllTuning{1,kk},2);
end
frllOccup = zeros(size(tuningbins,1), size(frllFR,2));
for kk = 1: size(frllFR,2)
   frllOccup(:,kk) = (hist(frllFR{1,kk}(:,2),tuningbins).*NewSr)';
end
frllOccup = sum(frllOccup,2);
frllTuningCurve = zeros(size(tuningbins,1), size(frllTuningSum,2));
for kk = 1: size(frllTuningSum,2)
    frllTuningCurve(:,kk) = frllTuningSum{1,kk}./frllOccup;
end
frllTuningCurve(isnan(frllTuningCurve)) = 0;
frllTuningCurveSm = zeros((size(frllTuningCurve,1)), size(frllTuningCurve,2));
for kk = 1: size(frllTuningCurve,2)
   frllTuningCurveSm(:,kk) = (smoothdata(frllTuningCurve(:,kk),'gaussian',6))+.01; % Adding 0.01Hz
end
close all
for kk = 1: size(frllTuningCurve,2)
   s(kk) = subplot(size(MUAvect,2),1,kk);
   hold on
   area(tuningbins, frllTuningCurveSm(:,kk),'FaceColor','g')
    line('XData', [.9 .9], 'YData', [0 55], 'LineStyle', '-.', ...
    'LineWidth', 2, 'Color',[0.5 0.5 0.5])
    line('XData', [.35 .35], 'YData', [0 55], 'LineStyle', '-.', ...
    'LineWidth', 2, 'Color',[0.5 0.5 0.5])
   hold off
   axis([0 1.3 0 inf])
end
xlabel('Position [m]')
ylabel('FR[Hz]')
title(s(1),'Forced Left to Left','FontSize',12)
txt1 = 'Start ';
text(.2,115,txt1,'HorizontalAlignment','Center','FontSize', 8)
txt2 = 'Center ';
text(.6,115,txt2,'HorizontalAlignment','Center','FontSize', 8)
txt3 = 'Reward ';
text(1.1,115,txt3,'HorizontalAlignment','Center','FontSize', 8)
txt4 = '\leftarrow';
text(1.23,120,txt4,'HorizontalAlignment','Center','FontSize', 40)
axis(s(1),[0 inf 0 25])
axis(s(2),[0 inf 0 30])
axis(s(3),[0 inf 0 53])
axis(s(4),[0 inf 0 40])
axis(s(5),[0 inf 0 10])
axis(s(6),[0 inf 0 25])
axis(s(7),[0 inf 0 20])
axis(s(8),[0 inf 0 10])
axis(s(9),[0 inf 0 50])
axis(s(10),[0 inf 0 45])
axis(s(11),[0 inf 0 20])
axis(s(12),[0 inf 0 30])
axis(s(13),[0 inf 0 45])
axis(s(14),[0 inf 0 5])

%%
frrlTuning = cell(1, size(MUAvectActive,2));
for kk = 1 : size(frrlFR,2)
   for jj = 1 : size(MUAvectActive,2)
       frrlTuning{1,jj}(:,kk) = (hist(frrlFR{1,kk}(frrlFRidx{1,kk}{1,jj},2),tuningbins))';
   end
end
frrlTuningSum = cell(1, size(frrlTuning,2));
for kk = 1: size(frrlTuning,2)
    frrlTuningSum{1,kk} = sum(frrlTuning{1,kk},2);
end
frrlOccup = zeros(size(tuningbins,1), size(frrlFR,2));
for kk = 1: size(frrlFR,2)
   frrlOccup(:,kk) = (hist(frrlFR{1,kk}(:,2),tuningbins).*NewSr)';
end
frrlOccup = sum(frrlOccup,2);
frrlTuningCurve = zeros(size(tuningbins,1), size(frrlTuningSum,2));
for kk = 1: size(frrlTuningSum,2)
    frrlTuningCurve(:,kk) = frrlTuningSum{1,kk}./frrlOccup;
end
frrlTuningCurve(isnan(frrlTuningCurve)) = 0;
frrlTuningCurveSm = zeros((size(frrlTuningCurve,1)), size(frrlTuningCurve,2));
for kk = 1: size(frrlTuningCurve,2)
   frrlTuningCurveSm(:,kk) = (smoothdata(frrlTuningCurve(:,kk),'gaussian',6))+.01; % Adding 0.01Hz
end
close all
for kk = 1: size(frrlTuningCurve,2)
   s(kk) = subplot(size(MUAvect,2),1,kk);
   hold on
   area(tuningbins, frrlTuningCurveSm(:,kk),'FaceColor','m')
    line('XData', [.9 .9], 'YData', [0 55], 'LineStyle', '-.', ...
    'LineWidth', 2, 'Color',[0.5 0.5 0.5])
    line('XData', [.35 .35], 'YData', [0 55], 'LineStyle', '-.', ...
    'LineWidth', 2, 'Color',[0.5 0.5 0.5])
   hold off
   axis([0 1.3 0 inf])
end
xlabel('Position [m]')
ylabel('FR[Hz]')
title(s(1),'Forced Right to Left','FontSize',12)
txt1 = 'Start ';
text(.2,115,txt1,'HorizontalAlignment','Center','FontSize', 8)
txt2 = 'Center ';
text(.6,115,txt2,'HorizontalAlignment','Center','FontSize', 8)
txt3 = 'Reward ';
text(1.1,115,txt3,'HorizontalAlignment','Center','FontSize', 8)
txt4 = '\leftarrow';
text(1.23,120,txt4,'HorizontalAlignment','Center','FontSize', 40)
axis(s(1),[0 inf 0 25])
axis(s(2),[0 inf 0 30])
axis(s(3),[0 inf 0 53])
axis(s(4),[0 inf 0 40])
axis(s(5),[0 inf 0 10])
axis(s(6),[0 inf 0 25])
axis(s(7),[0 inf 0 20])
axis(s(8),[0 inf 0 10])
axis(s(9),[0 inf 0 50])
axis(s(10),[0 inf 0 45])
axis(s(11),[0 inf 0 20])
axis(s(12),[0 inf 0 30])
axis(s(13),[0 inf 0 45])
axis(s(14),[0 inf 0 5])

%%
frlrTuning = cell(1, size(MUAvectActive,2));
for kk = 1 : size(frlrFR,2)
   for jj = 1 : size(MUAvectActive,2)
       frlrTuning{1,jj}(:,kk) = (hist(frlrFR{1,kk}(frlrFRidx{1,kk}{1,jj},2),tuningbins))';
   end
end
frlrTuningSum = cell(1, size(frlrTuning,2));
for kk = 1: size(frlrTuning,2)
    frlrTuningSum{1,kk} = sum(frlrTuning{1,kk},2);
end
frlrOccup = zeros(size(tuningbins,1), size(frlrFR,2));
for kk = 1: size(frlrFR,2)
   frlrOccup(:,kk) = (hist(frlrFR{1,kk}(:,2),tuningbins).*NewSr)';
end
frlrOccup = sum(frlrOccup,2);
frlrTuningCurve = zeros(size(tuningbins,1), size(frlrTuningSum,2));
for kk = 1: size(frlrTuningSum,2)
    frlrTuningCurve(:,kk) = frlrTuningSum{1,kk}./frlrOccup;
end
frlrTuningCurve(isnan(frlrTuningCurve)) = 0;
frlrTuningCurveSm = zeros((size(frlrTuningCurve,1)), size(frlrTuningCurve,2));
for kk = 1: size(frlrTuningCurve,2)
   frlrTuningCurveSm(:,kk) = (smoothdata(frlrTuningCurve(:,kk),'gaussian',6))+.01; % Adding 0.01Hz
end

for kk = 1: size(frlrTuningCurve,2)
   s(kk) = subplot(size(MUAvect,2),1,kk);
   hold on
   area(tuningbins, frlrTuningCurveSm(:,kk),'FaceColor','c')
    line('XData', [.9 .9], 'YData', [0 55], 'LineStyle', '-.', ...
    'LineWidth', 2, 'Color',[0.5 0.5 0.5])
    line('XData', [.35 .35], 'YData', [0 55], 'LineStyle', '-.', ...
    'LineWidth', 2, 'Color',[0.5 0.5 0.5])
   hold off
   axis([0 1.3 0 inf])
end
xlabel('Position [m]')
ylabel('FR[Hz]')
title(s(1),'Forced Left to Right','FontSize',12)
txt1 = 'Start ';
text(.2,115,txt1,'HorizontalAlignment','Center','FontSize', 8)
txt2 = 'Center ';
text(.6,115,txt2,'HorizontalAlignment','Center','FontSize', 8)
txt3 = 'Reward ';
text(1.1,115,txt3,'HorizontalAlignment','Center','FontSize', 8)
txt4 = '\leftarrow';
text(1.23,120,txt4,'HorizontalAlignment','Center','FontSize', 40)
axis(s(1),[0 inf 0 25])
axis(s(2),[0 inf 0 30])
axis(s(3),[0 inf 0 53])
axis(s(4),[0 inf 0 40])
axis(s(5),[0 inf 0 10])
axis(s(6),[0 inf 0 25])
axis(s(7),[0 inf 0 20])
axis(s(8),[0 inf 0 10])
axis(s(9),[0 inf 0 50])
axis(s(10),[0 inf 0 45])
axis(s(11),[0 inf 0 20])
axis(s(12),[0 inf 0 30])
axis(s(13),[0 inf 0 45])
axis(s(14),[0 inf 0 5])

%%
save('ForcedRRTuningCurves.mat','frrrFR', 'frrrFRidx', 'tuningbins', 'frrrTuningCurve', ...
    'frrrTuningCurveSm') 
save('ForcedLLTuningCurves.mat','frllFR', 'frllFRidx', 'tuningbins', 'frllTuningCurve', ...
    'frllTuningCurveSm')
save('ForcedRLTuningCurves.mat','frrlFR', 'frrlFRidx', 'tuningbins', 'frrlTuningCurve', ...
    'frrlTuningCurveSm') 
save('ForcedLRTuningCurves.mat','frlrFR', 'frlrFRidx', 'tuningbins', 'frlrTuningCurve', ...
    'frlrTuningCurveSm')
%%
close all
for i=1:size(frrrFR,2)
    figure(i)
    subplot(411)
    plot(frrrFR{1,i}(:,1), frrrFR{1,i}(:,2),'k', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('Position [m]')
    axis tight
    title('Forced Right to Right')
    
    subplot(412)
    plot(frrrFR{1,i}(:,1), frrrFR{1,i}(:,2)*(size(frrrFR{1,i},2)-8),'LineWidth',2,'Color', [0.7 0.7 0.7])
    for kk = 1: size(frrrFRidx{1,1},2)
        hold on
         plot(frrrFR{1,i}(frrrFRidx{1,i}{1,kk},1), frrrFR{1,i}(frrrFRidx{1,i}{1,kk},kk+5)+kk-1, '.K');
         hold off
    end
    axis([frrrFR{1,i}(1) frrrFR{1,i}(end,1) 0 kk+2])
    xlabel('Time [sec]')
    ylabel('Units #')
    
    subplot(413)
    plot(frrrFR{1,i}(:,1), frrrFR{1,i}(:,3),'r', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('m/sec')
    axis tight
    
    subplot(414)
    hold on
    plot(frrrFR{1,i}(:,1), frrrFR{1,i}(:,4),'b')
    hold on
    plot(frrrFR{1,i}(:,1), frrrFR{1,i}(:,5),'k', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('LPF[uV]')
    axis tight
end
%%
close all
for i=1:size(frllFR,2)
    figure(i)
    subplot(411)
    plot(frllFR{1,i}(:,1), frllFR{1,i}(:,2),'k', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('Position [m]')
    axis tight
    title('Forced Left to Left')
    
    subplot(412)
    plot(frllFR{1,i}(:,1), frllFR{1,i}(:,2)*(size(frllFR{1,i},2)-8),'LineWidth',2,'Color', [0.7 0.7 0.7])
    for kk = 1: size(frllFRidx{1,1},2)
        hold on
         plot(frllFR{1,i}(frllFRidx{1,i}{1,kk},1), frllFR{1,i}(frllFRidx{1,i}{1,kk},kk+5)+kk-1, '.K');
         hold off
    end
    axis([frllFR{1,i}(1) frllFR{1,i}(end,1) 0 kk+2])
    xlabel('Time [sec]')
    ylabel('Units #')
    
    subplot(413)
    plot(frllFR{1,i}(:,1), frllFR{1,i}(:,3),'r', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('m/sec')
    axis tight
    
    subplot(414)
    hold on
    plot(frllFR{1,i}(:,1), frllFR{1,i}(:,4),'b')
    hold on
    plot(frllFR{1,i}(:,1), frllFR{1,i}(:,5),'k', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('LPF[uV]')
    axis tight
end
%%
close all
for i=1:size(frrlFR,2)
    figure(i)
    subplot(411)
    plot(frrlFR{1,i}(:,1), frrlFR{1,i}(:,2),'k', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('Position [m]')
    axis tight
    title('Forced Right to Left')
    
    subplot(412)
    plot(frrlFR{1,i}(:,1), frrlFR{1,i}(:,2)*(size(frrlFR{1,i},2)-8),'LineWidth',2,'Color', [0.7 0.7 0.7])
    for kk = 1: size(frrlFRidx{1,1},2)
        hold on
         plot(frrlFR{1,i}(frrlFRidx{1,i}{1,kk},1), frrlFR{1,i}(frrlFRidx{1,i}{1,kk},kk+5)+kk-1, '.K');
         hold off
    end
    axis([frrlFR{1,i}(1) frrlFR{1,i}(end,1) 0 kk+2])
    xlabel('Time [sec]')
    ylabel('Units #')
    
    subplot(413)
    plot(frrlFR{1,i}(:,1), frrlFR{1,i}(:,3),'r', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('m/sec')
    axis tight
    
    subplot(414)
    hold on
    plot(frrlFR{1,i}(:,1), frrlFR{1,i}(:,4),'b')
    hold on
    plot(frrlFR{1,i}(:,1), frrlFR{1,i}(:,5),'k', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('LPF[uV]')
    axis tight
end

%%
close all
for i=1:size(frlrFR,2)
    figure(i)
    subplot(411)
    plot(frlrFR{1,i}(:,1), frlrFR{1,i}(:,2),'k', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('Position [m]')
    axis tight
    title('Forced Left to Right')
    
    subplot(412)
    plot(frlrFR{1,i}(:,1), frlrFR{1,i}(:,2)*(size(frlrFR{1,i},2)-8),'LineWidth',2,'Color', [0.7 0.7 0.7])
    for kk = 1: size(frlrFRidx{1,1},2)
        hold on
         plot(frlrFR{1,i}(frlrFRidx{1,i}{1,kk},1), frlrFR{1,i}(frlrFRidx{1,i}{1,kk},kk+5)+kk-1, '.K');
         hold off
    end
    axis([frlrFR{1,i}(1) frlrFR{1,i}(end,1) 0 kk+2])
    xlabel('Time [sec]')
    ylabel('Units #')
    
    subplot(413)
    plot(frlrFR{1,i}(:,1), frlrFR{1,i}(:,3),'r', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('m/sec')
    axis tight
    
    subplot(414)
    hold on
    plot(frlrFR{1,i}(:,1), frlrFR{1,i}(:,4),'b')
    hold on
    plot(frlrFR{1,i}(:,1), frlrFR{1,i}(:,5),'k', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('LPF[uV]')
    axis tight
end














