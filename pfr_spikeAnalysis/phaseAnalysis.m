NewSr = 1/2000; % 2Khz sampling rate after downsampling the LFP signal
% Spiking Activity with velocity larger than 
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
%% Choice Trials 
rrFRLFP = cell(1,size(ch_run,1)); 
llFRLFP = cell(1,size(ch_run,1));
for jj = 1:size(ch_run,1)
    if ch_run{jj,3} == 'R' && ch_run{jj,4} == 'R'
        ind1 = ch_run{jj,1} * (1/NewSr);
        ind2 = ch_run{jj,2} * (1/NewSr);
        rrFRLFP{1,jj}(:,1) = linTime(ind1: ind2); % Time
        rrFRLFP{1,jj}(:,2) = linPos(ind1: ind2); % Position
        rrFRLFP{1,jj}(:,3) = linVel(ind1: ind2); % Velocity
        rrFRLFP{1,jj}(:,4) = LFP7(ind1: ind2); %LFP
        rrFRLFP{1,jj}(:,5) = thetaLFP7(ind1: ind2); %Theta filtered LFP
        rrFRLFP{1,jj}(:,6) = slowGammaLFP7(ind1: ind2); %Slow Gamma filtered LFP
        rrFRLFP{1,jj}(:,7) = fastGammaLFP7(ind1: ind2); %Fast Gamma filtered LFP
        for kk = 1: size(MUAvectActive,2) % Spikes(0 or 1)
            rrFRLFP{1,jj}(:,7+kk) = MUAvectActive(ind1:ind2,kk);
        end
    elseif ch_run{jj,3} == 'L' && ch_run{jj,4} == 'L'
        ind1 = ch_run{jj,1} * (1/NewSr);
        ind2 = ch_run{jj,2} * (1/NewSr);
        llFRLFP{1,jj}(:,1) = linTime(ind1: ind2);
        llFRLFP{1,jj}(:,2) = linPos(ind1: ind2);
        llFRLFP{1,jj}(:,3) = linVel(ind1: ind2);
        llFRLFP{1,jj}(:,4) = LFP7(ind1: ind2);
        llFRLFP{1,jj}(:,5) = thetaLFP7(ind1: ind2);
        llFRLFP{1,jj}(:,6) = slowGammaLFP7(ind1: ind2);
        llFRLFP{1,jj}(:,7) = fastGammaLFP7(ind1: ind2);
        for kk = 1: size(MUAvectActive,2)
            llFRLFP{1,jj}(:,7+kk) = MUAvectActive(ind1:ind2,kk);
        end
    end
end
% Removing [] cell arrays
rrFRLFP = rrFRLFP(~cellfun(@isempty,rrFRLFP)); % Right to Right Choice
llFRLFP = llFRLFP(~cellfun(@isempty,llFRLFP)); % Left to Left Choice
% Finding Spikes
rrFRLFPidx = cell(1, size(rrFRLFP,2));
for jj = 1 : size(rrFRLFP,2)
    for kk = 1 : size(MUAvectActive,2)
        rrFRLFPidx{1,jj}{:,kk} = find(rrFRLFP{1,jj}(:, kk+7));
    end
end
llFRLFPidx = cell(1, size(llFRLFP,2));
for jj = 1 : size(llFRLFP,2)
    for kk = 1 : size(MUAvectActive,2)
        llFRLFPidx{1,jj}{:,kk} = find(llFRLFP{1,jj}(:, kk+7));
    end
end
%% Phase Calculation
rrAngleRad = cell(1, size(rrFRLFP,2));
rrAngle2Pi = cell(1, size(rrFRLFP,2));
rrAngle = cell(1, size(rrFRLFP,2));
rrAmpGammaSlow = cell(1, size(rrFRLFP,2));
rrAmpGammaFast = cell(1, size(rrFRLFP,2));
llAngleRad = cell(1, size(llFRLFP,2));
llAngle2Pi = cell(1, size(llFRLFP,2));
llAngle = cell(1, size(llFRLFP,2));
llAmpGammaSlow = cell(1, size(llFRLFP,2));
llAmpGammaFast = cell(1, size(llFRLFP,2));
for kk = 1: size(rrFRLFP,2)
    rrAngleRad{1,kk} = angle(hilbert(rrFRLFP{1,kk}(:,5)));
    rrAngle2Pi{1,kk} = wrapTo2Pi(rrAngleRad{1,kk});
%     rrAngle{1,kk} = wrapTo360(rad2deg(rrAngleRad{1,kk}));
    rrAngle{1,kk} = rad2deg(rrAngleRad{1,kk});
    rrAmpGammaSlow{1,kk} = abs(hilbert(rrFRLFP{1,kk}(:,6)));
    rrAmpGammaFast{1,kk} = abs(hilbert(rrFRLFP{1,kk}(:,7)));
end
for kk = 1: size(llFRLFP,2)
    llAngleRad{1,kk} = angle(hilbert(llFRLFP{1,kk}(:,5)));
    llAngle2Pi{1,kk} = wrapTo2Pi(llAngleRad{1,kk});
%     llAngle{1,kk} = wrapTo360(rad2deg(llAngleRad{1,kk}));
    llAngle{1,kk} = rad2deg(llAngleRad{1,kk});
    llAmpGammaSlow{1,kk} = abs(hilbert(llFRLFP{1,kk}(:,6)));
    llAmpGammaFast{1,kk} = abs(hilbert(llFRLFP{1,kk}(:,7)));
end
%%
% save('ChoiceSpikesLFPAngle.mat','rrFRLFP','llFRLFP','rrAngle','llAngle',...
%     'rrAmpGammaSlow','rrAmpGammaFast','llAmpGammaSlow','llAmpGammaFast')
%% Visual Inspection
close all
for i=1:size(rrFRLFP,2)
    figure(i)
    subplot(711)
    line('XData', [rrFRLFP{1,i}(1,1) rrFRLFP{1,i}(end,1)], 'YData', ...
        [.9 .9], 'LineStyle', '-.','LineWidth', 2, 'Color',[0.5 0.5 0.5])
    line('XData', [rrFRLFP{1,i}(1,1) rrFRLFP{1,i}(end,1)], 'YData', ...
        [.35 .35], 'LineStyle', '-.','LineWidth', 2, 'Color',[0.5 0.5 0.5])
    hold on
    plot(rrFRLFP{1,i}(:,1), rrFRLFP{1,i}(:,2),'k', 'LineWidth', 2)
    hold off
    xlabel('Time [sec]')
    ylabel('Position [m]')
    title('Forced/Choice Right to Right')
    axis tight
    
    subplot(712)
    plot(rrFRLFP{1,i}(:,1), rrFRLFP{1,i}(:,2)*(size(rrFRLFP{1,i},2)-10), ...
        'LineWidth',2,'Color', [0.7 0.7 0.7])
    hold on
    plot(rrFRLFP{1,i}(rrFRLFPidx{1,i}{1,1},1), ... 
        rrFRLFP{1,i}(rrFRLFPidx{1,i}{1,1},7+1)+1-1, '.k');
    for kk = 2: size(rrFRLFPidx{1,1},2)
         plot(rrFRLFP{1,i}(rrFRLFPidx{1,i}{1,kk},1), ...
             rrFRLFP{1,i}(rrFRLFPidx{1,i}{1,kk},kk+7)+kk-1, '.k');
    end
             hold off
    axis([rrFRLFP{1,i}(1) rrFRLFP{1,i}(end,1) 0 kk+2])
    xlabel('Time [sec]')
    ylabel('Units #')
    
    subplot(713)
    plot(rrFRLFP{1,i}(:,1), rrFRLFP{1,i}(:,3),'r', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('m/sec')
    axis tight
    
    subplot(714)
    hold on
    plot(rrFRLFP{1,i}(:,1), rrFRLFP{1,i}(:,4),'b')
    plot(rrFRLFP{1,i}(:,1), rrFRLFP{1,i}(:,5),'k', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('LPF[uV]')
    axis tight
    
    subplot(715)
    hold on
    plot(rrFRLFP{1,i}(:,1), rrAngle{1,i},'k', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('Phase[deg]')
    axis tight
    
    subplot(716)
    hold on
    plot(rrFRLFP{1,i}(:,1), rrFRLFP{1,i}(:,5),'b')
    plot(rrFRLFP{1,i}(:,1), rrFRLFP{1,i}(:,6),'k', 'LineWidth', 1)
    plot(rrFRLFP{1,i}(:,1), rrAmpGammaSlow{1,i},'r', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('Slow \gamma [uV]')
    axis tight
    
    subplot(717)
    hold on
    plot(rrFRLFP{1,i}(:,1), rrFRLFP{1,i}(:,5),'b')
    plot(rrFRLFP{1,i}(:,1), rrFRLFP{1,i}(:,7),'k', 'LineWidth', 1)
    plot(rrFRLFP{1,i}(:,1), rrAmpGammaFast{1,i},'r', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('Fast \gamma [uV]')
    axis tight
end
%%
close all
for i=1:size(llFRLFP,2)
    figure(i)
    subplot(711)
    line('XData', [llFRLFP{1,i}(1,1) llFRLFP{1,i}(end,1)], 'YData', ...
        [.9 .9], 'LineStyle', '-.','LineWidth', 2, 'Color',[0.5 0.5 0.5])
    line('XData', [llFRLFP{1,i}(1,1) llFRLFP{1,i}(end,1)], 'YData', ...
        [.35 .35], 'LineStyle', '-.','LineWidth', 2, 'Color',[0.5 0.5 0.5])
    hold on
    plot(llFRLFP{1,i}(:,1), llFRLFP{1,i}(:,2),'k', 'LineWidth', 2)
    hold off
    xlabel('Time [sec]')
    ylabel('Position [m]')
    title('Forced/Choice Right to Right')
    axis tight
    
    subplot(712)
    plot(llFRLFP{1,i}(:,1), llFRLFP{1,i}(:,2)*(size(llFRLFP{1,i},2)-10) ...
        ,'LineWidth',2,'Color', [0.7 0.7 0.7])
    hold on
    plot(llFRLFP{1,i}(llFRLFPidx{1,i}{1,1},1), ... 
        llFRLFP{1,i}(llFRLFPidx{1,i}{1,1},7+1)+1-1, '.k');
    for kk = 2: size(llFRLFPidx{1,1},2)
         plot(llFRLFP{1,i}(llFRLFPidx{1,i}{1,kk},1), ...
             llFRLFP{1,i}(llFRLFPidx{1,i}{1,kk},kk+7)+kk-1, '.k');
    end
    hold off
    axis([llFRLFP{1,i}(1) llFRLFP{1,i}(end,1) 0 kk+2])
    xlabel('Time [sec]')
    ylabel('Units #')
    
    subplot(713)
    plot(llFRLFP{1,i}(:,1), llFRLFP{1,i}(:,3),'r', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('m/sec')
    axis tight
    
    subplot(714)
    hold on
    plot(llFRLFP{1,i}(:,1), llFRLFP{1,i}(:,4),'b')
    plot(llFRLFP{1,i}(:,1), llFRLFP{1,i}(:,5),'k', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('LPF[uV]')
    axis tight
    
    subplot(715)
    hold on
    plot(llFRLFP{1,i}(:,1), llAngle{1,i},'k', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('Phase[deg]')
    axis([-inf inf -inf inf])
    
    subplot(716)
    hold on
    plot(llFRLFP{1,i}(:,1), llFRLFP{1,i}(:,5),'b')
    plot(llFRLFP{1,i}(:,1), llFRLFP{1,i}(:,6),'k', 'LineWidth', 1)
    plot(llFRLFP{1,i}(:,1), llAmpGammaSlow{1,i},'r', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('Slow Gamma [uV]')
    axis tight
    
    subplot(717)
    hold on
    plot(llFRLFP{1,i}(:,1), llFRLFP{1,i}(:,5),'b')
    plot(llFRLFP{1,i}(:,1), llFRLFP{1,i}(:,7),'k', 'LineWidth', 1)
    plot(llFRLFP{1,i}(:,1), llAmpGammaFast{1,i},'r', 'LineWidth', 1)
    hold off
    xlabel('Time [sec]')
    ylabel('Fast Gamma[uV]')
    axis tight
end

%% Plotting spikes over position
close all
NumCells = 14;
for kk = 6%: NumCells
    figure
    for ii = 1: size(rrFRLFPidx,2)
        subplot(211)
        hold on
        scatter(rrFRLFP{1,ii}(rrFRLFPidx{1,ii}{1,kk},2), ... 
        rrAngle{1,ii}(rrFRLFPidx{1,ii}{1,kk},1), 'ko','filled');
        title(['Unit/MUA #', num2str(kk),'Choice/Forced Right to Right'])
        hold off
        xlabel('Position [m]')
        ylabel('Firing Phase[deg]')
        axis([.7 1.1 -inf inf])
%         subplot(312)
%         hold on
%         scatter(rrFRLFP{1,ii}(rrFRLFPidx{1,ii}{1,kk},2), ... 
%         rrAngle2Pi{1,ii}(rrFRLFPidx{1,ii}{1,kk},1), 'ko','filled');
%         title(['Unit/MUA #', num2str(kk),'Choice/Forced Right to Right'])
%         hold off
%         xlabel('Position [m]')
%         ylabel('Phase[radians]')
%         axis([.7 1.1 -inf inf])
        subplot(212)
         area(tuningbins, rrTuningCurveSm(:,kk),'FaceColor','k')
         xlabel('Position [m]')
        ylabel('FR [spikes/sec]')
         axis([.7 1.1 0 inf])
    end
end

%%
close all
for kk = 6%: NumCells
    figure
    for ii = 1: size(llFRLFPidx,2)
        subplot(211)
        hold on
        scatter(llFRLFP{1,ii}(llFRLFPidx{1,ii}{1,kk},2), ... 
        llAngle{1,ii}(llFRLFPidx{1,ii}{1,kk},1), 'ro','filled');
        title(['Unit/MUA #', num2str(kk),'Choice/Forced Right to Right'])
        hold off
        xlabel('Position [m]')
        ylabel('Firing Phase[deg]')
        axis([.7 1.1 -inf inf])
%         subplot(312)
%         hold on
%         scatter(llFRLFP{1,ii}(llFRLFPidx{1,ii}{1,kk},2), ... 
%         llAngle2Pi{1,ii}(llFRLFPidx{1,ii}{1,kk},1), 'ro','filled');
%         title(['Unit/MUA #', num2str(kk),'Choice/Forced Right to Right'])
%         hold off
%         xlabel('Position [m]')
%         ylabel('Phase[radians]')
%         axis([.6 1.2 -inf inf])
        subplot(212)
         area(tuningbins, llTuningCurveSm(:,kk),'FaceColor','r')
         xlabel('Position [m]')
        ylabel('FR [spikes/sec]')
         axis([.7 1.1 0 inf])
    end
end
%%
[pg, t, pbc]=phaseampogram(rrFRLFP{1,4}(:,4), 2000,[6 12],[60 100],[1 .1], 1, 1);
























