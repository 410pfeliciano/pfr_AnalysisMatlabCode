BinTime = 0.2;
Timebins = cell(1,size(rrFR,2));
TimeTrialbins = cell(1,size(rrFR,2));
for kk = 1: size(rrFRidx,2)
    % Timebins to bin Spikes over time
    Timebins{1,kk} = (0: BinTime: (rrFR{1,kk}(end,1) - rrFR{1,kk}(1,1)))'; 
    % Timebins for bin Spikes over time/Trials
    TimeTrialbins{1,kk} = (rrFR{1,kk}(1,1): 0.2 : rrFR{1,kk}(end,1)); 
end
%% Histogram Spikes over time/trials
rrTimeHist = cell(1, size(rrFRidx,2));
for jj = 1: size(rrFRidx,2)
    for kk = 1: size(rrFRidx{1,1},2)
     rrTimeHist{1,jj}(:,kk) =(hist(rrFR{1,jj}(rrFRidx{1,jj}{1,kk},1), ...
         TimeTrialbins{1,jj}))';
    end
end

%% The exponetial multiplication and Sum of Units/MUA tuning curves 
rrSumFRExp = exp(-(BinTime).* sum(rrTuningCurveSm,2));

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
    colormap(flipud(bone))
    % caxis auto
    axis([TimeTrialbins{1,kk}(1) TimeTrialbins{1,kk}(end) 0 1.3])
    caxis([0 .6])
    xlabel('Time [sec]')
    ylabel('Linear Position [m]')
    zlabel('Probability') 
    legend(s,'Actual Position','Location','northwest')
    box on
    hold off
end