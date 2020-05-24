BinTime = 0.2;
Timebins = cell(1,size(frrrFR,2));
TimeTrialbins = cell(1,size(frrrFR,2));
for kk = 1: size(frrrFRidx,2)
    % Timebins to bin Spikes over time
    Timebins{1,kk} = (0: BinTime: (frrrFR{1,kk}(end,1) - frrrFR{1,kk}(1,1)))'; 
    % Timebins for bin Spikes over time/Trials
    TimeTrialbins{1,kk} = (frrrFR{1,kk}(1,1): 0.2 : frrrFR{1,kk}(end,1)); 
end
%% Histogram Spikes over time/trials
frrrTimeHist = cell(1, size(frrrFRidx,2));
for jj = 1: size(frrrFRidx,2)
    for kk = 1: size(frrrFRidx{1,1},2)
     frrrTimeHist{1,jj}(:,kk) =(hist(frrrFR{1,jj}(frrrFRidx{1,jj}{1,kk},1), ...
         TimeTrialbins{1,jj}))';
    end
end

%% The exponetial multiplication and Sum of Units/MUA tuning curves 
frrrSumFRExp = exp(-(BinTime).* sum(frrrTuningCurveSm,2));

%% Firing rate elevated by time spikes
frrrFR_eSpikes = cell(1, size(frrrFRidx,2));
for jj = 1: size(frrrFRidx,2)
    for tt = 1: length(frrrTimeHist{1,jj}(:,1))
        frrrFR_eSpikes{1, jj}(tt,:) = (prod(frrrTuningCurveSm(:,:).^ ...
            frrrTimeHist{1,jj}(tt,:),2))';
    end
end

%%
frrrProb = cell(1, size(frrrFRidx,2));
SumfrrrProb = cell(1, size(frrrFRidx,2));
frrrLikelihood = cell(1, size(frrrFRidx,2));
for kk = 1: size(frrrFRidx,2)
    frrrFR_eSpikes{1,kk} = (frrrFR_eSpikes{1,kk})';
    frrrProb{1,kk} = frrrSumFRExp.* frrrFR_eSpikes{1,kk};
    SumfrrrProb{1,kk} = sum(frrrProb{1,kk},1);
    frrrLikelihood{1,kk} = frrrProb{1,kk}./SumfrrrProb{1,kk};
end

%% Visual Inspection
close all
for kk = 1: size(frrrFRidx,2)
    figure(kk)
    hold on
    s = plot3(frrrFR{1,kk}(:,1), frrrFR{1,kk}(:,2), ...
        ones(size(frrrFR{1,kk}(:,2),1),1),'-.r','LineWidth', 0.5);
    alpha(s,0.1)
    s1 = surf(TimeTrialbins{1,kk}, tuningbins,frrrLikelihood{1,kk});
    s1.EdgeColor = 'none';
    c = colorbar;
    c.Label.String = 'Probability';
    colormap(flipud(bone))
    % caxis auto
    axis([TimeTrialbins{1,kk}(1) TimeTrialbins{1,kk}(end) 0.025 1.3])
    caxis([0 .4])
    xlabel('Time [sec]')
    ylabel('Linear Position [m]')
    zlabel('Probability') 
    legend(s,'Actual Position','Location','northeast')
    box on
    hold off
end