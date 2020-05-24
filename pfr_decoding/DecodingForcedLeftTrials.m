BinTime = 0.2;
Timebins = cell(1,size(frllFR,2));
TimeTrialbins = cell(1,size(frllFR,2));
for kk = 1: size(frllFRidx,2)
    % Timebins to bin Spikes over time
    Timebins{1,kk} = (0: BinTime: (frllFR{1,kk}(end,1) - frllFR{1,kk}(1,1)))'; 
    % Timebins for bin Spikes over time/Trials
    TimeTrialbins{1,kk} = (frllFR{1,kk}(1,1): 0.2 : frllFR{1,kk}(end,1)); 
end
%% Histogram Spikes over time/trials
frllTimeHist = cell(1, size(frllFRidx,2));
for jj = 1: size(frllFRidx,2)
    for kk = 1: size(frllFRidx{1,1},2)
     frllTimeHist{1,jj}(:,kk) =(hist(frllFR{1,jj}(frllFRidx{1,jj}{1,kk},1), ...
         TimeTrialbins{1,jj}))';
    end
end

%% The exponetial multiplication and Sum of Units/MUA tuning curves 
frllSumFRExp = exp(-(BinTime).* sum(frllTuningCurveSm,2));

%% Firing rate elevated by time spikes
frllFR_eSpikes = cell(1, size(frllFRidx,2));
for jj = 1: size(frllFRidx,2)
    for tt = 1: length(frllTimeHist{1,jj}(:,1))
        frllFR_eSpikes{1, jj}(tt,:) = (prod(frllTuningCurveSm(:,:).^ ...
            frllTimeHist{1,jj}(tt,:),2))';
    end
end

%%
frllProb = cell(1, size(frllFRidx,2));
SumfrllProb = cell(1, size(frllFRidx,2));
frllLikelihood = cell(1, size(frllFRidx,2));
for kk = 1: size(frllFRidx,2)
    frllFR_eSpikes{1,kk} = (frllFR_eSpikes{1,kk})';
    frllProb{1,kk} = frllSumFRExp.* frllFR_eSpikes{1,kk};
    SumfrllProb{1,kk} = sum(frllProb{1,kk},1);
    frllLikelihood{1,kk} = frllProb{1,kk}./SumfrllProb{1,kk};
end

%% Visual Inspection
close all
for kk = 1: size(frllFRidx,2)
    figure(kk)
    hold on
    s = plot3(frllFR{1,kk}(:,1), frllFR{1,kk}(:,2), ...
        ones(size(frllFR{1,kk}(:,2),1),1),'-.r','LineWidth', 0.5);
    alpha(s,0.1)
    s1 = surf(TimeTrialbins{1,kk}, tuningbins,frllLikelihood{1,kk});
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