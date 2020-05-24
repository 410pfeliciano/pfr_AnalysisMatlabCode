BinTime = 0.2;
Timebins = cell(1,size(frrlFR,2));
TimeTrialbins = cell(1,size(frrlFR,2));
for kk = 1: size(frrlFRidx,2)
    % Timebins to bin Spikes over time
    Timebins{1,kk} = (0: BinTime: (frrlFR{1,kk}(end,1) - frrlFR{1,kk}(1,1)))'; 
    % Timebins for bin Spikes over time/Trials
    TimeTrialbins{1,kk} = (frrlFR{1,kk}(1,1): 0.2 : frrlFR{1,kk}(end,1)); 
end
%% Histogram Spikes over time/trials
frrlTimeHist = cell(1, size(frrlFRidx,2));
for jj = 1: size(frrlFRidx,2)
    for kk = 1: size(frrlFRidx{1,1},2)
     frrlTimeHist{1,jj}(:,kk) =(hist(frrlFR{1,jj}(frrlFRidx{1,jj}{1,kk},1), ...
         TimeTrialbins{1,jj}))';
    end
end

%% The exponetial multiplication and Sum of Units/MUA tuning curves 
frrlSumFRExp = exp(-(BinTime).* sum(frrlTuningCurveSm,2));

%% Firing rate elevated by time spikes
frrlFR_eSpikes = cell(1, size(frrlFRidx,2));
for jj = 1: size(frrlFRidx,2)
    for tt = 1: length(frrlTimeHist{1,jj}(:,1))
        frrlFR_eSpikes{1, jj}(tt,:) = (prod(frrlTuningCurveSm(:,:).^ ...
            frrlTimeHist{1,jj}(tt,:),2))';
    end
end

%%
frrlProb = cell(1, size(frrlFRidx,2));
SumfrrlProb = cell(1, size(frrlFRidx,2));
frrlLikelihood = cell(1, size(frrlFRidx,2));
for kk = 1: size(frrlFRidx,2)
    frrlFR_eSpikes{1,kk} = (frrlFR_eSpikes{1,kk})';
    frrlProb{1,kk} = frrlSumFRExp.* frrlFR_eSpikes{1,kk};
    SumfrrlProb{1,kk} = sum(frrlProb{1,kk},1);
    frrlLikelihood{1,kk} = frrlProb{1,kk}./SumfrrlProb{1,kk};
end

%% Visual Inspection
close all
for kk = 1: size(frrlFRidx,2)
    figure(kk)
    hold on
    s = plot3(frrlFR{1,kk}(:,1), frrlFR{1,kk}(:,2), ...
        ones(size(frrlFR{1,kk}(:,2),1),1),'-.r','LineWidth', 0.5);
    alpha(s,0.1)
    s1 = surf(TimeTrialbins{1,kk}, tuningbins,frrlLikelihood{1,kk});
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