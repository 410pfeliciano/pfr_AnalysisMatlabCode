BinTime = 0.2;
Timebins = cell(1,size(frlrFR,2));
TimeTrialbins = cell(1,size(frlrFR,2));
for kk = 1: size(frlrFRidx,2)
    % Timebins to bin Spikes over time
    Timebins{1,kk} = (0: BinTime: (frlrFR{1,kk}(end,1) - frlrFR{1,kk}(1,1)))'; 
    % Timebins for bin Spikes over time/Trials
    TimeTrialbins{1,kk} = (frlrFR{1,kk}(1,1): 0.2 : frlrFR{1,kk}(end,1)); 
end
%% Histogram Spikes over time/trials
frlrTimeHist = cell(1, size(frlrFRidx,2));
for jj = 1: size(frlrFRidx,2)
    for kk = 1: size(frlrFRidx{1,1},2)
     frlrTimeHist{1,jj}(:,kk) =(hist(frlrFR{1,jj}(frlrFRidx{1,jj}{1,kk},1), ...
         TimeTrialbins{1,jj}))';
    end
end

%% The exponetial multiplication and Sum of Units/MUA tuning curves 
frlrSumFRExp = exp(-(BinTime).* sum(frlrTuningCurveSm,2));

%% Firing rate elevated by time spikes
frlrFR_eSpikes = cell(1, size(frlrFRidx,2));
for jj = 1: size(frlrFRidx,2)
    for tt = 1: length(frlrTimeHist{1,jj}(:,1))
        frlrFR_eSpikes{1, jj}(tt,:) = (prod(frlrTuningCurveSm(:,:).^ ...
            frlrTimeHist{1,jj}(tt,:),2))';
    end
end

%%
frlrProb = cell(1, size(frlrFRidx,2));
SumfrlrProb = cell(1, size(frlrFRidx,2));
frlrLikelihood = cell(1, size(frlrFRidx,2));
for kk = 1: size(frlrFRidx,2)
    frlrFR_eSpikes{1,kk} = (frlrFR_eSpikes{1,kk})';
    frlrProb{1,kk} = frlrSumFRExp.* frlrFR_eSpikes{1,kk};
    SumfrlrProb{1,kk} = sum(frlrProb{1,kk},1);
    frlrLikelihood{1,kk} = frlrProb{1,kk}./SumfrlrProb{1,kk};
end

%% Visual Inspection
close all
for kk = 1: size(frlrFRidx,2)
    figure(kk)
    hold on
    s = plot3(frlrFR{1,kk}(:,1), frlrFR{1,kk}(:,2), ...
        ones(size(frlrFR{1,kk}(:,2),1),1),'-.r','LineWidth', 0.5);
    alpha(s,0.1)
    s1 = surf(TimeTrialbins{1,kk}, tuningbins,frlrLikelihood{1,kk});
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