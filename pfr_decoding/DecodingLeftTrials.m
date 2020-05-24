BinTime = 0.2;
Timebins = cell(1,size(llFR,2));
TimeTrialbins = cell(1,size(llFR,2));
for kk = 1: size(llFRidx,2)
    % Timebins to bin Spikes over time
    Timebins{1,kk} = (0: BinTime: (llFR{1,kk}(end,1) - llFR{1,kk}(1,1)))'; 
    % Timebins for bin Spikes over time/Trials
    TimeTrialbins{1,kk} = (llFR{1,kk}(1,1): 0.2 : llFR{1,kk}(end,1)); 
end
%% Histogram Spikes over time/trials
llTimeHist = cell(1, size(llFRidx,2));
for jj = 1: size(llFRidx,2)
    for kk = 1: size(llFRidx{1,1},2)
     llTimeHist{1,jj}(:,kk) =(hist(llFR{1,jj}(llFRidx{1,jj}{1,kk},1), ...
         TimeTrialbins{1,jj}))';
    end
end

%% The exponetial multiplication and Sum of Units/MUA tuning curves 
llSumFRExp = exp(-(BinTime).* sum(llTuningCurveSm,2));

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
    axis([TimeTrialbins{1,kk}(1) TimeTrialbins{1,kk}(end) 0.025 1.3])
    caxis([0 .4])
    xlabel('Time [sec]')
    ylabel('Linear Position [m]')
    zlabel('Probability') 
    legend(s,'Actual Position','Location','northwest')
    box on
    hold off
end