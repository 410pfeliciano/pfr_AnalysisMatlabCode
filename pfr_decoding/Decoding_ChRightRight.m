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
%% Position Time histogram
Timebins = (0: 0.2: (rrFR{1,7}(end,1) - rrFR{1,7}(1,1)))'; % Timebins to bin position over time
Timeb = (rrFR{1,7}(1,1): 0.2 : rrFR{1,7}(end,1)); % Timebins for making the
% PosTimeHist = zeros(size(Timebins,1),1); % Memory Preallocation
% [a ,b] = hist(rrFR{1,7}(1:Timebins(1+1)*2000, 2),1); % Fist bin hist
% PosTimeHist(1,1) = b;
% for kk = 2: length(Timebins)-1
%     [a ,b] = hist(rrFR{1,7}(Timebins(kk)*2000:Timebins(kk+1)*2000,2),1);
%     PosTimeHist(kk,1) = b;
% end
% [a ,b] = hist(rrFR{1,7}(Timebins(50)*2000:Timebins(end)*2000, 2),1); % Fist bin hist
% PosTimeHist(51,1) = b;
%% Spike matrix 300ms
rrTimeHist = zeros(size(Timebins,1), size(rrFRidx{1,7},2));
for kk = 1: size(rrTimeHist,2)
    rrTimeHist(:,kk) =(hist(rrFR{1,7}(rrFRidx{1,7}{1,kk},1),Timeb))';%Hist positions @ spikes.
end

%% Summing the firing rates per position
rrSumFRate = (sum(rrTuningCurveSm,2));
rrSumFRate2 = (exp(-0.2.* rrSumFRate));
% rrSumFRate2 = -0.2.*(exp(rrSumFRate));
% Visual Inspection
figure
plot(tuningbins, rrSumFRate)
figure
plot(tuningbins, rrSumFRate2)
%% Firing rate elevated by time spikes
rrPosTime = zeros(length(rrTimeHist(:,1)), size(rrTuningCurveSm,1));
for tt = 1: length(rrTimeHist(:,1))
    rrPosTime(tt,:) = (prod(rrTuningCurveSm(:,:).^rrTimeHist(tt,:),2))';
end
%%
rrPosTime = rrPosTime';
%%
rrProb = rrSumFRate2.* rrPosTime;
%%
rrSumProb = sum(rrProb,1);
%%
rrProb3 = rrProb./rrSumProb;
%%
close all
figure
s = plot3(rrFR{1,7}(:,1), rrFR{1,7}(:,2),ones(size(rrFR{1,7}(:,2),1),1),'--r','LineWidth', 1);
hold on
s1 = surf(Timeb, tuningbins,rrProb3);
s1.EdgeColor = 'none';
c = colorbar;
c.Label.String = 'Probability';
colormap(flipud(bone))
% caxis auto
axis([Timeb(1) Timeb(end) 0.025 1.3])
caxis([0 .4])
xlabel('Time [sec]')
ylabel('Linear Position [m]')
zlabel('Probability') 
legend('Actual Position')
box on

