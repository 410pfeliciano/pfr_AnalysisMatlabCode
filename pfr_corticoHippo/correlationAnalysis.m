%{
You need in the workspace units.mat, CA1rippleClasification.mat, CA1Data.mat and 
ctxData.mat files.
%}
% Duration Distribution
Fs = 1000;
sint = 1/Fs;

%% Calculating Ripple Duration
close all
ctxHipIntData.ripDur.upst.dur = CA1ripple.indices.upstate.all.end.*sint - ...
    CA1ripple.indices.upstate.all.start.*sint;
ctxHipIntData.ripDur.dwst.dur = CA1ripple.indices.downstate.all.end.*sint - ...
    CA1ripple.indices.downstate.all.start.*sint;

 % Box plot 
 group = [    ones(size(ctxHipIntData.ripDur.upst.dur));
         2 * ones(size(ctxHipIntData.ripDur.dwst.dur));];
figure
boxplot([ctxHipIntData.ripDur.upst.dur; ctxHipIntData.ripDur.dwst.dur],group, 'PlotStyle','traditional', ...
    'Colors',['b' 'r' ], 'Notch', 'on','Symbol', '.')
set(gca,'XTickLabel',{'Up- State','Down- State'},'FontSize',12)
fix_xticklabels(gca,1,{'FontSize',12});
% xtickangle(0)
ylabel('SWR Duration [sec]')
ylim([0 0.2])
box off
%  hUpstate = adtest(ctxHipIntData.ripDur.upst) % test normality
%  hDownstate = adtest(ctxHipIntData.ripDur.dwst) % test normality
[h,p,ci,stats] = ttest2(ctxHipIntData.ripDur.upst.dur, ctxHipIntData.ripDur.dwst.dur) % parametric t-test
[pp,ph,pstats] = ranksum(ctxHipIntData.ripDur.upst.dur, ctxHipIntData.ripDur.dwst.dur) % non parametric Wilcoxon rank sum test
STATS=mwwtest(ctxHipIntData.ripDur.upst.dur', ctxHipIntData.ripDur.dwst.dur') % non parametric Mann-Whitney-Wilcoxon Test
... https://github.com/dnafinder/mwwtest/blob/master/mwwtest.m
%  upstateR = iqr(ctxHipIntData.ripDur.upst)
%  downstateR = iqr(ctxHipIntData.ripDur.dwst)

figure % Cumulative distribution plot
hold on
[hh1,statsh1] = cdfplot(ctxHipIntData.ripDur.upst.dur);
[hh2,statsh2] = cdfplot(ctxHipIntData.ripDur.dwst.dur);
hh1.Color = 'b';
hh1.LineWidth = 2;
hh2.Color = 'r';
hh2.LineWidth = 2;
hold off
legend([hh1 hh2],{'Up-State','Down-State'}, 'Location', 'southeast')
grid off
legend boxoff
box off
xlim([-inf 0.15])
xlabel('SWR Duration (sec)')
ylabel('Cumulative Probability')
set(gca,'FontSize',12)
ctxHipIntData.ripDur.upst.Cum.Stat = statsh1;
ctxHipIntData.ripDur.dwst.Cum.Stat = statsh2;
[h,p,ks2stat] = kstest2(ctxHipIntData.ripDur.upst.dur,ctxHipIntData.ripDur.dwst.dur) 
h = lillietest(ctxHipIntData.ripDur.dwst.dur) % testing nortmality, 1=non-normal, 0=normal

% Hitogram Data
[ctxHipIntData.ripDur.upst.Cum.data, ctxHipIntData.ripDur.upst.Cum.cumEdges]  = ...
    histcounts(ctxHipIntData.ripDur.upst.dur,'BinWidth',0.005,'BinEdges',0:0.005:0.2,...
    'Normalization','cdf');
[ctxHipIntData.ripDur.dwst.Cum.data, ctxHipIntData.ripDur.dwst.Cum.cumEdges]  = ...
    histcounts(ctxHipIntData.ripDur.dwst.dur,'BinWidth',0.005,'BinEdges',0:0.005:0.2,...
    'Normalization','cdf');

figure
plot(ctxHipIntData.ripDur.upst.Cum.cumEdges(2:end),ctxHipIntData.ripDur.upst.Cum.data,'b', ...
    'LineWidth',2)
hold on
plot(ctxHipIntData.ripDur.dwst.Cum.cumEdges(2:end),ctxHipIntData.ripDur.dwst.Cum.data, 'r', ...
    'LineWidth',2)
hold off
xlabel('SWR Duration (sec)')
ylabel('Cumulative Probability')
set(gca,'FontSize',12)


%% Calculating area of spike density
areaHippSpkRateUp = zeros(size(CA1ripple.indices.upstate.all.end,1),1);
areaHippSpkRateDown = zeros(size(CA1ripple.indices.downstate.all.end,1),1);
for kk = 1: size(CA1ripple.indices.upstate.all.end,1)
areaHippSpkRateUp(kk) = trapz(CA1Data.spikes.density(CA1ripple.indices.upstate.all.start(kk):...
    CA1ripple.indices.upstate.all.end(kk)));
end
for kk = 1: size(CA1ripple.indices.downstate.all.end,1)
    areaHippSpkRateDown(kk) = trapz(CA1Data.spikes.density(CA1ripple.indices.downstate.all.start(kk):...
    CA1ripple.indices.downstate.all.end(kk)));
end
ctxHipIntData.FR.upst.area.data = areaHippSpkRateUp;
ctxHipIntData.FR.dwst.area.data = areaHippSpkRateDown;
close all
% Box plot 
figure
 group = [    ones(size(ctxHipIntData.FR.upst.area.data));
         2 * ones(size(ctxHipIntData.FR.dwst.area.data));];
boxplot([ctxHipIntData.FR.upst.area.data; ctxHipIntData.FR.dwst.area.data],group, 'PlotStyle','traditional', ...
    'Colors',['b' 'r' ], 'Notch', 'on','Symbol', '.')
set(gca,'XTickLabel',{'Up- State','Down- State'},'FontSize',12)
fix_xticklabels(gca,1,{'FontSize',12});
% xtickangle(25)
ylabel('SWR MUA Spk Rate Area')
ylim([0 inf])
box off
set(gca,'FontSize',12)
 hUpstate = adtest(ctxHipIntData.FR.upst.area.data) % test normality
 hDownstate = adtest(ctxHipIntData.FR.dwst.area.data) % test normality
[h,p,ci,stats] = ttest2(ctxHipIntData.FR.upst.area.data, ctxHipIntData.FR.dwst.area.data) % parametric t-test
[pp,ph,pstats] = ranksum(ctxHipIntData.FR.upst.area.data, ctxHipIntData.FR.dwst.area.data) % non parametric Wilcoxon rank sum test
STATS=mwwtest(ctxHipIntData.FR.upst.area.data', ctxHipIntData.FR.dwst.area.data') % non parametric Mann-Whitney-Wilcoxon Test
%

figure
hold on
[hh1,statsh1] = cdfplot(ctxHipIntData.FR.upst.area.data);
[hh2,statsh2] = cdfplot(ctxHipIntData.FR.dwst.area.data);
hh1.Color = 'b';
hh1.LineWidth = 2;
hh2.Color = 'r';
hh2.LineWidth = 2;
hold off
legend([hh1 hh2],{'Up-State','Down-State'}, 'Location', 'southeast')
grid off
legend boxoff
box off
xlabel('SWR MUA Spk Rate Area')
ylabel('Cumulative Probability')
set(gca,'FontSize',12)
xlim([0 1.5*10.^4])

ctxHipIntData.FR.upst.area.Cum.Stat = statsh1;
ctxHipIntData.FR.dwst.area.Cum.Stat = statsh2;
[h,p,ks2stat] = kstest2(ctxHipIntData.FR.upst.area.data, ctxHipIntData.FR.dwst.area.data) 
h = lillietest(ctxHipIntData.FR.dwst.area.data) % testing nortmality, 1=non-normal, 0=normal

% Hitogram Data
[ctxHipIntData.FR.upst.area.Cum.cumData, ctxHipIntData.FR.upst.area.Cum.cumEdges]  = ...
    histcounts(ctxHipIntData.FR.upst.area.data,'BinWidth',100,'BinEdges',0:100:2*10.^4,...
    'Normalization','cdf');
[ctxHipIntData.FR.dwst.area.Cum.cumData, ctxHipIntData.FR.dwst.area.Cum.cumEdges]  = ...
    histcounts(ctxHipIntData.FR.dwst.area.data,'BinWidth',100,'BinEdges',0:100:2*10.^4,...
    'Normalization','cdf');

figure
plot(ctxHipIntData.FR.upst.area.Cum.cumEdges(2:end),ctxHipIntData.FR.upst.area.Cum.cumData,'b')
hold on
plot(ctxHipIntData.FR.dwst.area.Cum.cumEdges(2:end),ctxHipIntData.FR.dwst.area.Cum.cumData, 'r')
hold off

%% Calculating ripple power for ripples during Up- and Down-states
rippPowerUp = zeros(size(CA1ripple.indices.upstate.all.end,1),1);
rippPowerDown = zeros(size(CA1ripple.indices.downstate.all.end,1),1);
for kk = 1: size(CA1ripple.indices.upstate.all.end,1)
rippPowerUp(kk) = trapz(CA1Data.lfps.rippleEnvSmooth(CA1ripple.indices.upstate.all.start(kk):...
    CA1ripple.indices.upstate.all.end(kk)));
end
for kk = 1: size(CA1ripple.indices.downstate.all.end,1)
    rippPowerDown(kk) = trapz(CA1Data.lfps.rippleEnvSmooth(CA1ripple.indices.downstate.all.start(kk):...
    CA1ripple.indices.downstate.all.end(kk)));
end
ctxHipIntData.rippPower.upst.area.data = rippPowerUp;
ctxHipIntData.rippPower.dwst.area.data = rippPowerDown;
close all
% Box plot 
figure
 group = [    ones(size(ctxHipIntData.rippPower.upst.area.data));
         2 * ones(size(ctxHipIntData.rippPower.dwst.area.data));];
boxplot([ctxHipIntData.rippPower.upst.area.data; ctxHipIntData.rippPower.dwst.area.data],group, 'PlotStyle','traditional', ...
    'Colors',['b' 'r' ], 'Notch', 'on','Symbol', '.')
set(gca,'XTickLabel',{'Up- State','Down- State'},'FontSize',12)
fix_xticklabels(gca,1,{'FontSize',12});
% xtickangle(25)
ylabel('SWR Envp. Power Area')
ylim([0 10.^7])
box off
set(gca,'FontSize',12)
 hUpstate = adtest(ctxHipIntData.rippPower.upst.area.data) % test normality
 hDownstate = adtest(ctxHipIntData.rippPower.dwst.area.data) % test normality
[h,p,ci,stats] = ttest2(ctxHipIntData.rippPower.upst.area.data, ctxHipIntData.rippPower.dwst.area.data) % parametric t-test
[pp,ph,pstats] = ranksum(ctxHipIntData.rippPower.upst.area.data, ctxHipIntData.rippPower.dwst.area.data) % non parametric Wilcoxon rank sum test
STATS=mwwtest(ctxHipIntData.rippPower.upst.area.data', ctxHipIntData.rippPower.dwst.area.data') % non parametric Mann-Whitney-Wilcoxon Test
%
figure
hold on
[hh1,statsh1] = cdfplot(ctxHipIntData.rippPower.upst.area.data);
[hh2,statsh2] = cdfplot(ctxHipIntData.rippPower.dwst.area.data);
hh1.Color = 'b';
hh1.LineWidth = 2;
hh2.Color = 'r';
hh2.LineWidth = 2;
hold off
legend([hh1 hh2],{'Up-State','Down-State'}, 'Location', 'southeast')
grid off
legend boxoff
box off
axis tight
xlim([0 10.^7])
xlabel('SWR Envp. Power Area')
ylabel('Cumulative Probability')
set(gca,'FontSize',12)

ctxHipIntData.rippPower.upst.area.Cum.Stat = statsh1;
ctxHipIntData.rippPower.dwst.area.Cum.Stat = statsh2;
[h,p,ks2stat] = kstest2(ctxHipIntData.rippPower.upst.area.data, ctxHipIntData.rippPower.dwst.area.data) 
h = lillietest(ctxHipIntData.rippPower.dwst.area.data) % testing nortmality, 1=non-normal, 0=normal

% Hitogram Data
[ctxHipIntData.rippPower.upst.area.Cum.cumData, ctxHipIntData.rippPower.upst.area.Cum.cumEdges]  = ...
    histcounts(ctxHipIntData.rippPower.upst.area.data,'BinWidth',10000,'BinEdges',0:1000:2*10.^7,...
    'Normalization','cdf');
[ctxHipIntData.rippPower.dwst.area.Cum.cumData, ctxHipIntData.rippPower.dwst.area.Cum.cumEdges]  = ...
    histcounts(ctxHipIntData.rippPower.dwst.area.data,'BinWidth',10000,'BinEdges',0:1000: 2*10.^7,...
    'Normalization','cdf');
figure
plot(ctxHipIntData.rippPower.upst.area.Cum.cumEdges(2:end),ctxHipIntData.rippPower.upst.area.Cum.cumData,'b')
hold on
plot(ctxHipIntData.rippPower.dwst.area.Cum.cumEdges(2:end),ctxHipIntData.rippPower.dwst.area.Cum.cumData, 'r')
hold off
 %% Hippocampal MUA during up- and down-states
 % Hippocampal Spike Density
 close all
Tet = CA1Data.lfps.raw; % You need the length of the recording to calculate thebins
bins = (0 : 1/Fs : length(Tet) * (1/Fs) - (1/Fs))';

if exist('spkClust','var') == 0 % spkClust var from kilosort analysis
    warning('spiking data spkClust is not in your workspace');
    return
end

TetCa1 = [1,2,3,4,8,9,10,11,12,15]; % choose the tetrode of interest

CA1Tetfind = cell(size(spkClust,2),1);
for kk = 1 : size(spkClust,2)
    CA1Tetfind{kk} = find(spkClust(kk).tetNum == TetCa1);
end
CA1TetIdx = find(~cellfun('isempty', CA1Tetfind));

neuronType = cell(size(spkClust,2),1);
for kk = 1 : size(spkClust,2)
    neuronType{kk} = spkClust(kk).neurontype;
end

neurontype = cell(size(CA1TetIdx,1),1);
CA1pyr = cell(size(CA1TetIdx,1),1);
CA1int = cell(size(CA1TetIdx,1),1);
for kk = 1: size(CA1TetIdx,1)
    neurontype{kk} = convertCharsToStrings(neuronType{CA1TetIdx(kk)});
    CA1pyr{kk} = strfind(neurontype{kk}, 'pyr');
    CA1int{kk} = strfind(neurontype{kk}, 'int');
end

CA1pyrIdx = find(~cellfun('isempty', CA1pyr));
CA1intIdx = find(~cellfun('isempty', CA1int));

%%
CA1spiketrain = zeros(size(CA1TetIdx,1),length(bins));
for kk = 1 : size(CA1TetIdx,1)
    CA1spiketime = (spkClust(CA1TetIdx(kk)).spkTime)';
    CA1spiketrain(kk,:) = hist(CA1spiketime,bins);
end

% Spike density estimation
CA1totalSpk = sum(CA1spiketrain,1);
sigma = 0.02;
edges = -3*sigma: 0.001: 3*sigma;
kernel = (normpdf(edges,0,sigma)).*(1/10);
% kernel = kernel*(1/Fs); % kernel muitiplied by bin size
CA1spkDensity = conv(CA1totalSpk, kernel);
center = ceil(length(edges)/2);
CA1spkDensity = CA1spkDensity(center:length(CA1totalSpk)+center-1);

 %% 50ms Before and After Up-state
 eventExt = 0.1/ sint;
CA1ripple.MUA.upstate.matrix  =  cell(size(CA1ripple.indices.upstate.all.end, 1 ),1);
for kk = 1 : length(CA1ripple.indices.upstate.all.start)
  CA1ripple.MUA.upstate.matrix {kk} = CA1spiketrain( :, CA1ripple.indices.upstate.all.start(kk) : ...
      CA1ripple.indices.upstate.all.end(kk));
end
 
% Firing rate calculation and participation
unitsNumber = size(CA1ripple.MUA.upstate.matrix{1},1);
CA1ripple.MUA.upstate.spks =  cell(size(CA1ripple.indices.upstate.all.peak, 1 ),1);
CA1ripple.MUA.upstate.firingRate =  cell(size(CA1ripple.indices.upstate.all.peak, 1 ),1);
CA1ripple.MUA.upstate.participation = zeros(length(CA1ripple.indices.upstate.all.peak),1);
for kk = 1 :  length(CA1ripple.MUA.upstate.matrix)
    CA1ripple.MUA.upstate.spks{kk} = sum(CA1ripple.MUA.upstate.matrix{kk},2);
    CA1ripple.MUA.upstate.firingRate{kk} = CA1ripple.MUA.upstate.spks{kk} / 0.1;
    CA1ripple.MUA.upstate.participation(kk) = length(find(CA1ripple.MUA.upstate.spks{kk}))/ ...
        unitsNumber;
end
%%
CA1ripple.MUA.downstate.matrix  =  cell(size(CA1ripple.indices.downstate.all.end, 1 ),1);
for kk = 1 : length(CA1ripple.indices.downstate.all.start)
  CA1ripple.MUA.downstate.matrix {kk} = CA1spiketrain( :, CA1ripple.indices.downstate.all.start(kk) : ...
      CA1ripple.indices.downstate.all.end(kk));
end
 
% Firing rate calculation and participation
unitsNumber = size(CA1ripple.MUA.downstate.matrix{1},1);
CA1ripple.MUA.downstate.spks =  cell(size(CA1ripple.indices.downstate.all.peak, 1 ),1);
CA1ripple.MUA.downstate.firingRate =  cell(size(CA1ripple.indices.downstate.all.peak, 1 ),1);
CA1ripple.MUA.downstate.participation = zeros(length(CA1ripple.indices.downstate.all.peak),1);
for kk = 1 :  length(CA1ripple.indices.downstate.all.peak)
    CA1ripple.MUA.downstate.spks{kk} = sum(CA1ripple.MUA.downstate.matrix{kk},2);
    CA1ripple.MUA.downstate.firingRate{kk} = CA1ripple.MUA.downstate.spks{kk} / 0.1;
    CA1ripple.MUA.downstate.participation(kk) = length(find(CA1ripple.MUA.downstate.spks{kk}))/ ...
        unitsNumber;
end
%% Overall MUA participation per SWR
% Box Plot CA1 MUA participation in SWR occuring during "Up- or Down-States"
close all
group = [    ones(size(CA1ripple.MUA.upstate.participation));
         2 * ones(size(CA1ripple.MUA.downstate.participation));];
figure
boxplot([CA1ripple.MUA.upstate.participation; ...
    CA1ripple.MUA.downstate.participation],group,'PlotStyle','traditional', ...
    'Colors',['b' 'r' ], 'Notch', 'on', 'Symbol','.')
set(gca,'XTickLabel',{'Up- State','Down- State'},'FontSize',12)
fix_xticklabels(gca,1,{'FontSize',12});
% xtickangle(25)
ylabel('MUA SWR Participation Fraction')
ylim([0 1])
box off
[h,p,ci,stats] = ttest2(CA1ripple.MUA.upstate.participation, CA1ripple.MUA.downstate.participation)
[pp,ph,pstats] = ranksum(CA1ripple.MUA.upstate.participation, ...
    CA1ripple.MUA.downstate.participation)
STATS=mwwtest(CA1ripple.MUA.upstate.participation', ...
    CA1ripple.MUA.downstate.participation') % non parametric Mann-Whitney-Wilcoxon Test

 ctxHipIntData.ripParticipation.upst.data = CA1ripple.MUA.upstate.participation;
 ctxHipIntData.ripParticipation.dwst.data = CA1ripple.MUA.downstate.participation;

% Hitogram Data
[ctxHipIntData.ripParticipation.upst.Cum.cumData, ctxHipIntData.ripParticipation.upst.Cum.cumEdges]  = ...
    histcounts(ctxHipIntData.ripParticipation.upst.data,'BinWidth',0.001,'BinEdges',0:0.001:1,...
    'Normalization','cdf');
[ctxHipIntData.ripParticipation.dwst.Cum.cumData, ctxHipIntData.ripParticipation.dwst.Cum.cumEdges]  = ...
    histcounts(ctxHipIntData.ripParticipation.dwst.data,'BinWidth',0.001,'BinEdges',0:0.001:1,...
    'Normalization','cdf');
figure
hold on
p1=plot(ctxHipIntData.ripParticipation.upst.Cum.cumEdges(2:end),...
    ctxHipIntData.ripParticipation.upst.Cum.cumData,'b', 'LineWidth',2,'DisplayName','Up-State')
p2 = plot(ctxHipIntData.ripParticipation.dwst.Cum.cumEdges(2:end),...
    ctxHipIntData.ripParticipation.dwst.Cum.cumData, 'r', 'LineWidth',2,'DisplayName','Down-State')
hold off
axis([0 1 0 1])
legend([p1 p2],{'Up-State','Down-State'}, 'Location', 'southeast')
legend boxoff
xlabel('SWR MUA Participation')
ylabel('Cumulative Probability')
box off
figure
hold on
[hh1,statsh1] = cdfplot(ctxHipIntData.ripParticipation.upst.data);
[hh2,statsh2] = cdfplot(ctxHipIntData.ripParticipation.dwst.data);
hh1.Color = 'b';
hh1.LineWidth = 2;
hh2.Color = 'r';
hh2.LineWidth = 2;
hold off
legend([hh1 hh2],{'Up-State','Down-State'}, 'Location', 'southeast')
grid off
legend boxoff
box off
axis tight
% xlim([0 10.^7])
xlabel('SWR MUA Participation')
ylabel('Cumulative Probability')
set(gca,'FontSize',12)
%% Individual Participation of MUA activity during Ridpples occuring during Up-states
close all
muaUpParticipation = cell(unitsNumber,1);
for kk = 1:unitsNumber% memory allocation
    muaUpParticipation{kk} = zeros(length(CA1ripple.MUA.upstate.spks),1);
end

  for kk = 1: unitsNumber
      for ll = 1: length(CA1ripple.MUA.upstate.spks)
          if CA1ripple.MUA.upstate.spks{ll}(kk) == 0
                muaUpParticipation{kk}(ll) = 0;
          else
              muaUpParticipation{kk}(ll) = 1;
          end
      end
  end
  muaUpMeanParticipation = zeros(unitsNumber,1);
  for kk = 1: unitsNumber
      muaUpMeanParticipation(kk) = sum(muaUpParticipation{kk}) /...
          length(muaUpParticipation{kk});
  end
  ctxHipIntData.MUA.participation.upst.part = muaUpMeanParticipation;
  %%
  muaDownParticipation = cell(unitsNumber,1);
for kk = 1:unitsNumber% memory allocation
    muaDownParticipation{kk} = zeros(length(CA1ripple.MUA.downstate.spks),1);
end

  for kk = 1: unitsNumber
      for ll = 1: length(CA1ripple.MUA.downstate.spks)
          if CA1ripple.MUA.downstate.spks{ll}(kk) == 0
                muaDownParticipation{kk}(ll) = 0;
          else
              muaDownParticipation{kk}(ll) = 1;
          end
      end
  end
  muaDownMeanParticipation = zeros(unitsNumber,1);
  for kk = 1: unitsNumber
      muaDownMeanParticipation(kk) = sum(muaDownParticipation{kk}) /...
          length(muaDownParticipation{kk});
  end
  ctxHipIntData.MUA.participation.dsst.part = muaDownMeanParticipation;
 %%
close all
% figure
% loglog(muaUpMeanParticipation, muaDownMeanParticipation,'o', 'LineWidth',1, 'MarkerSize',5, ...
%     'MarkerEdgeColor', [104/255, 107/255, 110/255], 'MarkerFaceColor', [104/255, 107/255, 110/255])
% hold on
% loglog([0.01 10], [0.01 10], 'k--')
% hold off
% box off
% axis([10.^-1 1 10.^-1 1])
% xlabel('SWR Unit Participation during Up-State')
% ylabel('SWR Unit Participation during Down-State')
% set(gca,'FontSize',11)
 
figure
p1=plot(muaUpMeanParticipation, muaDownMeanParticipation,'o', ...
 'LineWidth',2,...
    'MarkerSize',5,...
    'MarkerEdgeColor', [104/255, 107/255, 110/255],...
    'MarkerFaceColor', [104/255, 107/255, 110/255])
hold on
p2=plot(muaUpMeanParticipation(CA1pyrIdx), muaDownMeanParticipation(CA1pyrIdx), ...
    'g.', 'MarkerSize', 20)
p3=plot(muaUpMeanParticipation(CA1intIdx), muaDownMeanParticipation(CA1intIdx), ...
    'r.', 'MarkerSize', 20)
plot([0 1], [0 1], 'k--')
hold off
legend([p1 p2 p3],{'MUA','Pyramidal', 'Interneuron'}, 'Location', 'SouthEast')
legend('boxoff')
box off
axis([0 1 0 1])
xlabel('SWR Unit Participation during Up-State')
ylabel('SWR Unit Participation during Down-State')
set(gca,'FontSize',11)

group = [ones(size(muaUpMeanParticipation));
         2 * ones(size(muaDownMeanParticipation));];
     figure
boxplot([muaUpMeanParticipation; muaDownMeanParticipation],group,'PlotStyle','traditional', ...
    'Colors',['b' 'r' ], 'Notch', 'on','Symbol', 'k.')
hold on
p = parallelcoords([muaUpMeanParticipation, muaDownMeanParticipation], 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
  'Marker', '.', 'MarkerSize', 10);
hold off
set(gca,'XTickLabel',{'Up- State','Down- State'},'FontSize',12)
fix_xticklabels(gca,1,{'FontSize',12});
% xtickangle(25)
ylabel('SWR Unit Participation Fraction')
box off
% ylim([0 1])
hUp = adtest(muaUpMeanParticipation) % Test Normality
hDown = adtest(muaDownMeanParticipation) % Test Normality
[h,p,ci,stats] = ttest(muaUpMeanParticipation,muaDownMeanParticipation) % Parametric paired T-Test
%  [p,h,stats] = signrank(muaUpMeanParticipation,muaDownMeanParticipation) % Non-parametric paired test
%% 
close all
meanParticipation = (muaUpMeanParticipation + muaDownMeanParticipation)./2;
ctxHipIntData.MUA.participation.meanPart = meanParticipation;
ratioParticipation = muaDownMeanParticipation./muaUpMeanParticipation; 
ctxHipIntData.MUA.participation.ratioPart = ratioParticipation;
hratioParticipation = adtest(ratioParticipation); % Test Normality
[rhoPart, pvalPart] = corr(meanParticipation,ratioParticipation, 'Type', 'Spearman','tail','both');
ctxHipIntData.MUA.participation.corr.coeff = rhoPart;
ctxHipIntData.MUA.participation.corr.pval = pvalPart;
figure
[bPart ,statsPart] = robustfit(meanParticipation, ratioParticipation,'bisquare');
ctxHipIntData.MUA.participation.fit.coeff = bPart;
ctxHipIntData.MUA.participation.fit.stat = statsPart;
xx = 0:0.01:1;
yy = bPart(1)+bPart(2).*xx;
% rsquare_robustfit = corr(ratioParticipation,b(1)+b(2)*meanParticipation)^2
hold on
% p = polyfit(meanParticipation, ratioParticipation,1); 
% f = polyval(p,meanParticipation); 
h = plot(meanParticipation, ratioParticipation,'.') ;
set(h, 'LineWidth',2, 'MarkerSize', 20, 'MarkerEdgeColor', [104/255, 107/255, 110/255],...
    'MarkerFaceColor', [104/255, 107/255, 110/255]);
hold on
hpyr = plot(meanParticipation(CA1pyrIdx), ratioParticipation(CA1pyrIdx),'g.','MarkerSize', 20, ...
    'DisplayName','Pyramidal');
hint = plot(meanParticipation(CA1intIdx), ratioParticipation(CA1intIdx),'r.','MarkerSize', 20, ...
    'DisplayName','Interneuron');
hrobust = plot(xx,yy,'-b', 'LineWidth',2,'DisplayName','Robust Regression');
hline = plot([0 1], [1 1], '--k');
box off
legend( [h, hrobust] , {'MUA', 'Robust Regression'}, 'Location', 'SouthEast' );
legend boxoff  

hold off
axis([ 0 1 0.2 1.2])
xlabel ('MUA SWR Participation')
ylabel ('Down-State/Up-State SWR Participation')
set(gca,'FontSize',11)
%% Firing Rate of MUA activity during Ridpple occuring during Up-states
muaFiringRate = cell(unitsNumber,1);
for kk = 1:unitsNumber% memory allocation
    muaFiringRate{kk} = zeros(length(CA1ripple.MUA.upstate.firingRate),1);
end

  for kk = 1: unitsNumber
      for ll = 1: length(CA1ripple.MUA.upstate.firingRate)
      muaFiringRate{kk}(ll) = CA1ripple.MUA.upstate.firingRate{ll}(kk);
      end
  end
  muaFiringRateUp = zeros(unitsNumber,1);
  for kk = 1: unitsNumber
      muaFiringRateUp(kk) = mean(muaFiringRate{kk});
  end
  %% Firing Rate of MUA activity during Ridpple occuring during Down-states
muaFiringRate2 = cell(unitsNumber,1);
for kk = 1:unitsNumber% memory allocation
    muaFiringRate2{kk} = zeros(length(CA1ripple.MUA.downstate.firingRate),1);
end

  for kk = 1: unitsNumber
      for ll = 1: length(CA1ripple.MUA.downstate.firingRate)
      muaFiringRate2{kk}(ll) = CA1ripple.MUA.downstate.firingRate{ll}(kk);
      end
  end
  muaFiringRateDown = zeros(unitsNumber,1);
  for kk = 1: unitsNumber
      muaFiringRateDown(kk) = mean(muaFiringRate2{kk});
  end
%%
close all
% figure
% loglog(muaFiringRateUp, muaFiringRateDown,'o', ...
%  'LineWidth',1,...
%     'MarkerSize',5,...
%     'MarkerEdgeColor', [104/255, 107/255, 110/255],...
%     'MarkerFaceColor', [104/255, 107/255, 110/255])
% hold on
% loglog([0.01 1000], [0.01 1000], 'k--')
% hold off
% box off
% axis([10.^-1 100 10.^-1 100])
% xlabel('SWR MUA spks/sec during Up-State')
% ylabel('SWR MUA spks/sec during Down-State')
% set(gca,'FontSize',11)
 ctxHipIntData.MUA.FR.upst.muaFR = muaFiringRateUp;
 ctxHipIntData.MUA.FR.dwst.muaFR = muaFiringRateDown;
figure
 p1=plot(muaFiringRateUp, muaFiringRateDown,'o', ...
 'LineWidth',2,...
    'MarkerSize',5,...
    'MarkerEdgeColor', [104/255, 107/255, 110/255],...
    'MarkerFaceColor', [104/255, 107/255, 110/255])
hold on
p2=plot(muaFiringRateUp(CA1pyrIdx), muaFiringRateDown(CA1pyrIdx), ...
    'g.', 'MarkerSize', 20)
p3=plot(muaFiringRateUp(CA1intIdx), muaFiringRateDown(CA1intIdx), ...
    'r.', 'MarkerSize', 20)
plot([0 1000], [0 1000], 'k--')
hold off
legend([p1 p2 p3],{'MUA','Pyramidal', 'Interneuron'}, 'Location', 'SouthEast')
legend('boxoff')
box off
axis([0 70 0 70])
xlabel('SWR MUA spks/sec during Up-State')
ylabel('SWR MUA spks/sec during Down-State')
set(gca,'FontSize',11)

% Plot
% close all
% figure();
% coordLineStyle = 'k.';
% group = [    ones(size(muaFiringRateUp));
%          2 * ones(size(muaFiringRateDown));];
% boxplot([muaFiringRateUp; muaFiringRateDown], group, 'Symbol', coordLineStyle); 
% hold on;
% p = parallelcoords([muaFiringRateUp, muaFiringRateDown], 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
%   'Marker', '.', 'MarkerSize', 10);
% hold off
%
group = [    ones(size(muaFiringRateUp));
         2 * ones(size(muaFiringRateDown));];
figure
boxplot([muaFiringRateUp; muaFiringRateDown],group,'PlotStyle','traditional', ...
    'Colors',['b' 'r' ], 'Notch', 'on','Symbol', 'k.')
hold on
p = parallelcoords([muaFiringRateUp, muaFiringRateDown], 'Color', 0.5*[1 1 1], 'LineStyle', '-',...
  'Marker', '.', 'MarkerSize', 10);
hold off
set(gca,'XTickLabel',{'Up- State','Down- State'},'FontSize',12)
fix_xticklabels(gca,1,{'FontSize',12});
% xtickangle(25)
ylabel('SWR MUA Firing Rate (spks/sec)')
box off
% ylim([0 1])
 [h,p,ci,stats] = ttest(muaFiringRateUp,muaFiringRateDown)
 % Correaltion analysis
 avgFRmua = (muaFiringRateUp + muaFiringRateDown)./2;
 ratioFRmua = muaFiringRateDown./muaFiringRateUp;
 ctxHipIntData.MUA.FR.muaFRavg = muaFiringRateUp;
 ctxHipIntData.MUA.FR.muaFRratio = muaFiringRateDown;

[rhoFR, pvalFR] = corr(avgFRmua, ratioFRmua, 'Type', 'Pearson')
ctxHipIntData.MUA.FR.corr.coeff = rhoFR;
ctxHipIntData.MUA.FR.corr.pval = pvalFR;
[bFR,statsFR] = robustfit(avgFRmua, ratioFRmua,'bisquare');
ctxHipIntData.MUA.FR.fit.coeff = bFR;
ctxHipIntData.MUA.FR.fit.stat = statsFR;
xx = 0:1:90;
yy = bFR(1)+bFR(2).*xx;
% rsquare_robustfit = corr(ratioFRmua,b(1)+b(2)*avgFRmua)^2
figure
hold on
% p = polyfit(avgFRmua, ratioFRmua,1); 
% f = polyval(p,avgFRmua); 
h = plot(avgFRmua, ratioFRmua,'.') ;
set(h, 'LineWidth',2, 'MarkerSize', 20, 'MarkerEdgeColor', [104/255, 107/255, 110/255],...
    'MarkerFaceColor', [104/255, 107/255, 110/255]);
% legend( h , 'MUA', 'Location', 'SouthEast' );
hold on
hhrobust = plot(xx,yy,'-b', 'LineWidth',2,'DisplayName','Robust Regression');
hpyr = plot(avgFRmua(CA1pyrIdx), ratioFRmua(CA1pyrIdx),'g.','MarkerSize', 20, ...
    'DisplayName','Pyramidal');
hint = plot(avgFRmua(CA1intIdx), ratioFRmua(CA1intIdx),'r.','MarkerSize', 20, ...
    'DisplayName','Interneuron');
hline = plot([0 1], [1 1], '--k');
hold off
box off
legend( [h, hhrobust] , {'MUA', 'Robust Regression'}, 'Location', 'SouthEast' );
legend boxoff  
axis([ 0 inf 0.2 1.2])
xlabel ('MUA SWR Firing Rate spks/sec')
ylabel ({'Ratio Down-/Up-State SWR', 'Firing Rate spks/sec'})
set(gca,'FontSize',11)
 
%%
% correlationData.participation.avg = meanParticipation;
% correlationData.participation.ratio = ratioParticipation;
% correlationData.SUAIdx.int = CA1intIdx;
% correlationData.SUAIdx.pyr = CA1pyrIdx;
% correlationData.FiringRate.avg = avgFRmua;
% correlationData.FiringRate.ratio =ratioFRmua;
%%
save('correlationData.mat', 'ctxHipIntData')
save('CA1rippleClasification.mat', 'CA1ripple')
 
 