% Correlation Participation Up- Down-States MUA/SUA during SWR
meanParticipation = [correlationData1.participation.avg; correlationData2.participation.avg; ...
    correlationData3.participation.avg];
ratioParticipation = [correlationData1.participation.ratio; correlationData2.participation.ratio; ...
    correlationData3.participation.ratio];
CA1pyrIdx =  [correlationData1.SUAIdx.pyr; correlationData2.SUAIdx.pyr+...
    size(correlationData1.participation.avg,1); correlationData3.SUAIdx.pyr+...
    size(correlationData1.participation.avg,1)+size(correlationData2.participation.avg,1)];
CA1intIdx = [correlationData1.SUAIdx.int; correlationData2.SUAIdx.int+...
    size(correlationData1.participation.avg,1); correlationData3.SUAIdx.int+...
    size(correlationData1.participation.avg,1)+size(correlationData2.participation.avg,1)];
%% Group ceration
pyr = cellfun(@(c)[c 'pyramidal'],cell(size(CA1pyrIdx,1),1),'uni',false);
int = cellfun(@(c)[c 'interneuron'],cell(size(CA1intIdx,1),1),'uni',false);
neurontype = [pyr; int];
neuronMeanPart = [meanParticipation(CA1pyrIdx); meanParticipation(CA1intIdx)];
neuronRatioPart = [ratioParticipation(CA1pyrIdx); ratioParticipation(CA1intIdx)];
%%
close all
brobust = fitlm(neuronMeanPart,neuronRatioPart,'RobustOpts','on'); % robust linear fit
bleast = fitlm(neuronMeanPart,neuronRatioPart); % least squares fit
xVal = 0:0.01:1;
yRobust = table2array(brobust.Coefficients(1,1)) + xVal.*table2array(brobust.Coefficients(2,1));
yLeast =table2array(bleast.Coefficients(1,1)) + xVal.*table2array(bleast.Coefficients(2,1));
scatterhist(neuronMeanPart,neuronRatioPart,'Group',neurontype,'Kernel','on','Location','SouthEast',...
    'Direction','out','Color','grr','LineStyle',{'-','-.',':'},...
    'LineWidth',[2,2,2],'Marker','..d','MarkerSize',[18,18,6]);
xlabel('SUA SWR Participation Fraction')
ylabel({'Down-/Up-State SWR' , 'Participation Fraction'})
hold on
plot(xVal, yLeast, '--b', 'LineWidth',2,'DisplayName','Least Squares')
plot(xVal, yRobust, '-b', 'LineWidth',2,'DisplayName','Robust Regres.')
plot([0 ,1],[1 1],'--k', 'LineWidth',1)
hold off
% box off
legend boxoff  
set(gca,'FontSize',11)
%%
figure
p = polyfit(meanParticipation, ratioParticipation,1); 
f = polyval(p,meanParticipation); 
h = plot(meanParticipation, ratioParticipation,'.',meanParticipation,f,'b-') 
set(h, 'LineWidth',2, 'MarkerSize', 20, 'MarkerEdgeColor', [104/255, 107/255, 110/255],...
    'MarkerFaceColor', [104/255, 107/255, 110/255])
legend( h , 'MUA', 'Linear Fit', 'Interneuron', 'Location', 'SouthEast' );
hold on
hpyr = plot(meanParticipation(CA1pyrIdx), ratioParticipation(CA1pyrIdx),'g.','MarkerSize', 20, ...
    'DisplayName','Pyramidal')
hint = plot(meanParticipation(CA1intIdx), ratioParticipation(CA1intIdx),'r.','MarkerSize', 20, ...
    'DisplayName','Interneuron')
plot([0 ,1],[1 1],'--k', 'LineWidth',1)
hold off
box off
legend boxoff  
axis([ 0 1 0.4 1.4])
xlabel ('avg MUA SWR Participation')
ylabel ('Down-State/Up-State SWR Participation')
set(gca,'FontSize',11)
[rho pval] = corr(meanParticipation,ratioParticipation, 'Type', 'Spearman')
%% Scatter plot with Histogram
x = meas(:,1);
y = meas(:,2);
scatterhist(x,y,'Group',species,'Kernel','on','Location','SouthEast',...
    'Direction','out','Color','kbr','LineStyle',{'-','-.',':'},...
    'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5,6]);
%% Firinf Rate Correaltion
avgFRmua = [correlationData1.FiringRate.avg; correlationData2.FiringRate.avg];
ratioFRmua = [correlationData1.FiringRate.ratio; correlationData2.FiringRate.ratio];
%%
figure
p = polyfit(avgFRmua, ratioFRmua,1); 
f = polyval(p,avgFRmua); 
h = plot(avgFRmua, ratioFRmua,'.',avgFRmua,f,'b-') ;
set(h, 'LineWidth',2, 'MarkerSize', 20, 'MarkerEdgeColor', [104/255, 107/255, 110/255],...
    'MarkerFaceColor', [104/255, 107/255, 110/255]);
legend( h , 'MUA', 'Linear Fit', 'Interneuron', 'Location', 'NorthEast' );
hold on
hpyr = plot(avgFRmua(CA1pyrIdx), ratioFRmua(CA1pyrIdx),'g.','MarkerSize', 20, ...
    'DisplayName','Pyramidal');
hint = plot(avgFRmua(CA1intIdx), ratioFRmua(CA1intIdx),'r.','MarkerSize', 20, ...
    'DisplayName','Interneuron');
hold off
box off
legend boxoff  
axis([ 0 inf 0 1.8])
xlabel('avg Firing Rate(spk/sec)')
ylabel ('Down-State/Up-State Firing Rate(spk/sec)')
set(gca,'FontSize',12)
[rho pval] = corr(avgFRmua, ratioFRmua, 'Type', 'Spearman')