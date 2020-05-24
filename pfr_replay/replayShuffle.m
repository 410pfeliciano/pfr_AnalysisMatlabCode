a = MUAreactLikelihoodRR{25};
numShuff = 3000;
b = cell(1, numShuff);
%%
for kk = 1: numShuff
    b{kk} = Shuffle(a,1);
end
%%
trial = 25;
sepVect = [10 18];
rsqShuff = zeros(1,numShuff); % R squared
rhoShuff = zeros(1,numShuff);
for kk = 1: numShuff
    ReplayCandShuff = b{kk}(sepVect(1):sepVect(2));
    [pShuff,Sshuff] = polyfit(TimeTrialbins{trial}(sepVect(1):sepVect(2)),ReplayCandShuff,1);
    [y_fitshuff,deltashuff] = polyval(pShuff,TimeTrialbins{trial}(sepVect(1):sepVect(2)),Sshuff); 
    yresidShuff = ReplayCandShuff - y_fitshuff;
    SSresidShuff = sum(yresidShuff.^2);
    SStotalShuff = (length(ReplayCandShuff)-1) * var(ReplayCandShuff);
    rsqShuff(1,kk) = 1 - SSresidShuff/SStotalShuff;
    
    rhoShuff(1,kk) = corr(TimeTrialbins{trial}(sepVect(1):sepVect(2))',ReplayCandShuff');
end
%%
close all
rsqShuffHist = histogram(rsqShuff);
axis([0 1 0 inf])
figure
bar(rsqShuffHist.BinEdges(1,2:end),rsqShuffHist.BinCounts)
line([rho rho], ylim,'LineWidth',2,'Color','red');
txt = {'Fit from', 'actual Event\rightarrow'};
text(.78,500,txt,'Color','red')
title('Column Shuffle Bootstrap Distribution')
xlabel('Line Fit R^2 from Shuffle Data')
ylabel('Counts')
axis([0 1 0 inf])
% pd = fitdist(rsqShuff','exponential');
% figure
% h = histfit(rsqShuff,100,'HalfNormal')
% axis([0 1 0 500])
% x_values = 0:.1:1;
% y = pdf(pd,x_values);
% figure
% plot(x_values,y,'LineWidth',2)
% axis([0 1 0 inf])

%%
rho = corr(TimeTrialbins{trial}(sepVect(1):sepVect(2))',ReplayCandShuff')