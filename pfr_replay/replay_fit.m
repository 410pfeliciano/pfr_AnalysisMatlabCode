llIndLikeli = cell(1, size(MUAreactLikelihoodRR,2));
for kk = 1: size(MUAreactLikelihoodRR,2)
    [MaxLikeli, llIndLikeli{kk}] = max(MUAreactLikelihoodRR{kk});
end
llMaxLikeliPos = cell(1, size(MUAreactLikelihoodRR,2));
for kk = 1:size(MUAreactLikelihoodRR,2)
    llMaxLikeliPos{kk} = tuningbins(llIndLikeli{kk});
end
%%
close all
for kk = 1:size(MUAreactLikelihoodRR,2)
    figure
    scatter(TimeTrialbins{kk}, llMaxLikeliPos{kk}, 'o', 'filled')
    axis([-inf inf 0 1.3])
end
%%
trial = 25;
sepVect = [10 18];
ReplayCand = llMaxLikeliPos{trial}(sepVect(1):sepVect(2));
%%
[p,S] = polyfit(TimeTrialbins{trial}(sepVect(1):sepVect(2)),ReplayCand',1); 
%%
close all
[y_fit,delta] = polyval(p,TimeTrialbins{trial}(sepVect(1):sepVect(2)),S); 
scatter(TimeTrialbins{trial}(sepVect(1):sepVect(2)),ReplayCand,'ko','filled')
hold on
plot(TimeTrialbins{trial}(sepVect(1):sepVect(2)),y_fit,'r-')
plot(TimeTrialbins{trial}(sepVect(1):sepVect(2)),y_fit+2*delta,'r--',...
    TimeTrialbins{trial}(sepVect(1):sepVect(2)),y_fit-2*delta,'r--')
title('Linear Fit of Data with 95% Prediction Interval')
legend('Data','Linear Fit','95% Prediction Interval')
xlabel('Time[sec]')
ylabel('Estimated Position[m]')
axis([-inf inf 0 1.3])
%%
yresid = ReplayCand - y_fit';
SSresid = sum(yresid.^2);
SStotal = (length(ReplayCand)-1) * var(ReplayCand);
rsq = 1 - SSresid/SStotal
rsq_adj = 1 - SSresid/SStotal * (length(ReplayCand)-1)/(length(ReplayCand)-length(p))
%%
rho = corr(TimeTrialbins{trial}(sepVect(1):sepVect(2))',ReplayCand);
