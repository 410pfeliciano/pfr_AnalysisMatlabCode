% Matrix Presleep
Fs = 1000;
sint = 1/Fs;
endrectime = 1609.6;
prebins = time;
prebins(prebins > 1609.6) = [];
presleepRippIdx= CA1pksRippleIdx;
presleepRippIdx(presleepRippIdx*sint > 1609.6) = [];
%%
presleepRippts = presleepRippIdx*sint;
rippletrain = hist(presleepRippts,prebins); 
% rippletrain(rippletrain ==1) = 2;
% rippletrain(rippletrain ==0) =1;
rippIdx = find(rippletrain);
%%
matTs = 0.1;% time to analyzed
tt = -(matTs):sint:0;
if ctxDownPreSlBelowIdx(end)+1 > endrectime
    ctxDownPreSlBelowIdx(end) = [];
elseif ctxDownPreSlBelowIdx(1)-length(tt) < prebins(1)
    ctxDownPreSlBelowIdx(1) = [];
else
end
rippMat = zeros(length(ctxDownPreSlBelowIdx),length(tt));
for kk = 1:length(ctxDownPreSlBelowIdx)
     rippMat(kk,1) = rippletrain(ctxDownPreSlBelowIdx(kk,1));
    for ll = 2: size(rippMat,2)
        rippMat(kk,ll) = rippletrain(1,ctxDownPreSlBelowIdx(kk)-ll+1);
    end
end 
%%
ttmat = 0:-sint:-(matTs);
figure
subplot(211)
plot(ttmat,rippMat(1,:),'k.')
hold on
for kk = 2:length(ctxDownPreSlBelowIdx)

    plot(ttmat,rippMat(kk,:)*kk,'k.')

end
hold off
axis([-inf inf 1 inf])

ylabel('Down-State #')
set(gca,'FontSize',11)
box off
subplot(212)
raterippmat = sum(rippMat)./length(ctxDownPreSlBelowIdx)./sint;
bar(ttmat,raterippmat,'k')

set(gca,'FontSize',11)
hold on
startPoints = [1 -0.03 0.001 0.001];
exclude = ttmat < -.06;
gaussEqn = 'a*exp(-((x-b)/c)^2)+d';
f = fit(ttmat',raterippmat',gaussEqn,'Start', startPoints, 'Exclude', exclude);
p = plot(f,'r');
p.LineWidth = 2;
ylabel('SWR rate ripple/sec')
xlabel('time rel. to Down-State(sec)')
box off
% legend off
hold off
%%
startPoints = [1 -0.03 0.01 0.1];
exclude = ttmat < -.06;
gaussEqn = 'a*exp(-((x-b)/c)^2)+d'
f=fit(ttmat',raterippmat',gaussEqn,'Start', startPoints, 'Exclude', exclude)
plot(f,ttmat',raterippmat')





