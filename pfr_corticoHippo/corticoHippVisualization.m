%% Visual Inspection of Cortical activity based on ripple detection
% close all
% eventExt = 0.1 / sint;
% for kk = 900:950 %: size(,1)
%    figure
%    subplot(611)
%    plot(time(CA1Data.rippleIdx.CA1RipSlBelowIdx(kk) - eventExt : CA1Data.rippleIdx.CA1RipSlAbovIdx(kk)+ eventExt), ...
%        CTXLFP(CA1Data.rippleIdx.CA1RipSlBelowIdx(kk) - eventExt : CA1Data.rippleIdx.CA1RipSlAbovIdx(kk)+ eventExt),'k')
%    axis([-inf inf, -inf inf])
%    
%    subplot(612)
%    plot(time(CA1Data.rippleIdx.CA1RipSlBelowIdx(kk) - eventExt : CA1Data.rippleIdx.CA1RipSlAbovIdx(kk)+ eventExt),...
%        ctxData.lfps.muaEnvp.smooth(CA1Data.rippleIdx.CA1RipSlBelowIdx(kk) - eventExt : CA1Data.rippleIdx.CA1RipSlAbovIdx(kk)+ eventExt),'k')
%    axis([-inf inf, 0 inf])
%    
%    subplot(613)
%    plot(time(CA1Data.rippleIdx.CA1RipSlBelowIdx(kk) - eventExt : CA1Data.rippleIdx.CA1RipSlAbovIdx(kk)+ eventExt), ...
%        CA1LFP(CA1Data.rippleIdx.CA1RipSlBelowIdx(kk) - eventExt : CA1Data.rippleIdx.CA1RipSlAbovIdx(kk)+ eventExt),'r')
%    axis([-inf inf, -inf inf])
%    
%    subplot(614)
%    plot(time(CA1Data.rippleIdx.CA1RipSlBelowIdx(kk) - eventExt : CA1Data.rippleIdx.CA1RipSlAbovIdx(kk)+ eventExt), ...
%        CA1LFPripple(CA1Data.rippleIdx.CA1RipSlBelowIdx(kk) - eventExt : CA1Data.rippleIdx.CA1RipSlAbovIdx(kk)+ eventExt),'r')
%    axis([-inf inf, -inf inf])
%    
%    subplot(615)
%    plot(time(CA1Data.rippleIdx.CA1RipSlBelowIdx(kk) - eventExt : CA1Data.rippleIdx.CA1RipSlAbovIdx(kk)+ eventExt),...
%        CA1Data.lfps.rippleEnvZscore(CA1Data.rippleIdx.CA1RipSlBelowIdx(kk) - eventExt : CA1Data.rippleIdx.CA1RipSlAbovIdx(kk)+ eventExt),'r')
%    axis([-inf inf, -inf inf])
%    
%    subplot(616)
%     plot(time(CA1Data.rippleIdx.CA1RipSlBelowIdx(kk) - eventExt : CA1Data.rippleIdx.CA1RipSlAbovIdx(kk)+ eventExt),... 
%         CA1Data.spikes.density(CA1Data.rippleIdx.CA1RipSlBelowIdx(kk) - eventExt : CA1Data.rippleIdx.CA1RipSlAbovIdx(kk)+ eventExt),... 
%         'r','LineWidth',2) 
%     xlabel('Time[sec]')
%     ylabel('spks/s')
%     axis tight
%    axis([-inf inf, 0 inf])
% end