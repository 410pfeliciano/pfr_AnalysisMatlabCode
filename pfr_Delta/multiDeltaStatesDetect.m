%Multi Delta-States detection
% multievents that occur withing 200ms
deltaCanddur = belowCorrAmp .*sint;
interEvents = zeros(length(deltaCanddur),1);
interEvents(2:end) = diff(deltaCanddur);
%%
multideltaEvents = interEvents <= 0.25 & interEvents >= 0.02;
difMulti = diff(multideltaEvents);
difMulti(end) = 0;
%% Detecting single event
startEventIdx = find(difMulti ==1) + 1;
endEventIdx = find(difMulti == -1) + 1;

%%
h = histogram(endEventIdx - startEventIdx, 'Normalization','probability');
h.EdgeColor = 'k';
h.FaceColor = 'w';
h.LineWidth = 2;
ylim([0 1])
xlabel('Down-State (Single/Bursts) Events')
ylabel('Probability')
legend('Down-States')
legend('boxoff') 
% xlim([0 1])