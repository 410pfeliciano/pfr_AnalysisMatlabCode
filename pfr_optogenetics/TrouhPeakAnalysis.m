% 
close all
figure
subplot(121)
h1 = polarhistogram(StimThetaAnglePK,20,'Normalization','probability',...
    'FaceColor',[210/255 105/255 30/255],'FaceAlpha',.4,'EdgeColor','none');
hold on
h2 = polarhistogram(StimThetaAngleTR,20,'Normalization','probability',...
    'FaceColor',[98/255 159/255 67/255],'FaceAlpha',.4,'EdgeColor','none');
hold off
% title('Theta Phase Angle Distribution')
subplot(122)
hold on
h3 = histogram(rad2deg(StimThetaAnglePK),20,'Normalization','probability',...
    'FaceColor',[210/255 105/255 30/255],'FaceAlpha',.4,'EdgeColor','none');
h4 = histogram(rad2deg(StimThetaAngleTR),20,'Normalization','probability',...
    'FaceColor',[98/255 159/255 67/255],'FaceAlpha',.4,'EdgeColor','none');
x = rad2deg(-pi:0.01:pi);
amplred = 0.3;
devzero = 0.02;
y =(amplred*cos(2*pi*.00135*x)).*(amplred*cos(2*pi*.00135*x))+devzero;
y(end) = NaN;
c = x;
patch(x,y,c,'LineWidth', 2,'EdgeColor','k', 'EdgeAlpha',0.3);
hold off
% title('Theta Phase Angle Distribution')
xlabel('Degrees')
ylabel('Probability')
axis([-181 180 0 inf])
legend('Peak', 'Trough','NumColumns',2, 'Location','northoutside')
legend('boxoff')
degreetick 'x'
sgtitle('Theta-Phase Stimulation Angle Distribution')
%% 
figure
subplot(211)
h1 = polarhistogram(StimThetaAnglePK,20,'Normalization','probability',...
    'FaceColor',[210/255 105/255 30/255],'FaceAlpha',.4,'EdgeColor','none');
hold on
h2 = polarhistogram(StimThetaAngleTR,20,'Normalization','probability',...
    'FaceColor',[98/255 159/255 67/255],'FaceAlpha',.4,'EdgeColor','none');
hold off
% title('Theta Phase Angle Distribution')
subplot(212)
hold on
h3 = histogram(rad2deg(StimThetaAnglePK),20,'Normalization','probability',...
    'FaceColor',[210/255 105/255 30/255],'FaceAlpha',.4,'EdgeColor','none');
h4 = histogram(rad2deg(StimThetaAngleTR),20,'Normalization','probability',...
    'FaceColor',[98/255 159/255 67/255],'FaceAlpha',.4,'EdgeColor','none');
x = rad2deg(-pi:0.01:pi);
amplred = 0.3;
devzero = 0.02;
y =(amplred*cos(2*pi*.00135*x)).*(amplred*cos(2*pi*.00135*x))+devzero;
y(end) = NaN;
c = x;
patch(x,y,c,'LineWidth', 3,'EdgeColor','k', 'EdgeAlpha',0.3);
hold off
% title('Theta Phase Angle Distribution')
xlabel('Degrees')
ylabel('Probability')
axis([-181 180 0 inf])
legend('Peak', 'Trough','NumColumns',2, 'Location','northoutside')
legend('boxoff')
degreetick 'x'
sgtitle('Theta-Phase Stimulation Angle Distribution')