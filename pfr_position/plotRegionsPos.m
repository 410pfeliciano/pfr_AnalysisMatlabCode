figure
subplot(311)
% plot(positiondata.poscm(:,1), positiondata.poscm(:,2),'ok')
hold on
col= {[1 0 0] [1 0.8 0] [0 0 1] [0 1 1] [0 0 0]};
for kk = 1:size(sites.regions,2)
plot(positiondata.poscm(sites.regions{kk},1),...
    positiondata.poscm(sites.regions{kk},2),'.','Color', col{kk})
end
hold off
xlabel(' x coordinates(cm)')
ylabel('y coordinates(cm)')
axis tight

subplot(312)
plot(positiondata.time,positiondata.linearpos, 'k','LineWidth',1)
for kk = 1:size(sites.regions,2)
    hold on
    plot(positiondata.time(sites.regions{kk}),...
        positiondata.linearpos(sites.regions{kk}),'.','Color', col{kk}, ...
        'MarkerSize',1)
end
xlabel('Time(sec)')
ylabel('Position(cm)')
hold off
axis tight 
subplot(313)
plot(positiondata.time,positiondata.linVel, 'k','LineWidth',1)
for kk = 1:size(sites.regions,2)
    hold on
    plot(positiondata.time(sites.regions{kk}),...
        positiondata.linVel(sites.regions{kk}),'.','Color', col{kk}, ...
        'MarkerSize',1)
end
xlabel('Time(sec)')
ylabel('cm/sec')
hold off
axis tight


