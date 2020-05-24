% Ploting the Results
figure
subplot(311)
plot(positiondata.position(:,1),positiondata.position(:,2),'k')
hold on
plot(positiondata.position(sites.regions{1,1},1),... 
    positiondata.position(sites.regions{1,1},2),'y.')
plot(positiondata.position(sites.regions{1,2},1),... 
    positiondata.position(sites.regions{1,2},2),'r.')
plot(positiondata.position(sites.regions{1,3},1),... 
    positiondata.position(sites.regions{1,3},2),'b.')
plot(positiondata.position(sites.regions{1,4},1),... 
    positiondata.position(sites.regions{1,4},2),'c.')
plot(positiondata.position(sites.regions{1,5},1),... 
    positiondata.position(sites.regions{1,5},2),'g.')
hold off
xlabel('x coordinates(m)')
ylabel('y coordinates(m)')
axis tight
subplot(312)
hold on
plot(positiondata.time, positiondata.linearpos, 'k')
plot(positiondata.time(sites.regions{1,1},1), ... 
    positiondata.linearpos(sites.regions{1,1},1), 'y.')
plot(positiondata.time(sites.regions{1,2},1), ... 
    positiondata.linearpos(sites.regions{1,2},1), 'r.')
plot(positiondata.time(sites.regions{1,3},1), ... 
    positiondata.linearpos(sites.regions{1,3},1), 'b.')
plot(positiondata.time(sites.regions{1,4},1), ... 
    positiondata.linearpos(sites.regions{1,4},1), 'c.')
plot(positiondata.time(sites.regions{1,5},1), ... 
    positiondata.linearpos(sites.regions{1,5},1), 'g.')
xlabel('Time(sec)')
ylabel('Linear Position(m)')
axis tight
subplot(313)
plot(positiondata.time, positiondata.linVel,'k')
xlabel('Time(sec)')
ylabel('Velocity(m/sec)')
axis tight