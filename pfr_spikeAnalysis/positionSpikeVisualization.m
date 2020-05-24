% Visual Inspection and increase the sampling of position data
videoSR = 1/15;
int_factor = 30;
linear_pos = interp(positiondata.linearpos,int_factor);
new_SR = videoSR / int_factor;
new_t = (0 : new_SR : length(linear_pos) * new_SR - new_SR)';
velocity = interp(positiondata.linVel,int_factor);
figure
subplot(2,1,1)
plot(new_t, linear_pos, 'r', 'LineWidth',1)
legend('Raw Position')
xlabel('Time [sec]')
ylabel('Position [cm]')
axis tight
subplot(2,1,2)
plot(new_t, velocity, 'r')
xlabel('Time [sec]')
ylabel('Speed [cm * sec^{-1}]')
legend('Interpolated Velocity')
axis tight
%% Load Spiking Information
close all
x_pos = positiondata.poscm(:,1);
new_x_pos = interp(x_pos,int_factor);
y_pos = positiondata.poscm(:,2);
new_y_pos = interp(y_pos,int_factor);
bins = (0 : new_SR : length(linear_pos) * new_SR - new_SR);
for kk = 44:54 % : size(spkClust,2)
spiketime = spkClust(kk).spkTime;
spiketrain = histcounts(spiketime,bins);
spiketrain2 = zeros(length(spiketrain),1);
    for qq = 1:length(spiketrain2)
        if velocity(qq) > 5
         
            spiketrain2(qq) = spiketrain(qq);
        else
        spiketrain2(qq) = 0;
        end
    end

%
spikeindex = find(spiketrain);
figure
subplot(3,1,1)
hold on
plot(positiondata.linearpos, positiondata.time, 'Color', [0.7,0.7,0.7])
plot(linear_pos(spikeindex),new_t(spikeindex),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 1)
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Single Unit Activity')
xlabel('Linear Position [cm]')
ylabel('Time [sec]')
axis([0 inf 1500 8000])

subplot(3,1,2)
spikeindex2 = find(spiketrain2);
hold on
plot(positiondata.linearpos, positiondata.time, 'Color', [0.7,0.7,0.7])
plot(linear_pos(spikeindex2),new_t(spikeindex2),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 2)
hold off
% legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Single Unit Activity During Movement')
xlabel('Linear Position [m]')
ylabel('Time [sec]')
axis([0 inf 1500 8000])

subplot(3,1,3)
plot(new_x_pos, new_y_pos,'k') 
axis tight
xlabel('x coordinate [cm]')
ylabel('y coordinate [cm]')
hold on;
%  Add dots where spikes occured. 
plot(new_x_pos(spikeindex2), new_y_pos(spikeindex2),'or', ...
    'MarkerSize', 1, 'MarkerFaceColor', 'r')
title('Position and Spikes')
hold off
end