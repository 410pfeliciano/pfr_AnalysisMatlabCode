% Cleaning the Workspace and Closing Figures
% close all
% clear
%% Filtering Positioning Data
filtpos = medfilt1(pos);
filledPos = fillmissing(filtpos,'nearest');
% filledPos = fillmissing(filledPos,'nearest');
figure
subplot(121)
plot(pos(:,1),pos(:,2))
legend('filledPos Raw')
subplot(122)
plot(filledPos(:,1),filledPos(:,2))
legend('Filling/Filt Pos')
% subplot(133)
% plot(pos(:,1),pos(:,2),'.')
% legend('Filling-Filt Pos')
%% Obtaining the Timing
close all 
Video_SR = 1/15; % PointGrey/Bonsai Video Sampling rate 15Hz
t = round((0 : Video_SR : length(filledPos) * Video_SR - Video_SR)',4); % Calculating time vector


%% Linearization,Interpolation and Filtering of Position Data
close all
linearPos = filledPos(:,1)+filledPos(:,2);
% Interpolation Filling NaN values in the linear position vector
% int_linear_centered_linear_pos = fillmissing(linearPos,'linear'); 
subplot(3,1,1)
plot(t, linearPos)
title('Estimated Linear Position')
xlabel('Time [sec]')
ylabel('Estimated Linear Position [pixels]')
axis tight
% Filtering and Plotting Linear Position
subplot(3,1,2)
plot(t,linearPos)
hold on
plot(t, linearPos)
linearPosFilt = movmean(linearPos,25);
plot(t, linearPosFilt, 'b', 'LineWidth', 2)
axis tight
title('Estimated Linear Position vs Filtered Estimated Position')
xlabel('Time [sec]')
ylabel('Estimated Linear Position [pixels]')
hold off
subplot(3,1,3)
plot(t, linearPosFilt, 'r')
title('Filtered Estimated Linear Position')
xlabel('Time [sec]')
ylabel('Estimated Linear Position [pixels]')
axis tight

%% Estimating Velocity
pospixels(:,1) = filledPos(:,1);
pospixels(:,2) = filledPos(:,2);
linvel = zeros(length(t),1);
linvel(2:end,:) = movmean(diff(linearPosFilt)./diff(t),25);
linvel(1) = linvel(2);

%% Putting the linear and Original Position in a Structure
positiondata.linearpos = linearPosFilt;
positiondata.time = t;
positiondata.pos = pospixels;
positiondata.linVel = linvel;

% save('PositionData.mat', 'positiondata');

%% Ploting the Results
figure
subplot(311)
plot(positiondata.pos(:,1),positiondata.pos(:,2),'.k')
xlabel('x coordinates(pixels)')
ylabel('y coordinates(pixels)')
axis tight
subplot(312)
plot(positiondata.time, positiondata.linearpos, 'k')
xlabel('Time(sec)')
ylabel('Linear Position(pixels)')
axis tight
subplot(313)
plot(positiondata.time, positiondata.linVel,'k')
xlabel('Time(sec)')
ylabel('Velocity(pixels/sec)')
axis tight
%%
save('PositionData.mat', 'positiondata');
%% Just in case you need to add more data to the file
save('pos.mat','positiondata','-append'); % use append to add variables to existing mat file
