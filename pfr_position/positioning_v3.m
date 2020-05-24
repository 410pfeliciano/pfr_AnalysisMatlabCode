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
% hold on
% plot(filledPos(:,1),filledPos(:,2),'.')
% hold off
legend('Filling/Filt Pos')
% subplot(133)
% plot(pos(:,1),pos(:,2),'.')
% legend('Filling-Filt Pos')
%% Obtaining the Timing
close all 
Video_SR = 1/15; % PointGrey/Bonsai Video Sampling rate 15Hz
t = round((0 : Video_SR : length(filledPos) * Video_SR - Video_SR)',4); % Calculating time vector

%% Choosing best candidates
close all
v = filledPos';
% choose a point which will be the center of rotation
x_center = (max(v(1,:))+min(v(1,:)))/2;
y_center = (max(v(2,:))+min(v(2,:)))/2;
% create a matrix which will be used later in calculations
center = repmat([x_center; y_center], 1, length(v));
% define a 60 degree counter-clockwise rotation matrix
deg = 265;
theta = deg * pi/180;       % 
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% do the rotation...
s = v - center;     % shift points in the plane so that the center of rotation is at the origin
so = R*s;           % apply the rotation about the origin
vo = so + center;   % shift again so the origin goes back to the desired center of rotation
% this can be done in one line as:
% vo = R*(v - center) + center
% pick out the vectors of rotated x- and y-data
x_rotated = vo(1,:);
y_rotated = vo(2,:);
% make a plot
plot(v(1,:), v(2,:), 'k.', x_rotated, y_rotated, 'r.');
axis equal
rotposition = [x_rotated;y_rotated]';
%% Centering the Maze
close all
y_min = min(rotposition(:,2));
x_min = min(rotposition(:,1));
centered_pos = zeros(size(filledPos));
centered_pos(:,1) = (rotposition(:,1) - x_min); 
centered_pos(:,2) =  (rotposition(:,2) - y_min);
figure
plot(centered_pos(:,1), centered_pos(:,2))
xmaxCentered = max(centered_pos(:,1));
ymaxCentered = max(centered_pos(:,2));
axis tight
%% OPTIONAL: Obtain center values 
% close all
% fig = figure;
% z = peaks;
% plot(centered_pos(:,1), centered_pos(:,2))
% [x_cent,y_cent] = ginput(1);
% close all
%% Rearrenge x coordinates
centered_pos2(:,1) = (centered_pos(:,1) - xmaxCentered/2); ...
    % subtracting the middle of the maze by their positional distribution
centered_pos2(:,2) = (centered_pos(:,2));
figure
plot(centered_pos2(:,1), centered_pos2(:,2))
title('After X coordinate Rearrange')
axis tight
%%
% Rearrenge x coordinates
% centered_pos2 = zeros(size(position));
% centered_pos2(:,1) = (centered_pos(:,1) - x_nbins(x_edges(x_Ind))); ...
%     % subtracting the middle of the maze by their positional distribution
% centered_pos2(:,2) = (centered_pos(:,2));
% figure
% plot(centered_pos2(:,1), centered_pos2(:,2))
% title('After X coordinate Rearrange')
% axis tight

%% Tranforming the H-Maze, adding the edges of the maze
% First Steps
close all
join_centered_pos = zeros(size(centered_pos2)); % 
join_centered_pos(:,2) = centered_pos2(:,2);
for ii = 1:length(centered_pos2)
    if centered_pos2(ii,1)  < 0
    join_centered_pos(ii,1) = centered_pos2(ii,1) * -1;
    elseif centered_pos2(ii,1) >= 0
        join_centered_pos(ii,1) = centered_pos2(ii,1);
    end
end
figure
plot(join_centered_pos(:,1), join_centered_pos(:,2), 'LineWidth',1)
title('Transformed Maze')
xlabel('Estimated Position [cm]')
ylabel('Estimated Position [cm]')
axis tight
% 2nd Steps
join_centered_pos2 = zeros(size(centered_pos2)); 
join_centered_pos2(:,2) = join_centered_pos(:,2);
for ii = 1:length(centered_pos2)
    if join_centered_pos(ii,2) < ymaxCentered/2
    join_centered_pos2(ii,1) = (join_centered_pos(ii,1) * -1)+ ymaxCentered/2;
    elseif join_centered_pos(ii,2) >= ymaxCentered/2
        join_centered_pos2(ii,1) = join_centered_pos(ii,1)+ ymaxCentered/2;
    end
end
figure
plot(join_centered_pos2(:,1), join_centered_pos2(:,2), 'LineWidth',1)
title('2nd Maze Transformation')
xlabel('Estimated Position [cm]')
ylabel('Estimated Position [cm]')
axis tight

%% Linearization,Interpolation and Filtering of Position Data
close all
linearPos = (join_centered_pos2(:,1) + join_centered_pos2(:,2));
% Interpolation Filling NaN values in the linear position vector
% int_linear_centered_linear_pos = fillmissing(linearPos,'linear'); 
subplot(3,1,1)
plot(t, linearPos)
title('Estimated Linear Position')
xlabel('Time [sec]')
ylabel('Estimated Linear Position [cm]')
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
ylabel('Estimated Linear Position [cm]')
hold off
subplot(3,1,3)
plot(t, linearPosFilt, 'r')
title('Filtered Estimated Linear Position')
xlabel('Time [sec]')
ylabel('Estimated Linear Position [cm]')
axis tight

%% Estimating Velocity and converting to meters
maxMazeLengthMeters = (0.6858) + 0.508;
maxMazeLengthCentimeters = maxMazeLengthMeters * 100;
linPosMax = max(linearPosFilt);
linPosM = linearPosFilt .* (maxMazeLengthCentimeters/linPosMax);
pos = zeros(size(centered_pos));
pos(:,1) = centered_pos(:,1) .* (0.6858*100/max(centered_pos(:,1)));
pos(:,2) = centered_pos(:,2) .* (0.508*100/max(centered_pos(:,2)));
linvel = zeros(length(t),1);
linvel(2:end,:) = movmean(diff(linPosM)./diff(t),25);
linvel(1) = linvel(2);

%% Putting the linear and Original Position in a Structure
positiondata.linearpos = linPosM;
positiondata.time = t;
positiondata.pos = pos;
positiondata.linVel = linvel;

% save('PositionData.mat', 'positiondata');

%% Ploting the Results
figure
subplot(311)
plot(positiondata.pos(:,1),positiondata.pos(:,2),'k.')
xlabel('x coordinates(cm)')
ylabel('y coordinates(cm)')
axis tight
subplot(312)
plot(positiondata.time, positiondata.linearpos, 'k')
xlabel('Time(sec)')
ylabel('Linear Position(cm)')
axis tight
subplot(313)
plot(positiondata.time, positiondata.linVel,'k')
xlabel('Time(sec)')
ylabel('Velocity(cm/sec)')
axis tight
%%
save('pos.mat','positiondata'); %'-append' use append to add variables to existing mat file
