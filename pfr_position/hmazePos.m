% Cleaning the Workspace and Closing Figures
% close all
% clear
% Filling missing values

position = fillmissing(position,'movmedian',100);
%% Filtering Positioning Data
filtpos = medfilt1(position, 11);
plot(position(:,1),position(:,2))
hold on
plot(filtpos(:,1),filtpos(:,2))
hold off
%% Obtaining the Timing
close all 
Video_SR = 1/15; % PointGrey/Bonsai Video Sampling rate 15Hz
t = round((0 : Video_SR : length(position) * Video_SR - Video_SR)',4); % Calculating time vector

%% Choosing best candidates
close all
v = filtpos';
% choose a point which will be the center of rotation
x_center = (max(v(1,:))+min(v(1,:)))/2;
y_center = (max(v(2,:))+min(v(2,:)))/2;
% create a matrix which will be used later in calculations
center = repmat([x_center; y_center], 1, length(v));
% define a 60 degree counter-clockwise rotation matrix
deg = 268.0;
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
centered_pos = zeros(size(position));
centered_pos(:,1) = (rotposition(:,1) - x_min); 
centered_pos(:,2) =  (rotposition(:,2) - y_min);
figure
plot(centered_pos(:,1), centered_pos(:,2),'.')
axis tight
%% X Distribution 
close all
% y_min = min(centered_pos(:,2)); 
% y_max = max(centered_pos(:,2));
% x_min = min(centered_pos(:,1));
% x_max = max(centered_pos(:,1));
% x_nbins = 0:2:x_max;
% x_edges = find(x_nbins>0 & x_nbins<300); % Determining the middle of the maze
% y_nbins = x_min:2:y_max;
% x_dist = hist(centered_pos(:,1), x_nbins);
% [Val,x_Ind] = max(x_dist(x_edges));
% y_dist = hist(centered_pos(:,2), y_nbins);
% figure
% bar(x_nbins(x_edges),x_dist(x_edges))
% title('Middle of the Maze X coordinate Position Distribution')
% figure
% bar(y_nbins, y_dist)
% title('Y coordinate Position Distribution')
figure
plot(centered_pos(:,1), centered_pos(:,2))
% Rearrenge x coordinates
centered_pos2 = zeros(size(position));
centered_pos2(:,1) = (centered_pos(:,1) - 232.8); ...
    % subtracting the middle of the maze by choosing the middle of the maze
centered_pos2(:,2) = (centered_pos(:,2));
hold on
plot(centered_pos2(:,1), centered_pos2(:,2))
title('Before and After X coordinate Rearrange')
axis tight
hold off

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
    if join_centered_pos(ii,2) < 35
    join_centered_pos2(ii,1) = (join_centered_pos(ii,1) * -1)+118;%Choose
    ...the largest value of x from the first step transformation
    elseif join_centered_pos(ii,2) >= 35
        join_centered_pos2(ii,1) = join_centered_pos(ii,1)+118;
    end
end
figure
plot(join_centered_pos2(:,1), join_centered_pos2(:,2), 'LineWidth',1)
title('2nd Maze Transformation')
xlabel('Estimated Position [cm]')
ylabel('Estimated Position [cm]')
axis([0 inf 0 inf])

%% Linearization,Interpolation and Filtering of Position Data
close all
linear_join_centered_pos = (join_centered_pos2(:,1) + join_centered_pos2(:,2));
% Interpolation Filling NaN values in the linear position vector
int_linear_centered_linear_pos = fillmissing(linear_join_centered_pos,'linear'); 
subplot(3,1,1)
plot(t, int_linear_centered_linear_pos)
title('Estimated Linear Position')
xlabel('Time [sec]')
ylabel('Estimated Linear Position [cm]')
axis tight
% Filtering and Plotting Linear Position
subplot(3,1,2)
plot(t,int_linear_centered_linear_pos)
hold on
plot(t, linear_join_centered_pos)
filt_centered_linear_pos = movmean(int_linear_centered_linear_pos,50);
plot(t, filt_centered_linear_pos, 'b', 'LineWidth', 2)
axis tight
title('Estimated Linear Position vs Filtered Estimated Position')
xlabel('Time [sec]')
ylabel('Estimated Linear Position [cm]')
hold off
subplot(3,1,3)
plot(t, filt_centered_linear_pos, 'r')
title('Filtered Estimated Linear Position')
xlabel('Time [sec]')
ylabel('Estimated Linear Position [cm]')
axis tight

%% Animated Graph
% close all
numpoints = length(t); 
x = t; 
y = filt_centered_linear_pos; 
x1 = centered_pos(:,1);
y1 = centered_pos(:,2);
x2 = join_centered_pos(:,1);
y2 = join_centered_pos(:,2);
x3 = join_centered_pos2(:,1);
y3 = join_centered_pos2(:,2);

subplot(2,2,2)  
h = animatedline; 
hold on
p = plot(x(1), y(1),'o','MarkerFaceColor','red','MarkerSize',8);
title('Estimated Linear Position')
xlabel('Time [sec]')
ylabel('Position [cm]')
hold off
axis([0 x(end) 0 420])

subplot(2,2,1)
h1 = animatedline;
axis([0 250 0 200])
hold on
p1 = plot(x1(1), y1(1),'o','MarkerFaceColor','red','MarkerSize',8);
title('H-Maze Position Data')
xlabel('X Coordinate Position [cm]')
ylabel('Y Coordinate Position [cm]')
hold off

subplot(2,2,3)
h2 = animatedline;
hold on
p2 = plot(x2(1), y2(1),'o','MarkerFaceColor','red','MarkerSize',8);
title('Transformed Position Data')
xlabel('X Coordinate Position [cm]')
ylabel('Y Coordinate Position [cm]')
axis([0 150 0 200])
hold off

subplot(2,2,4)
h3 = animatedline;
hold on
p3 = plot(x3(1), y3(1),'o','MarkerFaceColor','red','MarkerSize',8);
title('Transformed Position Data')
xlabel('X Coordinate Position [cm]')
ylabel('Y Coordinate Position [cm]')
axis([0 250 0 200])
hold off
  
 for k = 1:numpoints 
    addpoints(h,x(k),y(k))
    p.XData = x(k);
    p.YData = y(k);
    addpoints(h1,x1(k), y1(k))
    p1.XData = x1(k);
    p1.YData = y1(k);
    addpoints(h2,x2(k), y2(k))
    p2.XData = x2(k);
    p2.YData = y2(k);
    addpoints(h3,x3(k), y3(k))
    p3.XData = x3(k);
    p3.YData = y3(k);
    drawnow 
    pause(0.01);
 end 
%% Estimating Velocity and converting to meters
maxMazeLengthMeters = (0.6858) + 0.508;
linPosMax = max(filt_centered_linear_pos);
linPosM = filt_centered_linear_pos .* (maxMazeLengthMeters/linPosMax);
pos = zeros(size(centered_pos));
pos(:,1) = centered_pos(:,1) .* (0.6858/max(centered_pos(:,1)));
pos(:,2) = centered_pos(:,2) .* (0.508/max(centered_pos(:,2)));
linvel = zeros(length(t),1);
linvel(2:end,:) = movmean(diff(linPosM)./diff(t),25);
linvel(1) = linvel(2);

%% Putting the linear and Original Position in a Structure
positiondata.linearpos = linPosM;
positiondata.time = t;
positiondata.position = pos;
positiondata.linVel = linvel;

save('PositionData.mat', 'positiondata');

% Ploting the Results
figure
subplot(311)
plot(positiondata.position(:,1),positiondata.position(:,2),'k')
xlabel('x coordinates(m)')
ylabel('y coordinates(m)')
axis tight
subplot(312)
plot(positiondata.time, positiondata.linearpos, 'k')
xlabel('Time(sec)')
ylabel('Linear Position(m)')
axis tight
subplot(313)
plot(positiondata.time, positiondata.linVel,'k')
xlabel('Time(sec)')
ylabel('Velocity(m/sec)')
axis tight