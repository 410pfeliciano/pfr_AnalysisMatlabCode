%{ 
position is a matrix with x and y coordinates in cm
intFactor is interpolation factor to increase the  sampling rate
VideoSR is the video sampling rate
rotAngle is the rotating angle to straightening the position matrix(95.91)
%}

function [linpos, linVel, lintime] = linearPosVel(position_cm, intFactor, VideoSR, rotAngle)
% VideoSR = 1/15; % PointGrey/Bonsai Video Sampling rate 15Hz
t = round((0 : VideoSR : length(position_cm) * VideoSR - VideoSR)',4); % Calculating time vector
figure
plot(position_cm(:,1), position_cm(:,2)) % Ploting position data
axis tight
%%
% rotAngle = 95.91; % rotating angle counterclockwise
Rot_Matrix = [cosd(rotAngle) -sind(rotAngle); sin(rotAngle) cos(rotAngle)];
rot_centered_pos = zeros(size(position_cm));
    for ii = 1:length(position_cm)
    rot_centered_pos(ii,1) = (position_cm(ii,1)  * Rot_Matrix(1,1)) - (position_cm(ii,2) * Rot_Matrix(1,2));
    rot_centered_pos(ii,2) = (position_cm(ii,1) * Rot_Matrix(2,1)) - (position_cm(ii,2) * Rot_Matrix(2,2));
    end
figure
plot(rot_centered_pos(:,1), rot_centered_pos(:,2)) 
axis tight
%% Filling missing values from tracking
filled_rot_centered_pos = (fillmissing(rot_centered_pos,'linear'))/100;
figure
plot(filled_rot_centered_pos(:,1), filled_rot_centered_pos(:,2)) 
axis tight
%%
y_min = min(filled_rot_centered_pos(:,2));
x_min = min(filled_rot_centered_pos(:,1));
centered_pos = zeros(size(position_cm));
centered_pos(:,1) = (filled_rot_centered_pos(:,1) - x_min); 
centered_pos(:,2) =  (filled_rot_centered_pos(:,2) - y_min);
figure
plot(centered_pos(:,1), centered_pos(:,2))
axis tight
%%
y_min = min(centered_pos(:,2)); % important min and max values of x and y are updated
y_max = max(centered_pos(:,2));
x_min = min(centered_pos(:,1));
x_max = max(centered_pos(:,1));
x_nbins = 0:0.005:x_max;
x_edges = find(x_nbins > 0.3 & x_nbins< 0.45);
y_nbins = x_min:0.005:y_max;
x_dist = hist(centered_pos(:,1), x_nbins);
[Val,x_Ind] = max(x_dist(x_edges));
y_dist = hist(centered_pos(:,2), y_nbins);

figure
bar(x_nbins(x_edges),x_dist(x_edges))
title('Middle of the Maze X coordinate Position Distribution')

figure
bar(y_nbins, y_dist)
title('Y coordinate Position Distribution')
figure
plot(centered_pos(:,1), centered_pos(:,2))
title('Before X coordinate Rearrange')
axis tight
% Rearrenge x coordinates
centered_pos2 = zeros(size(position_cm));
centered_pos2(:,1) = (centered_pos(:,1) - x_nbins(x_edges(x_Ind))); ...
    % subtracting the middle of the maze by their positional distribution
centered_pos2(:,2) = (centered_pos(:,2) - y_min);

figure
plot(centered_pos2(:,1), centered_pos2(:,2))
title('After X coordinate Rearrange')
axis tight
%% Tranforming the H-Maze, adding the edges of the maze
% First Step
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
xlabel('Estimated Position [m]')
ylabel('Estimated Position [m]')
axis tight
% 2nd Step
join_centered_pos2 = zeros(size(centered_pos2)); 
join_centered_pos2(:,2) = join_centered_pos(:,2);
for ii = 1:length(centered_pos2)
    if join_centered_pos(ii,2) < 0.35
    join_centered_pos2(ii,1) = (join_centered_pos(ii,1) * - 1)+0.355;
    elseif join_centered_pos(ii,2) >= 0.35
        join_centered_pos2(ii,1) = join_centered_pos(ii,1)+0.355;
    end
end
figure
plot(join_centered_pos2(:,1), join_centered_pos2(:,2), 'LineWidth',1)
title('2nd Maze Transformation')
xlabel('Estimated Position [m]')
ylabel('Estimated Position [m]')
axis tight
%% Linearization,Interpolation and Filtering of Position Data
close all
linear_join_centered_pos = (join_centered_pos2(:,1) + join_centered_pos2(:,2));
subplot(3,1,1)
plot(t, linear_join_centered_pos)
title('Linear Position')
xlabel('Time [sec]')
ylabel('Linear Pos. [m]')
axis tight
% Filtering and Plotting Linear Position
subplot(3,1,2)
plot(t,linear_join_centered_pos)
hold on
plot(t, linear_join_centered_pos)
filt_centered_linear_pos = movmean(linear_join_centered_pos,25);
plot(t, filt_centered_linear_pos, 'b', 'LineWidth', 2)
axis tight
title('Linear Position vs Filtered Linear Position')
xlabel('Time [sec]')
ylabel('Linear Pos. [m]')
hold off
subplot(3,1,3)
plot(t, filt_centered_linear_pos, 'r')
title('Filtered Linear Position')
xlabel('Time [sec]')
ylabel('Linear Pos. [m]')
axis tight
%% Estimating velocity from linear position
close all
diff_time = diff(t); % First derivative of time
diff_pos = diff(filt_centered_linear_pos); % First derivative of position
velocity = abs(diff_pos./diff_time);
smooth_vel = movmean(velocity,20); % Filtering the velocity
figure
subplot(4,1,1)
plot(t(1:length(t)-1), velocity, 'k')
xlabel('Time [sec]')
ylabel('Speed [m * sec^{-1}]')
legend('Raw Speed')
axis tight
subplot(4,1,2)
hold on
plot(t(1:length(t)-1), velocity, 'k')
plot(t(1:length(t)-1), smooth_vel,'r', 'LineWidth',2)
legend('Raw Speed','Filtered Velocity')
xlabel('Time [sec]')
ylabel('Speed [m * sec^{-1}]')
axis tight
hold off
subplot(4,1,3)
plot(t(1:length(t)-1), smooth_vel,'r', 'LineWidth',1)
xlabel('Time [sec]')
ylabel('Speed [m * sec^{-1}]')
legend('Filtered Velocity')
axis tight
subplot(4,1,4)
plot(t, filt_centered_linear_pos,'k', 'LineWidth',1)
xlabel('Time [sec]')
ylabel('Position [m]')
legend('Position')
axis tight
%% Bin Spiking Data to Visualize Spiking in time and Position
close all
% Increase the sample rate of position and the input signal by a factor of r.
% intFactor = 30;
linpos = interp(filt_centered_linear_pos,intFactor);
new_SR = VideoSR / intFactor;
lintime = (0 : new_SR : length(linpos) * new_SR - new_SR)';
new_vel = zeros(length(smooth_vel)+1,1);
new_vel(1:length(smooth_vel))= smooth_vel;
new_vel(end) = smooth_vel(end);
linVel = interp(new_vel,intFactor);
figure
subplot(2,1,1)
plot(lintime, linpos, 'r', 'LineWidth', 2)
hold on
plot(t,filt_centered_linear_pos, 'k')
legend('Raw Position','Interpolated Position')
xlabel('Time [sec]')
ylabel('Position [m]')
axis tight
hold off
subplot(2,1,2)
plot(lintime, linVel, 'r')
xlabel('Time [sec]')
ylabel('Speed [m * sec^{-1}]')
legend('Interpolated Velocity')
axis tight
end
