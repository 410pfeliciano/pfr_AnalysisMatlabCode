% Cleaning the Workspace and Closing Figures
close all
clear
%% Identifying Outliers
plot(Drive_Pos(:,1),Drive_Pos(:,2))
%% Filtering Outliers From DeepLabCut
FiltDrivePos = medfilt1(Drive_Pos, 11);
figure
scatter(Drive_Pos(:,1),Drive_Pos(:,2),'.')
hold on
plot(FiltDrivePos(:,1),FiltDrivePos(:,2),'.')
hold off
FiltBodyPos = medfilt1(Body_Pos, 11);
figure
scatter(Body_Pos(:,1),Body_Pos(:,2),'.')
hold on
plot(FiltBodyPos(:,1),FiltBodyPos(:,2),'.')
hold off
FiltTailPos = medfilt1(Tail_Pos, 11);
figure
scatter(Tail_Pos(:,1),Tail_Pos(:,2),'.')
hold on
plot(FiltTailPos(:,1),FiltTailPos(:,2),'.')
hold off
figure
hold on
plot(FiltDrivePos(:,1),FiltDrivePos(:,2),'b.')
plot(FiltBodyPos(:,1),FiltBodyPos(:,2),'g.')
plot(FiltTailPos(:,1),FiltTailPos(:,2),'r.')
hold off
%% 
% Obtaining Positional Data
close all 
Video_SR = 1/15; % PointGrey/Bonsai Video Sampling rate 15Hz
t = round((0 : Video_SR : length(position_cm) * Video_SR - Video_SR)',4); % Calculating time vector
figure
plot(position_cm(:,1), position_cm(:,2)) % Ploting position data
axis tight

%% rotating the positional data by Creating a Rotating Matrix
...https://en.wikipedia.org/wiki/Rotation_matrix#In_two_dimensions
close all
for yy = 91:0.1:100
rot_angle = yy; % rotating angle counterclockwise
Rot_Matrix = [cosd(rot_angle) -sind(rot_angle); sin(rot_angle) cos(rot_angle)];
rot_centered_pos = zeros(size(position_cm));
    for ii = 1:length(position_cm)
    rot_centered_pos(ii,1) = (position_cm(ii,1) * Rot_Matrix(1,1)) - (position_cm(ii,2) * Rot_Matrix(1,2));
    rot_centered_pos(ii,2) = (position_cm(ii,1) * Rot_Matrix(2,1)) - (position_cm(ii,2) * Rot_Matrix(2,2));
    end
figure
plot(rot_centered_pos(:,1), rot_centered_pos(:,2)) 
axis tight
end
%% Choosing best candidates
close all
rot_angle = 95.91; % rotating angle counterclockwise
Rot_Matrix = [cosd(rot_angle) -sind(rot_angle); sin(rot_angle) cos(rot_angle)];
rot_centered_pos = zeros(size(position_cm));
    for ii = 1:length(position_cm)
    rot_centered_pos(ii,1) = (position_cm(ii,1)  * Rot_Matrix(1,1)) - (position_cm(ii,2) * Rot_Matrix(1,2));
    rot_centered_pos(ii,2) = (position_cm(ii,1) * Rot_Matrix(2,1)) - (position_cm(ii,2) * Rot_Matrix(2,2));
    end
figure
plot(rot_centered_pos(:,1), rot_centered_pos(:,2)) 
axis tight

%% Centering the Maze
close all
y_min = min(rot_centered_pos(:,2));
x_min = min(rot_centered_pos(:,1));
centered_pos = zeros(size(position_cm));
centered_pos(:,1) = (rot_centered_pos(:,1) - x_min); 
centered_pos(:,2) =  (rot_centered_pos(:,2) - y_min);
figure
plot(centered_pos(:,1), centered_pos(:,2))
axis tight
%% X Distribution 
close all
y_min = min(centered_pos(:,2)); % important min and max values of x and y are updated
y_max = max(centered_pos(:,2));
x_min = min(centered_pos(:,1));
x_max = max(centered_pos(:,1));
x_nbins = 0:0.5:x_max;
x_edges = find(x_nbins>25 & x_nbins<45);
y_nbins = x_min:0.5:y_max;
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
centered_pos2(:,2) = (centered_pos(:,2));
figure
plot(centered_pos2(:,1), centered_pos2(:,2))
title('After X coordinate Rearrange')
axis tight

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
    join_centered_pos2(ii,1) = (join_centered_pos(ii,1) * -1)+35.5;
    elseif join_centered_pos(ii,2) >= 35
        join_centered_pos2(ii,1) = join_centered_pos(ii,1)+35.5;
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
filt_centered_linear_pos = movmean(int_linear_centered_linear_pos,25);
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
p = plot(x(1), y(1),'o','MarkerFaceColor','red','MarkerSize',5);
title('Estimated Linear Position')
xlabel('Time [sec]')
ylabel('Position [cm]')
hold off
axis([0 1010 0 133])

subplot(2,2,1)
h1 = animatedline;
axis([0 70 0 70])
hold on
p1 = plot(x1(1), y1(1),'o','MarkerFaceColor','red','MarkerSize',12);
title('H-Maze Position Data')
xlabel('X Coordinate Position [cm]')
ylabel('Y Coordinate Position [cm]')
hold off

subplot(2,2,3)
h2 = animatedline;
hold on
p2 = plot(x2(1), y2(1),'o','MarkerFaceColor','red','MarkerSize',12);
title('Transformed Position Data')
xlabel('X Coordinate Position [cm]')
ylabel('Y Coordinate Position [cm]')
axis([0 35 0 70])
hold off

subplot(2,2,4)
h3 = animatedline;
hold on
p3 = plot(x3(1), y3(1),'o','MarkerFaceColor','red','MarkerSize',12);
title('Transformed Position Data')
xlabel('X Coordinate Position [cm]')
ylabel('Y Coordinate Position [cm]')
axis([0 70 0 70])
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

%% Bin Spiking Data to Visualize Spiking in time and Position
% close all
bins = (0 : Video_SR : length(position_cm) * Video_SR - Video_SR);
spiketime = tet2_wave_units.unit2.ts;
spiketrain = (hist(spiketime,bins))';
spikeindex = find(spiketrain);
figure
hold on
plot(filt_centered_linear_pos, t, 'Color', [0.7,0.7,0.7], 'LineWidth', 2)
plot(filt_centered_linear_pos(spikeindex),t(spikeindex),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 5)
hold off
legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Single Unit Activity')
xlabel('Linear Position [cm]')
ylabel('Time [sec]')
axis tight

figure
hold on
plot( t, filt_centered_linear_pos, 'Color', [0.7,0.7,0.7],'LineWidth',2)
plot(t(spikeindex),filt_centered_linear_pos(spikeindex),'o',...
    'MarkerEdgeColor','k','MarkerFaceColor',[1 0.2 0.2],'MarkerSize', 5)
hold off
legend('Linear Position','Single Unit Spikes','location','southoutside')
title('Single Unit Activity')
xlabel('Time [sec]')
ylabel('Linear Position [cm]')
axis tight


%% 
pos_bins = 0:10:120;


%% Separating trials
% Need to import 'right_trial' and 'left_trial'
right_trials = round((right_trial .* Video_SR),4); % converting frames to seconds
left_trials = round((left_trial .* Video_SR) ,4); % converting frames to seconds
right_ind = arrayfun(@(x)find(t==x,1),right_trials); % Find the indices of t
left_ind = arrayfun(@(x)find(t==x,1),left_trials);
right_choice_trials = cell(1,length(right_trials));
left_choice_trials = cell(1,length(left_trials));
for kk = 1:length(right_trials)
    right_choice_trials{kk} = t(right_ind(kk,1) : right_ind(kk,2));
end
for kk = 1:length(left_trials)
    left_choice_trials{kk} = (t(left_ind(kk,1)) : Video_SR : t(left_ind(kk,2)))';
end

%% Obtaining velocity from linear position
close all
diff_time = diff(t); % First derivative of time
diff_pos = diff(filt_centered_linear_pos); % First derivative of position
velocity = abs(diff_pos./diff_time);
smooth_vel = movmean(velocity,10);
figure
subplot(4,1,1)
plot(t(1:length(t)-1), velocity, 'k')
xlabel('Time [sec]')
ylabel('Speed [cm * sec^{-1}]')
legend('Raw Speed')
axis tight
subplot(4,1,2)
hold on
plot(t(1:length(t)-1), velocity, 'k')
plot(t(1:length(t)-1), smooth_vel,'r', 'LineWidth',2)
legend('Raw Speed','Filtered Velocity')
xlabel('Time [sec]')
ylabel('Speed [cm * sec^{-1}]')
axis tight
hold off
subplot(4,1,3)
plot(t(1:length(t)-1), smooth_vel,'r', 'LineWidth',1)
xlabel('Time [sec]')
ylabel('Speed [cm * sec^{-1}]')
legend('Filtered Velocity')
axis tight
subplot(4,1,4)
plot(t, filt_centered_linear_pos,'k', 'LineWidth',1)
xlabel('Time [sec]')
ylabel('Position [cm]')
legend('Position')
axis tight
%% Video play to separate trials

implay('vGAT_CreCh2_4_ForcedD32018-04-13T18_09_04.avi')

%%