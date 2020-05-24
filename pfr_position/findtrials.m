function [timeTrialsVal,timeTrialsIdx] = findtrials(positiondata, sites, combinationNumber)
%{
1-The combination number for the H-Maze is 8, For the T-Maze is 4
2-You need the position data obtained from the positioning_v3.m file
the file is stored in matlab drive
3-You need the sites data obtained from the MazeDivision.m function stored
in the matlab drive
4-First use the zoom to reach the data of interest
5-Press the space button from the keyboard to start selecting the trials
6-Once all the trials for a specific combination were choosen, use
the left click mouse to start with the following combination.
Possibles trial rutes/combinations in the H-Maze:
    1-Start Left to Reward Left(Choice Correct trial)
    2-Start Left to Reward right(Choice Incorrect Trial)
    3-Start Right to Reward Right(Choice Correct trial)
    4-Start Right to Reward Left (Choice Incorrect Trial)

    5-Reward Left to Start Left (Force)
    6-Reward Left to Start Right (Force)
    7-Reward Right to Start Right (Force)
    8-Reward Right to Start Left (Force)
%}
figure
subplot(211)
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

subplot(212)
plot(positiondata.time,positiondata.linearpos, 'k','LineWidth',1)
for kk = 1:size(sites.regions,2)
    hold on
    plot(positiondata.time(sites.regions{kk}),...
        positiondata.linearpos(sites.regions{kk}),'.','Color', col{kk}, ...
        'MarkerSize',5)
end
xlabel('Time(sec)')
ylabel('Position(cm)')
hold off
axis tight 
% subplot(313)
% plot(positiondata.time,positiondata.linVel, 'k','LineWidth',1)
% for kk = 1:size(sites.regions,2)
%     hold on
%     plot(positiondata.time(sites.regions{kk}),...
%         positiondata.linVel(sites.regions{kk}),'.','Color', col{kk}, ...
%         'MarkerSize',1)
% end
% xlabel('Time(sec)')
% ylabel('cm/sec')
% hold off
% axis tight
%%
zoom on; % Zoom in
pause() % When your image is okay, you press any key
zoom off; % Escape the zoom mode
trials{combinationNumber, 1} = [];
for kk = 1:combinationNumber
    button = 1;
    [x y] = deal([]);
    while button ==1  
        pan on; % to change the plot position
        pause() % you can change plot position with your mouse and ... 
        ... when your image is okay, you press any key. e.e space bar
        pan off; % to escape the pan mode
        [xt, yt, button] = ginput(2);
        x = [x , xt]; 
        y = [y yt];
    end
    xTrials{kk,1} = x;
    yTrials{kk,1} = y;
end
%% Calculate the closest value to time from ginput
timeTrialsIdx{combinationNumber,1} = [];
for ii = 1:combinationNumber
    for aa =1: size(xTrials{ii},2)
        for ll = 1:2
            d = sqrt((xTrials{ii}(ll,aa)-positiondata.time).^2 + (yTrials{ii}(ll,aa)-positiondata.linearpos).^2); 
            [~,minIdx] = min(d); 
            timeTrialsIdx{ii}(ll,aa) = minIdx;
        end
    end
end
%% Obtain the time points trials
timeTrialsVal{combinationNumber,1} = [];
for nn = 1:combinationNumber
   timeTrialsVal{nn,1}= positiondata.time(timeTrialsIdx{nn});
end
%% Eliminate Trials that are too short to be trials
diffTimeTrials{combinationNumber,1} = [];
for mm = 1:combinationNumber
   diffTimeTrials{mm} = timeTrialsVal{mm}(2,:) - timeTrialsVal{mm}(1,:); 
end
%%
badTrials{combinationNumber,1} = [];
for mm = 1:combinationNumber
   badTrials{mm,1} = find(diffTimeTrials{mm} < 2); 
end
for mm = 1:combinationNumber
    timeTrialsIdx{mm}(:,badTrials{mm}) = [];
    timeTrialsVal{mm}(:,badTrials{mm}) = [];
end
%%
save('MazeTrials.mat','timeTrialsIdx', 'timeTrialsVal')
end