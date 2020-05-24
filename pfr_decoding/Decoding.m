%% MUA activity not for decoding just for detecting regions of high activity
NewSr = 1/2000;
MUAvect2 = sum(MUAvect, 2);
sigma = 0.02; % Standard deviation of the kernel 15ms
edges = -3*sigma: NewSr : 3*sigma; % Evaluate ranges form -3*sd to 3*sd dev
kernel = normpdf(edges,0,sigma); % Evaluate the Gaussian kernel
MUAdensity = conv(MUAvect2, kernel); % Convolved MUA with the kernel
new_tt = (0 : NewSr : length(MUAdensity) * NewSr - NewSr);

%% Visual inspection of MUA activity
figure
subplot(411)
plot(linTime, linPos,'k', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('Position [m]')
axis([0 20 0 1.4])
subplot(412)
plot(linTime, linPos*12,'Color',[.7 .7 .7], 'LineWidth', 2)
hold on
plot(linTime(MUAidx{1}), MUAvect(MUAidx{1},1), '.K');
for kk = 2: size(MUAvect,2)
    plot(linTime(MUAidx{kk}), MUAvect(MUAidx{kk},kk)+kk-1, '.K');
end
hold off
xlabel('Time [sec]')
ylabel('Unit #')
axis([0 20 0 size(MUAvect,2)+3])
subplot(413)
plot(new_tt, MUAdensity,'b', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('MUA [Hz]')
axis([0 20 0 inf])
subplot(414)
plot(linTime, linVel,'k', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('m/sec')
axis([0 20 0 inf])
%% Finding MUA "reactivation"
threshold = 4 * mean(median(MUAdensity ./ .6745), 2);
MUAmax = MUAdensity >= threshold;
MUAmax = find(MUAmax == 1);
difMUAmax = find(diff(MUAmax) >= (0.02/NewSr)); 
z = zeros(length(difMUAmax)+1,1); 
z(1:end-1) = difMUAmax;
z(end) = length(MUAmax); 
st = difMUAmax+1; 
zst = zeros(length(difMUAmax)+1,1);
zst(2:end) = st;
zst(1)=1; 
MUAsteps = [MUAmax(zst), MUAmax(z)];  
duration = MUAsteps(:,2)- MUAsteps(:,1);
MUAsteps1 = MUAsteps((~(duration<=1)),:);
duration2 = (MUAsteps1(:,2) - MUAsteps1(:,1))*NewSr; %duration in sec.
MUAsidx = MUAsteps1((find(duration2<=0.550 & ...
    duration2 >= 0.01)),:); %Events < 150ms and > 15ms (Indices)
MUAsec = MUAsidx*NewSr;
below = MUAsidx(:,1)-(0.25/NewSr);
abov = MUAsidx(:,2)+(0.25/NewSr);
below(below < 0) = 1; % Substitute start negative value with 1
abov(abov > length(MUAdensity)) = length(MUAdensity); % Making MUA candidates finishing at the end of the LFP
belowidx = find(linVel(below) < 0.025);
below = below(belowidx);
abov = abov(belowidx);
% MUAregions = [below abov];
% MUAregionsSec = (MUAregions(:,2) .*NewSr) - (MUAregions(:,1) .*NewSr);
%% Visual Inspection of MUA "reactivation" detection
subplot(211)
plot(linTime,linPos)
axis([0 100 0 inf])
subplot(212)
plot(new_tt, MUAdensity)
hold on
for kk = 1: length(below)
    plot(new_tt(below(kk):abov(kk)),MUAdensity(below(kk):abov(kk)), 'r', 'LineWidth',2) 
end
% plot(MUAdensityt, filtMUAdensity,'k','LineWidth',1)
% plot(MUAdensityt, smsqfiltMUAdensity./80,'m', 'LineWidth',2)
hold off
axis([0 100 0 inf])
%% Creating a matrix with spikes, poistion, velocity and LFP based on MUA
MUAreact = cell(1,size(below,1));
for jj = 1:size(below,1)
        ind1 = below(jj,1);
        ind2 = abov(jj,1);
        MUAreact{1,jj}(:,1) = linTime(ind1: ind2); % Time
        MUAreact{1,jj}(:,2) = linPos(ind1: ind2); % Position
        MUAreact{1,jj}(:,3) = linVel(ind1: ind2); % Velocity
        MUAreact{1,jj}(:,4) = LFP15(ind1: ind2); %LFP
        MUAreact{1,jj}(:,5) = swrLFP15(ind1: ind2); %Theta filtered LFP
        MUAreact{1,jj}(:,6) = LFP2(ind1: ind2); %LFP
        MUAreact{1,jj}(:,7) = swrLFP2(ind1: ind2); %Theta filtered LFP
        for kk = 1: size(MUAvect,2) % Spikes(0 or 1)
            MUAreact{1,jj}(:,7+kk) = MUAvect(ind1:ind2,kk);
        end
end
%% Finding Spikes in the MUA
MUAreactIdx = cell (1, length(below));
for kk = 1: length(below)
    for jj = 1: size(MUAvect,2)
    MUAreactIdx{1, kk}{1,jj} = find(MUAreact{1,kk}(:, 7+jj));
    end
end
%% Visual Inspection
close all
for kk = 1: length(below)
    figure(kk)
    subplot(511)
    plot(MUAreact{1,kk}(:,1), MUAreact{1,kk}(:,2),'k', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('Position [m]')
    axis([-inf inf 0 1.3])
    subplot(512)
    plot(MUAreact{1,kk}(:,1), MUAreact{1,kk}(:,2)*12,'Color',[.7 .7 .7], 'LineWidth', 2)
    hold on
        for jj = 1: size(MUAvect,2)
            plot(MUAreact{1,kk}(MUAreactIdx{1,kk}{1,jj},1), ...
                MUAreact{1,kk}(MUAreactIdx{1,kk}{1,jj},jj+7)+jj-1, '.K');
        end
    hold off
    xlabel('Time [sec]')
    ylabel('Unit #')
    axis([-inf inf 0 size(MUAvect,2)+3])
    subplot(513)
    plot(new_tt(below(kk):abov(kk)), MUAdensity(below(kk):abov(kk)),'b', 'LineWidth', 2)
    xlabel('Time [sec]')
    ylabel('MUA [Hz]')
    axis([-inf inf 0 300])
    subplot(514)
    plot(linTime(below(kk):abov(kk)), LFP15(below(kk):abov(kk)),'b')
    hold on
    plot(linTime(below(kk):abov(kk)), swrLFP15(below(kk):abov(kk))-300,'k')
    hold off
    xlabel('Time [sec]')
    ylabel('uV')
    axis([-inf inf -inf inf])
    subplot(515)
    plot(linTime(below(kk):abov(kk)), LFP14(below(kk):abov(kk)),'b')
    hold on
    plot(linTime(below(kk):abov(kk)), swrLFP14(below(kk):abov(kk))-300,'k')
    hold off
    xlabel('Time [sec]')
    ylabel('uV')
    axis([-inf inf -inf inf])
end
%%
BinTime = 0.025;
Timebins = (0: BinTime:linTime(end))'; 
TimeHist = zeros(size(Timebins,1), size(MUAidx,2));
for kk = 1: size(TimeHist,2)
    TimeHist(:,kk) =(hist(linTime(MUAidx{1,kk}),Timebins))';
end

%%
SumFRExp = exp(-(BinTime).* sum(rrTuningCurveSm,2)); % Using the ChRight tuning curve
SumFRExp1 = exp(-(BinTime).* sum(llTuningCurveSm,2));
SumFRExp2 = exp(-(BinTime).* sum(frrrTuningCurveSm,2));
SumFRExp3 = exp(-(BinTime).* sum(frllTuningCurveSm,2));
SumFRExp4 = exp(-(BinTime).* sum(frrlTuningCurveSm,2));
SumFRExp5 = exp(-(BinTime).* sum(frlrTuningCurveSm,2));
%%
FRspikes = zeros(length(TimeHist(:,1)), size(rrTuningCurveSm,1));
for tt = 1: length(TimeHist(:,1))
    FRspikes(tt,:) = (prod(rrTuningCurveSm(:,:).^TimeHist(tt,:),2))';
end
FRspikes1 = zeros(length(TimeHist(:,1)), size(llTuningCurveSm,1));
for tt = 1: length(TimeHist(:,1))
    FRspikes1(tt,:) = (prod(llTuningCurveSm(:,:).^TimeHist(tt,:),2))';
end
FRspikes2 = zeros(length(TimeHist(:,1)), size(frrrTuningCurveSm,1));
for tt = 1: length(TimeHist(:,1))
    FRspikes2(tt,:) = (prod(frrrTuningCurveSm(:,:).^TimeHist(tt,:),2))';
end
FRspikes3 = zeros(length(TimeHist(:,1)), size(frllTuningCurveSm,1));
for tt = 1: length(TimeHist(:,1))
    FRspikes3(tt,:) = (prod(frllTuningCurveSm(:,:).^TimeHist(tt,:),2))';
end
FRspikes4 = zeros(length(TimeHist(:,1)), size(frrlTuningCurveSm,1));
for tt = 1: length(TimeHist(:,1))
    FRspikes4(tt,:) = (prod(frrlTuningCurveSm(:,:).^TimeHist(tt,:),2))';
end
FRspikes5 = zeros(length(TimeHist(:,1)), size(frlrTuningCurveSm,1));
for tt = 1: length(TimeHist(:,1))
    FRspikes5(tt,:) = (prod(frlrTuningCurveSm(:,:).^TimeHist(tt,:),2))';
end
%%
FRspikes = FRspikes';
FRspikes1 = FRspikes1';
FRspikes2 = FRspikes2';
FRspikes3 = FRspikes3';
FRspikes4 = FRspikes4';
FRspikes5 = FRspikes5';
%%
Prob = SumFRExp .* FRspikes;
SumProb = sum(Prob,1);
Likelihood = Prob./SumProb;
%
Prob1 = SumFRExp1 .* FRspikes1;
SumProb1 = sum(Prob1,1);
Likelihood1 = Prob1./SumProb1;
%
Prob2 = SumFRExp2 .* FRspikes2;
SumProb2 = sum(Prob2,1);
Likelihood2 = Prob2./SumProb2;
%
Prob3 = SumFRExp3 .* FRspikes3;
SumProb3 = sum(Prob3,1);
Likelihood3 = Prob3./SumProb3;
%
Prob4 = SumFRExp4 .* FRspikes4;
SumProb4 = sum(Prob4,1);
Likelihood4 = Prob4./SumProb4;
%
Prob5 = SumFRExp5 .* FRspikes5;
SumProb5 = sum(Prob5,1);
Likelihood5 = Prob5./SumProb5;
%% Summing Likelihoods for diferent routes
GenLikelihood = Likelihood + Likelihood1 + Likelihood2 + Likelihood3 + ...
    Likelihood4 + Likelihood5;
% GenLikelihood = Likelihood5+ Likelihood+ Likelihood1;
MaxGenLikeli = sum(GenLikelihood,1);
GenLikelihood2 = GenLikelihood./MaxGenLikeli;
%% Visual Inspection
% close all
figure
subplot(411)
plot(linTime, linPos,'k', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('Position [m]')
axis([19 21 0 1.4])
subplot(412)
plot(linTime(MUAidx{1}), MUAvect(MUAidx{1},1), '.K');
hold on
for kk = 2: size(MUAvect,2)
    plot(linTime(MUAidx{kk}), MUAvect(MUAidx{kk},kk)+kk-1, '.K');
end
hold off
xlabel('Time [sec]')
ylabel('Unit #')
axis([19 21 0 size(MUAvect,2)+3])
subplot(413)
plot(new_tt, MUAdensity,'b', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('MUA [Hz]')
axis([19 21 0 inf])
subplot(414)
hold on
s = plot3(linTime, linPos,ones(size(linPos,1),1),'--r','LineWidth', 1);
alpha(s,0.1)
s1 = surf(Timebins, tuningbins,GenLikelihood2);
s1.EdgeColor = 'none';
c = colorbar;
c.Label.String = 'Probability';
colormap(flipud(bone))
% caxis auto
axis([19 21 0 1.4])
caxis([0 .5])
xlabel('Time [sec]')
ylabel('Linear Position [m]')
zlabel('Probability') 
legend('Actual Position')
box on
hold off

% figure
% subplot(411)
% plot(linTime, linPos,'k', 'LineWidth', 2)
% xlabel('Time [sec]')
% ylabel('Position [m]')
% axis([19 21 0 1.4])
% subplot(412)
% plot(linTime(MUAidx{1}), MUAvect(MUAidx{1},1), '.K');
% hold on
% for kk = 2: size(MUAvect,2)
%     plot(linTime(MUAidx{kk}), MUAvect(MUAidx{kk},kk)+kk-1, '.K');
% end
% hold off
% xlabel('Time [sec]')
% ylabel('Unit #')
% axis([19 21 0 size(MUAvect,2)+3])
% subplot(413)
% plot(new_tt, MUAdensity,'b', 'LineWidth', 2)
% xlabel('Time [sec]')
% ylabel('MUA [Hz]')
% axis([19 21 0 inf])
% subplot(414)
% hold on
% s = plot3(linTime, linPos,ones(size(linPos,1),1),'--r','LineWidth', 1);
% alpha(s,0.1)
% s1 = surf(Timebins, tuningbins,Likelihood);
% s1.EdgeColor = 'none';
% c = colorbar;
% c.Label.String = 'Probability';
% colormap(flipud(bone))
% % caxis auto
% axis([19 21 0 1.4])
% caxis([0 .5])
% xlabel('Time [sec]')
% ylabel('Linear Position [m]')
% zlabel('Probability') 
% legend('Actual Position')
% box on
% hold off
% 
% figure
% subplot(411)
% plot(linTime, linPos,'k', 'LineWidth', 2)
% xlabel('Time [sec]')
% ylabel('Position [m]')
% axis([19 21 0 1.4])
% subplot(412)
% plot(linTime(MUAidx{1}), MUAvect(MUAidx{1},1), '.K');
% hold on
% for kk = 2: size(MUAvect,2)
%     plot(linTime(MUAidx{kk}), MUAvect(MUAidx{kk},kk)+kk-1, '.K');
% end
% hold off
% xlabel('Time [sec]')
% ylabel('Unit #')
% axis([19 21 0 size(MUAvect,2)+3])
% subplot(413)
% plot(new_tt, MUAdensity,'b', 'LineWidth', 2)
% xlabel('Time [sec]')
% ylabel('MUA [Hz]')
% axis([19 21 0 inf])
% subplot(414)
% hold on
% s = plot3(linTime, linPos,ones(size(linPos,1),1),'--r','LineWidth', 1);
% alpha(s,0.1)
% s1 = surf(Timebins, tuningbins,Likelihood1);
% s1.EdgeColor = 'none';
% c = colorbar;
% c.Label.String = 'Probability';
% colormap(flipud(bone))
% % caxis auto
% axis([19 21 0 1.4])
% caxis([0 .5])
% xlabel('Time [sec]')
% ylabel('Linear Position [m]')
% zlabel('Probability') 
% legend('Actual Position')
% box on
% hold off
% 
% figure
% subplot(411)
% plot(linTime, linPos,'k', 'LineWidth', 2)
% xlabel('Time [sec]')
% ylabel('Position [m]')
% axis([19 21 0 1.4])
% subplot(412)
% plot(linTime(MUAidx{1}), MUAvect(MUAidx{1},1), '.K');
% hold on
% for kk = 2: size(MUAvect,2)
%     plot(linTime(MUAidx{kk}), MUAvect(MUAidx{kk},kk)+kk-1, '.K');
% end
% hold off
% xlabel('Time [sec]')
% ylabel('Unit #')
% axis([19 21 0 size(MUAvect,2)+3])
% subplot(413)
% plot(new_tt, MUAdensity,'b', 'LineWidth', 2)
% xlabel('Time [sec]')
% ylabel('MUA [Hz]')
% axis([19 21 0 inf])
% subplot(414)
% hold on
% s = plot3(linTime, linPos,ones(size(linPos,1),1),'--r','LineWidth', 1);
% alpha(s,0.1)
% s1 = surf(Timebins, tuningbins,Likelihood2);
% s1.EdgeColor = 'none';
% c = colorbar;
% c.Label.String = 'Probability';
% colormap(flipud(bone))
% % caxis auto
% axis([19 21 0 1.4])
% caxis([0 .5])
% xlabel('Time [sec]')
% ylabel('Linear Position [m]')
% zlabel('Probability') 
% legend('Actual Position')
% box on
% hold off
% 
% figure
% subplot(411)
% plot(linTime, linPos,'k', 'LineWidth', 2)
% xlabel('Time [sec]')
% ylabel('Position [m]')
% axis([19 21 0 1.4])
% subplot(412)
% plot(linTime(MUAidx{1}), MUAvect(MUAidx{1},1), '.K');
% hold on
% for kk = 2: size(MUAvect,2)
%     plot(linTime(MUAidx{kk}), MUAvect(MUAidx{kk},kk)+kk-1, '.K');
% end
% hold off
% xlabel('Time [sec]')
% ylabel('Unit #')
% axis([19 21 0 size(MUAvect,2)+3])
% subplot(413)
% plot(new_tt, MUAdensity,'b', 'LineWidth', 2)
% xlabel('Time [sec]')
% ylabel('MUA [Hz]')
% axis([19 21 0 inf])
% subplot(414)
% hold on
% s = plot3(linTime, linPos,ones(size(linPos,1),1),'--r','LineWidth', 1);
% alpha(s,0.1)
% s1 = surf(Timebins, tuningbins,Likelihood3);
% s1.EdgeColor = 'none';
% c = colorbar;
% c.Label.String = 'Probability';
% colormap(flipud(bone))
% % caxis auto
% axis([19 21 0 1.4])
% caxis([0 .5])
% xlabel('Time [sec]')
% ylabel('Linear Position [m]')
% zlabel('Probability') 
% legend('Actual Position')
% box on
% hold off
% 
% figure
% subplot(411)
% plot(linTime, linPos,'k', 'LineWidth', 2)
% xlabel('Time [sec]')
% ylabel('Position [m]')
% axis([19 21 0 1.4])
% subplot(412)
% plot(linTime(MUAidx{1}), MUAvect(MUAidx{1},1), '.K');
% hold on
% for kk = 2: size(MUAvect,2)
%     plot(linTime(MUAidx{kk}), MUAvect(MUAidx{kk},kk)+kk-1, '.K');
% end
% hold off
% xlabel('Time [sec]')
% ylabel('Unit #')
% axis([19 21 0 size(MUAvect,2)+3])
% subplot(413)
% plot(new_tt, MUAdensity,'b', 'LineWidth', 2)
% xlabel('Time [sec]')
% ylabel('MUA [Hz]')
% axis([19 21 0 inf])
% subplot(414)
% hold on
% s = plot3(linTime, linPos,ones(size(linPos,1),1),'--r','LineWidth', 1);
% alpha(s,0.1)
% s1 = surf(Timebins, tuningbins,Likelihood4);
% s1.EdgeColor = 'none';
% c = colorbar;
% c.Label.String = 'Probability';
% colormap(flipud(bone))
% % caxis auto
% axis([19 21 0 1.4])
% caxis([0 .5])
% xlabel('Time [sec]')
% ylabel('Linear Position [m]')
% zlabel('Probability') 
% legend('Actual Position')
% box on
% hold off
% 
% figure
% subplot(411)
% plot(linTime, linPos,'k', 'LineWidth', 2)
% xlabel('Time [sec]')
% ylabel('Position [m]')
% axis([19 21 0 1.4])
% subplot(412)
% plot(linTime(MUAidx{1}), MUAvect(MUAidx{1},1), '.K');
% hold on
% for kk = 2: size(MUAvect,2)
%     plot(linTime(MUAidx{kk}), MUAvect(MUAidx{kk},kk)+kk-1, '.K');
% end
% hold off
% xlabel('Time [sec]')
% ylabel('Unit #')
% axis([19 21 0 size(MUAvect,2)+3])
% subplot(413)
% plot(new_tt, MUAdensity,'b', 'LineWidth', 2)
% xlabel('Time [sec]')
% ylabel('MUA [Hz]')
% axis([19 21 0 inf])
% subplot(414)
% hold on
% s = plot3(linTime, linPos,ones(size(linPos,1),1),'--r','LineWidth', 1);
% alpha(s,0.1)
% s1 = surf(Timebins, tuningbins,Likelihood5);
% s1.EdgeColor = 'none';
% c = colorbar;
% c.Label.String = 'Probability';
% colormap(flipud(bone))
% % caxis auto
% axis([19 21 0 1.4])
% caxis([0 .5])
% xlabel('Time [sec]')
% ylabel('Linear Position [m]')
% zlabel('Probability') 
% legend('Actual Position')
% box on
% hold off



















