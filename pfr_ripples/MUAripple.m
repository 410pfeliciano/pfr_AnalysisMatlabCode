NewSr = 1/2000;
linPos = interp1(new_t, lpos, LFPt);
linVel = interp1(new_t, new_vel, LFPt);
linTime = 0:NewSr: length(linPos)* NewSr - NewSr;

%%
figure
plot(linTime, linPos)
% hold on
% plot(new_t,linear_pos)
% hold off
figure
plot(linTime, linVel)
% hold on
% plot(new_t, new_vel)
%%
bins = (0 : NewSr : length(linPos) * NewSr - NewSr); % bin width
Timebins = 0:0.02: LFPt(end);
% Spike times
semiMUAtimes{1} = semiMUATet7.unit1.ts;
semiMUAtimes{2} = semiMUATet7.unit2.ts;
semiMUAtimes{3} = semiMUATet7.unit3.ts;
semiMUAtimes{4} = semiMUATet7.unit4.ts;
semiMUAtimes{5} = semiMUATet8.unit1.ts;
semiMUAtimes{6} = semiMUATet8.unit2.ts;
semiMUAtimes{7} = semiMUATet8.unit3.ts;
semiMUAtimes{8} = semiMUATet8.unit4.ts;
semiMUAtimes{9} = semiMUATet8.unit5.ts;
semiMUAtimes{10} = semiMUATet8.unit6.ts;
semiMUAtimes{11} = semiMUATet8.unit7.ts;
semiMUAtimes{12} = semiMUATet8.unit8.ts;
semiMUAtimes{13} = semiMUATet9.unit1.ts;
semiMUAtimes{14} = semiMUATet9.unit2.ts;
semiMUAtimes{15} = semiMUATet9.unit3.ts;
semiMUAtimes{16} = semiMUATet9.unit4.ts;
semiMUAtimes{17} = semiMUATet10.unit1.ts;
semiMUAtimes{18} = semiMUATet10.unit2.ts;
semiMUAtimes{19} = semiMUATet10.unit3.ts;
semiMUAtimes{20} = semiMUATet10.unit4.ts;
semiMUAtimes{21} = semiMUATet10.unit5.ts;
semiMUAtimes{22} = semiMUATet10.unit6.ts;
semiMUAtimes{23} = semiMUATet10.unit7.ts;
semiMUAtimes{24} = semiMUATet10.unit8.ts;
semiMUAtimes{25} = semiMUATet10.unit9.ts;
semiMUAtimes{26} = semiMUATet10.unit10.ts;
semiMUAtimes{27} = semiMUATet10.unit11.ts;
semiMUAtimes{28} = semiMUATet11.unit1.ts;
semiMUAtimes{29} = semiMUATet11.unit2.ts;
semiMUAtimes{30} = semiMUATet11.unit3.ts;
semiMUAtimes{31} = semiMUATet11.unit4.ts;
semiMUAtimes{32} = semiMUATet15.unit1.ts;
semiMUAtimes{33} = semiMUATet15.unit2.ts;
semiMUAtimes{34} = semiMUATet15.unit3.ts;
semiMUAtimes{35} = semiMUATet15.unit4.ts;

%% Binary Spiking Vector
MUAvect = zeros(length(LFPt), size(semiMUAtimes,2));
for kk = 1 : size(semiMUAtimes,2)
MUAvect(:,kk) = (hist(semiMUAtimes{kk},bins)); 
end

%%
MUAvectActive = zeros(size(MUAvect));
for kk = 1 : size(MUAvect,2)
    for qq = 1:size(MUAvect,1)
        if linVel(qq) > 0.03
            MUAvectActive(qq,kk) = MUAvect(qq,kk);
        else
            MUAvectActive(qq,kk) = 0;
        end
    end
end
%%
MUAidx = cell(1,size(MUAvect,2));
for kk = 1: size(MUAvect,2)
    MUAidx{1,kk} = find(MUAvect(:,kk));
end

MUAidx2 = cell(1,size(MUAvectActive,2));
for kk = 1: size(MUAvectActive,2)
    MUAidx2{1,kk} = find(MUAvectActive(:,kk));
end
%%
save('ForcedD3MUAvect.mat','MUAvect', 'MUAidx');
save('ForcedD3MUAvectActive_3cm_sec.mat','MUAvectActive', 'MUAidx2');
%% Position dependent Firing rates
Posbins = (0: 0.025 : 1.3)';
Occup =(hist(linPos,Posbins)* NewSr)';	%Convert occupancy to seconds.
MUAhist = zeros(size(Posbins,1), size(MUAvect,2));
for kk = 1: size(MUAidx,2)
    MUAv =(hist(linPos(MUAidx{1,kk}),Posbins))';%Hist positions @ spikes.
    MUAsps = MUAv./Occup; % MUA spike per second 
    MUAsps(isnan(MUAsps)) = 0;
    MUAhist(:,kk) = MUAsps;
end
%%
MUAhist2 = zeros(size(Posbins,1), size(MUAvect,2));
for kk = 1: size(MUAidx2,2)
    MUAv =(hist(linPos(MUAidx2{1,kk}),Posbins))';%Hist positions @ spikes.
    MUAsps = MUAv./Occup; % MUA spike per second 
    MUAsps(isnan(MUAsps)) = 0;
    MUAhist2(:,kk) = MUAsps;
end
%% Visual Inspection
% close all
% for kk = 1: size(MUAhist,2)
%    figure
%    bar(Posbins, MUAhist(:,kk),'b') 
%    hold on
%    bar(Posbins, MUAhist2(:,kk),'r')
%    hold off
% end
%%
Possigma = 0.05; % Standard deviation of the kernel 15ms
Posedges = -3*Possigma: Possigma : 3*Possigma; % Evaluate ranges form -3*sd to 3*sd dev
Poskernel = normpdf(Posedges,0,Possigma); % Evaluate the Gaussian kernel
Poskernel = Poskernel * Possigma;
MUAhistSm = zeros((size(MUAhist,1))+6, size(MUAhist2,2));
for kk = 1: size(MUAhist2,2)
   MUAhistSm(:,kk) = (conv(MUAhist2(:,kk), Poskernel))+.01; % Adding 0.01Hz
end
%%
close all
nSR = 1.3/(length(MUAhistSm(:,1))-1);
nPosbins = (0 : nSR : 1.3)';
C = {'k','b','r','g','y','m','c','k','b','r','g','y','m'}; % Cell array of colros.
for kk = 1: size(MUAhist2,2)
%    figure
%    bar(Posbins, MUAhist2(:,kk)) 
   subplot(size(MUAvect,2),1,kk)
%    plot(nPosbins, MUAhistSm(:,kk), 'LineWidth',2)
   area(nPosbins, MUAhistSm(:,kk),'FaceColor','flat')
   axis([0 1.3 0 inf])
%    hold off
end
xlabel('Position [m]')
ylabel('FR[Hz]')
%% Gaussian Kernel over time
close all
MUAvect2 = sum(MUAvect, 2);
sigma = 0.01; % Standard deviation of the kernel 15ms
edges = -3*sigma: NewSr : 3*sigma; % Evaluate ranges form -3*sd to 3*sd dev
kernel = normpdf(edges,0,sigma); % Evaluate the Gaussian kernel
MUAdensity = conv(MUAvect2, kernel); % Convolved MUA with the kernel
new_tt = (0 : NewSr : length(MUAdensity) * NewSr - NewSr);

% Gaussian kernel for individual cells
sigma2 = 0.1; % Standard deviation of the kernel 15ms
edges2 = -3*sigma2: NewSr : 3*sigma2; % Evaluate ranges form -3*sd to 3*sd dev
kernel2 = (normpdf(edges2,0,sigma2)); % Evaluate the Gaussian kernel
MUAseq = cell(1,size(MUAvect,2));
for kk = 1 : size(MUAvect,2)
    MUAseq{:,kk} = conv(MUAvect(:,1), kernel2);
end
MUAseq = cell2mat(MUAseq)';
seqtt = (0 : NewSr : length(MUAseq(1,:)) * NewSr - NewSr);
%% MUA Theshhold
% MUAhist = histogram(MUAdensity, 500);
% axis([0 200 0 3*10.^4])
% MUAmedian = mean(MUAdensity);
% MUAmad = std(MUAdensity);
% bar(MUAhist)
%% Visualization 
figure
subplot(511)
plot(linTime, linPos,'k', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('Position [m]')
axis([ 519 523 0 1.4])
subplot(512)
plot(new_tt, MUAdensity,'b', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('MUA [Hz]')
axis([ 519 523 0 inf])
subplot(513)
plot(linTime, linVel,'k', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('m/sec')
axis([ 519 523 -inf inf])
subplot(514)
% plot(LFPt, filtLFP,'k','LineWidth',2)
% hold on
% plot(LFPt, LFP, 'b')
% for kk = 1: length(below)
%     plot(LFPt(below(kk):abov(kk)),LFP(below(kk):abov(kk)), 'r', 'LineWidth',2) 
% end
% plot(LFPt, filtLFP*2,'k','LineWidth',1)
% hold off
% axis([519 523 -500 500])
% subplot(515)
plot(linTime(MUAidx{1}), MUAvect(MUAidx{1},1), '.K');
hold on
for kk = 2: size(MUAvect,2)
    plot(linTime(MUAidx{kk}), MUAvect(MUAidx{kk},kk)+kk-1, '.K');
end
hold off
xlabel('Time [sec]')
ylabel('Unit #')
axis([ 519 523 0 size(MUAvect,2)+3])
%%
figure
subplot(611)
plot(linTime, linPos,'k', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('Position [m]')
axis([76 90 0 1.4])
subplot(612)
plot(linTime(MUAidx2{1}), MUAvect(MUAidx2{1},1), '.k');
hold on
plot(linTime, linPos*10,'LineWidth',2,'Color', [0.7 0.7 0.7])
for kk = 2: size(MUAvect,2)
    plot(linTime(MUAidx2{kk}), MUAvect(MUAidx2{kk},kk)+kk-1, '.k');
end
hold off
xlabel('Time [sec]')
ylabel('Unit #')
axis([76 90 0 size(MUAvect,2)+3])

subplot(613)
plot(new_tt,MUAdensity,'b', 'LineWidth', 1)
xlabel('Time [sec]')
ylabel('MUA [Hz]')
axis([76 90 0 inf])

subplot(614)
plot(linTime, linVel,'k', 'LineWidth', 2)
xlabel('Time [sec]')
ylabel('m/sec')
axis([76 90 0 .4])

subplot(615)
hold on
plot(LFPt, LFP2, 'b')
for kk = 1: length(below)
    plot(LFPt(below(kk):abov(kk)),LFP2(below(kk):abov(kk)), 'r', 'LineWidth',2) 
end
plot(LFPt, (filtLFP-200)*4,'k','LineWidth',0.5)
xlabel('Time [sec]')
ylabel('LPF[uV]')
hold off
axis([76 90 -1500 500])

subplot(616)
hold on
plot(LFPt, LFP2, 'b')
plot(LFPt, (swrLFP2-500)*1.25,'k','LineWidth',0.5)
hold off
xlabel('Time [sec]')
ylabel('LPF[uV]')
axis([76 90 -1000 500])

% subplot(717)
% hold on
% plot(LFPt, LFP7, 'b')
% plot(LFPt, (swrLFP7-500)*1.25,'k','LineWidth',0.5)
% hold off
% xlabel('Time [sec]')
% ylabel('LPF[uV]')
% axis([76 90 -1000 700])