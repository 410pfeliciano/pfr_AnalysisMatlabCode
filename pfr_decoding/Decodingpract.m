NewSr = 1/2000;
linPos = interp1(new_t, lpos, LFPt);
linVel = interp1(new_t, new_vel, LFPt);
linTime = 0:NewSr: length(linPos)* NewSr - NewSr;

%%
figure
plot(linTime, linPos)
hold on
plot(new_t,linear_pos)
hold off
figure
plot(linTime, linVel)
hold on
plot(new_t, new_vel)
%%
bins = (0 : NewSr : length(linPos) * NewSr - NewSr); % bin width
Timebins = 0:0.02: LFPt(end);
spiketime1 = sorting_tet2_v2.unit1.ts; % Obtaining the spikes
spiketime2 = sorting_tet2_v2.unit2.ts; % Obtaining the spikes
spiketime3 = sorting_tet2_v2.unit3.ts; % Obtaining the spikes
spiketime4 = sorting_tet2_v2.unit4.ts; % Obtaining the spikes
spiketime5 = sorting_tet8_v2.unit1.ts; % Obtaining the spikes
spiketime6 = sorting_tet8_v2.unit2.ts; % Obtaining the spikes
spiketime7 = sorting_tet8_v2.unit4.ts; % Obtaining the spikes
spiketime8 = sorting_tet8_v2.unit5.ts; % Obtaining the spikes
spiketime9 = sorting_tet9_w2.unit1.ts; % Obtaining the spikes
spiketime10 = sorting_tet9_w2.unit2.ts; % Obtaining the spikes
spiketime11 = sorting_tet9_w2.unit3.ts; % Obtaining the spikes
spiketime12 = sorting_tet10_w1.unit1.ts; % Obtaining the spikes
spiketime13 = sorting_tet10_w1.unit2.ts; % Obtaining the spikes
spiketime14 = sorting_tet10_w1.unit3.ts; % Obtaining the spikes
spiketime15 = sorting_tet10_w1.unit4.ts; % Obtaining the spikes
spiketime16 = sorting_tet11_w1.unit1.ts; % Obtaining the spikes
spiketime17 = sorting_tet11_w1.unit2.ts; % Obtaining the spikes
spiketime18 = sorting_tet11_w1.unit3.ts; % Obtaining the spikes
spiketime19 = sorting_tet15_v2.unit1.ts; % Obtaining the spikes
spiketime20 = sorting_tet15_v2.unit2.ts; % Obtaining the spikes
spiketime21 = sorting_tet15_v2.unit3.ts; % Obtaining the spikes
% spiketime22 = sorting_tet15_v2.unit4.ts; % Obtaining the spikes

spiketrain1 = (hist(spiketime1,bins))'; % MUA Binary vector
spiketrain2 = (hist(spiketime2,bins))'; % MUA Binary vector
spiketrain3 = (hist(spiketime3,bins))'; % MUA Binary vector
spiketrain4 = (hist(spiketime4,bins))'; % MUA Binary vector
spiketrain5 = (hist(spiketime5,bins))'; % MUA Binary vector
spiketrain6 = (hist(spiketime6,bins))'; % MUA Binary vector
spiketrain7 = (hist(spiketime7,bins))'; % MUA Binary vector
spiketrain8 = (hist(spiketime8,bins))'; % MUA Binary vector
spiketrain9 = (hist(spiketime9,bins))'; % MUA Binary vector
spiketrain10 = (hist(spiketime10,bins))'; % MUA Binary vector
spiketrain11 = (hist(spiketime11,bins))'; % MUA Binary vector
spiketrain12 = (hist(spiketime12,bins))'; % MUA Binary vector
spiketrain13 = (hist(spiketime13,bins))'; % MUA Binary vector
spiketrain14 = (hist(spiketime14,bins))'; % MUA Binary vector
spiketrain15 = (hist(spiketime15,bins))'; % MUA Binary vector
spiketrain16 = (hist(spiketime16,bins))'; % MUA Binary vector
spiketrain17 = (hist(spiketime17,bins))'; % MUA Binary vector
spiketrain18 = (hist(spiketime18,bins))'; % MUA Binary vector
spiketrain19 = (hist(spiketime19,bins))'; % MUA Binary vector
spiketrain20 = (hist(spiketime20,bins))'; % MUA Binary vector
spiketrain21 = (hist(spiketime21,bins))'; % MUA Binary vector
% spiketrain22 = (hist(spiketime22,bins))'; % MUA Binary vector


MUAvect = [spiketrain1 spiketrain2 spiketrain3 spiketrain4 spiketrain5 ...
    spiketrain6 spiketrain7 spiketrain8 spiketrain10 ...
    spiketrain12 spiketrain13 ...
    spiketrain16 spiketrain17 spiketrain18 spiketrain19 spiketrain20 ...
    ];
%%
MUAvectActive = zeros(size(MUAvect));
for kk = 1 : size(MUAvect,2)
    for qq = 1:size(MUAvect,1)
        if linVel(qq) > 0.05
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
%% Position dependent Firing rates
Posbins = (0: 0.025 : 1.4)';
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
close all
for kk = 1: size(MUAhist,2)
   figure
   bar(Posbins, MUAhist(:,kk),'b') 
   hold on
   bar(Posbins, MUAhist2(:,kk),'r')
   hold off
end
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
nSR = 1.4/(length(MUAhistSm(:,1))-1);
nPosbins = (0 : nSR : 1.4)';
for kk = 1: size(MUAhist2,2)
   figure
   bar(Posbins, MUAhist2(:,kk)) 
   hold on
   plot(nPosbins, MUAhistSm(:,kk),'r', 'LineWidth',2)
   hold off
end
%% Pos/time histogram
Timebins = (0 : 0.2 : LFPt(end))'; % Time bins
PostimeHist = zeros(size(Timebins,1),1); % Memory Preallocation
[a ,b] = hist(linPos(1:Timebins(1+1)*2000), 1); % Fist bin hist
PostimeHist(1,1) = b;
for kk = 2: length(Timebins)-1
    [a ,b] = hist(linPos(Timebins(kk)*2000:Timebins(kk+1)*2000), 1);
    PostimeHist(kk,1) = b;
end
%% Spike matrix 300ms
MUAtimehist = zeros(size(Timebins,1), size(rrFRidx{1,7},2));
for kk = 1: size(MUAtimehist,2)
    MUAtimehist(:,kk) =(hist(linTime(MUAidx{1,kk}),Timebins))';%Hist positions @ spikes.
end

%% Summing the firing rates per position
MUASumFRate = (sum (MUAhistSm,2))';
MUASumFRate2 = (exp(-0.3.* MUASumFRate))';
%% Firing rate elevated by time spikes
MUAposTime = zeros(length(MUAtimehist(:,1)), size(MUAhistSm,1));
for tt = 1: length(MUAtimehist(:,1))
    MUAposTime(tt,:) = (prod(MUAhistSm(:,:).^MUAtimehist(tt,:),2))';
end
%%
MUAposTime = MUAposTime';
%%
MUAprob = MUASumFRate2.* MUAposTime;
%%
MUAprob2 = max(MUAprob);
%%
MUAprob3 = MUAprob./MUAprob2;