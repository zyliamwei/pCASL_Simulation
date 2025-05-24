clear all
clc
close all

% Simulating the adibatic inversion process
% ZWei, January 21, 2025.

%% Loading the hanning pulse
hanning=load('hanning.txt');
hanning=hanning(:,1).';
hanarea=sum(hanning);
hanning=hanning./hanarea;
figure;plot(hanning,'r-');
xlabel('Time point');
ylabel('Normalized amplitude');
W_Plot;
title('Hanning Pulse');
% Loading in the Fourier form of Hanning window
load('hanningwin.mat');

%% Gradient and pulse parameters
% Duration of the labeling pulse; unit: ms
labdur=0.4;         
% Inter-labeling-pulse delay, comprising the length of labeling pulse and a
% blank period; unit: ms
labdelay=0.8;       
% Gradient strength accompanying the labeling pulse; unit: gauss/cm; 
% Bruker uses percentage gradient, namely, percentage of the maximal value
% and 74 Gauss/cm is the maxial gradient strength at the 11.7T scanner
fovsatgrad=0.14*74; 
% Gradient strength accompanying the blank period; unit: gauss/cm; 
fovsatgradreph=-0.16*74;
% T1 relaxationt time of blood; unit: ms 
% Referencing Magn Reson Mater Phy 2012;25:245-249.
T1b=2813; 
% T2 relaxation time of blood; unit: ms
% Referencing Magn Reson Imaging 2017;38:234-249.
T2b=54.05;
% Gyromagnetic ratio; unit: Hz/gauss
gama=4260;  
% Spatial range for simulation; labeling pulse is assumed to label the
% spins at -0.5-0.5 mm, covering a slice thickness of 1.0 mm; unit: mm
z=[-5 5]; 
% Flip angle of the labeling pulse; unit: degree
flip=40; 
% Ramp time for gradient; unit: ms
ramp=0.048; 
ramppoint=fix(ramp./labdur.*length(hanning));
% Flow velocity; unit: cm/s
gvel=0.5:0.1:20;
gvel=20.1:0.1:30;
for mi=1:length(gvel)
    vel=gvel(mi);
% Number of pulses seen during the flowing through the designated range
labnum=ceil(abs(diff(z))/10./vel*1000./(2*labdur+3*ramp));
% temporal resolution of the simulation; unit: ms
timeres=labdur./length(hanning);

%% Basic labeling module
labmod_rf=[zeros(1,ramppoint) hanning zeros(1,ramppoint) ...
    zeros(size(hanning)) zeros(1,ramppoint)];
figure;subplot(2,1,1);
plot(labmod_rf,'r');title('RF labeling module');
labmod_grad=[linspace(0,fovsatgrad,ramppoint) fovsatgrad*ones(size(hanning)) ...
    linspace(fovsatgrad,fovsatgradreph,ramppoint) fovsatgradreph*ones(size(hanning))...
    linspace(fovsatgradreph,0,ramppoint)];
subplot(2,1,2);plot(labmod_grad,'b');title('Gradient labeling module');


%% Adibatic flow inversion
label_rf=[];
label_grad=[];
for ni=1:labnum
    label_rf=[label_rf labmod_rf];
    label_grad=[label_grad labmod_grad];
end
figure;subplot(2,1,1);
plot(label_rf,'r');title('Labeling pulse train');
ylabel('Normalized amplitude');
subplot(2,1,2);plot(label_grad,'b');title('Labeling gradient train');
ylabel('Gradient (G/cm)');

disp([num2str(mi) ' / ' num2str(length(vel))]);
% Starting magnetization 
Minit=[0 0 1].';
% Total number of time points
label_count=length(label_rf);
% Dynamic signals during the labeling
sig=zeros(3,label_count);
% Real-time Z position; unit: cm (convert from mm to cm)
label_z=linspace(z(1),z(2),label_count)/10;
% Off-resonance frequency depending on the Z-position; freq=gama*grad*z
label_offres=gama.*label_grad.*label_z;
% Real-time RF flipping angle; unit: degree
label_angle=label_rf*40;

% Display the off-resonance frequency
figure;subplot(2,1,1);
plot(label_angle,'r');
xlabel('Labeling process');ylabel('Angle (degree)');
subplot(2,1,2);plot(label_offres,'r');
xlabel('Labeling process');ylabel('Frequency (Hz)');

for ni=1:length(label_rf)
    % Initialize signal
    if ni==1
        sig(:,ni)=Minit;
    else
        sig(:,ni)=sig(:,ni-1);
    end
    [tempA tempB]=freeprecess(timeres/2,T1b,T2b,label_offres(ni));
    sig(:,ni)=tempA*sig(:,ni)+tempB;
    sig(:,ni)=ZRot(label_angle(ni),'y')*sig(:,ni);
    sig(:,ni)=tempA*sig(:,ni)+tempB; 
end

% figure;subplot(3,1,1);
% plot(sig(1,:),'b');
% xlabel('Labeling process');
% ylabel('X magnetization');
% subplot(3,1,2);
% plot(sig(2,:),'b');
% xlabel('Labeling process');
% ylabel('Y magnetization');
% subplot(3,1,3);
% plot(sig(3,:),'r');
% xlabel('Labeling process');
% ylabel('Z magnetization');

% Calculate
inveff(mi)=(1-mean(sig(3,end-9:end)))/2;

close all
end
    





