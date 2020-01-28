% Newmark's Method with constant average acceleration
% Jared Rivera 2019

clear all; close all; clc;

%Initialization
dispdata=csvread('DispEdit.csv');
acceldata=csvread('AccelEdit.csv');
numsteps=length(acceldata);

uddg=acceldata(:,3);

ud0=0;
u0=0;

%Decay estimation
decaytime=dispdata(13000:23000,1);
decaydisp=dispdata(13000:23000,3);

%Experimental parameters
k=1.079;
m=3.465;

zeta=0.135; %Damping ratio
c=zeta*2*(m*k)^0.5;

gamma=0.5;
beta=0.25;

p=-m*uddg;
p0=p(1);

%Initial calculations
udd0=(p0-c*ud0-k*u0)/m;
dt=0.1*(k/m)^0.5;
khat=k+(gamma/beta/dt)*c+(1/beta/dt^2)*m;
khatinv=khat^(-1);

%Updates per timestep
u=zeros(1,numsteps);
u(1)=u0;
ud=zeros(1,numsteps);
ud(1)=ud0;
udd=zeros(1,numsteps);
udd(1)=udd0;
for i=2:numsteps 
    dp=p(i)-p(i-1);
    
    dphat=dp+(m/(beta*dt)+gamma*c/beta)*ud(i-1)+(m/(2*beta)+dt*(1-gamma/(2*beta))*c)*udd(i-1);
    du=khatinv*dphat;
    dud=(gamma/(beta*dt))*du-(gamma/beta)*ud(i-1)+(dt*(1-gamma/(2*beta)))*udd(i-1);
    dudd=(1/(beta*dt^2))*du-(1/(beta*dt))*ud(i-1)-(1/(2*beta))*udd(i-1);
    
    u(i)=u(i-1)+du;
    ud(i)=ud(i-1)+dud;
    udd(i)=udd(i-1)+dudd;  
end

time=linspace(0,max(acceldata(:,1)),numsteps);

timeoffset=21.35;

figure(1)
plot(time,u,dispdata(:,1)-timeoffset,dispdata(:,3))
title('Displacement Comparison')
xlabel('Time (s)')
ylabel('Displacement (in)')
legend('Newmark','Experiment')
xlim([10,290])

figure(2)
plot(time,-udd,acceldata(:,1),acceldata(:,3))
title('Acceleration Comparison')
xlabel('Time (s)')
ylabel('Acceleration (g)')
legend('Newmark','Experiment')
ylim([-0.5,0.5])
xlim([10,290])
