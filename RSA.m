% Newmark's Method with constant average acceleration
% Applied to Response Spectrum Analysis
% Jared Rivera 2019

clear all; close all; clc;

%Initialization
dispdata=csvread('DispEdit.csv');
acceldata=csvread('AccelEdit.csv');
numsteps=length(acceldata);

uddg=acceldata(:,3);

ud0=0;
u0=0;

T=linspace(0.1,3,3/0.1);
%c=0.5221; % experimental case
c=0.1934; % 5 percent case
m=3.465;
k=(2*pi./T)*m;

gamma=0.5;
beta=0.25;

p=-m*uddg;
p0=p(1);

response=zeros(length(T));

for j=1:length(T)
    
    %Initial calculations
    udd0=(p0-c*ud0-k(j)*u0)/m;
    dt=0.1;
    khat=k(j)+(gamma/beta/dt)*c+(1/beta/dt^2)*m;
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
    
    response(j)=max(udd);
    
end

plot(T,response-0.225)
title('RSA with 5% damping')
xlabel('Structural Period')
ylabel('Spectral Response Acceleration')
ylim([0,1.3])


