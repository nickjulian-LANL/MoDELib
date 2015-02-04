clear all
close all
clc

F=load('F_0.txt');

runID=F(:,1);
dt=F(:,2);
ddLength=F(:,3);
u3=F(:,4);
f3=F(:,5);

L=1000; % cube side length
V=L^3;
A=L^2;

e33=u3/L;
s33=f3/A;

figure(1)
plot(-e33,-s33)
grid on

figure(5)
plot(runID,u3)
grid on


figure(2)
plot(runID,-s33)
grid on

figure(3)
subplot(2,1,1)
plot(runID,F(:,6:end)/V)
v=axis;
axis([0 max(runID) v(3) v(4)])
xlabel('runID')
grid on

subplot(2,1,2)
plot(runID,dt)
axis([0 max(runID) 0 1.5*max(dt)])
grid on
xlabel('runID')
ylabel('dt')

figure(4)
plot(runID,cumsum(F(:,6:end).*repmat(dt,1,size(F(:,6:end),2))/V))

