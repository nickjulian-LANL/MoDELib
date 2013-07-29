close all
clear all
clc

Nx=1001;
Ny=101;

Sn=load('S/S_0.txt');
X=reshape(Sn(:,1),Nx,Ny);
Y=reshape(Sn(:,2),Nx,Ny);
s11n=reshape(Sn(:,4),Nx,Ny);
s12n=reshape(Sn(:,5),Nx,Ny);
s22n=reshape(Sn(:,8),Nx,Ny);

%% Plot
figure(1)
%plot3(X,Y,s11n)
surf(X,Y,s11n,'EdgeAlpha',0.05)

figure(2)
plot3(X,Y,s12n)

figure(3)
plot3(X,Y,s22n)