clear all
close all
clc

T=load('T/T_0.txt');
N=load('N/N_0.txt');

figure(1)
plot3(N(:,2),N(:,3),N(:,4),'.k')
axis equal

%% edges

hold on
id=[3 923]+1;

figure(1)
plot3(N(id,2),N(id,3),N(id,4),'-r','Linewidth',2)