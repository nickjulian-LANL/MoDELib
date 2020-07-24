clear all
close all
clc
fontSize=16;

material='Zr';
%materialFile='../../tutorials/DislocationDynamics/MaterialsLibrary/W.txt'
%materialFile='../../tutorials/DislocationDynamics/MaterialsLibrary/Cu.txt'
materialFile=['../../tutorials/DislocationDynamics/MaterialsLibrary/' material '.txt'];

system(['./mobility ' materialFile])

nT=101;

dataS=load('velocityS.txt');
plotMobility(dataS,nT,fontSize)
print(gcf,[material '_prismScreMobility'], '-dpng', '-r300');

dataE=load('velocityE.txt');
plotMobility(dataE,nT,fontSize)
print(gcf,[material '_prismEdgewMobility'], '-dpng', '-r300');


function plotMobility(data,nT,fontSize)
datasize=size(data,1);

T=reshape(data(:,2),datasize/nT,nT);
S=reshape(data(:,1),datasize/nT,nT);
V=reshape(data(:,3),datasize/nT,nT);

figure
clf
hold on
surf(T,log10(S),log10(V+1e-10),'edgecolor','none','FaceAlpha',0.2)

isolevels=[1e-3];
for e=[-7:-1]
    for k=[10]
        isolevels=[isolevels k*10^e];
    end
end
%isolevels=[BrunnerV/cs isolevels];
%isolevels=[0.05:0.05:1]*cs
[C1,h1]=contour(T,log10(S),V,isolevels,'k','Linewidth',1);
clabel(C1,h1,'FontSize',fontSize,'Color','k','labelspacing', 700)
xlabel('T/Tm','FontSize',fontSize)
ylabel('log_{10}(\tau/\mu)','FontSize',fontSize)
grid on
end
