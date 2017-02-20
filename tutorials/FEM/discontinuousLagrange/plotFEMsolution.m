clc
clear all
close all

%U=load('U/U_50.txt');
S=load('D/D_0.txt');

%% Plot function
comp={'x_1','x_2','value'};
figure(2)
hold on
%axis equal

clrCol=3;
sMax=max(S(:,clrCol))
sMin=min(S(:,clrCol))
clrID=round((S(:,clrCol)-sMin)/(sMax-sMin)*size(colormap,1));
%caxis([sMin sMax])

%h=colorbar;
for f=1:3:size(S,1)
    %fill(S(f:f+2,1),S(f:f+2,2),clrID(f:f+2))
    fill3(S(f:f+2,1),S(f:f+2,2),S(f:f+2,clrCol),clrID(f:f+2))

end
grid on
xlabel('X_1','FontSize',14)
ylabel('X_2','FontSize',14)
%zlabel('X_3','FontSize',14)

%set(h, 'ylim', [sMin sMax])
title(comp{clrCol},'FontSize',14)

