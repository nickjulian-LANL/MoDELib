clear all
close all
clc

system('rm input.txt');
system('rm output.txt');

%% Data points for gamma surface


%structure='gamma';
%structure='gammaPrime';
%structure='hcpBasal';
structure='hcpPrismatic';


fid=fopen('input.txt','w');

switch structure
    
    case 'gamma'
        
        A=[1 0.5;
            0 sqrt(3)/2];
        
        N=[2 2];
        D=[1 1];
        
        f=[0 0 0
            0.5 sqrt(3)/6 42;    % ISF
            0.25 sqrt(3)/12 182; % USF
            0.75 sqrt(3)/12 182; % USF
            0.5  sqrt(3)/3 182; % USF
            ];
        
        df=[0.25 sqrt(3)/12 -0.5 sqrt(3)/2 0; % USF
            0.75 sqrt(3)/12 0.5 sqrt(3)/2 0; % USF
            %            0.5  sqrt(3)/3 1 0 0; % USF
            ];
        
        printMatrixToFile(fid,N,'N');
        printMatrixToFile(fid,D,'D');
        readWaveVectors=0;
        
        cutDirs=[A(:,1) A(:,2) A(:,1)+A(:,2)];
        cutDirsLabels={'a_1', 'a_2', 'a_1+a_2'};
        cutDirsStyles={'k-', 'm--', 'b-'};

    case 'gammaPrime'
        
        A=[1 0.5;
            0 sqrt(3)/2];
        
        N=[3 3];
        D=[2 2];
        
        
        waveVec_ctr=0;
        
        for i=0:N(1)-1
            for j=0:N(2)-1
                % if( ~(i==(N(1)-1) & j==(N(2)-1) ) )
                waveVec_ctr= waveVec_ctr + 1;
                waveVec(waveVec_ctr,:)=[i/D(1) j/D(2)];
                % end
            end
        end
        
        %extra_points of waveVec
        
        % % N1=N2=3 case
        waveVec_ctr=waveVec_ctr+1;
        waveVec(waveVec_ctr,:)= [1/D(1) -1/D(2)];
        
        APB=175;
        SISF=10;
        CESF=270;
        CISF=230;
        SESF=75;
        
        
        
        
        f=[0 0 0;
            1 sqrt(3)/3 SISF;
            1 0 APB;
            0.5 sqrt(3)/2 APB;
            1.5 sqrt(3)/2 APB;
            0.5 sqrt(3)/6 CESF;
            1 sqrt(3)*2/3 CESF;
            1.5 sqrt(3)/6 CESF; % ---- till here cell 1
            2 2*sqrt(3)/3 SESF;
            2 sqrt(3)/3 CISF;
            1.5 5*sqrt(3)/6 CISF;
            2.5 5*sqrt(3)/6 CISF;
            ]; %  12 conditions
        
        
        df=[0 0 1 0 0; %PC
            0 0 0 1 0; %PC
            1 sqrt(3)/3 1 0 0; %SISF
            1 sqrt(3)/3 0 1 0; %SISF
            1 0 1 0 0; %APB
            0.5 sqrt(3)/2 0.5 sqrt(3)/2 0; %APB
            1.5 sqrt(3)/2 -0.5 sqrt(3)/2 0 %APB
            ]; % 7 conditions
        
        printMatrixToFile(fid,waveVec,'waveVecInt');
        readWaveVectors=1;
        
        cutDirs=2*[A(:,1) A(:,2) A(:,1)+A(:,2)];
        cutDirsLabels={'a_1', 'a_2', 'a_1+a_2'};
        cutDirsStyles={'k-', 'm--', 'b-'};

        
    case 'hcpBasal'
        
        A=[1 0.5;
            0 sqrt(3)/2];
        
        D=[1 1];
        
        ISF=213;
        USF=261;
        
        waveVec=[0  0
            1  0
            0  1
            1  1
            1  2
            2  1
            1  -1
            ];
        
        
        f=[0 0 0
            0.5 sqrt(3)/6 ISF;    % ISF
            0.25 sqrt(3)/12 USF; % USF
            0.75 sqrt(3)/12 USF; % USF
            0.5  sqrt(3)/3 USF; % USF
            1 sqrt(3)/3 653
            ];
        
        df=[
            0 0 1 0 0
            0 0 0 1 0
            1 sqrt(3)/3 1 0 0
            1 sqrt(3)/3 0 1 0
            0.5 0 1 0 0
            0.25 sqrt(3)/4 0.5 sqrt(3)/2 0
            0.75 sqrt(3)/4 -0.5 sqrt(3)/2 0
            ];
        
        printMatrixToFile(fid,waveVec,'waveVecInt');
        readWaveVectors=1;
        
        figure(3)
        clf
        xx = [0.05:0.05:0.95];
        yy = [30.79 93.21 174.8 242.99 256.51 225.92 216.43 272.143 351.82 448.77 ...
            540.92 616.76 653.3 648.6 585.4 461.8 303.7 155.1 43];
        plot(xx,yy,'o')
        
        cutDirs=[A(:,1) A(:,2) A(:,1)+A(:,2)];
        cutDirsLabels={'a_1', 'a_2', 'a_1+a_2'};
        cutDirsStyles={'k-', 'm--', 'b-'};

    case 'hcpPrismatic'
        c=sqrt(8/3);
        A=[1 0;
            0 c];
        
        D=[1 1];
        
        waveVec=[
            0  0
            1  0
            0  1
            1 1
            1 -1
            2 0];
        
        
        h=c/4; % cannot be 0 or c/2
        hv=645;
        
        f=[0 0 0;
            1/2 0 211
            0 c/2 1350
            1/2 c/2 700
            1/4 h hv
            3/4 h hv
            ];
% 
% f=[ 
%             0 0 0;
%             0.5 0 211;
%             0 c/2 1350;
%             0.25 h 645;
%             0.75 h 645;
%             0.5 c/2 700;
%             ];.
        
        df=[0 0 1 0 0
            0 0 0 1 0
            1/2 0 1 0 0
            1/2 c/2 1 0 0
            1/2 c/2 0 1 0
            ];
        
        printMatrixToFile(fid,waveVec,'waveVecInt');
        readWaveVectors=1;
        
        
        
        figure(3)
        clf
        xx = [0.05:0.05:0.95];
        
        yy = [12.23 48.38 93.21 139.05 173.68 205.26 220.54 220.07 215.45 210.866 ...
            216.47 222.07 221.56 203.22 174.19 138.54 92.7 48.9 12.22];
        
        plot(xx,yy,'ko')
        
        cutDirs=[A(:,1) A(:,2) A(:,1)+A(:,2)];
        cutDirsLabels={'a', 'c', 'a+c'}
        cutDirsStyles={'k-', 'm-', 'b-'};
end

B=2*pi*inv(A');

printMatrixToFile(fid,A,'A');
printMatrixToFile(fid,f,'f');
printMatrixToFile(fid,df,'df');
fclose(fid)

%% run code
system(['./test ' num2str(readWaveVectors)])

%% Load data and plot
data=load('output.txt')

%% Plot Wave vectors
figure(1)
hold on;
grid on
axis equal
plot(data(:,1)/2/pi,data(:,2)/2/pi,'bo','Linewidth',2)
plot(-data(:,1)/2/pi,-data(:,2)/2/pi,'mx','Linewidth',2)
plot([-max(data(:,1)) max(data(:,1))]/2/pi,[0 0],'k')
plot([0 0],[-max(data(:,2)) max(data(:,2))]/2/pi,'k')
title('wave vectors / 2\pi')
set(gca,'Fontsize',16)
print([structure '_waveVectors'], '-dpng', '-r300');


%% Construct GammaSurface and plot it
np=400;
%[X,Y] = meshgrid([0:np]/np*2*D(1),[0:np]/np*sqrt(3)/2*2*D(2));
[X,Y] = meshgrid([0:np]/np*sum(cutDirs(1,[1 2])),[0:np]/np*sum(cutDirs(2,[1 2])));
fs=zeros(size(X));
for i=1:size(data,1)
    k=data(i,[1 2]);
    S=data(i,3);
    C=data(i,4);
    fs=fs+S*sin(k(1)*X+k(2)*Y)+C*cos(k(1)*X+k(2)*Y);
end
fMax=max(max(fs));

figure(2)
hold on
surf(X,Y,fs,'edgecolor','none')
plot3([0 D(1)*A(1,1)],[0 D(1)*A(2,1)],[fMax fMax]+1,'m','Linewidth',2)
plot3([0 D(2)*A(1,2)],[0 D(2)*A(2,2)],[fMax fMax]+1,'m','Linewidth',2)
plot3([D(2)*A(1,2) D(1)*A(1,1)],[D(2)*A(2,2) D(1)*A(2,1)],[fMax fMax]+1,'m','Linewidth',2)
plot3([D(2)*A(1,2) D(2)*A(1,2)+D(1)*A(1,1)],[D(2)*A(2,2)  D(2)*A(2,2)],[fMax fMax]+1,'m','Linewidth',2)
plot3([D(1)*A(1,1) D(2)*A(1,2)+D(1)*A(1,1)],[D(1)*A(2,1)  D(2)*A(2,2)],[fMax fMax]+1,'m','Linewidth',2)
plot3(f(:,1),f(:,2),f(:,1)*0+fMax+1,'ok','Linewidth',2)
%plot3(df(:,1),df(:,2),df(:,1)*0+fMax+1,'xg','Linewidth',2)
quiver3(df(:,1),df(:,2),df(:,1)*0+fMax+1,df(:,3),df(:,4),df(:,1)*0,0.2,'k','Linewidth',2)
grid on
xlabel('x')
ylabel('y')
h = get(gca,'DataAspectRatio');
if h(3)==1
    set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
else
    set(gca,'DataAspectRatio',[1 1 h(3)])
end
colormap jet
h=colorbar
h.Label.String = 'mJ/m^2';
axis image
print([structure '_gammaSurface'], '-dpng', '-r300');

%% Plot cuts of GammaSurface
figure(3)
hold on
cutPlots=[];
for c=1:size(cutDirs,2)
    fcut=functionCut(data,cutDirs(:,c),np);
    cutPlots=[cutPlots plot([0:(np-1)]/np,fcut,cutDirsStyles{c},'Linewidth',2)];
end
grid on
xlabel('reaction coordinate')
ylabel('\gamma-surface [mJ/m^2]')
legend(cutPlots,cutDirsLabels)
set(gca,'Fontsize',16)
print([structure '_gammaSurfaceCuts'], '-dpng', '-r300');


