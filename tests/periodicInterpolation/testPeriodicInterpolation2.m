clear all
close all
clc

system('rm input.txt');

%% Data points for gamma surface
A=[2 1;
    0 sqrt(3)];

N=[3 3];

APB=175;
SISF=10;
CESF=270;
CSF=230;

f=[0 0 0;
    1 sqrt(3)/3 SISF;
    1 0 APB;
    0.5 sqrt(3)/2 APB;
    1.5 sqrt(3)/2 APB;
    0.5 sqrt(3)/6 CESF;
    1 sqrt(3)*2/3 CESF;
    1.5 sqrt(3)/6 CESF;
    ];

df=[1 sqrt(3)/3 1 0 0; % SISF
    1 sqrt(3)/3 0 1 0; % SISF
    1 0 1 0 0; % APB
    0.5 sqrt(3)/2  0.5 sqrt(3)/2 0; % APB
    0.5 sqrt(3)/6 1 0 0;
    0.5 sqrt(3)/6 0 1 0;
    1 sqrt(3)*2/3 1 0 0;
    1 sqrt(3)*2/3 0 1 0;
    1.5 sqrt(3)/6 1 0 0;
    1.5 sqrt(3)/6 0 1 0;
    ];



%% Write input file
fid=fopen('input.txt','w')
printMatrixToFile(fid,A,'A');
printMatrixToFile(fid,N,'N');
printMatrixToFile(fid,f,'f');
printMatrixToFile(fid,df,'df');
fclose(fid)

%% Write input file
fid=fopen('input.txt','w')
printMatrixToFile(fid,A,'A');
printMatrixToFile(fid,N,'N');
printMatrixToFile(fid,f,'f');
printMatrixToFile(fid,df,'df');
fclose(fid)

%% run code
system('./test')

%% Load data and plot
data=load('output.txt')

[X,Y] = meshgrid([0:0.01:2],[0:0.01:sqrt(3)*2]);

f=zeros(size(X));
for i=1:size(data,1)
    k=data(i,[1 2]);
    S=data(i,3);
    C=data(i,4);
    f=f+S*sin(k(1)*X+k(2)*Y)+C*cos(k(1)*X+k(2)*Y);
end


fMax=max(max(f));

figure
clf
hold on
surf(X,Y,f)
plot3([0 A(1,1)],[0 A(2,1)],[fMax fMax]+1,'m','Linewidth',2)
plot3([0 A(1,2)],[0 A(2,2)],[fMax fMax]+1,'m','Linewidth',2)
grid on
xlabel('x')
ylabel('y')
%axis equal
colormap jet
colorbar

function printMatrixToFile(fid,M,label)
fprintf(fid,[label '=']);

format='';
for(c=1:size(M,2))
    format=[format '%1.15e '];
end

for(k=1:size(M,1))
    if k<size(M,1)
        fprintf(fid,[format '\n'],M(k,:));
    else
        fprintf(fid,[format ';\n'],M(k,:));
    end
end
end