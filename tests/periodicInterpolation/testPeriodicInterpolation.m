clear all
close all
clc

system('rm input.txt');

%% Data points for gamma surface
A=[1 0.5;  
  0 8.660254037844386e-01];

N=[2 2];

f=[0 0 0
  0.5 2.886751345948129e-01 42;
  0.25 1.443375672974064e-01 182];

df=[0.5 2.886751345948129e-01  1 0 0;
    0.5 2.886751345948129e-01 0 1 0;
    0.25 1.443375672974064e-01 1 0 0;
    0.25 1.443375672974064e-01 0 1 0];

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

[X,Y] = meshgrid([0:0.01:1],[0:0.01:sqrt(3)]);

f=zeros(size(X));
for i=1:size(data,1)
k=data(i,[1 2]);
S=data(i,3);
C=data(i,4);
    f=f+S*sin(k(1)*X+k(2)*Y)+C*cos(k(1)*X+k(2)*Y);
end

figure
clf
hold on
surf(X,Y,f)
grid on
xlabel('x')
ylabel('y')

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