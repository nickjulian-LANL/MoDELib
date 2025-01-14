clear all
close all
clc
fontSize=16;

C=[0 0 0];
a=1;
eT=eye(3)*0.01;
nu=1/3;
mu=1;

Lx=3;
Ly=3;
Lz=3;

nx=100;
ny=100;
nz=100;

x=[-nx:nx]/nx*Lx;
y=[-ny:ny]/ny*Ly;
z=[-nz:nz]/nz*Lz;

%% Write input file
fid=fopen('inputFile.txt','w');
printMatrixToFile(fid,C,'C');
printMatrixToFile(fid,a,'a');
printMatrixToFile(fid,eT,'eT');
printMatrixToFile(fid,nu,'nu');
printMatrixToFile(fid,mu,'mu');
printMatrixToFile(fid,x,'x');
printMatrixToFile(fid,y,'y');
printMatrixToFile(fid,z,'z');
fclose(fid)

%% run code
system('./test')

%% Load output
data=load('inclusion.txt');
Nx=2*nx+1;
Ny=2*ny+1;
Nz=2*nz+1;
xx=reshape(data(:,1),Ny,Nx,Nz);
yy=reshape(data(:,2),Ny,Nx,Nz);
zz=reshape(data(:,3),Ny,Nx,Nz);


%% Plot
colLabels={'x','y','z','s11','s12','s13','s21','s22','s23','s31','s32','s33'};
plotCols=[4 5 6 8 9 12];

for col=plotCols

curData=reshape(data(:,col),Ny,Nx,Nz);

maxData=max(data(:,col));
minData=min(data(:,col));

nIso=5;
dIso=(maxData-minData)/(nIso-1);
isoLevs=[minData:dIso:maxData];

figure(1)
clf
hold on
caxis([minData maxData]);
cmap = colormap;
m = length(cmap);
for isoLev=isoLevs
p = patch(isosurface(xx,yy,zz,curData,isoLev),'FaceAlpha',0.2);
isonormals(xx,yy,zz,curData,p);
cIndex = fix((isoLev-minData)/(maxData-minData)*m)+1; %A
RGB = ind2rgb(cIndex,cmap);
p.FaceColor = RGB;
p.EdgeColor = 'none';
end
daspect([1 1 1])
view(3); 
axis image
camlight 
lighting gouraud
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')
set(gca,'FontSize',fontSize)
colorbar
grid on
print(gcf,'-depsc',['stress_' colLabels{col}]);
end
