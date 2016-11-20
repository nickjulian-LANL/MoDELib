clc
close all
clear all

MODEL_DIR='../../../..';

%% Define output file name
meshID=0; % creates ../N/N_1.txt and ../T/T_1.txt
targetElements=1e4;

filename='tricrystal'; % this creates file tricrystal.poly
t=200;      	% Film thickness (of columnar grains)
L=400;      	% Side length of cube containing 3 crystals
X=0.5*L*tan(acos(-1/3)-pi/2);
V=L*L*t;	%Volume of the triple crystal

%t=t/2;
%L=L/2;



P(1,:)=[0,-t/2,0];
P(2,:)=[L/2,-t/2,-X];
P(3,:)=[-L/2,-t/2,-X];
P(4,:)=[L/2,-t/2,-L/2];
P(5,:)=[-L/2,-t/2,-L/2];

P(6,:)=[L/2,-t/2,0];
P(7,:)=[L/2,-t/2,L/2];
P(8,:)=[0,-t/2,L/2];
P(9,:)=[-L/2,-t/2,L/2];
P(10,:)=[-L/2,-t/2,0];


for k=1:10
P(10+k,:)=P(k,:)+[0,t,0];
end


figure(1)
clf
plot3(P(:,1),P(:,2),P(:,3),'ro','Linewidth',2)
hold on
text(P(:,1),P(:,2),P(:,3),num2str([1:size(P,1)]'),'FontSize',16)
axis equal
xlabel('x')
ylabel('y')
grid on


%% FRONT BOTTOM 	~~~~~~~~~~~~ pentagon grain face points:
a=1;
b=2;
c=3;
d=4;
e=5;

k=1;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
FT(k+2,:)=[c,d,e];
k=k+3;
%% FRONT RIGHT 		~~~~~~~~~~~~ pentagon grain face points:
a=1;
b=2;
c=8;
d=6;
e=7;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
FT(k+2,:)=[c,d,e];
k=k+3;
%%  FRONT LEFT	 	~~~~~~~~~~~~ pentagon grain face points:
a=1;
b=3;
c=8;
d=10;
e=9;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
FT(k+2,:)=[c,d,e];
k=k+3;



%% BACK BOTTOM 	~~~~~~~~~~~~ pentagon grain face points:
a=1+10;
b=2+10;
c=3+10;
d=4+10;
e=5+10;

FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
FT(k+2,:)=[c,d,e];
k=k+3;
%% BACK RIGHT 		~~~~~~~~~~~~ pentagon grain face points:
a=1+10;
b=2+10;
c=8+10;
d=6+10;
e=7+10;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
FT(k+2,:)=[c,d,e];
k=k+3;
%%  BACK LEFT	 	~~~~~~~~~~~~ pentagon grain face points:
a=1+10;
b=3+10;
c=8+10;
d=10+10;
e=9+10;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
FT(k+2,:)=[c,d,e];
k=k+3;



%%  $$$$$$$$$ INTERNAL MIDDLE LEFT GRAIN BOUNDARY	$$$$$$$$$$$$$$$$$$$$$$$ 
a=1;
b=3;
c=11;
d=13;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
k=k+2;
%%  $$$$$$$$$ INTERNAL MIDDLE RIGHT GRAIN BOUNDARY	$$$$$$$$$$$$$$$$$$$$$$$ 
a=1;
b=2;
c=11;
d=12;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
k=k+2;
%%  $$$$$$$$$ INTERNAL MIDDLE RIGHT GRAIN BOUNDARY	$$$$$$$$$$$$$$$$$$$$$$$ 
a=1;
b=8;
c=11;
d=18;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
k=k+2;

%%  $$$$$$$$$ BOTTOM CUBE FACE
a=4;
b=5;
c=14;
d=15;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
k=k+2;

%%  $$$$$$$$$ RIGHT-BOTTOM CUBE FACE
a=2;
b=4;
c=a+10;
d=b+10;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
k=k+2;

%%  $$$$$$$$$ RIGHT-MID CUBE FACE
a=2;
b=6;
c=a+10;
d=b+10;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
k=k+2;

%%  $$$$$$$$$ RIGHT-TOP CUBE FACE
a=7;
b=6;
c=a+10;
d=b+10;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
k=k+2;

%%  $$$$$$$$$ TOP-RIGHT CUBE FACE
a=7;
b=8;
c=a+10;
d=b+10;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
k=k+2;

%%  $$$$$$$$$ TOP-LEFT CUBE FACE
a=8;
b=9;
c=a+10;
d=b+10;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
k=k+2;

%%  $$$$$$$$$ LEFT-TOP CUBE FACE
a=9;
b=10;
c=a+10;
d=b+10;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
k=k+2;

%%  $$$$$$$$$ LEFT-MID CUBE FACE
a=3;
b=10;
c=a+10;
d=b+10;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
k=k+2;


%%  $$$$$$$$$ LEFT-BOTTOM CUBE FACE
a=3;
b=5;
c=a+10;
d=b+10;
FT(k,:)=[a,b,c];
FT(k+1,:)=[b,c,d];
k=k+2;

FT
for f=1:size(FT,1)
    vID=FT(f,:);
    fill3(P(vID,1),P(vID,2),P(vID,3),'g','Facealpha',0.2)
    drawnow
    %pause(0.1)
end

%% Define Regions
% Each row in the following matrix correponds to a region.
% Each row contains the vertices that bound that region.
RP=[	[1  2  3  4  5  11 12 13 14 15]
	[1  2  6  7  8  11 12 16 17 18]
	[1  8  10 9  3  11 18 20 19 13];];

size(RP);
% The rows of the matrix Xm are the barycenters of each region
for r=1:size(RP,1)
Xm(r,:)=mean(P(RP(r,:),:));
end
plot3(Xm(:,1),Xm(:,2),Xm(:,3),'bx','Linewidth',2)
text(Xm(:,1),Xm(:,2),Xm(:,3),num2str([1:size(Xm,1)]'),'FontSize',16)

%% Write .poly file
polyFile = fopen([filename '.poly'],'w');

% Part 1- the node list.
pointFormat='%i %.15e %.15e %.15e \n';
fprintf(polyFile,'# Part 1 - the node list.\n');
fprintf(polyFile,'%i 3 0 0 \n',size(P,1));  % number of nodes
for k=1:size(P,1)
fprintf(polyFile,pointFormat,[k P(k,:)]);
end
fprintf(polyFile,'\n');
fprintf(polyFile,'\n');

% Part 2- the facet list.
fprintf(polyFile,'# Part 2 - the facet list.\n');
fprintf(polyFile,'%i 0 \n',size(FT,1)); % number of facets
for r=1:size(FT,1)
fprintf(polyFile,'1 \n'); 
fprintf(polyFile,'3   %d %d %d \n',FT(r,:));
end

% Part 3- the hole list.
fprintf(polyFile,'# Part 3 - the hole list.\n');
fprintf(polyFile,'0 \n\n');

% Part 4- the region list.
fprintf(polyFile,'# Part 4 - the region list.\n');
fprintf(polyFile,'%i \n',size(Xm,1));

meshSize=ones(size(Xm,1),1)*V/targetElements;
for r=1:size(Xm,1)
fprintf(polyFile,'%i %.15e %.15e %.15e %i %.15e \n',[r Xm(r,:) r meshSize(2)]);
end

fclose(polyFile);

%% Run Tetgen 
system([MODEL_DIR '/scripts/tetgenPOLY.sh ' filename]);

%% Create T and N files and clean tetgent output
system([MODEL_DIR '/scripts/tetgen2TN.sh ' filename ' ' num2str(meshID)]);


%% FT

% C2G1 =  [[-5 2 5]/sqrt(54);
%         [1  0  1]/sqrt(2);
%         [1 5 -1]/sqrt(27)]
% 
%     C2G2 = [[1 1  -1]/sqrt(3);
%        [1  0  1]/sqrt(2);
%        [1 -2 -1]/sqrt(6)]
       
%C2G3 = [[1 -1 -1]/sqrt(3);
%       [1  0  1]/sqrt(2);
%       [-1 -2  1]/sqrt(6)]
    
C2G2=[[1 1 -1]/sqrt(3);[1 0 1]/sqrt(2);[1 -2 -1]/sqrt(6)]

R12=angleAxis([1 0 1]',acos(-1/3));
C2G1=C2G2*R12

R13=angleAxis([1 0 1]',-acos(-1/3));
C2G3=C2G2*R13