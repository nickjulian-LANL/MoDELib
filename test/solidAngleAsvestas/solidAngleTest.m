clear all
close all
clc

L1=10;
L2=10;

nx=50;
ny=50;

X=(rand(3,1)-[1;1;1]*0.5)*30

P=[0  0  0;
   L1 0  0;
   L1 L2 0;
   0  L2 0];

figure(1)
axis equal
hold on
plot3([P(:,1);P(1,1)],[P(:,2);P(1,2)],[P(:,3);P(1,3)])
plot3(X(1),X(2),X(3),'ro')

dx=L1/nx;
dy=L2/ny;
dA=dx*dy;
S=0;
for i=1:nx
for j=1:ny
    X1=[-dx/2+i*dx;-dy/2+j*dy;0];
%    plot3(X1(1),X1(2),X1(3),'.')
    r=X1-X;
    S=S+r(3)/norm(r)^3*dA;
end
end

S

I=0;
P=P';
q=[0 0 -1];
for i=1:4
    i1=i+1;
if i1==5
    i1=1;
end

    Rn=P(:,i)-X;
    Rm=P(:,i1)-X;
    Rn=Rn/norm(Rn);
    Rm=Rm/norm(Rm);
    
    RdR=dot(Rn,Rm);
    RcR=cross(Rn,Rm);
    RcR=RcR/norm(RcR);
    A=sqrt((1-RdR)/(1+RdR));
    
    num=dot(q,RcR)*A;
    den=1-dot(q,Rn)-dot(q,cross(Rn,RcR))*A;
    
    I=I+2*atan(num/den);
    
end

I