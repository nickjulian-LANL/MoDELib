close all
clear all
clc

s0=rand(6,1);
e0=rand(6,1);
C=rand(6,6);
C=C'*C;
eP=rand(6,1);

%alpha=[1e10 0e10 1e10 1e10 1e10 1e10];
alpha=[0 1 1 0 0 0]*1e3;
%alpha=[0 0 0 0 0 0];

strategy=2;

switch strategy
    case 1 % OLD WRONG STRATEGY
        A=diag(alpha);
        A1=diag(1./(1+alpha));
        A2=diag(alpha./(1+alpha));
        s=A1*s0+A2*C*(e0-eP);
    case 2 % NEW WORKING STRATEGY
        D=diag(alpha)*diag(diag(C));
        eye(6)+D*inv(C);
        A=inv(eye(6)+D*inv(C));
        B=A*D;
        C+D
        inv(C+D)
        inv(C+D)*D
        %inv(C)*A
        %inv(C)*B
        %s=inv(eye(6)+D*inv(C))*(s0+A*(e0-eP));
        s=A*s0+B*(e0-eP);
end
e=inv(C)*s+eP;

[e-e0 s-s0]
