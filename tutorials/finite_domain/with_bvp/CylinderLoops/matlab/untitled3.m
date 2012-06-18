clc
clear all
close all

N=60;
for k=0:N

A=load(['A/A_' num2str(k) '.txt']);

A3(1:size(A,1),k+1,:)=A(:,2:end);

end

%%
u1=A3(:,:,1);
u2=A3(:,:,2);
u3=A3(:,:,3);

u1inf=A3(:,:,4);
u2inf=A3(:,:,5);
u3inf=A3(:,:,6);

figure(1)
plot(u1inf')

figure(2)
plot(u2inf')

figure(3)
plot(u3inf')

figure(4)
plot(u1')

figure(5)
plot(u2')

figure(6)
plot(u3')