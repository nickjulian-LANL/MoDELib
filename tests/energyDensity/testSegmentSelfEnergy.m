clear all
close all
clc

system('rm input.txt')
system('rm points.txt')

computeSelfEnergiesOnly=0;
nGP=100;
Nseg=100;

%% Write input file
fid=fopen('input.txt','w');
fprintf(fid,['nGP=' num2str(nGP) ';\n']);
segformat='%1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e %1.15e \n';
    fprintf(fid, "S=");
    segments=zeros(Nseg,9);
for s=1:Nseg
    P0=(rand(1,3)-0.5)*10^(3*rand(1));
    P1=(rand(1,3)-0.5)*10^(3*rand(1));
    b=(rand(1,3)-0.5); b=b/norm(b);

    fprintf(fid,segformat, [P0 P1 b]);
segments(s,:)=[P0 P1 b];
end
fprintf(fid, ";");
fclose(fid)

%% Run code
system(['./test ' num2str(computeSelfEnergiesOnly)])

%% Read output file
points=load('points.txt');

%% Plots
figure(1)
clf
hold on
for s=1:size(segments,1)
plot3([segments(s,1) segments(s,4)],[segments(s,2) segments(s,5)],[segments(s,3) segments(s,6)])
end

for p=1:size(points,1)
plot3(points(p,3),points(p,4),points(p,5),'k.')
end

figure(2)
clf
hold on
segIDs=unique(points(:,1));
for k=1:length(segIDs)
    s=segIDs(k);
pointIDs=find(points(:,1)==s);
E(k)=sum(points(pointIDs,6).*points(pointIDs,7)); % integral e*dL
W(k)=mean(points(pointIDs,8));
plot(points(pointIDs,6))
end

figure(3)
clf
yyaxis left
plot(segIDs,E,'bo')
hold on
plot(segIDs,W,'rx')
yyaxis right
plot(segIDs,abs(E-W)./W)
grid on
xlabel('segment ID')