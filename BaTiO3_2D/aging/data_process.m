clc;clear;close all;
step=[200,500,1029,1045,1083,1200,1500];
n=length(step);
hyloop_ini=cell(n,1);
l=0;
for i=step
    l=l+1;
    load(['step',num2str(step(l)),'.mat']);
    hyloop_ini{l}=[E_hyloop(40:200)',Px_mean(40:200)];
end
% plot hyloop
for i=1:n
    plot(hyloop_ini{i}(:,1),hyloop_ini{i}(:,2));
    hold on;
end
legend('200','500','1029','1045','1083','1200','1500');
% 
% 
hyloop=cell(n,1);
hyloop_neg=cell(n,1);
hyloop_forward=cell(n,1);
energy=zeros(n,1);
energy2=zeros(n,1);
for i=1:2
    hyloop{i,1}=[hyloop_ini{i}(1:41,1),hyloop_ini{i}(1:41,2)];
    energy(i,1)=intg_Newton(flipud(hyloop{i,1}(:,2)),flipud(hyloop{i,1}(:,1)))/1e6;
    hyloop_forward{i,1}=[[1.5631e7;hyloop_ini{i}(132:161,1)],[0;hyloop_ini{i}(132:161,2)]];
    energy2(i,1)=intg_Newton(hyloop_forward{i,1}(:,2),hyloop_forward{i,1}(:,1))/1e6;
end
for i=3:7
    hyloop{i,1}=[hyloop_ini{i}(1:41,1),hyloop_ini{i}(1:41,2)];
    energy(i,1)=intg_Newton(flipud(hyloop{i,1}(:,2)),flipud(hyloop{i,1}(:,1)))/1e6;
    hyloop_forward{i,1}=[hyloop_ini{i}(121:161,1),hyloop_ini{i}(121:161,2)];
    energy2(i,1)=intg_Newton(hyloop_forward{i,1}(:,2),hyloop_forward{i,1}(:,1))/1e6;
end

figure;
plot(step,energy(:,1),'ro');


eff=energy./energy2;
figure;
plot(step,eff,'ro');






