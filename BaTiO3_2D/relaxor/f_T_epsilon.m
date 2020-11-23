clc;clear;close all;
filen=2;
P=cell(filen,1);
load('f500_T_relaxorre2.mat');
P{1,1}=Px_mean;
load('f_T_relaxorre2.mat');
P{2,1}=Px_mean;
Emax=1e6;
epsilon1=zeros(length(T_series),filen);
delta=zeros(length(T_series),filen);
sampling=100;
t=(2*pi/sampling:(2*pi/sampling):2*pi);
for i=1:filen
    for j=1:length(T_series)
        epsilon1(j,i)=(max(P{i,1}(j,:))-min(P{i,1}(j,:)))/(2*Emax);
        [fitresult,error]=sinfit(t,P{i,1}(j,:));
        delta(j,i)=-fitresult.c;
    end
end
epsilon2=tan(delta).*epsilon1;
for i=1:filen
subplot(1,3,1);
plot(T_series,epsilon1(:,i));
hold on;
subplot(1,3,2);
plot(T_series,epsilon2(:,i));
hold on;
subplot(1,3,3);
plot(T_series,tan(delta(:,i)));
hold on;
end
legend;
% T_maxepsilon1=zeros(filen,1);
% T_maxepsilon2=zeros(filen,1);
% for i=1:filen
%     [max_index1,temp]=find(epsilon1==max(epsilon1(:,i)));
%     T_maxepsilon1(i)=T_series(max_index1);
%     [max_index2,temp]=find(epsilon2==max(epsilon2(:,i)));
%     T_maxepsilon2(i)=T_series(max_index2);
% end

% figure;
% plot(t,Px_mean(321,:)*30-0.2);
% hold on;
% plot(t,acE/1e7);

% figure;
% plot(delta);