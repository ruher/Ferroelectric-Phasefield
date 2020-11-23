clear;clc;
load('ferro373K.mat','Px2','Py2');
Px=Px2; Py=Py2;
scale=size(Px);
scale=scale(1);
deltat=2e-18;
iter_num=1500;
rho=5e6;
alpha2=1; alpha1=-2*alpha2*rho^2; alpha3=2.5*alpha2; mu=5e22;
kb=1.38e-23;
T=298;
GAMMA=4e4;
rhox=zeros(scale); rhoy=zeros(scale);
rhohistory=zeros(iter_num,1);
interfield_history=cell(8,3);
n=1;
steps=zeros(4,1);
for i=1:iter_num
    
    rhointer_i0=mean2(sqrt(rhox.^2+rhoy.^2));
    drhox=2*alpha1*rhox+4*alpha2*rhox.^3+2*alpha3*rhox.*rhoy.^2;
    drhoy=2*alpha1*rhoy+4*alpha2*rhoy.^3+2*alpha3*rhoy.*rhox.^2;
    intereactionx=-mu*Px;
    intereactiony=-mu*Py;
    ditem_rhox=-deltat*(drhox+intereactionx);
    ditem_rhoy=-deltat*(drhoy+intereactiony);

    rhox=rhox+ditem_rhox;
    rhoy=rhoy+ditem_rhoy;
    
    rhointer_i1=mean2(sqrt(rhox.^2+rhoy.^2));
    rhohistory(i)=rhointer_i1;
    
    disp(rhointer_i1);
    
    if rhointer_i0<n*1e6&&rhointer_i1>n*1e6
        interfield_history{n,1}=rhox;
        interfield_history{n,2}=rhoy;
        interfield_history{n,3}=i;
        n=n+1;
    end

   
end
% [x,y]=meshgrid(1:128,1:128);
% quiver(x,y,rhox,rhoy);
plot(rhohistory);