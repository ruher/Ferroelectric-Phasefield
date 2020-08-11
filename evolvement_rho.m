clear;clc;
load('298K_ela025.mat','Px2','Py2');
Px=Px2; Py=Py2;
scale=size(Px);
scale=scale(1);
deltat=8e-18;
iter_num=2000;
rho=4.1e7;
alpha2=1; alpha1=-2*alpha2*rho^2; alpha3=2*alpha2; mu=1e4;
kb=1.38e-23;
T=298;
GAMMA=4e4;
rhox=zeros(scale); rhoy=zeros(scale);
rhohistory=zeros(iter_num,1);
interfield_history=cell(8,2);
n=0;
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
    if i==200
        n=n+1;
        interfield_history{n,1}=rhox;
        interfield_history{n,2}=rhoy;
        steps(n)=i;
    end
%     
    if i==500
        n=n+1;
        interfield_history{n,1}=rhox;
        interfield_history{n,2}=rhoy;
        steps(n)=i;
    end
%     
    if rhointer_i0<1e7&&rhointer_i1>1e7
        n=n+1;
        interfield_history{n,1}=rhox;
        interfield_history{n,2}=rhoy;
        steps(n)=i;
    end
%     
    if rhointer_i0<2e7&&rhointer_i1>2e7
        n=n+1;
        interfield_history{n,1}=rhox;
        interfield_history{n,2}=rhoy;
        steps(n)=i;
    end
%     
    if rhointer_i0<3e7&&rhointer_i1>3e7
        n=n+1;
        interfield_history{n,1}=rhox;
        interfield_history{n,2}=rhoy;
        steps(n)=i;
    end

%     
    if rhointer_i0<4e7&&rhointer_i1>4e7
        n=n+1;
        interfield_history{n,1}=rhox;
        interfield_history{n,2}=rhoy;
        steps(n)=i;
    end
%     
    if i==1200
        n=n+1;
        interfield_history{n,1}=rhox;
        interfield_history{n,2}=rhoy;
        steps(n)=i;
    end
%     
    if i==1500
        n=n+1;
        interfield_history{n,1}=rhox;
        interfield_history{n,2}=rhoy;
        steps(n)=i;
    end
   
end
[x,y]=meshgrid(1:128,1:128);
plot(rhohistory);