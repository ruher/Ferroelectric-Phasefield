clear;clc;close all;
scale=0.3;
d=0.01;
[Px,Py]=meshgrid(-scale:d:scale,scale:-d:-scale);
[n,m]=size(Px);


T=260;
alpha1=3.34e5*(T-108-273);
alpha11=4.69e6*(T-120-273)-2.02e8;
alpha12=3.23e8;
alpha111=-5.52e7*(T-120-273)+2.76e9;
alpha112=4.47e9;
Q11=0.1104; Q12=-0.0452;

E1=0;
E2=0;

Ex=E1*ones(n);
Ey=E2*ones(n);

GL_alpha1=alpha1*Px.^2+alpha1*Py.^2;
GL_alphabar11=alpha11*Px.^4+alpha11*Py.^4;
GL_alphabar12=alpha12*Px.^2.*Py.^2;
GL_alpha111=alpha111*Px.^6+alpha111*Py.^6;
GL_alpha112=alpha112*(Px.^4.*Py.^2+Px.^2.*Py.^4);


e1=0;

G_GL=GL_alpha1+GL_alphabar11+GL_alphabar12+GL_alpha111+GL_alpha112-Ex.*Px-Ey.*Py-e1*Q11*Px.^2;
% plot(Px((n+1)/2,:),G_GL((n+1)/2,:));
% plot(Px(4,:),G_GL(4,:));
surf(G_GL);