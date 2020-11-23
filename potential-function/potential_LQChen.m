clear;clc;close all;
scale=0.3;
d=0.01;
[Px,Py]=meshgrid(-scale:d:scale,scale:-d:-scale);
[n,m]=size(Px);
% T1=[280,320,360,400];
% for T=T1

Tc=115+273;
T=350;
alpha1=4.124*1e5*(T-Tc);
alpha11=-20.97*1e7;alpha12=7.974e8;
alpha111=1.294*1e9; alpha112=-1.950e9;
alpha1111=3.863*1e10; alpha1112=2.529e10; alpha1122=1.637e10;
Q11=0.1104; Q12=-0.0452; Q44=0.0289;

E1=0;
E2=0;

Ex=E1*ones(n);
Ey=E2*ones(n);

GL_alpha1=alpha1*Px.^2+alpha1*Py.^2;
GL_alphabar11=alpha11*Px.^4+alpha11*Py.^4;
GL_alphabar12=alpha12*Px.^2.*Py.^2;
GL_alpha111=alpha111*Px.^6+alpha111*Py.^6;
GL_alpha112=alpha112*(Px.^4.*Py.^2+Px.^2.*Py.^4);
GL_alpha1111=alpha1111*Px.^8+alpha1111*Py.^8;
GL_alpha1112=alpha1112*(Px.^6.*Py.^2+Px.^2.*Py.^6);
GL_alpha1122=alpha1122*Px.^4.*Py.^4;

e1=1e8;
e2=2e8;
e3=1e9;

G_GL=GL_alpha1+GL_alphabar11+GL_alphabar12+GL_alpha111+GL_alpha112+GL_alpha1111+GL_alpha1112+GL_alpha1122-Ex.*Px-Ey.*Py ...
-e1*(Q11*Px.^2+Q12*Py.^2)-e2*(Q12*Px.^2+Q11*Py.^2)-e3*Q44*Px.*Py;
% plot(Px((n+1)/2,:),G_GL((n+1)/2,:));
% hold on;
% end
% legend;
% plot(Px(4,:),G_GL(4,:));
% surf(G_GL);
% axis off;
% figure;
contour(Px,Py,G_GL,25);
% axis off;

% dG_GL_k_Px_alpha1=2.*alpha1.*Px;
% dG_GL_k_Px_alpha11=4.*alpha11.*Px.^3;
% dG_GL_k_Px_alpha12=2.*alpha12.*Px.*Py.^2;
% dG_GL_k_Px_alpha111=6.*alpha111.*Px.^5;
% dG_GL_k_Px_alpha112=alpha112.*(4.*Px.^3.*Py.^2+2.*Px.*Py.^4);
% dG_GL_k_Px_alpha1111=alpha1111.*8.*Px.^7;
% dG_GL_k_Px_alpha1112=alpha1112.*(6.*Px.^5.*Py.^2+2.*Px.*Py.^6);
% dG_GL_k_Px_alpha1122=alpha1122.*4.*Px.^3.*Py.^4;
% dG_GL_k_Px=dG_GL_k_Px_alpha1+dG_GL_k_Px_alpha11+dG_GL_k_Px_alpha12+dG_GL_k_Px_alpha111+dG_GL_k_Px_alpha112+dG_GL_k_Px_alpha1111+dG_GL_k_Px_alpha1112+dG_GL_k_Px_alpha1122;
% 
% dG_GL_k_Py_alpha1=2.*alpha1.*Py;
% dG_GL_k_Py_alpha11=4.*alpha11.*Py.^3;
% dG_GL_k_Py_alpha12=2.*alpha12.*Py.*Px.^2;
% dG_GL_k_Py_alpha111=6.*alpha111.*Py.^5;
% dG_GL_k_Py_alpha112=alpha112.*(4.*Py.^3.*Px.^2+2.*Py.*Px.^4);
% dG_GL_k_Py_alpha1111=alpha1111.*8.*Py.^7;
% dG_GL_k_Py_alpha1112=alpha1112.*(6.*Py.^5.*Px.^2+2.*Py.*Px.^6);
% dG_GL_k_Py_alpha1122=alpha1122.*4.*Py.^3.*Px.^4;
% dG_GL_k_Py=dG_GL_k_Py_alpha1+dG_GL_k_Py_alpha11+dG_GL_k_Py_alpha12+dG_GL_k_Py_alpha111+dG_GL_k_Py_alpha112+dG_GL_k_Py_alpha1111+dG_GL_k_Py_alpha1112+dG_GL_k_Py_alpha1122;
% 
% figure;
% quiver(Px,Py,dG_GL_k_Px,dG_GL_k_Py,2);
