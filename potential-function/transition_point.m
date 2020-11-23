clear;clc;close all;
scale=0.3;
d=0.001;
[Px,Py]=meshgrid(-scale:d:scale,scale:-d:-scale);
[n,m]=size(Px);

Tc=115+273;
alpha11=-20.97*1e7;alpha12=7.974e8;
alpha111=1.294*1e9; alpha112=-1.950e9;
alpha1111=3.863*1e10; alpha1112=2.529e10; alpha1122=1.637e10;
E1=0;
E2=0;
Ex=E1*ones(n);
Ey=E2*ones(n);

for i=1:200
    T=200+i;
    alpha1=4.124*1e5*(T-Tc);
    GL_alpha1=alpha1*Px.^2+alpha1*Py.^2;
    GL_alphabar11=alpha11*Px.^4+alpha11*Py.^4;
    GL_alphabar12=alpha12*Px.^2.*Py.^2;
    GL_alpha111=alpha111*Px.^6+alpha111*Py.^6;
    GL_alpha112=alpha112*(Px.^4.*Py.^2+Px.^2.*Py.^4);
    GL_alpha1111=alpha1111*Px.^8+alpha1111*Py.^8;
    GL_alpha1112=alpha1112*(Px.^6.*Py.^2+Px.^2.*Py.^6);
    GL_alpha1122=alpha1122*Px.^4.*Py.^4;
    G_GL=GL_alpha1+GL_alphabar11+GL_alphabar12+GL_alpha111+GL_alpha112+GL_alpha1111+GL_alpha1112+GL_alpha1122-Ex.*Px-Ey.*Py;
    G_min=min(min(G_GL));
    a=find(G_GL==G_min);
    y=mod(a(1),n);
    x=floor(a(1)/n);
    P(i)=sqrt(Px(x,y)^2+Py(x,y)^2);
end
T=(201:400);
plot(T,P);
