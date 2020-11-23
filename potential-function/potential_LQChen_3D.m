clear;clc;close all;
L=0.2;
d=0.01;
[Px,Py,Pz]=meshgrid(-L:d:L,L:-d:-L,L:-d:-L);
[n,m,l]=size(Px);

T=403;

alpha12=7.974e8;alpha11=-2.097e8;
alpha111=1.294e9; alpha112=-1.950e9; alpha123=-2.5e9;
alpha1111=3.863e10; alpha1112=2.529e10; alpha1122=1.637e10; alpha1123=1.367e10;
E1=0;
E2=0;
E3=0;
Ex=E1*ones(n);
Ey=E2*ones(n);
Ez=E3*ones(n);

for i=4:9
alpha1=(i+0.124)*1e5*(T-388);
GL_alpha1=alpha1*(Px.^2+Py.^2+Pz.^2);
GL_alphabar11=alpha11*(Px.^4+Py.^4+Pz.^4);
GL_alphabar12=alpha12*(Px.^2.*Py.^2+Py.^2.*Pz.^2+Px.^2.*Pz.^2);
GL_alpha111=alpha111*(Px.^6+Py.^6+Pz.^6);
GL_alpha112=alpha112*(Px.^2.*(Py.^4+Pz.^4)+Py.^2.*(Px.^4+Pz.^4)+Pz.^2.*(Py.^4+Px.^4));
GL_alpha123=alpha123*Px.^2.*Py.^2.*Pz.^2;
GL_alpha1111=alpha1111*(Px.^8+Py.^8+Pz.^8);
GL_alpha1112=alpha1112*(Px.^6.*(Py.^2+Pz.^2)+Py.^6.*(Px.^2+Pz.^2)+Pz.^6.*(Py.^2+Px.^2));
GL_alpha1122=alpha1122*(Px.^4.*Py.^4+Pz.^4.*Py.^4+Px.^4.*Pz.^4);
GL_alpha1123=(Px.^4.*Py.^2.*Pz.^2+Py.^4.*Pz.^2.*Px.^2+Pz.^4.*Px.^2.*Py.^2);
G_GL=GL_alpha1+GL_alphabar11+GL_alphabar12+GL_alpha111+GL_alpha112+GL_alpha123+GL_alpha1111+GL_alpha1112+GL_alpha1122+GL_alpha1123-Ex.*Px-Ey.*Py;
plot(Px((n+1)/2,:,(n+1)/2),G_GL((n+1)/2,:,(n+1)/2));
hold on;
end
legend('4','5','6','7','8','9');
% hold on;
% plot(Py(:,(n+1)/2,(n+1)/2),G_GL(:,(n+1)/2,(n+1)/2));
% hold on;
% plot(reshape(Pz((n+1)/2,(n+1)/2,:),n,1),reshape(G_GL((n+1)/2,(n+1)/2,:),n,1));
