clear;clc;
% scale=128;
% Px=0.26*ones(scale,scale);Py=zeros(scale,scale);
load('domain_xy','Px','Py');


[x,y]=meshgrid(1:128,128:-1:1);


Px=flipud(Px);Py=flipud(Py);
% Px(abs(Px)<0.25)=0;Px(abs(Px)>0.25)=-0.26;
% Py(Px==0)=0.26;Py(Py~=0.26)=0;

Px_pad=Padding(Px);Py_pad=Padding(Py);



L=0.125e-6; scale=128; h=L/scale; d=h^2;
epsilon0=8.854187e-12;epsilon=2000.*epsilon0;
% C0=6.02e3.*0.005/0.2332; c0=C0*ones(scale,scale); c=c0;%concentration
% Z=2;e=1.6e-19;


d_Px_x=Px_pad(:,3:(scale+2))-Px_pad(:,1:scale);
d_Px_x=d_Px_x(2:(scale+1),:)/(2*h);

d_Py_y=Py_pad(1:scale,:)-Py_pad(3:(scale+2),:);
d_Py_y=d_Py_y(:,2:(scale+1))/(2*h);

div_P=d_Px_x+d_Py_y;


figure;
subplot(2,2,1);
imagesc(Px_pad);
subplot(2,2,2);
surf(d_Px_x);
subplot(2,2,3);
imagesc(Py_pad);
subplot(2,2,4);
surf(d_Py_y);
figure;
subplot(1,2,1);
quiver(x,y,Px,Py);
subplot(1,2,2);
surf(div_P);

phi_dcre=zeros(scale^2,scale^2);
for i=1:scale
    for j=1:scale
        p=(i-1)*scale+j;
        phi_dcre(p,(f(i,scale)-1)*scale+j)=1;
        phi_dcre(p,(b(i,scale)-1)*scale+j)=1;
        phi_dcre(p,(i-1)*scale+f(j,scale))=1;
        phi_dcre(p,(i-1)*scale+b(j,scale))=1;
        phi_dcre(p,(i-1)*scale+j)=-4;
    end
end

laplace_phi=div_P/epsilon;
% laplace_phi=-(e*Z*c-div_P)/epsilon;
f_unfold=laplace_phi';
f_unfold=f_unfold(:);
b=d*f_unfold(:);
phi=phi_dcre\b;
phi=reshape(phi,scale,scale);
phi=phi';
surf(phi);
% % % % 
phi_pad=Padding(phi);
d_phi_x=phi_pad(:,3:(scale+2))-phi_pad(:,1:scale);
Ex=-d_phi_x(2:(scale+1),:)/(2*h);
d_phi_y=phi_pad(1:scale,:)-phi_pad(3:(scale+2),:);
Ey=-d_phi_y(:,2:(scale+1))/(2*h);
% % % % 
% % % dd_phi_x=phi_pad(:,3:(scale+2))+phi_pad(:,1:scale)-2*phi_pad(:,2:(scale+1));
% % % dd_phi_x=dd_phi_x(2:(scale+1),:)/d;
% % % dd_phi_y=phi_pad(1:scale,:)+phi_pad(3:(scale+2),:)-2*phi_pad(2:(scale+1),:);
% % % dd_phi_y=dd_phi_y(:,2:(scale+1))/d;
% % % laplace_phi=dd_phi_x+dd_phi_y;
% % % % figure;
% % % % surf(laplace_phi*epsilon);
% % % % 
figure;
quiver(x,y,Ex,Ey);
% figure;
% quiver(x,y,Px,Py);
% % 
