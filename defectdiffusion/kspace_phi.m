clear;clc;
load('domain_298K','Px','Py');

% Px(abs(Px)<0.26)=0;
% Py(Px~=0)=0;
% Py(Py~=0)=-0.2602;

L=0.125e-6; scale=128; h=L/scale; d=h^2;
epsilon0=8.854187e-12;epsilon1=3600*epsilon0;epsilon3=188*epsilon0;

K=2*pi*scale/L;dk=K/scale;
[kx,ky]=meshgrid(-K/2:dk:(K/2-dk),(-K/2):dk:(K/2-dk));
kx=fftshift(kx);
ky=fftshift(ky);

k_Px=d*fft2(fftshift(Px));k_Py=d*fft2(fftshift(Py));

phi_k=-1i*(kx.*k_Px+ky.*k_Py)./(epsilon1*kx.^2+epsilon3*ky.^2);
phi_k(1,1)=0;
phi=real(fftshift(ifft2(phi_k/d)));
% phi=real(ifft2(phi_k/d));
                                                                                                                                                                                               

grad_phi_x_k=1i.*kx.*phi_k;
grad_phi_y_k=1i.*ky.*phi_k;
laplace_phi_k=-(kx.^2+ky.^2).*phi_k;
laplace_phi=real(fftshift(ifft2(laplace_phi_k/d)));
d_phi_x=real(fftshift(ifft2(grad_phi_x_k)/d));
Ex=-d_phi_x;
d_phi_y=real(fftshift(ifft2(grad_phi_y_k)/d));
Ey=-d_phi_y;

[x,y]=meshgrid(-L/2:h:L/2-h,-L/2:h:L/2-h);


% figure;
% quiver(x,y,Px,Py);
% figure;
% quiver(x,y,Ex,Ey);
% figure;
% subplot(1,2,1);
% imagesc(Ex);
% subplot(1,2,2);
% imagesc(Ey);
figure;
surf(phi);


