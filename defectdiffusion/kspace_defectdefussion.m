clear;clc;
load('domain_298K','Px','Py');
L=0.0625e-6; scale=128; h=L/scale; d=h^2; [x,y]=meshgrid(-L/2:h:L/2-h,-L/2:h:L/2-h);%meshgrid

K=2*pi*scale/L;dk=K/scale;
[kx,ky]=meshgrid(-K/2:dk:(K/2-dk),(-K/2):dk:(K/2-dk));
kx=fftshift(kx);
ky=fftshift(ky);
k_Px=d*fft2(ifftshift(Px));k_Py=d*fft2(ifftshift(Py));

Oxygen_messper=0.005;
N_A=6.02e23;
C0=1e6*233.1922*Oxygen_messper/(233.1922*16/6.017);
c_O=normrnd(C0,4e2,scale,scale); c_K=2*c_O;%concentration
Z=2;e=1.6e-19;%elementary charge
k=1.3806e-23;%Boltzman contant
T=25+273;%temperature
mu=1.73e-20; D=mu*k*T/(Z*e);%diffusion constant
epsilon0=8.854187e-12;epsilon1=3600*epsilon0;epsilon3=188*epsilon0;%dielectric constant
L_d=(2000*epsilon0*k*T/(4*e^2*C0*N_A))^0.5;
t_d=2000*epsilon0/(e*C0*N_A*mu);
%iteration item
deltat=5e-1;% delta t
iter_epoch=20; iter_batch=1000;
iter_num=iter_epoch*iter_batch+1;
%record
compositionfield=cell(1,20);
c_history=zeros(1,iter_num);
c_k_K=d*fft2(ifftshift(c_K));
for i=1:iter_num
    c_O_k=d*fft2(ifftshift(c_O));
    %phi
    phi_k=((-e*Z*c_O_k+e*c_k_K)*N_A-1i*(kx.*k_Px+ky.*k_Py))./(epsilon1*kx.^2+epsilon3*ky.^2);
%     phi_k=(-e*Z*c_O_k*N_A-1i*(kx.*k_Px+ky.*k_Py))./(epsilon1*kx.^2+epsilon3*ky.^2);
    phi_k(1,1)=0;
    phi=real(fftshift(ifft2(phi_k/d)));

    grad_phi_x_k=1i.*kx.*phi_k;
    grad_phi_y_k=1i.*ky.*phi_k;
    laplace_phi_k=-(kx.^2+ky.^2).*phi_k;
    d_phi_x=real(fftshift(ifft2(grad_phi_x_k)/d));
    d_phi_y=real(fftshift(ifft2(grad_phi_y_k)/d));
    laplace_phi=real(fftshift(ifft2(laplace_phi_k/d)));
    
    d_c_x_k=1i*kx.*c_O_k;
    d_c_y_k=1i*ky.*c_O_k;
    dd_c_x_k=-kx.^2.*c_O_k;
    dd_c_y_k=-ky.^2.*c_O_k;
    d_c_x=real(fftshift(ifft2(d_c_x_k)/d));
    d_c_y=real(fftshift(ifft2(d_c_y_k)/d));
    dd_c_x=real(fftshift(ifft2(dd_c_x_k)/d));
    dd_c_y=real(fftshift(ifft2(dd_c_y_k)/d));
    
    d_item=deltat*(D*(dd_c_x+dd_c_y)+mu*(-c_O.*laplace_phi-(d_phi_x.*d_c_x+d_phi_y.*d_c_y)));
    disp(c_O(50,50));
%     d_item=deltat*(D*(dd_c_x+dd_c_y));
    c_O=c_O+d_item;
    c_history(i)=c_O(10,10);
    if rem(i,iter_batch)==0
       compositionfield{i/iter_batch}=c_O;
    end
end
% set(gcf,'Position',[0, 0, 1300, 1000]);
for p=1:iter_epoch
    subplot(4,5,p);
    imagesc(compositionfield{p});
    axis off;
end
figure;
% quiver(x,y,Px,Py);
% figure;
quiver(x,y,-d_phi_x,-d_phi_y);
surf(phi);