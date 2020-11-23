clear;clc;
% Px=0.26*ones(scale,scale);Py=zeros(scale,scale);
load('domain_xy','Px','Py');
Px=flipud(Px);Py=flipud(Py);
load('phi_dcre_matrix.mat');
L=0.125e-6; scale=128; h=L/scale; d=h^2; [x,y]=meshgrid(-L/2:h:L/2-h,-L/2:h:L/2-h);%meshgrid
C0=6.02e3.*0.005/0.2332; c0=normrnd(C0,30,scale,scale); c=c0;%concentration
Z=2;e=1.6e-19;%elementary charge
k=1.3806e-23;%Boltzman contant
T=45+273;%temperature
mu=8.4e-22; D=mu.*k.*T;%diffusion constant
epsilon0=8.854187e-12;epsilon=2000.*epsilon0;%dielectric constant
%iteration item
deltat=2e22;% delta t
iter_epoch=20; iter_batch=1000;
iter_num=iter_epoch*iter_batch+1;
%record
compositionfield=cell(1,20);
c_history=zeros(1,iter_num);
for i=1:iter_num
    c_pad=Padding(c);
    d_c_x=c_pad(:,3:(scale+2))-c_pad(:,1:scale);
    d_c_x=d_c_x(2:(scale+1),:)/(2*h);
    d_c_y=c_pad(1:scale,:)-c_pad(3:(scale+2),:);
    d_c_y=d_c_y(:,2:(scale+1))/(2*h);
    
    dd_c_x=c_pad(:,3:(scale+2))+c_pad(:,1:scale)-2*c_pad(:,2:(scale+1));
    dd_c_x=dd_c_x(2:(scale+1),:)/d;
    dd_c_y=c_pad(3:(scale+2),:)+c_pad(1:scale,:)-2*c_pad(2:(scale+1),:);
    dd_c_y=dd_c_y(:,2:(scale+1))/d;
    
    
    Px_pad=Padding(Px);Py_pad=Padding(Py);
    d_Px_x=Px_pad(:,3:(scale+2))-Px_pad(:,1:scale);
    d_Px_x=d_Px_x(2:(scale+1),:)/(2*h);

    d_Py_y=Py_pad(1:scale,:)-Py_pad(3:(scale+2),:);
    d_Py_y=d_Py_y(:,2:(scale+1))/(2*h);

    div_P=d_Px_x+d_Py_y;
    
    laplace_phi=-(e*Z*c-div_P)/epsilon;
%     
    if rem(i-1,2000)==0
        f_unfold=laplace_phi';
        f_unfold=f_unfold(:);
        b=d*f_unfold(:);
        phi=phi_dcre\b;
        phi=reshape(phi,scale,scale);
        phi=phi';
    end
    
    
    phi_pad=Padding(phi);
    d_phi_x=phi_pad(:,3:(scale+2))-phi_pad(:,1:scale);
    d_phi_x=d_phi_x(2:(scale+1),:)/(2*h);
    d_phi_y=phi_pad(1:scale,:)-phi_pad(3:(scale+2),:);
    d_phi_y=d_phi_y(:,2:(scale+1))/(2*h);

    d_item=deltat*(D*(dd_c_x+dd_c_y)+e*Z*mu*(c.*laplace_phi+d_phi_x.*d_c_x+d_phi_y.*d_c_y));
    disp(c(10,10));
%     d_item=deltat*(D*(dd_c_x+dd_c_y));
    c=c+d_item;
    c_history(i)=c(10,10);
    if rem(i,iter_batch)==0
       compositionfield{i/iter_batch}=c;
    end
end
for p=1:iter_epoch
    subplot(4,5,p);
    imagesc(compositionfield{p});
end