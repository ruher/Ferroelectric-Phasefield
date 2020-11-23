clear;clc;
scale=128;
L=0.125e-6;
kb=1.38e-23;
T=600;
d=(L/scale)^2;
GAMMA=4e4; deltat=1e-14; delta_tao=GAMMA*deltat;
flucstd=2*kb*T*GAMMA;
Px0=normrnd(0,flucstd,[scale,scale]);
Py0=normrnd(0,flucstd,[scale,scale]);
% Px0=flucstd*ones(scale);
% Py0=0*ones(scale);
%%------------alpha------------%%
Tc=492+273;
alpha1=3.6574e5*(T-Tc);
alpha11=-7.253e7;
alpha12=7.5e8;
alpha111=2.606e8;
alpha112=6.1e8;

epsilon0=8.854187e-12;epsilon11=93.2*epsilon0;
%%--------------g--------------%%
g11=5e-10;
g44=2.7e-11;
% g11=0;
% g44=0;
%%--------------kspace--------------%%
K=2*pi*scale/L;dk=K/scale;
[kx,ky]=meshgrid(-K/2:dk:(K/2-dk),(-K/2):dk:(K/2-dk));
kx=fftshift(kx);
ky=fftshift(ky);
kx_n=kx./sqrt(kx.^2+ky.^2);ky_n=ky./sqrt(kx.^2+ky.^2);
kx_n(1,1)=0;ky_n(1,1)=0;
%%------------electric mechanical coefficient and kernal------------%%
C11=1.796e11; C12=7.937e10; C44=1.111e11;
Q11=0.089; Q12=-0.026; Q44=0.0675;
q11=C11*Q11+2*C12*Q12; q12=C12*Q11+(C11+C12)*Q12; q44=2*C44*Q44;

q11_c=q11+2*q12; q22_c=q11-q12;
Xi=(C11-C12-2*C44)/C44;
dx=C44*(1+Xi*kx_n.^2);dy=C44*(1+Xi*ky_n.^2);
Chi=kx_n.^2./dx+ky_n.^2./dy;
D=(C12+C44)./(1+(C12+C44)*Chi);
%----------Beta----------%
beta11=kx_n.^2./dx; beta12=0; beta13=0; beta14=0; beta15=0; beta16=kx_n.*ky_n./dx;
beta21=0; beta22=ky_n.^2./dy; beta23=0; beta24=0; beta25=0; beta26=kx_n.*ky_n./dy;
beta31=0; beta32=0; beta33=0; beta34=0; beta35=0; beta36=0;
beta41=0; beta42=0; beta43=0; beta44=ky_n.^2/C44; beta45=kx_n.*ky_n/C44; beta46=0;
beta51=0; beta52=0; beta53=0; beta54=kx_n.*ky_n/C44; beta55=kx_n.^2/C44; beta56=0;
beta61=kx_n.*ky_n./dx; beta62=kx_n.*ky_n./dy; beta63=0; beta64=0; beta65=0; beta66=kx_n.^2./dy+ky_n.^2./dx;
%----------Theta----------%
theta1=kx_n.^2./dx; theta2=ky_n.^2./dy; theta3=0; theta4=0; theta5=0; theta6=kx_n.*ky_n.*(1./dx+1./dy);
%----------B----------%
for i=1:6
    for j=1:6
        eval(['B',num2str(i),num2str(j),'=beta',num2str(i),num2str(j),'-D.*theta',num2str(i),'.*','theta',num2str(j),';']);
    end
end
%----------A----------%
A11=q22_c^2*B11+2*q12*q22_c*(B11+B12+B13)+q12^2*(B11+B12+B13+B21+B22+B23+B31+B32+B33);
A12=q22_c^2*B12-q12*q22_c*(B13+B23+B33)+q11*q12*(B11+B12+B13+B21+B22+B23+B31+B32+B33);
A16=q22_c*q44*B16+q12*q44*(B16+B26+B36);
A21=A12;
A22=q22_c^2*B22+2*q12*q22_c*(B21+B22+B23)+q12^2*(B11+B12+B13+B21+B22+B23+B31+B32+B33);
A26=q22_c*q44*B26+q12*q44*(B16+B26+B36);
A61=A16; A62=A26; A66=q44^2*B66;
%%----------gradient term kernal----------%
grad_kernal1=g11*kx.^2+g44*ky.^2;
grad_kernal2=g11*ky.^2+g44*kx.^2;
%%----------electrostatic term kernal----------%
elesta_kernal=4*pi*(epsilon11*kx.^2+epsilon11*ky.^2);
%%--------------history--------------%%
iter_batch=800; iter_epoch=20;
iter_num=iter_epoch*iter_batch;
domain_history=cell(1,iter_epoch);
Px_history=zeros(1,iter_num);
[x,y]=meshgrid(1:scale,1:scale);
ela_reduc_para=0.25;


k_Px0=d*fft2(fftshift(Px0)); k_Py0=d*fft2(fftshift(Py0));
[dG_k_Px0,dG_k_Py0]=getPZTF_2D(Px0,Py0,scale,ela_reduc_para);
k_Px1=(k_Px0-delta_tao*dG_k_Px0)./(1+delta_tao*grad_kernal1);
k_Py1=(k_Py0-delta_tao*dG_k_Py0)./(1+delta_tao*grad_kernal2);
Px1=real(fftshift(ifft2(k_Px1/d))); Py1=real(fftshift(ifft2(k_Py1/d)));

[dG_k_Px1,dG_k_Py1]=getPZTF_2D(Px1,Py1,scale,ela_reduc_para);
k_Px2=(4*k_Px1-k_Px0-delta_tao*(3*dG_k_Px1-dG_k_Px0))./(3+2*delta_tao*grad_kernal1);
k_Py2=(4*k_Py1-k_Py0-delta_tao*(3*dG_k_Py1-dG_k_Py0))./(3+2*delta_tao*grad_kernal2);
Px2=real(fftshift(ifft2(k_Px2/d))); Py2=real(fftshift(ifft2(k_Py2/d)));
for i=1:iter_num
    %%-----------------derivative of Ginzburg-Landau energy-----------------%%
    dG_GL_k_Px2_alpha1=2.*alpha1.*Px2;
    dG_GL_k_Px2_alpha11=4.*alpha11.*Px2.^3;
    dG_GL_k_Px2_alpha12=2.*alpha12.*Px2.*Py2.^2;
    dG_GL_k_Px2_alpha111=6.*alpha111.*Px2.^5;
    dG_GL_k_Px2_alpha112=alpha112.*(4.*Px2.^3.*Py2.^2+2.*Px2.*Py2.^4);
    
    dG_GL_k_Py2_alpha1=2.*alpha1.*Py2;
    dG_GL_k_Py2_alpha11=4.*alpha11.*Py2.^3;
    dG_GL_k_Py2_alpha12=2.*alpha12.*Py2.*Px2.^2;
    dG_GL_k_Py2_alpha111=6.*alpha111.*Py2.^5;
    dG_GL_k_Py2_alpha112=alpha112.*(4.*Py2.^3.*Px2.^2+2.*Py2.*Px2.^4);

    dG_GL_k_Px2=d*fft2(ifftshift(dG_GL_k_Px2_alpha1+dG_GL_k_Px2_alpha11+dG_GL_k_Px2_alpha12+dG_GL_k_Px2_alpha111+dG_GL_k_Px2_alpha112));
    dG_GL_k_Py2=d*fft2(ifftshift(dG_GL_k_Py2_alpha1+dG_GL_k_Py2_alpha11+dG_GL_k_Py2_alpha12+dG_GL_k_Py2_alpha111+dG_GL_k_Py2_alpha112));
    
    %%----------------------------gradient energy----------------------------%%
    k_Px2=d*fft2(fftshift(Px2));
    k_Py2=d*fft2(fftshift(Py2));
    %%----------------------------gradient energy----------------------------%%
    %%-------------------derivative of electrostatic energy------------------%%
    temp1=kx.*k_Px2+ky.*k_Py2;
    dG_k_es_Px2=temp1.*kx./elesta_kernal;
    dG_k_es_Py2=temp1.*ky./elesta_kernal;
    dG_k_es_Px2(1,1)=0;
    dG_k_es_Py2(1,1)=0;
    %------------------derivative of electromechanic energy-----------------%%
    Y1=Px2.^2; Y2=Py2.^2; Y6=Px2.*Py2;
    Y1_k=d*fft2(fftshift(Y1)); Y2_k=d*fft2(fftshift(Y2)); Y6_k=d*fft2(fftshift(Y6));
    kernal_sum1=real(fftshift(ifft2((A11.*Y1_k+A12.*Y2_k+A16.*Y6_k)/d)));
    kernal_sum2=real(fftshift(ifft2((A21.*Y1_k+A22.*Y2_k+A26.*Y6_k)/d)));
    kernal_sum6=real(fftshift(ifft2((A61.*Y1_k+A62.*Y2_k+A66.*Y6_k)/d)));
    dG_k_em_Px2=-0.5*d*fft2(fftshift(4*Px2.*kernal_sum1+2*Py2.*kernal_sum6))*ela_reduc_para;
    dG_k_em_Py2=-0.5*d*fft2(fftshift(4*Py2.*kernal_sum2+2*Px2.*kernal_sum6))*ela_reduc_para;

    dG_k_em_Px2(1,1)=0;
    dG_k_em_Py2(1,1)=0;
    %--------------------------------evolvetion------------------------------%%
%     dG_k_Px2=dG_GL_k_Px2+dG_k_em_Px2;
%     dG_k_Py2=dG_GL_k_Py2+dG_k_em_Py2;
%     dG_k_Px2=dG_GL_k_Px2+dG_k_es_Px2;
%     dG_k_Py2=dG_GL_k_Py2+dG_k_es_Py2;
    dG_k_Px2=dG_GL_k_Px2+dG_k_es_Px2+dG_k_em_Px2;
    dG_k_Py2=dG_GL_k_Py2+dG_k_es_Py2+dG_k_em_Py2;

    k_Px3=(18*k_Px2-9*k_Px1+2*k_Px0-delta_tao*(11.5*dG_k_Px2-8*dG_k_Px1+2.5*dG_k_Px0))./(11+6*delta_tao*grad_kernal1);

    k_Py3=(18*k_Py2-9*k_Py1+2*k_Py0-delta_tao*(11.5*dG_k_Py2-8*dG_k_Py1+2.5*dG_k_Py0))./(11+6*delta_tao*grad_kernal2);
    Px3=real(fftshift(ifft2(k_Px3/d))); Py3=real(fftshift(ifft2(k_Py3/d)));
    
    if rem(i,iter_batch)==0
        domain_history{i/iter_batch}=Px2;
    end
    
    Px_history(i)=Px3(100,100);

    disp(Px3(100,100));
    
    k_Px0=k_Px1; k_Py0=k_Py1;
    k_Px1=k_Px2; k_Py1=k_Py2;
    k_Px2=k_Px3; k_Py2=k_Py3;
    
    Px2=Px3; Py2=Py3;
    
    dG_k_Px0=dG_k_Px1; dG_k_Py0=dG_k_Py1;
    dG_k_Px1=dG_k_Px2; dG_k_Py1=dG_k_Py2;
end
for p=1:iter_epoch
    subplot(4,5,p);
    imagesc(domain_history{p});
    axis off;
end

figure;
quiver(x,y,Px2,Py2);
mean2(sqrt(Px2.^2+Py2.^2))