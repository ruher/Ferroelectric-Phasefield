clear;clc;
scale=128;
scalex=scale;scaley=scale;scalez=scale/8;
Lx=0.125e-6;Ly=Lx;Lz=Lx/8;
kb=1.38e-23;
T=298;
d=(Lx/scalex)^3;
GAMMA=4e4; deltat=1e-14; delta_tao=GAMMA*deltat;
Px0=normrnd(0,2*kb*T*GAMMA,[scalex,scaley,scalez]);
Py0=normrnd(0,2*kb*T*GAMMA,[scalex,scaley,scalez]);
Pz0=normrnd(0,2*kb*T*GAMMA,[scalex,scaley,scalez]);
%%------------alpha------------%%
alpha1=4.124e5*(T-388); 
alpha11=-2.097e8; alpha12=7.974e8;
alpha111=1.294e9; alpha112=-1.950e9; alpha123=-2.5e9;
alpha1111=3.863e10; alpha1112=2.529e10; alpha1122=1.637e10; alpha1123=1.367e10;
epsilon0=8.854187e-12; epsilon11=7.35*epsilon0; epsilon33=7.35*epsilon0;
%%--------------g--------------%%
g11=5.1e-10;g44=2e-11;
%%--------------kspace--------------%%
Kx=2*pi*scalex/Lx;dkx=Kx/scalex;
Ky=2*pi*scaley/Ly;dky=Ky/scaley;
Kz=2*pi*scalez/Lz;dkz=Kz/scalez;
[kx,ky,kz]=meshgrid(-Kx/2:dkx:(Kx/2-dkx),(-Ky/2):dky:(Ky/2-dky),(-Kz/2):dkz:(Kz/2-dkz));
kx=fftshift(kx);
ky=fftshift(ky);
kz=fftshift(kz);
kx_n=kx./sqrt((kx.^2+ky.^2+kz.^2));ky_n=ky./sqrt((kx.^2+ky.^2+kz.^2));kz_n=kz./sqrt((kx.^2+ky.^2+kz.^2));
kx_n(1,1,1)=0;ky_n(1,1,1)=0;kz_n(1,1,1)=0;
%%------------electric mechanical coefficient and kernal------------%%
C11=27.5e10; C12=17.9e10; C44=5.43e10;
q11=14.20e9; q12=-0.74e9; q44=1.57e9;
q11_c=q11+2*q12; q22_c=q11-q12;
Xi=(C11-C12-2*C44)/C44;
dx=C44*(1+Xi*kx_n.^2);dy=C44*(1+Xi*ky_n.^2);dz=C44*(1+Xi*kz_n.^2);
Chi=kx_n.^2./dx+ky_n.^2./dy+kz_n.^2./dz;
D=(C12+C44)./(1+(C12+C44)*Chi);
%----------Beta----------%
beta11=kx_n.^2./dx; beta12=0; beta13=0; beta14=0; beta15=kx_n.*kz_n./dx; beta16=kx_n.*ky_n./dx;
beta21=0; beta22=ky_n.^2./dy; beta23=0; beta24=ky_n.*kz_n./dy; beta25=0; beta26=kx_n.*ky_n./dy;
beta31=0; beta32=0; beta33=kz_n.^2./dz; beta34=ky_n.*kz_n./dz; beta35=kx_n.*kz_n./dz; beta36=0;
beta41=0; beta42=ky_n.*kz_n./dy; beta43=ky_n.*kz_n./dz; beta44=ky_n.^2./dz+kz_n.^2./dy; beta45=kx_n.*ky_n./dz; beta46=kx_n.*kz_n./dy;
beta51=kx_n.*kz_n./dx; beta52=0; beta53=kx_n.*kz_n./dz; beta54=kx_n.*ky_n./dz; beta55=kx_n.^2./dz+kz_n.^2./dx; beta56=ky_n.*kz_n./dx;
beta61=kx_n.*ky_n./dx; beta62=kx_n.*ky_n./dy; beta63=0; beta64=kx_n.*kz_n./dy; beta65=ky_n.*kz_n./dx; beta66=kx_n.^2./dy+ky_n.^2./dx;
%----------Theta----------%
theta1=kx_n.^2./dx; theta2=ky_n.^2./dy; theta3=kz_n.^2./dz; theta4=ky_n.*kz_n.*(1./dy+1./dz); theta5=kx_n.*kz_n.*(1./dx+1./dz); theta6=kx_n.*ky_n.*(1./dx+1./dy);
%----------B----------%
for i=1:6
    for j=1:6
        eval(['B',num2str(i),num2str(j),'=beta',num2str(i),num2str(j),'-D.*theta',num2str(i),'.*','theta',num2str(j),';']);
    end
end
%----------A----------%
A11=q22_c^2*B11+2*q12*q22_c*(B11+B12+B13)+q12^2*(B11+B12+B13+B21+B22+B23+B31+B32+B33);
A12=q22_c^2*B12-q12*q22_c*(B13+B23+B33)+q11*q12*(B11+B12+B13+B21+B22+B23+B31+B32+B33);
A13=q22_c^2*B13-q12*q22_c*(B12+B22+B32)+q11*q12*(B11+B12+B13+B21+B22+B23+B31+B32+B33);
A14=q22_c*q44*B14+q12*q44*(B14+B24+B34);
A15=q22_c*q44*B15+q12*q44*(B15+B25+B35);
A16=q22_c*q44*B16+q12*q44*(B16+B26+B36);
A21=A12;
A22=q22_c^2*B22+2*q12*q22_c*(B21+B22+B23)+q12^2*(B11+B12+B13+B21+B22+B23+B31+B32+B33);
A23=q22_c^2*B23-q12*q22_c*(B11+B21+B31)+q11*q12*(B11+B12+B13+B21+B22+B23+B31+B32+B33);
A24=q22_c*q44*B24+q12*q44*(B14+B24+B34);
A25=q22_c*q44*B25+q12*q44*(B15+B25+B35);
A26=q22_c*q44*B26+q12*q44*(B16+B26+B36);
A31=A13;
A32=A23;
A33=q22_c^2*B33+2*q12*q22_c*(B31+B32+B33)+q12^2*(B11+B12+B13+B21+B22+B23+B31+B32+B33);
A34=q22_c*q44*B34+q12*q44*(B14+B24+B34);
A35=q22_c*q44*B35+q12*q44*(B15+B25+B35);
A36=q22_c*q44*B36+q12*q44*(B16+B26+B36);
A41=A14; A42=A24; A43=A34; A44=q44^2*B44; A45=q44^2*B45; A46=q44^2*B46;
A51=A15; A52=A25; A53=A35; A54=q44^2*B54; A55=q44^2*B55; A56=q44^2*B56;
A61=A16; A62=A26; A63=A36; A64=q44^2*B64; A65=q44^2*B65; A66=q44^2*B66;
%%----------gradient term kernal----------%
grad_kernal1=g11*kx.^2+g44*(ky.^2+kz.^2);
grad_kernal2=g11*ky.^2+g44*(kx.^2+kz.^2);
grad_kernal3=g11*kz.^2+g44*(ky.^2+kx.^2);
%%----------electrostatic term kernal----------%
elesta_kernal=4*pi*(epsilon11*kx.^2+epsilon11*ky.^2+epsilon33*kz.^2);
%%--------------history--------------%%
iter_batch=1200; iter_epoch=20;
iter_num=iter_epoch*iter_batch;
domain_history=cell(1,iter_epoch);
Px_history1=zeros(1,iter_num);
Px_history2=zeros(1,iter_num);
[x,y,z]=meshgrid(-scalex/2:scalex/2-1,-scaley/2:scaley/2-1,-scalez/2:scalez/2-1);


k_Px0=d*fftn(fftshift(Px0)); k_Py0=d*fftn(fftshift(Py0)); k_Pz0=d*fftn(fftshift(Pz0));
[dG_k_Px0,dG_k_Py0,dG_k_Pz0]=getF_3D(Px0,Py0,Pz0,scalex,scaley,scalez);
k_Px1=(k_Px0-delta_tao*dG_k_Px0)./(1+delta_tao*grad_kernal1);
k_Py1=(k_Py0-delta_tao*dG_k_Py0)./(1+delta_tao*grad_kernal2);
k_Pz1=(k_Pz0-delta_tao*dG_k_Pz0)./(1+delta_tao*grad_kernal3);
Px1=real(fftshift(ifftn(k_Px1/d))); Py1=real(fftshift(ifftn(k_Py1/d))); Pz1=real(fftshift(ifftn(k_Pz1/d)));

[dG_k_Px1,dG_k_Py1,dG_k_Pz1]=getF_3D(Px1,Py1,Pz1,scalex,scaley,scalez);
k_Px2=(4*k_Px1-k_Px0-delta_tao*(3*dG_k_Px1-dG_k_Px0))./(3+2*delta_tao*grad_kernal1);
k_Py2=(4*k_Py1-k_Py0-delta_tao*(3*dG_k_Py1-dG_k_Py0))./(3+2*delta_tao*grad_kernal2);
k_Pz2=(4*k_Pz1-k_Pz0-delta_tao*(3*dG_k_Pz1-dG_k_Pz0))./(3+2*delta_tao*grad_kernal3);
Px2=real(fftshift(ifftn(k_Px2/d))); Py2=real(fftshift(ifftn(k_Py2/d))); Pz2=real(fftshift(ifftn(k_Pz2/d)));

for i=1:iter_num
    %%-----------------derivative of Ginzburg-Landau energy-----------------%%
    dG_GL_k_Px2_alpha1=2.*alpha1.*Px2;
    dG_GL_k_Px2_alpha11=4.*alpha11.*Px2.^3;
    dG_GL_k_Px2_alpha12=2.*alpha12.*Px2.*(Py2.^2+Pz2.^2);
    dG_GL_k_Px2_alpha111=6.*alpha111.*Px2.^5;
    dG_GL_k_Px2_alpha112=alpha112.*(4.*Px2.^3.*(Py2.^2+Pz2.^2)+2.*Px2.*(Py2.^4+Pz2.^4));
    dG_GL_k_Px2_alpha123=alpha123*2*Px2.*Py2.^2.*Pz2.^2;
    dG_GL_k_Px2_alpha1111=alpha1111.*8.*Px2.^7;
    dG_GL_k_Px2_alpha1112=alpha1112.*(6.*Px2.^5.*(Py2.^2+Pz2.^2)+2.*Px2.*(Py2.^6+Pz2.^6));
    dG_GL_k_Px2_alpha1122=alpha1122.*4.*Px2.^3.*(Py2.^4+Pz2.^4);
    dG_GL_k_Px2_alpha1123=alpha1123*(4*Px2.^3.*Py2.^2.*Pz2.^2+2*Px2.*(Py2.^4.*Pz2.^2+Pz2.^4.*Py2.^2));

    dG_GL_k_Py2_alpha1=2.*alpha1.*Py2;
    dG_GL_k_Py2_alpha11=4.*alpha11.*Py2.^3;
    dG_GL_k_Py2_alpha12=2.*alpha12.*Py2.*(Px2.^2+Pz2.^2);
    dG_GL_k_Py2_alpha111=6.*alpha111.*Py2.^5;
    dG_GL_k_Py2_alpha112=alpha112.*(4.*Py2.^3.*(Px2.^2+Pz2.^2)+2.*Py2.*(Px2.^4+Pz2.^4));
    dG_GL_k_Py2_alpha123=alpha123*2*Py2.*Px2.^2.*Pz2.^2;
    dG_GL_k_Py2_alpha1111=alpha1111.*8.*Py2.^7;
    dG_GL_k_Py2_alpha1112=alpha1112.*(6.*Py2.^5.*(Px2.^2+Pz2.^2)+2.*Py2.*(Px2.^6+Pz2.^6));
    dG_GL_k_Py2_alpha1122=alpha1122.*4.*Py2.^3.*(Px2.^4+Pz2.^4);
    dG_GL_k_Py2_alpha1123=alpha1123*(4*Py2.^3.*Px2.^2.*Pz2.^2+2*Py2.*(Px2.^4.*Pz2.^2+Pz2.^4.*Px2.^2));

    dG_GL_k_Pz2_alpha1=2.*alpha1.*Pz2;
    dG_GL_k_Pz2_alpha11=4.*alpha11.*Pz2.^3;
    dG_GL_k_Pz2_alpha12=2.*alpha12.*Pz2.*(Px2.^2+Py2.^2);
    dG_GL_k_Pz2_alpha111=6.*alpha111.*Pz2.^5;
    dG_GL_k_Pz2_alpha112=alpha112.*(4.*Pz2.^3.*(Px2.^2+Py2.^2)+2.*Pz2.*(Px2.^4+Py2.^4));
    dG_GL_k_Pz2_alpha123=alpha123*2*Pz2.*Px2.^2.*Py2.^2;
    dG_GL_k_Pz2_alpha1111=alpha1111.*8.*Pz2.^7;
    dG_GL_k_Pz2_alpha1112=alpha1112.*(6.*Pz2.^5.*(Px2.^2+Py2.^2)+2.*Pz2.*(Px2.^6+Py2.^6));
    dG_GL_k_Pz2_alpha1122=alpha1122.*4.*Pz2.^3.*(Px2.^4+Py2.^4);
    dG_GL_k_Pz2_alpha1123=alpha1123*(4*Pz2.^3.*Px2.^2.*Py2.^2+2*Pz2.*(Py2.^4.*Px2.^2+Px2.^4.*Py2.^2));

    dG_GL_k_Px2=d*fftn(ifftshift(dG_GL_k_Px2_alpha1+dG_GL_k_Px2_alpha11+dG_GL_k_Px2_alpha12+dG_GL_k_Px2_alpha111+dG_GL_k_Px2_alpha112+dG_GL_k_Px2_alpha123+dG_GL_k_Px2_alpha1111+dG_GL_k_Px2_alpha1112+dG_GL_k_Px2_alpha1122+dG_GL_k_Px2_alpha1123));
    dG_GL_k_Py2=d*fftn(ifftshift(dG_GL_k_Py2_alpha1+dG_GL_k_Py2_alpha11+dG_GL_k_Py2_alpha12+dG_GL_k_Py2_alpha111+dG_GL_k_Py2_alpha112+dG_GL_k_Py2_alpha123+dG_GL_k_Py2_alpha1111+dG_GL_k_Py2_alpha1112+dG_GL_k_Py2_alpha1122+dG_GL_k_Py2_alpha1123));
    dG_GL_k_Pz2=d*fftn(ifftshift(dG_GL_k_Pz2_alpha1+dG_GL_k_Pz2_alpha11+dG_GL_k_Pz2_alpha12+dG_GL_k_Pz2_alpha111+dG_GL_k_Pz2_alpha112+dG_GL_k_Pz2_alpha123+dG_GL_k_Pz2_alpha1111+dG_GL_k_Pz2_alpha1112+dG_GL_k_Pz2_alpha1122+dG_GL_k_Pz2_alpha1123));

    
    %%----------------------------gradient energy----------------------------%%
    %%-------------------derivative of electrostatic energy------------------%%
    temp1=kx.*k_Px2+ky.*k_Py2+kz.*k_Pz2;
    dG_k_es_Px2=temp1.*kx./elesta_kernal;
    dG_k_es_Py2=temp1.*ky./elesta_kernal;
    dG_k_es_Pz2=temp1.*kz./elesta_kernal;
    dG_k_es_Px2(1,1,1)=0;
    dG_k_es_Py2(1,1,1)=0;
    dG_k_es_Pz2(1,1,1)=0;
    %------------------derivative of electromechanic energy-----------------%%
    Y1=Px2.^2; Y2=Py2.^2; Y3=Pz2.^2; Y4=Py2.*Pz2; Y5=Px2.*Pz2; Y6=Px2.*Py2;
    Y1_k=d*fftn(fftshift(Y1)); Y2_k=d*fftn(fftshift(Y2)); Y3_k=d*fftn(fftshift(Y3)); Y4_k=d*fftn(fftshift(Y4)); Y5_k=d*fftn(fftshift(Y5)); Y6_k=d*fftn(fftshift(Y6));
    kernal_sum1=real(fftshift(ifftn((A11.*Y1_k+A12.*Y2_k+A13.*Y3_k+A14.*Y4_k+A15.*Y5_k+A16.*Y6_k)/d)));
    kernal_sum2=real(fftshift(ifftn((A21.*Y1_k+A22.*Y2_k+A23.*Y3_k+A24.*Y4_k+A25.*Y5_k+A26.*Y6_k)/d)));
    kernal_sum3=real(fftshift(ifftn((A31.*Y1_k+A32.*Y2_k+A33.*Y3_k+A34.*Y4_k+A35.*Y5_k+A36.*Y6_k)/d)));
    kernal_sum4=real(fftshift(ifftn((A41.*Y1_k+A42.*Y2_k+A43.*Y3_k+A44.*Y4_k+A45.*Y5_k+A46.*Y6_k)/d)));
    kernal_sum5=real(fftshift(ifftn((A51.*Y1_k+A52.*Y2_k+A53.*Y3_k+A54.*Y4_k+A55.*Y5_k+A56.*Y6_k)/d)));
    kernal_sum6=real(fftshift(ifftn((A61.*Y1_k+A62.*Y2_k+A63.*Y3_k+A64.*Y4_k+A65.*Y5_k+A66.*Y6_k)/d)));

    dG_k_em_Px2=-0.5*d*fftn(fftshift(4*Px2.*kernal_sum1+2*Py2.*kernal_sum6+2*Pz2.*kernal_sum5));
    dG_k_em_Py2=-0.5*d*fftn(fftshift(4*Py2.*kernal_sum2+2*Pz2.*kernal_sum4+2*Px2.*kernal_sum6));
    dG_k_em_Pz2=-0.5*d*fftn(fftshift(4*Pz2.*kernal_sum3+2*Px2.*kernal_sum5+2*Py2.*kernal_sum4));
    dG_k_em_Px2(1,1,1)=0;
    dG_k_em_Py2(1,1,1)=0;
    dG_k_em_Pz2(1,1,1)=0;
    %------------------derivative of electromechanic energy-----------------%%
    %--------------------------------evolvetion------------------------------%%
    dG_k_Px2=dG_GL_k_Px2+dG_k_es_Px2+dG_k_em_Px2;
    dG_k_Py2=dG_GL_k_Py2+dG_k_es_Py2+dG_k_em_Py2;
    dG_k_Pz2=dG_GL_k_Pz2+dG_k_es_Pz2+dG_k_em_Pz2;
    %%%%%
    k_Px3=(18*k_Px2-9*k_Px1+2*k_Px0-delta_tao*(11.5*dG_k_Px2-8*dG_k_Px1+2.5*dG_k_Px0))./(11+6*delta_tao*grad_kernal1);
    k_Py3=(18*k_Py2-9*k_Py1+2*k_Py0-delta_tao*(11.5*dG_k_Py2-8*dG_k_Py1+2.5*dG_k_Py0))./(11+6*delta_tao*grad_kernal2);
    k_Pz3=(18*k_Pz2-9*k_Pz1+2*k_Pz0-delta_tao*(11.5*dG_k_Pz2-8*dG_k_Pz1+2.5*dG_k_Pz0))./(11+6*delta_tao*grad_kernal3);
    Px3=real(fftshift(ifftn(k_Px3/d))); Py3=real(fftshift(ifftn(k_Py3/d))); Pz3=real(fftshift(ifftn(k_Pz3/d)));
    
    if rem(i,iter_batch)==0
        domain_history{i/iter_batch}=Px2;
    end
    
    Px_history1(i)=Px3(60,60,4);
    Px_history2(i)=Px3(64,64,2);
    disp(Px3(60,60,4));
    
    k_Px0=k_Px1; k_Py0=k_Py1; k_Pz0=k_Pz1;
    k_Px1=k_Px2; k_Py1=k_Py2; k_Pz1=k_Pz2;
    k_Px2=k_Px3; k_Py2=k_Py3; k_Pz2=k_Pz3;
    
    Px2=Px3; Py2=Py3; Pz2=Pz3;
    
    dG_k_Px0=dG_k_Px1; dG_k_Py0=dG_k_Py1; dG_k_Pz0=dG_k_Pz1;
    dG_k_Px1=dG_k_Px2; dG_k_Py1=dG_k_Py2; dG_k_Pz1=dG_k_Pz2;
    


end
figure;
quiver3(x,y,z,Px3,Py3,Pz3);
figure;
for p=1:iter_epoch
    subplot(5,4,p);
    imagesc(domain_history{p}(:,:,1));
end
figure;
xslice=[-scalex/2,scalex/2-1];
yslice=[-scalex/2,scalex/2-1];
zslice=scalez/2-1;
slice(x,y,z,Px3,xslice,yslice,zslice);
axis off;