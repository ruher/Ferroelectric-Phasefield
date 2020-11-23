clear;clc;
load('3Ddomain_298K.mat');
%%--------------generate Temperature gradient--------------%%
Tmin=299;Tmax=500;
Tgrad=generateTgrad(Tmin,Tmax,Tmax-Tmin,true);
%%----------iteration setup----------%%
relaxtime=100;
iter_num=length(Tgrad)*relaxtime;
P_mean=zeros(length(Tgrad),1);
Px_history=zeros(iter_num,1);
[x,y,z]=meshgrid(-scalex/2:scalex/2-1,-scaley/2:scaley/2-1,-scalez/2:scalez/2-1);
for i=1:iter_num
    T=Tgrad(ceil(i/relaxtime));
    alpha1=4.124e5*(T-388);
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
    
    disp(Px3(60,60,4));
    Px_history(i)=Px3(60,60,4);
    
    k_Px0=k_Px1; k_Py0=k_Py1; k_Pz0=k_Pz1;
    k_Px1=k_Px2; k_Py1=k_Py2; k_Pz1=k_Pz2;
    k_Px2=k_Px3; k_Py2=k_Py3; k_Pz2=k_Pz3;
    
    Px2=Px3; Py2=Py3; Pz2=Pz3;
    
    dG_k_Px0=dG_k_Px1; dG_k_Py0=dG_k_Py1; dG_k_Pz0=dG_k_Pz1;
    dG_k_Px1=dG_k_Px2; dG_k_Py1=dG_k_Py2; dG_k_Pz1=dG_k_Pz2;

    if rem(i,relaxtime)==0
        P_mean(i/relaxtime)=mean2(sqrt(Px2.^2+Py2.^2+Pz2.^2));
    end

end
plot(Tgrad,P_mean);
figure;
quiver3(x,y,z,Px3,Py3,Pz3);
