function [dG_k_Px,dG_k_Py]=getF_2D(Px,Py,scale,ela_reduc_para)
    L=0.125e-6;
    T=298;
    d=(L/scale)^2;
    %%------------alpha------------%%
    alpha1=4.124e5*(T-388); 
    alpha11=-2.097e8; alpha12=7.974e8;
    alpha111=1.294e9; alpha112=-1.950e9;
    alpha1111=3.863e10; alpha1112=2.529e10; alpha1122=1.637e10;
    epsilon0=8.854187e-12;epsilon11=7.35*epsilon0;
    %%--------------kspace--------------%%
    K=2*pi*scale/L;dk=K/scale;
    [kx,ky]=meshgrid(-K/2:dk:(K/2-dk),(-K/2):dk:(K/2-dk));
    kx=fftshift(kx);
    ky=fftshift(ky);
    kx_n=kx./sqrt(kx.^2+ky.^2);ky_n=ky./sqrt(kx.^2+ky.^2);
    kx_n(1,1)=0;ky_n(1,1)=0;
    %%------------electric mechanical coefficient and kernal------------%%
    C11=27.5e10; C12=17.9e10; C44=5.43e10;
    q11=14.2e9; q12=-0.74e9; q44=1.57e9;
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
    %%----------electrostatic term kernal----------%
    elesta_kernal=4*pi*(epsilon11*kx.^2+epsilon11*ky.^2);

    dG_GL_k_Px_alpha1=2.*alpha1.*Px;
    dG_GL_k_Px_alpha11=4.*alpha11.*Px.^3;
    dG_GL_k_Px_alpha12=2.*alpha12.*Px.*Py.^2;
    dG_GL_k_Px_alpha111=6.*alpha111.*Px.^5;
    dG_GL_k_Px_alpha112=alpha112.*(4.*Px.^3.*Py.^2+2.*Px.*Py.^4);
    dG_GL_k_Px_alpha1111=alpha1111.*8.*Px.^7;
    dG_GL_k_Px_alpha1112=alpha1112.*(6.*Px.^5.*Py.^2+2.*Px.*Py.^6);
    dG_GL_k_Px_alpha1122=alpha1122.*4.*Px.^3.*Py.^4;

    dG_GL_k_Py_alpha1=2.*alpha1.*Py;
    dG_GL_k_Py_alpha11=4.*alpha11.*Py.^3;
    dG_GL_k_Py_alpha12=2.*alpha12.*Py.*Px.^2;
    dG_GL_k_Py_alpha111=6.*alpha111.*Py.^5;
    dG_GL_k_Py_alpha112=alpha112.*(4.*Py.^3.*Px.^2+2.*Py.*Px.^4);
    dG_GL_k_Py_alpha1111=alpha1111.*8.*Py.^7;
    dG_GL_k_Py_alpha1112=alpha1112.*(6.*Py.^5.*Px.^2+2.*Py.*Px.^6);
    dG_GL_k_Py_alpha1122=alpha1122.*4.*Py.^3.*Px.^4;

    dG_GL_k_Px=d*fft2(ifftshift(dG_GL_k_Px_alpha1+dG_GL_k_Px_alpha11+dG_GL_k_Px_alpha12+dG_GL_k_Px_alpha111+dG_GL_k_Px_alpha112+dG_GL_k_Px_alpha1111+dG_GL_k_Px_alpha1112+dG_GL_k_Px_alpha1122));
    dG_GL_k_Py=d*fft2(ifftshift(dG_GL_k_Py_alpha1+dG_GL_k_Py_alpha11+dG_GL_k_Py_alpha12+dG_GL_k_Py_alpha111+dG_GL_k_Py_alpha112+dG_GL_k_Py_alpha1111+dG_GL_k_Py_alpha1112+dG_GL_k_Py_alpha1122));

    %%----------------------------gradient energy----------------------------%%
    k_Px=d*fft2(fftshift(Px));
    k_Py=d*fft2(fftshift(Py));
    %%----------------------------gradient energy----------------------------%%
    %%-------------------derivative of electrostatic energy------------------%%
    temp1=kx.*k_Px+ky.*k_Py;
    dG_k_es_Px=temp1.*kx./elesta_kernal;
    dG_k_es_Py=temp1.*ky./elesta_kernal;
    dG_k_es_Px(1,1)=0;
    dG_k_es_Py(1,1)=0;
    %------------------derivative of electromechanic energy-----------------%%
    Y1=Px.^2; Y2=Py.^2; Y6=Px.*Py;
    Y1_k=d*fft2(fftshift(Y1)); Y2_k=d*fft2(fftshift(Y2)); Y6_k=d*fft2(fftshift(Y6));
    kernal_sum1=real(fftshift(ifft2((A11.*Y1_k+A12.*Y2_k+A16.*Y6_k)/d)));
    kernal_sum2=real(fftshift(ifft2((A21.*Y1_k+A22.*Y2_k+A26.*Y6_k)/d)));
    kernal_sum6=real(fftshift(ifft2((A61.*Y1_k+A62.*Y2_k+A66.*Y6_k)/d)));

    dG_k_em_Px=-0.5*d*fft2(fftshift(4*Px.*kernal_sum1+2*Py.*kernal_sum6))*ela_reduc_para;
    dG_k_em_Py=-0.5*d*fft2(fftshift(4*Py.*kernal_sum2+2*Px.*kernal_sum6))*ela_reduc_para;
    dG_k_em_Px(1,1)=0;
    dG_k_em_Py(1,1)=0;


    dG_k_Px=dG_GL_k_Px+dG_k_es_Px+dG_k_em_Px;
    dG_k_Py=dG_GL_k_Py+dG_k_es_Py+dG_k_em_Py;
end