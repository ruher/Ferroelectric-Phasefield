GL_alpha1=alpha1*Px2.^2+alpha1*Py2.^2;
GL_alphabar11=alpha11*Px2.^4+alpha11*Py2.^4;
GL_alphabar12=alpha12*Px2.^2.*Py2.^2;
GL_alpha111=alpha111*Px2.^6+alpha111*Py2.^6;
GL_alpha112=alpha112*(Px2.^4.*Py2.^2+Px2.^2.*Py2.^4);
GL_alpha1111=alpha1111*Px2.^8+alpha1111*Py2.^8;
GL_alpha1112=alpha1112*(Px2.^6.*Py2.^2+Px2.^2.*Py2.^6);
GL_alpha1122=alpha1122*Px2.^4.*Py2.^4;
G_GL=GL_alpha1+GL_alphabar11+GL_alphabar12+GL_alpha111+GL_alpha112+GL_alpha1111+GL_alpha1112+GL_alpha1122-Ex.*Px2-Ey.*Py2-exx.*(Q11*Px2.^2+Q12*Py2.^2)-eyy.*(Q11*Py2.^2+Q12*Px2.^2)-0.5*Q44*exy.*Px2.*Py2;

G_grad=0.5*g11*(-kx.^2.*k_Px2.^2-ky.^2.*k_Py2.^2)+0.5*g44*(-ky.^2.*k_Px2.^2-kx.^2.*k_Py2.^2);

G_es=-0.5*(k_Px2.*kx+k_Py2.*Py)./(epsilon11*(kx.^2+ky.^2));

Y1=Px2.^2; Y2=Py2.^2; Y6=Px2.*Py2;
Y1_k=d*fft2(fftshift(Y1)); Y2_k=d*fft2(fftshift(Y2)); Y6_k=d*fft2(fftshift(Y6));
G_em=A11.*Y1.*Y1+A12.*Y1.*Y2+A21.*Y2.*Y1+A22.*Y2.*Y2;

G(j)=sum(sum(G_GL+G_grad+G_es+G_em*ela_reduc_para));

G_temp=[G;0;0]+[0;0;G]-2*[0;G;0];
Cv=G_temp(2:(length(T_series)-1)).*T_series(2:(length(T_series)-1));

