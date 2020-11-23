clear;clc;
load('ferro298K.mat')
load('interfield298K.mat');

relaxtime1=1000;
wave_amount=1;%generate ac electric field
sampling=100;%generate ac electric field
Emax=1e6;%generate ac electric field
acE=acEfield('sin',wave_amount,Emax,sampling);%generate ac electric field
f=500;relaxtime2=round(1/(deltat*f*sampling));%set polarize iteration parameter
Px_mean=zeros(length(acE),15);
for k=1:15
    j=1;rhox=interfield_history{k,1}; rhoy=interfield_history{k,2};
    for i=1:relaxtime1
        dG_GL_k_Px2_alpha1=2.*alpha1.*Px2;
        dG_GL_k_Px2_alpha11=4.*alpha11.*Px2.^3;
        dG_GL_k_Px2_alpha12=2.*alpha12.*Px2.*Py2.^2;
        dG_GL_k_Px2_alpha111=6.*alpha111.*Px2.^5;
        dG_GL_k_Px2_alpha112=alpha112.*(4.*Px2.^3.*Py2.^2+2.*Px2.*Py2.^4);
        dG_GL_k_Px2_alpha1111=alpha1111.*8.*Px2.^7;
        dG_GL_k_Px2_alpha1112=alpha1112.*(6.*Px2.^5.*Py2.^2+2.*Px2.*Py2.^6);
        dG_GL_k_Px2_alpha1122=alpha1122.*4.*Px2.^3.*Py2.^4;

        dG_GL_k_Py2_alpha1=2.*alpha1.*Py2;
        dG_GL_k_Py2_alpha11=4.*alpha11.*Py2.^3;
        dG_GL_k_Py2_alpha12=2.*alpha12.*Py2.*Px2.^2;
        dG_GL_k_Py2_alpha111=6.*alpha111.*Py2.^5;
        dG_GL_k_Py2_alpha112=alpha112.*(4.*Py2.^3.*Px2.^2+2.*Py2.*Px2.^4);
        dG_GL_k_Py2_alpha1111=alpha1111.*8.*Py2.^7;
        dG_GL_k_Py2_alpha1112=alpha1112.*(6.*Py2.^5.*Px2.^2+2.*Py2.*Px2.^6);
        dG_GL_k_Py2_alpha1122=alpha1122.*4.*Py2.^3.*Px2.^4;

        dG_GL_k_Px2=d*fft2(ifftshift(dG_GL_k_Px2_alpha1+dG_GL_k_Px2_alpha11+dG_GL_k_Px2_alpha12+dG_GL_k_Px2_alpha111+dG_GL_k_Px2_alpha112+dG_GL_k_Px2_alpha1111+dG_GL_k_Px2_alpha1112+dG_GL_k_Px2_alpha1122-rhox));
        dG_GL_k_Py2=d*fft2(ifftshift(dG_GL_k_Py2_alpha1+dG_GL_k_Py2_alpha11+dG_GL_k_Py2_alpha12+dG_GL_k_Py2_alpha111+dG_GL_k_Py2_alpha112+dG_GL_k_Py2_alpha1111+dG_GL_k_Py2_alpha1112+dG_GL_k_Py2_alpha1122-rhoy));

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
        dG_k_Px2=dG_GL_k_Px2+dG_k_es_Px2+dG_k_em_Px2;
        dG_k_Py2=dG_GL_k_Py2+dG_k_es_Py2+dG_k_em_Py2;

        k_Px3=(18*k_Px2-9*k_Px1+2*k_Px0-delta_tao*(11.5*dG_k_Px2-8*dG_k_Px1+2.5*dG_k_Px0))./(11+6*delta_tao*grad_kernal1);
        k_Py3=(18*k_Py2-9*k_Py1+2*k_Py0-delta_tao*(11.5*dG_k_Py2-8*dG_k_Py1+2.5*dG_k_Py0))./(11+6*delta_tao*grad_kernal2);
        Px3=real(fftshift(ifft2(k_Px3/d))); Py3=real(fftshift(ifft2(k_Py3/d)));

        disp(Px3(60,60));

        k_Px0=k_Px1; k_Py0=k_Py1;
        k_Px1=k_Px2; k_Py1=k_Py2;
        k_Px2=k_Px3; k_Py2=k_Py3;

        Px2=Px3; Py2=Py3;

        dG_k_Px0=dG_k_Px1; dG_k_Py0=dG_k_Py1;
        dG_k_Px1=dG_k_Px2; dG_k_Py1=dG_k_Py2;
    end


    for E1=acE
        Ex=E1*ones(scale);
        for i=1:relaxtime2
            %%-----------------derivative of Ginzburg-Landau energy-----------------%%
            dG_GL_k_Px2_alpha1=2.*alpha1.*Px2;
            dG_GL_k_Px2_alpha11=4.*alpha11.*Px2.^3;
            dG_GL_k_Px2_alpha12=2.*alpha12.*Px2.*Py2.^2;
            dG_GL_k_Px2_alpha111=6.*alpha111.*Px2.^5;
            dG_GL_k_Px2_alpha112=alpha112.*(4.*Px2.^3.*Py2.^2+2.*Px2.*Py2.^4);
            dG_GL_k_Px2_alpha1111=alpha1111.*8.*Px2.^7;
            dG_GL_k_Px2_alpha1112=alpha1112.*(6.*Px2.^5.*Py2.^2+2.*Px2.*Py2.^6);
            dG_GL_k_Px2_alpha1122=alpha1122.*4.*Px2.^3.*Py2.^4;

            dG_GL_k_Py2_alpha1=2.*alpha1.*Py2;
            dG_GL_k_Py2_alpha11=4.*alpha11.*Py2.^3;
            dG_GL_k_Py2_alpha12=2.*alpha12.*Py2.*Px2.^2;
            dG_GL_k_Py2_alpha111=6.*alpha111.*Py2.^5;
            dG_GL_k_Py2_alpha112=alpha112.*(4.*Py2.^3.*Px2.^2+2.*Py2.*Px2.^4);
            dG_GL_k_Py2_alpha1111=alpha1111.*8.*Py2.^7;
            dG_GL_k_Py2_alpha1112=alpha1112.*(6.*Py2.^5.*Px2.^2+2.*Py2.*Px2.^6);
            dG_GL_k_Py2_alpha1122=alpha1122.*4.*Py2.^3.*Px2.^4;

            dG_GL_k_Px2=d*fft2(ifftshift(dG_GL_k_Px2_alpha1+dG_GL_k_Px2_alpha11+dG_GL_k_Px2_alpha12+dG_GL_k_Px2_alpha111+dG_GL_k_Px2_alpha112+dG_GL_k_Px2_alpha1111+dG_GL_k_Px2_alpha1112+dG_GL_k_Px2_alpha1122-rhox-Ex));
            dG_GL_k_Py2=d*fft2(ifftshift(dG_GL_k_Py2_alpha1+dG_GL_k_Py2_alpha11+dG_GL_k_Py2_alpha12+dG_GL_k_Py2_alpha111+dG_GL_k_Py2_alpha112+dG_GL_k_Py2_alpha1111+dG_GL_k_Py2_alpha1112+dG_GL_k_Py2_alpha1122-rhoy));

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
            dG_k_Px2=dG_GL_k_Px2+dG_k_es_Px2+dG_k_em_Px2;
            dG_k_Py2=dG_GL_k_Py2+dG_k_es_Py2+dG_k_em_Py2;

            k_Px3=(18*k_Px2-9*k_Px1+2*k_Px0-delta_tao*(11.5*dG_k_Px2-8*dG_k_Px1+2.5*dG_k_Px0))./(11+6*delta_tao*grad_kernal1);

            k_Py3=(18*k_Py2-9*k_Py1+2*k_Py0-delta_tao*(11.5*dG_k_Py2-8*dG_k_Py1+2.5*dG_k_Py0))./(11+6*delta_tao*grad_kernal2);
            Px3=real(fftshift(ifft2(k_Px3/d))); Py3=real(fftshift(ifft2(k_Py3/d)));

            disp(Px3(100,100));

            k_Px0=k_Px1; k_Py0=k_Py1;
            k_Px1=k_Px2; k_Py1=k_Py2;
            k_Px2=k_Px3; k_Py2=k_Py3;

            Px2=Px3; Py2=Py3;

            dG_k_Px0=dG_k_Px1; dG_k_Py0=dG_k_Py1;
            dG_k_Px1=dG_k_Px2; dG_k_Py1=dG_k_Py2;
        end
        Px_mean(j,k)=mean2(Px2);
        j=j+1;
    end
end

% t=(2*pi/sampling:(2*pi/sampling):2*pi);
epsilon1=(max(Px_mean,[],1)-min(Px_mean,[],1))/(2*Emax);
plot(epsilon1);
