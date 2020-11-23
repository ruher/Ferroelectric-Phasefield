clear;clc;
load('relaxorS.mat')
%%--------------cooling----------------%%
relaxtime1=20000;
Px_history=zeros(relaxtime1,1);
Tmin=220;
for i=1:relaxtime1
    alpha1=4.124e5*(Tmin-388);
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

    dG_GL_k_Px2=d*fft2(ifftshift(dG_GL_k_Px2_alpha1+dG_GL_k_Px2_alpha11+dG_GL_k_Px2_alpha12+dG_GL_k_Px2_alpha111+dG_GL_k_Px2_alpha112+dG_GL_k_Px2_alpha1111+dG_GL_k_Px2_alpha1112+dG_GL_k_Px2_alpha1122-2*Px2.*(Q11*exx+Q12*eyy)-0.5*Q44*exy.*Py2-Ernd_x));
    dG_GL_k_Py2=d*fft2(ifftshift(dG_GL_k_Py2_alpha1+dG_GL_k_Py2_alpha11+dG_GL_k_Py2_alpha12+dG_GL_k_Py2_alpha111+dG_GL_k_Py2_alpha112+dG_GL_k_Py2_alpha1111+dG_GL_k_Py2_alpha1112+dG_GL_k_Py2_alpha1122-2*Py2.*(Q11*eyy+Q12*exx)-0.5*Q44*exy.*Px2-Ernd_y));

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
    
    Px_history(i)=Px3(100,100);
    disp(Px3(100,100));

    k_Px0=k_Px1; k_Py0=k_Py1;
    k_Px1=k_Px2; k_Py1=k_Py2;
    k_Px2=k_Px3; k_Py2=k_Py3;

    Px2=Px3; Py2=Py3;

    dG_k_Px0=dG_k_Px1; dG_k_Py0=dG_k_Py1;
    dG_k_Px1=dG_k_Px2; dG_k_Py1=dG_k_Py2;
end
% pause;
%%--------------calculate conductivity versus temperature--------------%%
Tmax=420;%generate temperature gradient
T_series=generateTgrad(Tmin,Tmax,Tmax-Tmin,false);%generate temperature gradient
j=1;k=0;l=0;%set heating iteration parameter
relaxtime2=500*ones(length(T_series),1);
for i=1:length(T_series)
    if T_series(i)>=250&&T_series(i)<=280
        relaxtime2(i)=500;
    elseif T_series(i)>=380&&T_series(i)<=410
        relaxtime2(i)=8000;
    end
end
wave_amount=1;%generate ac electric field
sampling=100;%generate ac electric field
Emax=1e6;%generate ac electric field
acE=acEfield('sin',wave_amount,Emax,sampling);%generate ac electric field
f=500;relaxtime3=round(1/(deltat*f*sampling));%set polarize iteration parameter
Px_mean=zeros(length(T_series),length(acE));
domain_history=cell(length(T_series),2);
for T=T_series
%%-------------heating--------------%%
    l=l+1;
    for i=1:relaxtime2(l)
        %%-----------------derivative of Ginzburg-Landau energy-----------------%%
        alpha1=4.124e5*(T-388);
        
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

        dG_GL_k_Px2=d*fft2(ifftshift(dG_GL_k_Px2_alpha1+dG_GL_k_Px2_alpha11+dG_GL_k_Px2_alpha12+dG_GL_k_Px2_alpha111+dG_GL_k_Px2_alpha112+dG_GL_k_Px2_alpha1111+dG_GL_k_Px2_alpha1112+dG_GL_k_Px2_alpha1122-2*Px2.*(Q11*exx+Q12*eyy)-0.5*Q44*exy.*Py2-Ernd_x));
        dG_GL_k_Py2=d*fft2(ifftshift(dG_GL_k_Py2_alpha1+dG_GL_k_Py2_alpha11+dG_GL_k_Py2_alpha12+dG_GL_k_Py2_alpha111+dG_GL_k_Py2_alpha112+dG_GL_k_Py2_alpha1111+dG_GL_k_Py2_alpha1112+dG_GL_k_Py2_alpha1122-2*Py2.*(Q11*eyy+Q12*exx)-0.5*Q44*exy.*Px2-Ernd_y));

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

        disp(Px2(60,60));

        k_Px0=k_Px1; k_Py0=k_Py1;
        k_Px1=k_Px2; k_Py1=k_Py2;
        k_Px2=k_Px3; k_Py2=k_Py3;

        Px2=Px3; Py2=Py3;

        dG_k_Px0=dG_k_Px1; dG_k_Py0=dG_k_Py1;
        dG_k_Px1=dG_k_Px2; dG_k_Py1=dG_k_Py2;
    end
%%--------------apply electric field and calculate dielectric constant------------%%  
    j=1;
    for E1=acE
        Ex=E1*ones(scale);
        for i=1:relaxtime3
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

            dG_GL_k_Px2=d*fft2(ifftshift(dG_GL_k_Px2_alpha1+dG_GL_k_Px2_alpha11+dG_GL_k_Px2_alpha12+dG_GL_k_Px2_alpha111+dG_GL_k_Px2_alpha112+dG_GL_k_Px2_alpha1111+dG_GL_k_Px2_alpha1112+dG_GL_k_Px2_alpha1122-2*Px2.*(Q11*exx+Q12*eyy)-0.5*Q44*exy.*Py2-Ernd_x-Ex));
            dG_GL_k_Py2=d*fft2(ifftshift(dG_GL_k_Py2_alpha1+dG_GL_k_Py2_alpha11+dG_GL_k_Py2_alpha12+dG_GL_k_Py2_alpha111+dG_GL_k_Py2_alpha112+dG_GL_k_Py2_alpha1111+dG_GL_k_Py2_alpha1112+dG_GL_k_Py2_alpha1122-2*Py2.*(Q11*eyy+Q12*exx)-0.5*Q44*exy.*Px2-Ernd_y));

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
    Px_mean(l,j)=mean2(Px2);
 
    j=j+1;
    end
    domain_history{l,1}=Px2; domain_history{l,2}=Py2; 
end
t=(2*pi/sampling:(2*pi/sampling):2*pi);


% epsilon=Px_mean(:,50)./Emax;
% plot(T_series,epsilon)
epsilon1=zeros(length(T_series),1);
delta=zeros(length(T_series),1);
for i=1:length(T_series)
epsilon1(i)=(max(Px_mean(i,:))-min(Px_mean(i,:)))/2*Emax;
[fitresult,error]=sinfit(t,Px_mean(i,:));
delta(i)=-fitresult.c;
end
epsilon2=tan(delta).*epsilon1;
figure;
subplot(1,3,1);
plot(T_series,epsilon1,'o');
subplot(1,3,2);
plot(T_series,epsilon2,'o');
subplot(1,3,3)
plot(T_series,delta,'o');

% in=1;
% in2=60;
% subplot(1,2,1)
% quiver(x,y,domain_history{in,1},domain_history{in,2})
% subplot(1,2,2)
% quiver(x,y,domain_history{in2,1},domain_history{in2,2})



% for i=101:102
% plot(acE,Px_mean(i,:));
% hold on;
% end
% 
% plot(t,acE/5e7,'o');
% hold on;
% for i=1:20
% plot(t,Px_mean(118+i,:));
% hold on;
% end
% legend;