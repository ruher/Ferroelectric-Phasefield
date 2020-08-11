clear;clc;
load('domain_298K_ela035_poled.mat');
%%----------aging----------%%
% load('intfield_aging.mat');%aging before poling
% rhox=1e4*ones(scale,scale);rhoy=zeros(scale);%aging after poling
rhox=zeros(scale);rhoy=zeros(scale);%no aging
relaxtion=100;
%%--------------temperature change--------------%%
Tmin=298;Tmax=450;
T_series=generateTgrad(Tmin,Tmax,Tmax-Tmin,true);

for i=1:relaxtion
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
%%----------change temperature after poling and aging----------%%
j=1;
iter_num=zeros(length(T_series),1);
for i=1:((length(T_series)+1)/2)
    if T_series(i)>=402&&T_series(i)<=410
        iter_num(i)=300;
    else
        iter_num(i)=50;
    end
end
for i=((length(T_series)+1)/2):length(T_series)
    if T_series(i)>=380&&T_series(i)<=390
        iter_num(i)=2000;
    elseif T_series(i)>390
        iter_num(i)=50;
    else
        iter_num(i)=100;
    end
end
Px_history=zeros(sum(iter_num),1);
P_mean=zeros(length(T_series),1);
P_ext=zeros(length(T_series),1);
P_ext_Px=zeros(length(T_series),1);
domain_history=cell(7,2);
k=0;
ela_reduc_para=0.5;
for T=T_series
    allsteps_present=sum(iter_num(1:j))-iter_num(j);
    flucstd=13e-5*exp(-(T-Tc)^2/4)+2*kb*T*GAMMA;
    for i=1:iter_num(j)
        alpha1=4.124e5*(T-Tc);
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
        
        noise=normrnd(0,flucstd,[scale,scale]);
        k_Px3=(18*k_Px2-9*k_Px1+2*k_Px0-delta_tao*(11.5*dG_k_Px2-8*dG_k_Px1+2.5*dG_k_Px0)+6*deltat*noise)./(11+6*delta_tao*grad_kernal1);%%+6*deltat*noise
        noise=normrnd(0,flucstd,[scale,scale]);
        k_Py3=(18*k_Py2-9*k_Py1+2*k_Py0-delta_tao*(11.5*dG_k_Py2-8*dG_k_Py1+2.5*dG_k_Py0)+6*deltat*noise)./(11+6*delta_tao*grad_kernal2);
        Px3=real(fftshift(ifft2(k_Px3/d))); Py3=real(fftshift(ifft2(k_Py3/d)));

        Px_history(allsteps_present+i)=Px3(60,60);

        disp(Px3(60,60));

        k_Px0=k_Px1; k_Py0=k_Py1;
        k_Px1=k_Px2; k_Py1=k_Py2;
        k_Px2=k_Px3; k_Py2=k_Py3;

        Px2=Px3; Py2=Py3;

        dG_k_Px0=dG_k_Px1; dG_k_Py0=dG_k_Py1;
        dG_k_Px1=dG_k_Px2; dG_k_Py1=dG_k_Py2;
    end
    P_mean(j)=mean2(sqrt(Px2.^2+Py2.^2));
    P_ext(j)=sqrt(mean2(Px2)^2+mean2(Py2)^2);
    P_ext_Px(j)=sqrt(mean2(Px2)^2);
    j=j+1;
    if T==298||T==350||T==386||T==450
        k=k+1;
        domain_history{k,1}=Px2;
        domain_history{k,2}=Py2;
    end
end
subplot(2,2,1);
plot(T_series,P_mean);
title('P_{mean}')
subplot(2,2,2);
plot(T_series,P_ext);
title('P_{ext}')
subplot(2,2,3);
plot(T_series,P_ext_Px);
title('P_{ext Px}')
figure;
quiver(x,y,Px2,Py2);