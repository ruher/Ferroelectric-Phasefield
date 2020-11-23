clear; clc; close all;
rho=2e7;
rho0=2.5; interval=0.1;
[rhox,rhoy]=meshgrid(-rho0:interval:rho0,-rho0:interval:rho0);
alpha2=1e-30; alpha1=-2*alpha2*rho^2; alpha3=2.5*alpha2; mu=1e4;
% G=exp(alpha1.*(rhox.^2+rhoy.^2)+alpha2.*(rhox.^4+rhoy.^4)+alpha3.*rhox.^2.*rhoy.^2);
% surf(G);
% min(min(G))
surf(-log(rhox.^2+rhoy.^2));