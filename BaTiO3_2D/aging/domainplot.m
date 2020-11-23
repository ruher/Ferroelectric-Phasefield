clc;clear;close all;
[x,y]=meshgrid(1:128,1:128);
load('298K_ela025.mat');
quiver(x,y,Px2,Py2)
color_quiver(Px2,Py2,7,'white');
axis off;

color1=[228 146 56]/256;
color2=[244 214 182]/256;
color3=[161 241 249]/256;
color4=[50 164 234]/256;

n=100;
map1=zeros(n,3);
for i=1:n
map1(i,:)=color1-(color1-color2)*(i-1)/n;
end
map2=zeros(n,3);
for i=1:n
map2(i,:)=color3-(color3-color4)*(i-1)/n;
end

map=[map1;map2];

colormap(map);