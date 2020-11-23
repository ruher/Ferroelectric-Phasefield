clear;clc;close all;
width=840;height=700;
load('oxygenvacancy.mat');
set(gcf,'Position',[0,0,width,height]);
for p=1:iter_epoch
    subplot(4,5,p);
    imagesc(compositionfield{p});
    axis off;
end
figure;
set(gcf,'Position',[0,0,width,height]);
imagesc(c_O);
axis off;
figure;
set(gcf,'Position',[0,0,width,height]);
surf(phi);
axis off;
figure;
set(gcf,'Position',[0,0,width,height]);
quiver(x,y,-d_phi_x,-d_phi_y);
axis off;
figure;
set(gcf,'Position',[0,0,width,height]);
imagesc((d_phi_x.^2+d_phi_y.^2).^0.5)
axis off;
