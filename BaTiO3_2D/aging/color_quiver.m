function color_quiver(Px,Py,scale,color)
%scale control the length of arrow in quiver plot
%color control the color of arrow
[m,n]=size(Px);

Pvx=zeros(m,n);Pvy=zeros(m,n);

arrowdensity=16;

for i=1:m
    for j=1:n
        if rem(i,arrowdensity)==1&&rem(j,arrowdensity)==1
            Pvx(i,j)=Px(i,j);
            Pvy(i,j)=Py(i,j);
        end
            
    end
end
Pvx=reshape(Pvx,m,n);Pvy=reshape(Pvy,m,n);
[x,y]=meshgrid(1:m,1:n);
imagesc(Px+Py*2,[-0.56,0.56]);
hold on;
quiver(x,y,Pvx,Pvy,(sqrt(mean2(Px.^2+Py.^2))/0.2)*scale,color,'LineWidth',2.5);
end