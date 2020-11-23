function [strenchstress,shearstress]=bendstress(scale,L,F)
d=L/scale;
F_d=F/(scale/2);
M=zeros(scale);
for i=1:scale
    M(i,:)=[(F_d:F_d:F),((F-F_d):-F_d:0)]'*L/4;
end
[x,y]=meshgrid(-L/2:d:(L/2-d),-L/2:d:(L/2-d));
strenchstress=M.*y/(L^4/12);
shearstress=6*(F/2)*(L^2/4-y.^2)/(L^4);
end