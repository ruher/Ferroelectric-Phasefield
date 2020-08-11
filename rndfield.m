function [stressfield,elecfieldx,elecfieldy]=rndfield(n,L,random_point_number,range)
d=L/n;
stressfield=zeros(n);
random_direction=rand(random_point_number,1);
direc_x=cos(random_direction);direc_y=sin(random_direction);
elecfieldx=zeros(n);
elecfieldy=zeros(n);
random_point=randi(128*128,random_point_number,1);
random_x=round(random_point/128);
random_y=mod(random_point,128);
[x,y]=meshgrid(1:128,1:128);
for i=1:random_point_number
    distance=d^2*((x-random_x(i)).^2+(y-random_y(i)).^2);
    stressfield=stressfield+5e8*exp(-distance/(range));
    elecfieldx=elecfieldx+1e6*exp(-distance/(range)).*direc_x(i);
    elecfieldy=elecfieldy+1e6*exp(-distance/(range)).*direc_y(i);
end
end