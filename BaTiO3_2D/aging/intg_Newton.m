function I=intg_Newton(x,y)
%{
this function used Newton method to calculate the numerical integral

x(vector): points on horizontal coordinate.X should be increasing;

y(vector): the function value.
%}
n=length(x);
h=x(2:n)-x(1:n-1);
I=0;
for i=1:n-1
    I=I+(h(i)/2)*(y(i)+y(i+1));
end
end