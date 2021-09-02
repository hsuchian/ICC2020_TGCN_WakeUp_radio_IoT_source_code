
n = 10;
xc = 0;
yc=0;
r=1;

theta = (0:n-1)*(2*pi/n);
x = xc + r*cos(theta);
y = yc + r*sin(theta);
P = polyshape(x,y);
plot(P);