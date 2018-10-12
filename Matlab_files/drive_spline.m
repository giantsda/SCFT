x=N33forread(2:end-1,2);
y=N33forread(2:end-1,3);
xp=linspace(0,0.7,100000);
Nxp=length(xp);
Nx=length(x);
yp=spline_chen(x,y,xp,8);
plot(x,y,'*-');
hold on;
plot(xp,yp,'.:')

