% x,y is the table waitting to be interpolated.
% xp is the x vector you want to know yp at.
% if not spcify m, use not-a-knot end condition
% m=0 for natural spline;otherwise m is the ddy at boundarys
x=[1 2 3 4 5];
yp=spline_chen(x,x,1.2345,0) 