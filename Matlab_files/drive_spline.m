% N33forread=table2array(N33forread);
% solutionyita1DN43=table2array(solutionyita1DN43);
% x=N33forread(2:end-1,2);
% y=N33forread(2:end-1,3);
% x_new=solutionyita1DN43(2:end-1,2);
% y_new=solutionyita1DN43(2:end-1,3);
% xp=linspace(0,0.7,100000);
% Nxp=length(xp);
% Nx=length(x);
% yp=spline_chen(x,y,xp);
% plot(x,y,'*-','MarkerSize',10);
% hold on;
% plot(xp,yp,'.:')
% xlim([0 0.7])
% plot(x_new,y_new,'o-','MarkerSize',10);








x=[1 2 3];
y=[-10.603773500914594, -16.99329984275693, -12.77544930644134];
xp=[ 0, 1, 0, 1, 0, 1, 0.5, 0.5, 0.5, 2, 2, 2, 1.5, 1.5, 1.5, 3, 3, 3, 2.5, 2.5, 2.5, 4, 4, 4, 3.5, 3.5, 3.5]

yp=spline_chen(x,y,xp)
plot(x,y,'*-','MarkerSize',10);
hold on;
plot(xp,yp,'o')
xlim([0 0.7])
plot(xp,yp,'o-','MarkerSize',10);
