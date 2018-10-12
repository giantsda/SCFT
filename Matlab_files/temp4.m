% solutionyita1DN141=table2array(solutionyita1DN141);
% solutionyita1DN225=table2array(solutionyita1DN225);
% solutionyita1DN33=table2array(solutionyita1DN33);
% solutionyita1DN333=table2array(solutionyita1DN333);
% solutionyita1DN43=table2array(solutionyita1DN43);
% solutionyita1DN522=table2array(solutionyita1DN522);
% solutionyita1DN59=table2array(solutionyita1DN59);
% solutionyita1DN91=table2array(solutionyita1DN91);
close all;
x=solutionyita1DN59(2:end-1,2);
y=solutionyita1DN59(2:end-1,3);
xp=linspace(0,3.7237,50000);
yp=spline_chen([0 1 3 2],[0 0 2 2],[1]);
yp=spline_chen(x,y,xp);
ypm = interp1(x,y,xp,'spline');
plot(xp,yp,xp,ypm);
legend('wo', 'mat');
xlim([0 0.2])


% plot(x,y,'*-');
% hold on;
% plot(xp,yp,'.:');
% plot(solutionyita1DN91(2:end-1,2),solutionyita1DN91(2:end-1,3),'o-');
% legend('old', 'spline', 'true');
% xlim([0 0.2])