
function drive_adm_chen
x_old=[1;2;3];
% [x_old, found] = AndMix (@myfun, x_old, 0.000001, 50000, 3)
x_old=adm_chen (@myfun, 3, x_old, 1e-12, 50000,0.99,30);
  

function out=myfun(in)
x=in(1);
y=in(2);
z=in(3);

out(1) = x*y*z-12.;
out(2) = x*x+y*y-8;
out(3) = x+y+z-511;
