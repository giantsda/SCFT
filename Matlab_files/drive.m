function drive
N=32+1;
folder_path='C:\Users\chenshen\Downloads' ;
file_name='Exp_m32_n2048_IE.res';
data=[];
filename1= [folder_path '\' file_name ] ;
disp(filename1);
delimiterIn = ' '; %read txt file
headerlinesIn = 9;
data_m = importdata(filename1,delimiterIn,headerlinesIn);
data=data_m.data;
yita_answer=data(:,3);
 
% out=myfun(yita_answer)
yita_answer(1)=[];
yita_answer(end)=[];
% yita_answer=yita_answer*0.18;
[Xr, found] = broyden_wq(@myfun, yita_answer, 1e-4, 1e-7, 20000, length(yita_answer));
dd=Xr(1)
Xr=Xr*1;
% [Xr, found] = AndMix (@myfun,yita_answer , 0.00001, 20, length(yita_answer))
tryd=[];
% found=0
% while (found==0)
%     [Xr, found] = AndMix (@myfun,Xr ,  1e-7, 20, length(yita_answer))
%     tryd=[tryd Xr];
% end

plot(Xr)
hold on;
plot(yita_answer)
legend('my','his')
Xr(1)

 

function out=myfun(yita)
f0=simple_FEM_1D_transient(yita);
tau=0.5302;
L= 3.72374;
N=32+1;
x=linspace(0,L,N);
x_left=x(1:ceil(N*tau/L));
f0_given_left=(exp(4*tau*x_left./(tau*tau-x_left.*x_left))-1).^2./((exp(4*tau*x_left./(tau*tau-x_left.*x_left))+1).^2);
f0_given=ones(1,N);
f0_given(1:length(f0_given_left))=f0_given_left;
f0_given(end-length(f0_given_left)+1:end)=fliplr(f0_given_left);
f0_given(isnan(f0_given)) = 1;
f0_given=f0_given.';
out=f0_given-f0;
out(1)=[];
out(end)=[];
% plot(f0);
% hold on;
% plot(f0_given)
% norm(out)
% pause()