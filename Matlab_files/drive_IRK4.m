global info;
info.tau=5.30252230020752e-01;
info.L= 3.72374357332160;
folder_path='./inputFiles';
file_name='Exp_m32_n2048_IE.res';
% file_name='solution_yita_1D_N= 43.txt';
% file_name='solution_yita_1D_N= 59.txt';
filename= [folder_path '/' file_name ] ;
%% run
read_input(filename)
set_f0_given(tau,L);
N=info.N;
yita_middle_1D=info.yita_middle_1D;
% yita_middle_1D(:)=0;
out=IRK4(yita_middle_1D)
% [x_old, found] = broyden_wq(@simple_FEM_1D_transient, yita_middle_1D,0.00001, 1e-7, 20000, N-2);
% x_old=adm_chen (@IRK4,N-2,yita_middle_1D, 1e-7, 50000,0.9,10); 