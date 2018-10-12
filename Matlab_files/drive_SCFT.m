% function drive

global info;
tau=0.5302;
L= 3.72374;
folder_path='./inputFiles';
file_name='N=33_for_read.txt';
% file_name='solution_yita_1D_N= 43.txt';
% file_name='solution_yita_1D_N= 59.txt';
filename= [folder_path '/' file_name ] ;
%% run
read_input(filename)
set_f0_given(tau,L);
N=info.N;
yita_middle_1D=info.yita_middle_1D;
% [x_old, found] = broyden_wq(@simple_FEM_1D_transient, yita_middle_1D,0.00001, 1e-7, 20000, N-2);
x_old=adm_chen (@simple_FEM_1D_transient,N-2,yita_middle_1D, 1e-7, 50000,0.99);


for i=1:10
    yita_middle_1D=refine_mesh(x_old);
    N=info.N;
    info.yita_middle_1D=yita_middle_1D;
    set_f0_given(tau,L);
%     [x_old, found] = broyden_wq(@simple_FEM_1D_transient, yita_middle_1D,0.00001, 1e-7, 20000, N-2);
    x_old=adm_chen (@simple_FEM_1D_transient,N-2,yita_middle_1D, 1e-7, 50000,0.59);
end

 