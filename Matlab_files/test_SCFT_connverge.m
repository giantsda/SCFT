% % myvars = who;
% % for i=1:length(myvars)
% %     v=char(myvars(i));
% %     eval([v '=table2array(' v ')']);
% % end
% % close all;
% clc;
% 
% global info;
% tau=0.5302;
% L= 3.72374;
% mesh=solutionmatlabN33(:,2);
% x_old=solutionmatlabN33(:,3);
% x_old=x_old(2:end-1);
% N=length(mesh);
% info.N=N;
% info.mesh=mesh;
% set_f0_given(tau,L);
% yita_middle_1D=refine_mesh(x_old);
% N=info.N;  
% info.yita_middle_1D=yita_middle_1D;
% xp=linspace(0,1,5000);
% yp=spline_chen(mesh(2:end-1),x_old,xp);
% set_f0_given(tau,L);
% 
% % plot(mesh(2:end-1),x_old,'*-','MarkerSize',10);
% % hold on;
% % plot(xp,yp,'LineWidth',2);
% % plot(info.mesh(2:end-1),yita_middle_1D,'o-','MarkerSize',10);
% % xlim([0 1])
% % legend('old','interpolated','interpolated');
% 
% out=adm_chen (@simple_FEM_1D_transient,N-2,yita_middle_1D, 1e-4, 50000,0.99,5);
% 
%  



% function drive

global info;
tau=0.5302;
L= 3.72374;
folder_path='./inputFiles';
% file_name='N=33_for_read.txt';
% file_name='solution_yita_1D_N= 43.txt';
file_name='solution_yita_1D_N= 59.txt';
filename= [folder_path '/' file_name ] ;
%% run
read_input(filename)
set_f0_given(tau,L);
N=info.N;
yita_middle_1D=info.yita_middle_1D;
% [x_old, found] = broyden_wq(@simple_FEM_1D_transient, yita_middle_1D,0.00001, 1e-7, 20000, N-2);
% x_old=adm_chen (@simple_FEM_1D_transient,N-2,yita_middle_1D, 1e-7, 50000,0.9,10);

yita_middle_1D=adm_chen (@simple_FEM_1D_transient,N-2,yita_middle_1D, 0.1, 500,0.99,2);
yita_middle_1D=adm_chen (@simple_FEM_1D_transient,N-2,yita_middle_1D, 1e-4, 50000,0.9,5);
yita_middle_1D=adm_chen (@simple_FEM_1D_transient,N-2,yita_middle_1D, 1e-7, 50000,0.9,15);

