% % function drive
% 
% 
% global info;
% tau=0.5302;
% L= 3.72374;
% folder_path='/home/chen/Desktop/project/SCFT/Matlab_files';
% file_name='N=33_for_read.txt';
% filename= [folder_path '/' file_name ] ;
% %% run
% read_input(filename)
% set_f0_given(tau,L);
% out=simple_FEM_1D_transient(info.yita_middle_1D)
% % [Xr, found] = broyden_wq(@myfun, yita_answer, 1e-4, 1e-7, 20000, length(yita_answer));
% N=info.N;
% yita_middle_1D=info.yita_middle_1D;
% x_old=adm_chen (@simple_FEM_1D_transient,N-2,yita_middle_1D, 1e-7, 50000,0.9);
% 
% 
% 
% mesh1=info.mesh;
% sol1=x_old;
% yita_middle_1D=refine_mesh(x_old);
% N=info.N;
% info.yita_middle_1D=yita_middle_1D;
% set_f0_given(tau,L);

% for i=1:10
%     yita_middle_1D=refine_mesh(x_old);
%     N=info.N;
%     info.yita_middle_1D=yita_middle_1D;
%     set_f0_given(tau,L);
%     x_old=adm_chen (@simple_FEM_1D_transient,N-2,yita_middle_1D, 1e-7, 50000,0.009);
% end

plot(mesh1(2:end-1),sol1,'x-');
hold on;
plot(info.mesh(2:end-1),yita_middle_1D,'*-')
plot(solutionyita1DN43(2:end-1,2),solutionyita1DN43(2:end-1,3),'r:o')
legend('solution of N=33','interpolated N=43', 'real N=43');