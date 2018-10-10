function new_solution=refine_mesh(x_old)
global info;
node=info.mesh;
N=info.N;
solution=[Inf;x_old;Inf];
for cell=1:N-1
    left=node(cell);
    right=node(cell+1);
    error(cell,2)=abs((solution(cell+1)-solution(cell))/(node(cell+1)-node(cell)));
end
error(:,1)=1:N-1;
refine=[];
threshold=median(error(:,2))*10;
for cell=1:N-1
    if error(cell,2)>=threshold
        error(cell,3)=1;
    end
end
%error (index  solution  flag(1 for refine))
%% generator new mesh
mesh_new=[];
for cell=1:N-1;
    mesh_new=[mesh_new node(cell)];
    if error(cell,3)==1
        mesh_new=[mesh_new (node(cell)+node(cell+1))/2];
    end
end
mesh_new=[mesh_new node(end)];
info.mesh=mesh_new; 
info.N=length(mesh_new);
%% interpolate new solution; 
% plot(node(2:end-1), x_old  ,'*-');
% hold on;
new_solution = interp1(node(2:end-1),x_old,mesh_new(2:end-1),'spline');
new_solution=new_solution.';
% plot( mesh_new(2:end-1),new_solution, ':.');

