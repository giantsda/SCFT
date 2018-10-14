function out=simple_FEM_1D_transient(in)
%% for homogenous and inhomogenous mesh
global info;
N=info.N;
L=info.L;
tau=info.tau;
mesh=info.mesh;
in=[0; in; 0];
h=L/(N-1);

xold=ones(N,1);
xold(1)=0;
xold(end)=0;
x_initial=xold;
X=zeros(N,1);
t=0;
time_step=2048;
dt= 1/(time_step-1);
solution_store=zeros(N+1,time_step);

%% allcate memory
% A=sparse(N,N);
% B=sparse(N,N);
% C=sparse(N,N);
% D=sparse(N,N);
% b=zeros(N,1);

A=zeros(N,N);
B=zeros(N,N);
C=zeros(N,N);
D=zeros(N,N);
b=zeros(N,1);

%% A=(f1,f2) come with the modified diffusion term
for i=2:length(mesh)-1
    a1=mesh(i)-mesh(i-1);
    a2=mesh(i+1)-mesh(i);
    A(i,i)=a1/3+a2/3;
    A(i,i+1)=a2/6;
    A(i,i-1)=a1/6;
end
a2=mesh(2)-mesh(1);
A(1,1)=a2/3;
A(1,2)=a2/6;
a1=mesh(end)-mesh(end-1);
A(end,end)=a1/3;
A(end,end-1)=a1/6;


%% B come with (df1, df2)
for i=2:length(mesh)-1
    a1=mesh(i)-mesh(i-1);
    a2=mesh(i+1)-mesh(i);
    B(i,i)=1/a1+1/a2;
    B(i,i+1)=-1/a2;
    B(i,i-1)=-1/a1;
end
a2=mesh(2)-mesh(1);
B(1,1)=1/a2;
B(1,2)=-1/a2;
a1=mesh(end)-mesh(end-1);
B(end,end)=1/a1;
B(end,end-1)=-1/a1;

%% C come with (yita_i*f1,f2)
C=A;
D=A;
for i=1:N
    C(i,:)=C(i,:)*in(i);
end

%% It is so important to choose implicit Euler method here,
%% explicit Euler results in a very small time step.
D=A+dt*(B+C);
%% apply BC
D(1,:)=0;
D(N,:)=0;
D(1,1)=1;
D(N,N)=1;
D_inv=inv(D);
for T=1:time_step
    b=A*xold;
    b(1)=0;
    b(N)=0;
    %% solving
    X=D_inv*b;
    t=t+dt;
    solution_store(:,T)=[t;X];
    xold=X;
end
solution_store=[[0;x_initial] solution_store];

%% integrate for f0

f0=zeros(N+1,1);
for j=1:time_step
    value_left=solution_store(:,j).*solution_store(:,time_step-j+2);
    value_right=solution_store(:,j+1).*solution_store(:,time_step-(j+1)+2);
    f0=f0+0.5*(value_left+value_right)*dt;
end
f0(1)=[];
out=info.f0_given-f0;
out(1)=[];
out(end)=[];