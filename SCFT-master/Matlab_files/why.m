% function f0=simple_FEM_1D_transient(yita)

N=32;
yita=ones(N,1)*1;
% yita=yita*-2000;
% yita=[0;yita;0];
% yita=[0;1;2;3;0]
% yita=yita_answer;
L=3.72374;
h=L/(N-1);
tau=0.5302;

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
for i=2:N-1
    A(i,i-1:i+1)=[h/6 2/3*h h/6];
end
A(1, 1:2)=[2/3*h h/6];
A(N,N-1:N)=[h/6 2/3*h];

%% B come with (df1, df2)
for i=2:N-1
    B(i,i-1:i+1)=[-1/h 2/h -1/h];
end
B(1, 1:2)=[2/h -1/h];
B(N,N-1:N)=[-1/h 2/h];

%% C come with (yita_i*f1,f2)
C=A;
D=A;
for i=1:N
    C(i,:)=C(i,:)*yita(i);
end

%% It is so important to choose implicit Euler method here,
%% explicit Euler results in a very small time step.
D=A+dt*(B+C);
    %% apply BC
    D(1,:)=0;
    D(N,:)=0;
    D(1,1)=1;
    D(N,N)=1;
for T=1:3
    b=A*xold;
    b(1)=0;
    b(N)=0;
    %% solving
    X=D\b;
    t=t+dt;
    solution_store(:,T)=[t;X];
    xold=X;
end
solution_store=[[0;x_initial] solution_store];

 