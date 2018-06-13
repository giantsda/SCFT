M=10;
h=1/(M-1);

%% A come with (df1, df2)
A=zeros(M);
for i=2:M-1
    A(i,i-1:i+1)=[-1/h 2/h -1/h];
end
A(1, 1:2)=[2/h -1/h];
A(M,M-1:M)=[-1/h 2/h];

%% B=(f1,f2) come with the modified diffusion term
B=zeros(M);
for i=2:M-1
    B(i,i-1:i+1)=[h/6 2/3*h h/6];
end
B(1, 1:2)=[2/3*h h/6];
B(M,M-1:M)=[h/6 2/3*h];
% A=A+B;

b=ones(M,1)*h;
% b=zeros(M,1);
b(1)=0;
b(M)=0;
A(1,:)=0;
% A(:,1)=0;
A(M,:)=0;
% A(:,M)=0;
A(1,1)=1;
A(M,M)=1;

x=A\b;
t=linspace(0,1,M)
plot(t,x)