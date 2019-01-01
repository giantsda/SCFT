% x,y is the table waitting to be interpolated.
% xp is the x vector you want to know yp at.
% if not spcify m, use not-a-knot end condition
% m=0 for natural spline;otherwise m is the ddy at boundarys
function yp=spline_chen(x,y,xp,m)
Nx=length(x);
Nxp=length(xp);
A=zeros(Nx);
b=zeros(Nx,1);
for i=2:Nx-1
    A(i,i-1:i+1)=[(x(i)-x(i-1))/6 (x(i+1)-x(i-1))/3 (x(i+1)-x(i))/6];
    b(i)=(y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1));
end

%% not-a-knot end condition
if (nargin==3)
    if Nx<=3
         error('spline_chen:To use not-a-knot end condition, more than 3 x-y should be provided!\n');
    end
    A(1,1:3)=[1/(x(2)-x(1)) -1/(x(2)-x(1))-1/(x(3)-x(2)) 1/(x(3)-x(2))];
    A(Nx,Nx-2:Nx)=[1/(x(Nx-1)-x(Nx-2)) -1/(x(Nx-1)-x(Nx-2))-1/(x(Nx)-x(Nx-1)) 1/(x(Nx)-x(Nx-1))];
    b(1)=0;
    b(Nx)=0;
elseif (nargin==4)
    A(1,1)=1;   
    A(Nx,Nx)=1;
    % A(1,2)=-1*m;
    % A(Nx-1,Nx)=-1*m;
    b(1)=m;
    b(Nx)=m;
else
    error('spline_chen: Too many input args!\n');
end






X=A\b; % X is the solution of ddy at each point of x;

yp=zeros(Nxp,1);
for i=1:Nxp 
    klo=1;
    khi=Nx;
    while ((khi-klo)>1)
        k=floor((khi+klo)/2);
        if (x(k)>xp(i))
            khi=k;
        else
            klo=k;
        end
    end
    h=x(khi)-x(klo);
    if (h==0)
        error('spline_chen: Bad x input. x(%d)=x(%d)=%f',khi,klo,x(khi));
    end
    a=(x(khi)-xp(i))/h;
    b=(xp(i)-x(klo))/h;
    yp(i)=a*y(klo)+b*y(khi)+((a*a*a-a)*X(klo)+(b*b*b-b)*X(khi))*h*h/6;
end

