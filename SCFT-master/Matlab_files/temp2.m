%% solve

% N=300;
% yita=ones(N,1)*1
% f0=simple_FEM_1D_transient(yita)
% plot(f0)
function temp2
 
% %   [Xr, found] = AndMix (@myfun, 1.5, 0.00001, 20, 1)
% %   % Solving n nonlinear equations F(Xr)=Xr using Anderson mixing. The 
% % % function handle F(Xr) needs to be created in the calling program. 
% % % ee is the convergence criterion. imax is the maximum number of iterations.
% % % Xr in the input is the initial guess, and in the output is the solution. 
% % % "found" returns the status: 0 -- error, 
% % %                             1 -- max|F(Xr)| < ee, 
% % %                             2 -- reached maximum number of iterations.
% 
%  [Xr, found] = broyden (@myfun, 1.5, 0.0001, 0.000001, 200, 1)
% % Solving n nonlinear equations F(Xr)=0 using Broyden's method. The 
% % function handle F(Xr) needs to be created in the calling program. dlt is
% % the small step size used in the first-order forward finite difference as 
% % the initial B0; the identity matrix is used as the initial B0 if dlt=0.
% % ee is the convergence criterion. imax is the maximum number of iterations.
% % Xr in the input is the initial guess, and in the output is the solution. 
% % "found" returns the status: 0 -- error, 
% %                             1 -- max|F(Xr)| < ee, 
% %                             2 -- reached maximum number of iterations.


N=300;
yita=ones(N,1)*-10;
re=myfun(yita)
plot(re)
% x = fsolve(@myfun,yita)
[x, ithist] = broyden( @myfun, yita )

function out=myfun(yita)
f0=simple_FEM_1D_transient(yita);
tau=0.2;
N=300;
x=linspace(0,1,N);
x_left=x(1:N*tau);
f0_given_left=(exp(4*tau*x_left./(tau*tau-x_left.*x_left))-1).^2./((exp(4*tau*x_left./(tau*tau-x_left.*x_left))+1).^2);
f0_given=ones(1,N);
f0_given(1:length(f0_given_left))=f0_given_left;
f0_given(end-length(f0_given_left)+1:end)=fliplr(f0_given_left);
f0_given(isnan(f0_given)) = 1;
f0_given=f0_given.';
out=f0_given-f0;
norm(out)
% pause()