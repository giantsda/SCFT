function set_f0_given(tau,L)
global info;
N=info.N;
mesh=info.mesh;
info.L=L;
info.tau=tau;
x_left=mesh(1:ceil(N*tau/L));
f0_given_left=(exp(4*tau*x_left./(tau*tau-x_left.*x_left))-1).^2./((exp(4*tau*x_left./(tau*tau-x_left.*x_left))+1).^2);
f0_given=ones(1,N);
f0_given(1:length(f0_given_left))=f0_given_left;
f0_given(end-length(f0_given_left)+1:end)=fliplr(f0_given_left);
f0_given(isnan(f0_given)) = 1;
f0_given=f0_given.';
info.f0_given=f0_given;