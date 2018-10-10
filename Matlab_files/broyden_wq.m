function [Xr, found] = broyden_wq(F, Xr, dlt, ee, imax, n)
% Solving n nonlinear equations F(Xr)=0 using Broyden's method. The 
% function handle F(Xr) needs to be created in the calling program. dlt is
% the small step size used in the first-order forward finite difference as 
% the initial B0; the identity matrix is used as the initial B0 if dlt=0.
% ee is the convergence criterion. imax is the maximum number of iterations.
% Xr in the input is the initial guess, and in the output is the solution. 
% "found" returns the status: 0 -- error, 
%                             1 -- max|F(Xr)| < ee, 
%                             2 -- reached maximum number of iterations.

% Output header
fprintf ('iter            x                     y            max|F(Xr)|\n');

% Check the initial guess for a solution
iter = 0;
found = 0;
Y = F(Xr);
err = 0;
for i = 1:n
    if (abs(Y(i)) > err)
        err = abs(Y(i));
    end
end
fprintf (' %2d ', iter);
for i = 1:n
    fprintf ('  %.14e', Xr(i));
end
fprintf ('  %e\n', err);
if (err < ee)
    found = 1;
else
    % Calculate the inverse of the initial B0, stored as B1
    if (dlt == 0)
        B1 = eye(n);
    else
        for j = 1:n
            Xr(j) = Xr(j) + dlt;
            Yn = F(Xr);
            for i = 1:n
                B0(i,j) = (Yn(i) - Y(i)) / dlt;
            end
            Xr(j) = Xr(j) - dlt;
        end
%  output2d (B0, '%.14e');
        [B1] = inverse (B0, 1e-14);
%         if err,  exit,  end
    end
end

while (found==0 && iter<=imax)
    % Calculate the next Xr
    dX = -B1 * Y;
    Xr = Xr + dX;
    iter = iter + 1;

    % Check for a solution
    Yn = F(Xr);
    err = 0;
    for i = 1:n
        if (abs(Yn(i)) > err)
            err = abs(Yn(i));
        end
    end
    fprintf (' %2d ', iter);
    for i = 1:n
        fprintf ('  %.14e', Xr(i));
    end
    fprintf ('  %e\n', err);
    if (err < ee)
        found = 1;
        fprintf("broydn_wq: solved successfully!\nThe solution is:\n");
        fprintf('%2.15f\n',Xr);
        break;
    elseif (iter == imax)
        found = 2;
        break;
    end        
    
    % Calculate the next B1 using the Shermann-Morrison formula
    dF = Yn - Y;
    B1dF = B1 * dF;
    dXB1dF = dX' * B1dF;
    if (abs(dXB1dF) > 1e-14)
        B1 = B1 + (dX-B1dF) * dX' * B1 / dXB1dF;
    else
        break;
    end
    Y = Yn;
%     pause()
end

 
