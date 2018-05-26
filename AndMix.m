function [Xr, found] = AndMix (F, Xr, ee, imax, n)
% Solving n nonlinear equations F(Xr)=Xr using Anderson mixing. The 
% function handle F(Xr) needs to be created in the calling program. 
% ee is the convergence criterion. imax is the maximum number of iterations.
% Xr in the input is the initial guess, and in the output is the solution. 
% "found" returns the status: 0 -- error, 
%                             1 -- max|F(Xr)| < ee, 
%                             2 -- reached maximum number of iterations.

% Output header
fprintf ('iter            x                     y            max|F(Xr)|\n');

lbd = 0.9;
lk = lbd;

if (n >= 30)
    nm = 30;
else
    nm = n;
end

k = 1;
X(:,k) = Xr;    % x^(0)
found = 0;
while (found==0 && k<=imax)
    % Calculate f^(k) and check for a solution
    Y(:,k) = F(X(:,k));
    err = 0;
    for i = 1:n
        if (abs(Y(i,k)) > err)
            err = abs(Y(i,k));
        end
    end
    fprintf (' %2d ', k-1);
    for i = 1:n
        fprintf ('  %.14e', X(i,k));
    end
    fprintf ('  %e\n', err);
    if (err < ee)
        found = 1;
        Xr = X(:,k);
        break;
    end
    
    % Calculate the matrix U and the column vector v
    if (k-1 <= nm)
        m = k-1;
    else
        m = nm;
    end
    for i = 1:m
        dltY = (Y(:,k) - Y(:,k-i))';    % a row vector
        v(i) = dltY * Y(:,k);
        for j = 1:m
            U(i,j) = dltY * (Y(:,k) - Y(:,k-j));
        end
    end

    % Calculate c = U^(-1) * v using LU decomposition
    if (m > 0)
        [U, o, err] = Doolittle (U, 1e-10);
        if ~err
            [c] = Doolsub (U, v, o, m);     % c is a row vector
        else
            break;
        end
    end
        
    % Calculate the next x^(k)
    for i = 1:n
        cx = 0;
        cd = 0;
        for j = 1:m
            cx = cx + c(j) * (X(i,k-j) - X(i,k));
            cd = cd + c(j) * (Y(i,k-j) - Y(i,k));
        end
        X(i,k+1) = X(i,k) + cx + (1-lk)*(Y(i,k)+cd);
    end
    
    k = k + 1;
    lk = lk * lbd;
end    

if (found==0 && k>imax)
    Xr = X(:,k);
    found = 2;
    dltY = F(Xr);
    err = 0;
    for i = 1:n
        if (abs(dltY(i)) > err)
            err = abs(dltY(i));
        end
    end
    fprintf (' %2d ', k-1);
    for i = 1:n
        fprintf ('  %.14e', Xr(i));
    end
    fprintf ('  %e\n', err);
end
Xr=X(:,end);
