function [a, o, err] = Doolittle (a, tol)
% Doolittle decompsition of square matrix A with implicit and in-place pivoting.
% tol specifies the smallest acceptable absolute value of scaled pivot elements.
% err=0 for successful, 1 for too small absolute value of scaled pivot element.

err = 0;
n = size(a,1);  % Get the size of A
for i = 1:n
    o(i) = i;   % initial value for the order vector used for in-place pivoting
    % Find s(i), the largest absolute value in row i of A (used for scaling)
    s(i) = abs(a(i,1));
    for j = 2:n
        if abs(a(i,j)) > s(i)
            s(i) = abs(a(i,j));
        end
    end
end

% Gauss elimination. Note that the 0's below the diagonal of A (with pivoting)
% after elimination are not calculated or stored.
for k = 1 : n-1         % the k^th row
    [o] = pivot (a, o, s, n, k);
    if abs(a(o(k),k))/s(o(k)) < tol
        err = 1;
        fprintf ('Scaled pivot %i is %e < %e !\n', k, abs(a(o(k),k))/s(o(k)), tol);
        break;
    end
    for i = k+1 : n         % eliminate row i of A
        factor = a(o(i),k) / a(o(k),k);
        a(o(i),k) = factor;     % Store the L element
        for j = k+1 : n
            a(o(i),j) = a(o(i),j) - factor*a(o(k),j);
        end
    end
end
if ~err
    if abs(a(o(n),n))/s(o(n)) < tol
        err = 1;
        fprintf ('Scaled pivot %i is %e < %e !\n', n, abs(a(o(n),n))/s(o(n)), tol);
    end
end
