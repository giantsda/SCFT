function [x, err] = gauss (A, b, tol)
% Gauss elimination solving Ax=b with implicit and in-place pivoting.
% tol specifies the smallest acceptable absolute value of scaled pivot elements.
% err=0 for successful, 1 for too small absolute value of scaled pivot element.

err = 0;
n = size(A,1);  % Get the size of A
for i = 1:n
    o(i) = i;   % Initialize the order vector used for in-place pivoting
    % Find s(i), the largest absolute value in row i of A (used for scaling)
    s(i) = abs(A(i,1));
    for j = 2:n
        if abs(A(i,j)) > s(i)
            s(i) = abs(A(i,j));
        end
    end
end

% Gauss elimination. Note that the 0's below the diagonal of A (with pivoting) 
% after elimination are not calculated or stored.
for k = 1 : n-1         % the k^th row
    [o] = pivot (A, o, s, n, k);
    if abs(A(o(k),k))/s(o(k)) < tol
        err = 1;
        fprintf ('Scaled pivot %i is %e < %e !\n', k, abs(A(o(k),k))/s(o(k)), tol);
        break;
    end
    for i = k+1 : n         % eliminate row i of A and b
        factor = A(o(i),k) / A(o(k),k);
        for j = k+1 : n
            A(o(i),j) = A(o(i),j) - factor*A(o(k),j);
        end
        b(o(i)) = b(o(i)) - factor*b(o(k));
    end
end
if ~err
    if abs(A(o(n),n))/s(o(n)) < tol
        err = 1;
        fprintf ('Scaled pivot %i is %e < %e !\n', n, abs(A(o(n),n))/s(o(n)), tol);
    end
end

if ~err     % Back substitution
    [x] = backsub (A, b, o, n);
end
