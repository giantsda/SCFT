function [A1] = inverse (A, tol)
% Calculate the inverse A1 of a square matrix A using Gauss elimination with pivoting.
% tol specifies the smallest acceptable absolute value of scaled pivot elements.

n = size(A,1);
A1 = eye(n);            % Initialize A1 as the identity matrix

err = 0;
for i = 1:n
    [A1(:,i), err] = gauss (A, A1(:,i), tol);
    if err
        fprintf ('Scaled pivot element < %e\n', tol);
        break;
    end
end
