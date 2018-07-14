function [x] = backsub (A, b, o, n)
% Back substitution with in-place pivoting.

for i = n : -1 : 1
    sum = 0;
    for j = i+1 : n
        sum = sum + A(o(i),j)*x(j);
    end
    x(i) = (b(o(i))-sum) / A(o(i),i);
end