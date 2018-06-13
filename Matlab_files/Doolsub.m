function [b] = Doolsub (a, b, o, n)
% Forward and back substitution after Doolittle decomposition with pivoting
% to calculate A^(-1)*b, which is returned in b.

% Calculate L^(-1)*b and store in b
for i = 2:n
    sum = b(o(i));
    for j = 1 : i-1
        sum = sum - a(o(i),j)*b(o(j));
    end
    b(o(i)) = sum;
end

% Calculate U^(-1)*b and store in b 
[b] = backsub (a, b, o, n);
