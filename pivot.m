function [o] = pivot (A, o, s, n, k)
% Implicit and in-place partial pivoting for row k of A.

% Find the row number p of the pivot element "big"
p = k;
big = abs(A(o(k),k)/s(o(k)));
for i = k+1 : n
    dummy = abs(A(o(i),k))/s(o(i));
    if dummy > big
        p = i;
        big = dummy;
    end
end

if p ~= k
    dummy = o(p);
    o(p) = o(k);
    o(k) = dummy;
end
