function output2d (a, fmt)
% Output the 2D array A(m,n) using the format fmt for each element.

[m,n] = size (a);   % Get m and n
for i = 1:m
    for j = 1:n
        fprintf (fmt, a(i,j));  fprintf (' ');
    end
    fprintf ('\n');
end
