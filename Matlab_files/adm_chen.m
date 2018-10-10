function [x_old,found] = adm_chen (F, n, x_old, tol, maxIteration,lmd_in)
found=0;
kGlobal = 0;
X_end = x_old;

Y = F(X_end);
err=max(abs(Y));
fprintf("adm_chen:n=%d Before solving, error=%2.15f \n",n,err);
% if err>1
    fprintf("adm_chen: Before solving, error>1 Continue? \n");
    pause();
% end

while (err>tol && kGlobal<=maxIteration)
    [X_end, kLocal]=mixing(F,n,X_end,tol,kGlobal,maxIteration,lmd_in);
    kGlobal=kGlobal+kLocal-1;
    Y = F(X_end);
    err=max(abs(Y));
end

if err<tol
    found=1;
else
    fprintf("And_chen failed after %d iterations :(  Try to increase max iteration allowed\n",maxIteration);
end
x_old=X_end;

function [X_end, kLocal]=mixing(F,n,x_old,tol,kGlobal,maxIteration,lmd_in)
lmd = lmd_in;
lk = lmd;
found=0;
nm=min(30,n);
kLocal = 1;
X(:,kLocal) = x_old;
err=Inf;
U=1;
while (err>tol && kLocal+kGlobal-1<=maxIteration)
    
    Y(:,kLocal) = F(X(:,kLocal));
    err=max(abs(Y(:,kLocal)));
    if kLocal>1
        fprintf('adm iteration: %d,n=%d, error: %.14e\n',kGlobal+kLocal-1,n,err);
    end
    
    if (err < tol)
        found = 1;
        X_end = X(:,kLocal);
        fprintf('*****And_chen: Solved equation successfully!*****\nThe solution is:\n');
        fprintf('%2.15f\n',X(:,kLocal))
        return;
    end
    
    if (err > 1000)
        explode=1
    end
    % Calculate the matrix U and the column vector v
    if (kLocal-1 <= nm)
        m = kLocal-1;
    else
        m = nm;
    end
    
    for i = 1:m
        V(i,1) = dot(Y(:,kLocal) - Y(:,kLocal-i),Y(:,kLocal));
        for j = 1:m
            U(i,j) =dot(Y(:,kLocal) - Y(:,kLocal-i),Y(:,kLocal) - Y(:,kLocal-j));
        end
    end
      
        % Calculate c = U^(-1) * v using LU decomposition
    if (m > 0)
        [U, o, err] = Doolittle (U, 1e-10);
        if ~err
            [c] = Doolsub (U, V, o, m);     % c is a row vector
        else
            fprintf("And_chen: Singular Matrix detected And_chen restarted!\n");
        X_end=X(:,kLocal);
        return;
        end
    end
    
    
    % Calculate the next x^(k)
    for i = 1:n
        cx = 0;
        cd = 0;
        for j = 1:m
            cx = cx + c(j) * (X(i,kLocal-j) - X(i,kLocal));
            cd = cd + c(j) * (Y(i,kLocal-j) - Y(i,kLocal));
        end
        X(i,kLocal+1) = X(i,kLocal) + cx + (1-lk)*(Y(i,kLocal)+cd);
    end
    
    kLocal = kLocal + 1;
    lk = lk * lmd;
end

X_end=X(:,kLocal);

