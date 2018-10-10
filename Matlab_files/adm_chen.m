function [x_old,found] = adm_chen (F, n, x_old, tol, maxIteration,lmd)
% lmd = 0.9;
lk = lmd;
found=0;
nm=min(30,n);

kGlobal = 0;
X_end = x_old;
err=Inf;

Y = F(X_end);
err=max(abs(Y));
fprintf("adm_chen: Before solving, error=%2.15f \n",err);
if err>1
    fprintf("adm_chen: Before solving, error>1 Continue? \n",err);
    pause();
end

while (err>tol && kGlobal<=maxIteration)
    [X_end, kLocal]=mixing(F,n,X_end,tol,kGlobal,maxIteration);
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

function [X_end, kLocal]=mixing(F,n,x_old,tol,kGlobal,maxIteration)
lmd = 0.9;
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
        fprintf('adm iteration: %d, error: %.14e\n',kGlobal+kLocal-1,err);
    end
    
    if (err < tol)
        found = 1;
        X_end = X(:,kLocal);
        fprintf('*****And_chen: Solved equation successfully!*****\nThe solution is:\n');
        fprintf('%2.15f\n',X(:,kLocal))
        return;
    end
    
    if (err > 1e20)
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
    
    if rcond(U)<1e-10
        fprintf("And_chen: Singular Matrix detected And_chen restarted!\n");
        X_end=X(:,kLocal);
        return;
    end
    
    if kLocal==1
        c=0;
    else
        c=U\V;
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

