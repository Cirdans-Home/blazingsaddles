function [x,it,r] = chebyshev(A, b, x0, maxit, tol, lMax, lMin)
%%CHEBYSHEV Iteration, useful for solving systems with the mass-matrix
% It can be modified to accomodate using a preconditioner.

d = (lMax + lMin) / 2;
c = (lMax - lMin) / 2;
x = x0;
r = b - A * x;
r0 = norm(r);
for it = 1:min(maxit,size(A,1))
    z = r; % If you want to implement a preconditioner you have to use it
    % here!
    if (it == 1)
        p = z;
        alpha = 1/d;
    elseif (it == 2)
        beta = (1/2) * (c * alpha)^2;
        alpha = 1/(d - beta / alpha);
        p = z + beta * p;
    else
        beta = (c * alpha / 2)^2;
        alpha = 1/(d - beta / alpha);
        p = z + beta * p;
    end
    
    x = x + alpha * p;
    r = b - A * x; % Or equivalently: r = r - alpha * A * p
    if (norm(r)/r0 < tol)
        % stop on tolerance with the relative residual
        break; 
    end 
end

end