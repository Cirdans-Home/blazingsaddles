function [x,it,r] = chebyshev(A, b, x0, maxit, tol, lMax, lMin)
%%CHEBYSHEV Iteration, useful for solving systems with the mass-matrix
% It can be modified to accomodate using a preconditioner.

theta = (lMax + lMin) / 2;
delta = (lMax - lMin) / 2;
x = x0;
r = b - A * x;
r0 = norm(r);
sigma = theta/delta;
rho = 1/sigma;
d = r/theta;
for it = 1:min(maxit,size(A,1))
    x = x + d;
    r = r - A*d;
    rhonew = 1/(2*sigma - rho);
    d = rhonew*rho*d + (2*rhonew/delta)*r;
    rho = rhonew;
    if (norm(r)/r0 < tol)
        % stop on tolerance with the relative residual
        break; 
    end 
end

end