%**********************************************************************;
% Project           : Iterative methods and preconditioning for large 
%			and sparse linear systems with applications
%
% Program name      : mgm_main.m
%
% Author            : F. Durastante -- fdurastante@uninsubria.it
%
% Date created      : 28 July 2017
%
% Purpose           : Multigrid calls
%
% Revision History  :
%
% Date        Author      	Ref    Revision (Date in DD/MM/YYYY format) 
% 28/07/2017  F. Durastante     1      File Created
%
%**********************************************************************;
function [ x,err ] ...
    = mgm_main( A,b,x,numiter,tol,a,v1,v2,omega1,omega2,gamma,xtrue )
%MGM_MAIN Solve Ax = b with V-cycle or with a W-cycle with gamma recursive
%calls, Damped Jacobi with omega as smoother
%   Ax = b linear system to be solved
%   numiter maximum number of iteration
%   tol tolerance
%   a projector mask
%   omega damping parameter for Jacobi
%   gamma number of recursive calls for the W-cycle
%   xtrue true solution



m = fix(log2(length(x)))-1; %Number of levels - 1,i.e., m = 1 TGM
%m = 1;
tol = tol*norm(b);

for i = 1:numiter
    x = mgm(A,b,x,a,omega1,omega2,gamma,v1,v2,0,m);
    r = b-A*x;
    err(i) = norm(x - xtrue);
    if norm(r)<tol
        break
    end
end

err = err/norm(xtrue); %relative residual: change if no xtrue available

end

