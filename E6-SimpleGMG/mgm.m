%**********************************************************************;
% Project           : Iterative methods and preconditioning for large 
%			and sparse linear systems with applications
%
% Program name      : mgm.m
%
% Author            : F. Durastante -- fdurastante@uninsubria.it
%
% Date created      : 28 July 2017
%
% Purpose           : One cycle of multigrid (either V-- or W--)
%
% Revision History  :
%
% Date        Author      	Ref    Revision (Date in DD/MM/YYYY format) 
% 28/07/2017  F. Durastante     1      File Created
%
%**********************************************************************;
function [x] = mgm(A,b,x,a,omega1,omega2,gamma,v1,v2,liv,m)
%MGM One Multigrid Cycle for Ax = b
%   a = mask of the projector
%   omega = damped jacobi parameter
%   v1 = iterations of the pre-smoother
%   v2 = iterations of the post-smoother
%   liv = actual level
%   m = overall number of levels - 1


if liv == m
    x = A\b;
    %[x,~] = bicgstab(A,b);
else
    n = length(b); % size at current level
    x = jacobi(A,b,x,omega1,v1); %Pre-smoother
    r = b - A*x; %Residuo
    P = projector( a,n ); % residual at current level
    r = P'*r; % projected residual
    B = P'*A*P;
    e = zeros(size(r)); % error initialization on smaller grid
    for i = 1:gamma
        e = mgm(B,r,e,a,omega1,omega2,gamma,v1,v2,liv+1,m);
    end
    x = x + P*e;
    x = jacobi(A,b,x,omega2,v2); %Post-smoother
end

