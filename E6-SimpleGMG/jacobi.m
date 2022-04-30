%**********************************************************************;
% Project           : Iterative methods and preconditioning for large 
%			and sparse linear systems with applications
%
% Program name      : jacobi.m
%
% Author            : F. Durastante -- fdurastante@uninsubria.it
%
% Date created      : 28 July 2017
%
% Purpose           : Jacobi smoother
%
% Revision History  :
%
% Date        Author      	Ref    Revision (Date in DD/MM/YYYY format) 
% 28/07/2017  F. Durastante     1      File Created
%
%**********************************************************************;
function [ x ] = jacobi( A,b,x,omega,numiter )
A = sparse(A);
D = diag(A);
for i=1:numiter
    x = x + omega*(b - A*x)./D;
end
end 

