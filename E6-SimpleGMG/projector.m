%**********************************************************************;
% Project           : Iterative methods and preconditioning for large 
%			and sparse linear systems with applications
%
% Program name      : projector.m
%
% Author            : F. Durastante -- fdurastante@uninsubria.it
%
% Date created      : 28 July 2017
%
% Purpose           : Projector for the MGM routine
%
% Revision History  :
%
% Date        Author      	Ref    Revision (Date in DD/MM/YYYY format) 
% 28/07/2017  F. Durastante     1      File Created
%
%**********************************************************************;
function [ P ] = projector( a,n )
%PROIETTORE P = T*K', where K is the downsampling matrices and T is
%the Toeplitz matrix built from a
%   a = mask for the projector
%   n = size of the projector
if mod(n,2) == 0
    error('n needs to be odd!')
end
t = zeros(n,1);
t(1:length(a))=a(:);
T = toeplitz(t);
P = sparse(T(:,2:2:end-1)); 
end

