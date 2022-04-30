%**********************************************************************;
% Project           : Iterative methods and preconditioning for large
%			and sparse linear systems with applications
%
% Program name      : ex_toepmultigrid.m
%
% Author            : F. Durastante -- fdurastante@uninsubria.it
%
% Date created      : 28 July 2017
%
% Purpose           : Example of Toeplitz multigrid for 1D problems
%
% Revision History  :
%
% Date        Author      	Ref    Revision (Date in DD/MM/YYYY format)
% 28/07/2017  F. Durastante     1      File Created
%
%**********************************************************************;
% Test for the Toeplitz Multigrid Method (1D)

clear;
clc;

fprintf('Select the problem:\n');
fprintf('1) f(x) = 2-2 cos(x)\n');
fprintf('2) f(x) = 6-6 cos(x)\n');
fprintf('3) f(x) =  x^4 + 1\n');
fprintf('4) f(x) = x^4 + x^2 + 1\n');
fprintf('5) f(x) = x^2\n');
fprintf('6) f(x) = x^4\n');
fprintf('7) f(x) = |x|^3\n');
problem = input('problem = ');
fprintf('Select the interpolation operator: \n');
fprintf('1) Linear interpolation\n')
fprintf('2) Bilinear interpolation\n')
liniterp = input('Interpolator (1/2) = ');
switch liniterp
    case 1
        a = [2,1];
    case 2
        a = [6 4 1];
    otherwise
        a = [2,1];
end
for k = 6:12
    n = 2^k - 1;    % Odd power at each level
    
    if problem == 1
        t = zeros(1,n);
        t(1:2) = [2 -1];
        maxf = 4;
    elseif problem == 2
        t = zeros(1,n);
        t(1:3) = [6 -4 1]/8;
        maxf = 3/2;
    elseif problem == 3
        t = zeros(n,1);
        t(1) = 1+(pi^4)/5;
        for j=2:n
            t(j) = (4*(-1)^(j-1))*(-6 + (j-1)^2*pi^2)/((j-1)^4);
        end
        maxf = pi^4+1;
    elseif problem == 4
        t = zeros(n,1);
        t(1) = 1 + pi^(2)/3 + pi^(4)/5;
        for j=2:n
            t(j) = (2*(-1)^(j-1))*(-12 + (j-1)^2*(1+2*pi^2))/((j-1)^4);
        end
        maxf = 1+pi^2+pi^4;
    elseif problem == 5
        t = toepgenerator(n,6);
        maxf = pi^2;
    elseif problem == 6
        t = toepgenerator(n,8);
        maxf = pi^4;
    elseif problem == 7
        t = toepgenerator(n,10);
        maxf = pi^3;
    else
        error('Problem unrecognized!');
    end
    
    % Building matrix and reference solution:
    tt = linspace(0,pi,n)';
    x = sin(4*tt) + 2*sin(50*tt) - 1; % Building frequences
    A = sparse(toeplitz(t));
    b = A*x;
    
    %The Multigrid Method
    
    xtrue = x;
    omega1 = t(1)/maxf;
    omega2 = 2*t(1)/maxf;
    v1 = 2;
    v2 = 2;
    
    gamma = 2;
    numiter = 1000;
    tol = 1e-6;
    x0 = zeros(size(x));
    tic;
    [x,err] = mgm_main(A,b,x0,numiter,tol,a,v1,v2,omega1,omega2,gamma,xtrue);
    time = toc;
    fprintf('Size %d Iteration %d Time %1.2e \n',n,length(err),time);
    figure(1)
    semilogy(err,'--ko');
    title(sprintf('Converged in %d it',length(err)));
    xlabel('Iteration');
    ylabel('Error');
end
