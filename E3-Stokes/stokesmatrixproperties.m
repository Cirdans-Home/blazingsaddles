%% Test of the properties of the Stokes matrix for a stable scheme
% This example uses iFEM, so we make sure of having it.
here = pwd;
cd('../ifem/');
setpath;
cd(here);

clear; clc; close all;

%% Setup

maxIt = 4;
lambdamin = zeros(maxIt,1);
lambdamax = zeros(maxIt,1);
h = zeros(maxIt,1);

%% Generate initial mesh
[node,elem] = squaremesh([0 1 0 1], 0.25);
bdFlag = setboundary(node,elem,'Dirichlet');

%% PDE and options
pde = Stokesdata1;
options.solver = 'none';

%% Computing spectra and h
for k = 1:maxIt
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    [soln,eqn] = StokesP2P1(node,elem,bdFlag,pde,options);
    A = eqn.A(eqn.ufreeDof,eqn.ufreeDof);
    lambdamin(k) = eigs(A,1,'smallestreal','MaxIterations',10000);
    lambdamax(k) = eigs(A,1,'largestreal','MaxIterations',10000);
    h(k) = 1./(sqrt(size(node,1))-1);
end

disp([h.^2,lambdamin,h.^2./lambdamin,lambdamax])

%% Inf-Sup Condition
% Build the Schur complement

[node,elem] = squaremesh([0 1 0 1], 0.50);
bdFlag = setboundary(node,elem,'Dirichlet');
maxIt = 4;
bound = zeros(maxIt,1);
maxval = zeros(maxIt,1);
for k = 1:maxIt
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    [soln,eqn] = StokesP2P1(node,elem,bdFlag,pde,options);
    % Build the Schur-complement
    A = eqn.A(eqn.ufreeDof,eqn.ufreeDof);
    Np = size(node,1);
    M = [eqn.A,eqn.B';eqn.B,sparse(Np,Np)];
    B = M(eqn.ufreeDof(end)+1:end,eqn.ufreeDof);
    S = B*(A\B');
    % Build the mass matrix
    [~,Q,area] = assemblematrix(node,elem);
    lval = eigs(S,Q,2,'smallestreal');
    maxval(k) = eigs(S,Q,1,'largestreal','MaxIterations',10000);
    h(k) = 1./(sqrt(size(node,1))-1);
    bound(k) = lval(2);
 end




