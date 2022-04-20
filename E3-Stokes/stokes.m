%% Stokes problem
% Some example with ifem for the Stokes problem.
% This example uses iFEM, so we make sure of having it.
here = pwd;
cd('../ifem/');
setpath;
cd(here);

clear; clc; close all;

%% Building the mesh
RefinementLevels = 2;
square = [-1,1,-1,1];
h = 0.25;
[node,elem] = squaremesh(square,h);
figure(1)
subplot(1,2,1);
showmesh(node,elem);
for i=1:RefinementLevels
    [node,elem] = uniformrefine(node,elem);
end
subplot(1,2,2);
showmesh(node,elem);
title(sprintf("Refinement Levels %d",RefinementLevels));

%% Building the test problem: colliding flows
bdFlag = setboundary(node,elem,'Dirichlet');
pde = Stokesdata1;
options.solver='direct';  % We just perform the build, don't solve
[soln,eqn,info] = StokesP2P1(node,elem,bdFlag,pde,options);

%% Solving phase
[elem2dof,edge,bdDof] = dofP2(elem);
N = size(node,1);  NT = size(elem,1);  Nu = N+size(edge,1);   Np = N;
A = eqn.A;
B = eqn.B;

freedof = [eqn.ufreeDof;2*Nu+eqn.pDof];
C = sparse(size(B,1),size(B,1));
M = [A,B';B,C];
b = [eqn.f;eqn.g];

x = [soln.u;soln.p];
x(freedof) = M(freedof,freedof)\b(freedof);
uh = x(1:2*Nu);
ph = x(2*Nu+1:end);

uh = reshape(uh,Nu,2);

figure(2)
subplot(1,3,1);
showsolution(node,elem,uh(:,1))
axis square
subplot(1,3,2);
showsolution(node,elem,uh(:,2))
axis square
subplot(1,3,3);
showsolution(node,elem,ph)
axis square

figure(2)
subplot(2,4,[1,2,5,6])
spy(M(freedof,freedof))
subplot(2,4,3)
spy(M(eqn.ufreeDof,eqn.ufreeDof))
subplot(2,4,4)
spy(M(eqn.ufreeDof,eqn.ufreeDof(end)+1:end))
subplot(2,4,7)
spy(M(eqn.ufreeDof(end)+1:end,eqn.ufreeDof))
subplot(2,4,8)
spy(M(eqn.ufreeDof(end)+1:end,eqn.ufreeDof(end)+1:end))

figure(3)
spy(M(freedof,freedof))