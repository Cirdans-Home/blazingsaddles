%% Build Mixed FEM Poisson
% This example uses iFEM, so we make sure of having it.
here = pwd;
cd('../ifem/');
setpath;
cd(here);
addpath('../bdm_mfem/');
clear; clc; close all;

%% Build the mesh
RefinementLevels = 2;
node = [-1 1; 0 1; 1 1; -0.5 0.5; 0.5 0.5; -1 0; 0 0; 1 0; ...
    -0.5 -0.5; 0.5 -0.5;-1.0 -1.0; 0.0 -1.0;1.0 -1.0];
elem =[4 2 1;4 1 6;4 6 7;4 7 2;5 3 2;5 2 7;5 7 8;5 8 3;9 7 6;...
    9 6 11;9 11 12;9 12 7;10 8 7;10 7 12;10 12 13;10 13 8];
bdEdge = [2 0 0;1 0 0; 0 0 0;0 0 0;2 0 0;0 0 0;0 0 0;1 0 0;...
    0 0 0;1 0 0;1 0 0;0 0 0;0 0 0;0 0 0;1 0 0;1 0 0];
for i=1:RefinementLevels
    [node,elem,bdEdge] = uniformrefine(node,elem,bdEdge);
end

[edge,elem2edge,signedge] = geomrelations(elem);
figure(1)
showmesh(node,elem);

%% Build the Saddle-Point system
NT = size(elem,1);      % Number of triangles
NE = size(edge,1);      % Number of edges
sol = zeros(2*NE+NT,1); % Space to store the solution

inva =1./exactalpha((node(elem(:,1))+node(elem(:,2))+node(elem(:,3)))/3);

%% Assemblying the matrix:
[a,b,area] = gradlambda(node,elem);
M = assemblebdm(NT,NE,a,b,area,elem2edge,signedge,inva);
[M,b,sol,freeDof] = rhside(node,elem,edge,bdEdge,area,M,sol,@f,@gD,@gN);

figure(2)
Mred = M(freeDof,freeDof);
A = Mred(1:2*NE,1:2*NE);
BT = Mred(1:2*NE,2*NE+1:end);
B = Mred(2*NE+1:end,1:2*NE);
C = Mred(2*NE+1:end,2*NE+1:end);
subplot(2,4,[1,2,5,6])
spy(Mred)
subplot(2,4,3)
spy(A)
subplot(2,4,4)
spy(BT)
subplot(2,4,7)
spy(B)
subplot(2,4,8)
spy(C)


%% Spectral Analysis
if length(freeDof) < 5000
    lambda = eig(M(freeDof,freeDof));
    mun = eigs(A,1,'smallestabs');
    mu1 = eigs(A,1,'largestabs');
    sigma1 = svds(BT,1,'largest');
    sigmam = svds(BT,1,'smallest');
    
    Iminus(1) = 0.5*(mun - sqrt(mun^2+4*sigma1^2));
    Iminus(2) = 0.5*(mu1 - sqrt(mu1^2+4*sigmam^2));
    Iplus(1) = mun;
    Iplus(2) = 0.5*(mu1 + sqrt(mu1^2 + 4*sigma1^2));
    
    figure(4)
    plot(1:length(freeDof),lambda,'x',...
        1:length(freeDof),Iminus(1)*ones(size(freeDof)),'k--',...
        1:length(freeDof),Iminus(2)*ones(size(freeDof)),'k--',...
        1:length(freeDof),Iplus(1)*ones(size(freeDof)),'r--',...
        1:length(freeDof),Iplus(2)*ones(size(freeDof)),'r--',...
        'LineWidth',2);
    xlabel('n')
    ylabel('\lambda_n')
    axis tight
    
end

%% Solving the linear system:
sol(freeDof) = M(freeDof,freeDof)\b(freeDof);
sigma = sol(1:2*NE);    u = sol(2*NE+1:end);

%% Visualizing the solution
figure(3)
showresult(node,elem,u);