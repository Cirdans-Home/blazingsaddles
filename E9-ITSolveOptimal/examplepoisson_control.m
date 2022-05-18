%% Optimal Control of Poisson Equation with unbounded constraints
% This example needs IFISS so we look for it.
clear; clc; close all;
try
    here = pwd;
    ifiss
    cd(here)
catch
    warning("This example needs IFISS to be installed to run, I'm adding it to the path");
    here = pwd;
    cd ../ifiss3.6
    setpath
    cd(here)
end

clear variables

% Declare global variables for scalar and vector problems
global amg_grid amg_smoother number_vcycles

%% Build the grid
grid_type = 1;
square_type = 2;
% ----- Grid parameters ------------------------------------------------- %
nc = 6;
alpha = 1e-5;
qmethod = 2; % Q1 = 1, Q2 = 2
%% Build the discrete operators
% How large is resulting stiffness or mass matrix?
np = (2^nc+1)^2; % entire matrix system is thus 3*np-by-3*np in size

% Compute matrices specifying location of nodes
[x,y,xy,bound,mv,mbound] = square_domain_x(nc,grid_type);

% Compute connectivity, stiffness and mass matrices, for Q1 or Q2 elements
if qmethod==1
    [ev,ebound] = q1grid(xy,mv,bound,mbound);
    [K,M] = femq1_diff(xy,ev);
elseif qmethod==2
    [x,y,xy] = q2grid(x,y,xy,mv,bound); % to generate plot of Q2 elements
    [K,M] = femq2_diff(xy,mv);
else
    error('illegal parameter choice ''qmethod'', try again.')
end
[yhat_vec,bc_nodes] = poissoncontrol_rhs_ex1(xy,bound,square_type);
%% Assemble the Saddle-Point matrix
% Initialise RHS vector corresponding to desired state
Myhat = M*yhat_vec;
% Enforce Dirichlet BCs on state, and zero Dirichlet BCs on adjoint
[K,d] = nonzerobc_input(K,zeros(np,1),xy,bound,bc_nodes);
[M,Myhat] = nonzerobc_input(M,Myhat,xy,bound,bc_nodes);
% Saddle-point:
Matrix = [M sparse(np,np) K; ...
    sparse(np,np) alpha*M -M; ...
    K -M sparse(np,np)];
rhs = [Myhat; zeros(np,1); d];

figure(1)
spy(Matrix);

%% Solve the system
x_it = Matrix\rhs;
sol_y = x_it(1:np); sol_u = x_it(np+1:2*np);

% Make contour and surface plots of solutions for state y and control u
solplot_poissoncontrol(sol_y,sol_u,xy,x,y,2)

%% Iterative Solution Algorithms
% We test here different solution strategies

maxit = 100;
tol = 1e-7;

%% No preconditioner
MA = 'mass_identity';
MS = 'schur_identity';
tic %%start timing
[x_it,resvec,iter,flag] = poissoncontrol_minres(M,K,alpha,Myhat,d,maxit,tol,MA,[],MS,[]);
etoc = toc;
results{1,1} = resvec;
results{1,2} = iter;
results{1,3} = 'No Preconditioner';
results{1,4} = etoc;
%whatdidiget

%% Diagonal Preconditioner
D=diag((diag(K)./diag(M)).*diag(K)+1/alpha*diag(M));
MA = 'mass_diagonal'; massparams = struct('M',M,'beta',alpha);
MS = 'schur_diagonal'; schurparams = struct('D',D);

tic %%start timing
[x_it,resvec,iter,flag] = poissoncontrol_minres(M,K,alpha,Myhat,d,maxit,tol,MA,massparams,MS,schurparams);
etoc = toc;
%whatdidiget
results{2,1} = resvec;
results{2,2} = iter;
results{2,3} = 'Diagonal Preconditioner';
results{2,4} = etoc;

%% Block preconditioner with Schur Approximate 1
% x_it = schurparams.L\(schurparams.M*(schurparams.L\x_it));
MA = 'mass_exact'; massparams = struct('M',M,'beta',alpha);
MS = 'schur_exact'; schurparams = struct('M',M,'L',K);

tic %%start timing
[x_it,resvec,iter,flag] = poissoncontrol_minres(M,K,alpha,Myhat,d,maxit,tol,MA,massparams,MS,schurparams);
etoc = toc;
%whatdidiget
results{3,1} = resvec;
results{3,2} = iter;
results{3,3} = '\begin{tabular}{c}Approximate Block Diagonal\\$M$,$K M^{-1} K$\end{tabular}';
results{3,4} = etoc;

%% Block preconditioner with Schur Approximate 2
% x_it = schurparams.L\(schurparams.M*(schurparams.L\x_it));
MA = 'mass_exact'; massparams = struct('M',M,'beta',alpha);
MS = 'schur_exact'; schurparams = struct('M',M,'L',K+1/sqrt(alpha)*M);

tic %%start timing
[x_it,resvec,iter,flag] = poissoncontrol_minres(M,K,alpha,Myhat,d,maxit,tol,MA,massparams,MS,schurparams);
etoc = toc;
%whatdidiget
results{4,1} = resvec;
results{4,2} = iter;
results{4,3} = '\begin{tabular}{c}Approximate Block Diagonal\\$M$,$(K + 1/\sqrt{\alpha} M) M^{-1} (K + 1/\sqrt{\alpha} M)$\end{tabular}';
results{4,4} = etoc;

%% Block preconditioner diagonal approximation of the Mass matrix, AMG for Schur

MA = 'mass_diagonal'; massparams = struct('M',M,'beta',alpha);
number_vcycles = default('number of AMG V-Cycles? (default 2)',2);
amg_grid = amg_grids_setup(K);
sm_ch = default('AMG smoother; point Jacobi (1), point Gauss-Seidel (2) or ILU (3)? (default 1)',1);
if sm_ch==1
    smoother_type = 'PDJ';
elseif sm_ch==2
    smoother_type = 'PGS';
elseif sm_ch==3
    smoother_type = 'ILU';
else
    error('illegal parameter choice ''smoother_type'', try again.')
end
no_sweeps = default('number of pre- and post-smoothing steps? (default 2)',2);
if no_sweeps>=0 && fix(no_sweeps)==no_sweeps
    smoother_params = amg_smoother_params(amg_grid,smoother_type,no_sweeps);
else
    error('illegal parameter choice ''no_sweeps'', try again.')
end
amg_smoother = amg_smoother_setup(amg_grid, smoother_params);
MS = 'schur_amg'; schurparams = struct('M',M,'A',K);

tic %%start timing
[x_it,resvec,iter,flag] = poissoncontrol_minres(M,K,alpha,Myhat,d,maxit,tol,MA,massparams,MS,schurparams);
etoc = toc;
%whatdidiget
results{5,1} = resvec;
results{5,2} = iter;
results{5,3} = '\begin{tabular}{c}Approximate Block Diagonal\\$D(M)$, AMG for $K M^{-1} K$\end{tabular}';
results{5,4} = etoc;

%% Block diagonal preconditioner diagonal mass matrix and AMG for Schur of type 2
MA = 'mass_diagonal'; massparams = struct('M',M,'beta',alpha);
number_vcycles = default('number of AMG V-Cycles? (default 2)',2);
amg_grid = amg_grids_setup(K+1/sqrt(alpha)*M);
sm_ch = default('AMG smoother; point Jacobi (1), point Gauss-Seidel (2) or ILU (3)? (default 1)',1);
if sm_ch==1
    smoother_type = 'PDJ';
elseif sm_ch==2
    smoother_type = 'PGS';
elseif sm_ch==3
    smoother_type = 'ILU';
else
    error('illegal parameter choice ''smoother_type'', try again.')
end
no_sweeps = default('number of pre- and post-smoothing steps? (default 2)',2);
if no_sweeps>=0 && fix(no_sweeps)==no_sweeps
    smoother_params = amg_smoother_params(amg_grid,smoother_type,no_sweeps);
else
    error('illegal parameter choice ''no_sweeps'', try again.')
end
amg_smoother = amg_smoother_setup(amg_grid, smoother_params);
MS = 'schur_amg'; schurparams = struct('M',M,'A',K+1/sqrt(alpha)*M);

tic %%start timing
[x_it,resvec,iter,flag] = poissoncontrol_minres(M,K,alpha,Myhat,d,maxit,tol,MA,massparams,MS,schurparams);
etoc = toc;
%whatdidiget
results{6,1} = resvec;
results{6,2} = iter;
results{6,3} = '\begin{tabular}{c}Approximate Block Diagonal\\$D(M)$, AMG for $(K + 1/\sqrt{\alpha} M) M^{-1} (K + 1/\sqrt{\alpha} M)$ \end{tabular}';
results{6,4} = etoc;

%% Block diagonal preconditioner Chebyshev for Mass, AMG for Schur version 1
fprintf('Chebyshev mass matrix approx., Schur comp. approx. 1, AMG ...\n')
cheb_its = default('number of Chebyshev iterations? (default 10)',10);
% Call mass matrix 'Q', not 'M', below to agree with structure of code
MA = 'mass_chebyshev';
massparams = struct('Q',M,'its',cheb_its,'qmethod',qmethod,'beta',alpha);
number_vcycles = default('number of AMG V-Cycles? (default 2)',2);
amg_grid = amg_grids_setup(K);
sm_ch = default('AMG smoother; point Jacobi (1), point Gauss-Seidel (2) or ILU (3)? (default 1)',1);
if sm_ch==1
    smoother_type = 'PDJ';
elseif sm_ch==2
    smoother_type = 'PGS';
elseif sm_ch==3
    smoother_type = 'ILU';
else
    error('illegal parameter choice ''smoother_type'', try again.')
end
no_sweeps = default('number of pre- and post-smoothing steps? (default 2)',2);
if no_sweeps>=0 && fix(no_sweeps)==no_sweeps
    smoother_params = amg_smoother_params(amg_grid,smoother_type,no_sweeps);
else
    error('illegal parameter choice ''no_sweeps'', try again.')
end
amg_smoother = amg_smoother_setup(amg_grid, smoother_params);
MS = 'schur_amg'; schurparams = struct('M',M,'A',K);

tic %%start timing
[x_it,resvec,iter,flag] = poissoncontrol_minres(M,K,alpha,Myhat,d,maxit,tol,MA,massparams,MS,schurparams);
etoc = toc;
%whatdidiget
results{7,1} = resvec;
results{7,2} = iter;
results{7,3} = '\begin{tabular}{c}Approximate Block Diagonal\\Chebyshev for $M$, AMG for $K M^{-1} K$\end{tabular}';
results{7,4} = etoc;

%% Block diagonal preconditioner Chebyshev for Mass matrix, AMG for Schur Version 2

fprintf('Chebyshev mass matrix approx., Schur comp. approx. 2, AMG ...\n')
cheb_its = default('number of Chebyshev iterations? (default 10)',10);
% Call mass matrix 'Q', not 'M', below to agree with structure of code
MA = 'mass_chebyshev';
massparams = struct('Q',M,'its',cheb_its,'qmethod',qmethod,'beta',alpha);
number_vcycles = default('number of AMG V-Cycles? (default 2)',2);
amg_grid = amg_grids_setup(K+1/sqrt(alpha)*M);
sm_ch = default('AMG smoother; point Jacobi (1), point Gauss-Seidel (2) or ILU (3)? (default 1)',1);
if sm_ch==1
    smoother_type = 'PDJ';
elseif sm_ch==2
    smoother_type = 'PGS';
elseif sm_ch==3
    smoother_type = 'ILU';
else
    error('illegal parameter choice ''smoother_type'', try again.')
end
no_sweeps = default('number of pre- and post-smoothing steps? (default 2)',2);
if no_sweeps>=0 && fix(no_sweeps)==no_sweeps
    smoother_params = amg_smoother_params(amg_grid,smoother_type,no_sweeps);
else
    error('illegal parameter choice ''no_sweeps'', try again.')
end
amg_smoother = amg_smoother_setup(amg_grid, smoother_params);
MS = 'schur_amg'; schurparams = struct('M',M,'A',K+1/sqrt(alpha)*M);

tic %%start timing
[x_it,resvec,iter,flag] = poissoncontrol_minres(M,K,alpha,Myhat,d,maxit,tol,MA,massparams,MS,schurparams);
etoc = toc;
%whatdidiget
results{8,1} = resvec;
results{8,2} = iter;
results{8,3} = '\begin{tabular}{c}Approximate Block Diagonal\\Chebyshev for $M$, AMG for $(K + 1/\sqrt{\alpha} M) M^{-1} (K + 1/\sqrt{\alpha} M)$ \end{tabular}';
results{8,4} = etoc;

%% Visualize performances
linetypes = {'-','--',':','-.','-','--',':','-.'};
markers = {'','o','+','x','h','s','d','v'};
figure(3)
clf
hold on
for i=1:size(results,1)
    semilogy(0:results{i,2},results{i,1},strcat(markers{i},linetypes{i}),...
        'LineWidth',2,...
        'Displayname',results{i,3});
end
hold off
set(gca,'YScale','log');
legend('Interpreter','latex','Location','westoutside');
title(sprintf('$\\alpha = %1.2e$ - ndof = %d',alpha,np),'Interpreter','latex');


figure(4)
Y = cell2mat(results(:,4)); % Times
bar(Y)
xtickangle(40)
xticklabels(results(:,3))
xaxisproperties= get(gca,'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
text(1:length(Y),Y,num2str(cell2mat(results(:,2))),'vert','bottom','horiz','center');
title(sprintf('$\\alpha = %1.2e$ - ndof = %d',alpha,np),'Interpreter','latex');