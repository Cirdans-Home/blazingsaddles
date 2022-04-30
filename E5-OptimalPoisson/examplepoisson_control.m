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

%% Build the grid
grid_type = 1;
square_type = 2;
% ----- Grid parameters ------------------------------------------------- %
nc = 4;
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