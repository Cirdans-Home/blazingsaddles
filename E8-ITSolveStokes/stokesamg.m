%% An AMG Example from IFISS
% This example needs IFISS 3.6, so we make sure of having it
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
clear; clc; close all;
%% ---------------------------------------------------------------------- %
%% Build the discrete problem
%% ---------------------------------------------------------------------- %
pde=3; domain=10; enclosed=0;
outbnd=default('channel domain: L? (default unity)',1);
fprintf('\n\nGrid generation for extended channel domain.\n')
if exist('outbnd','var')==0,
    fprintf('setting channel length to default\n'), outbnd=1; end
nc=default('grid parameter? (default is 4)',4); nc=nc-1; %%% compatibility
if nc<1, error('illegal parameter choice, try again.'), end
grid_type=1;
ny=2^nc; npy=ny/2; nqy=ny/4;                   %%% modified from channel_domain
nx=ny*(outbnd/2);  npx=nx/2;  nqx=nx/4;        %%% to handle rectangular domain
nx = min(ny*(outbnd/2),ny*100/2);   npx = nx/2;   nqx = nx/4;
% compute (x,y) coordinates of vertices
yy=(1/npy:1/npy:1);
ypos=[0,yy]; yneg=-yy(length(yy):-1:1);
yyy=[yneg,ypos];
xxx=(-1:(outbnd/2)/npx:outbnd);
[xy,mv,left,right,bottom,top,mbound,x,y] = grid_xblock(xxx,yyy);
outbc=0;
% exclude outflow boundary data
bound=unique([left;top;bottom]);
kk=find(mbound(:,2)~=2); mbound=mbound(kk,:);
%% specify boundary information for graphics
% bndxy: (x,y)-coordinates of vertices that define the domain and
%         obstacle(s)
% bnde: boundary edges (node1 node2 1(for dirichlet)/0(for neumann))
% obs: obstacles (node1 node2 node3 node4)
% sbnde: boundary edges near which stretching is needed (edge1 edge2 ...)
% 'obs' and/or 'sbnde' can be absent if there is no obstacle in the problem
% and/or only uniform grid is needed
bndxy = [-1,-1; outbnd,-1; outbnd,1; -1,1];
bnde = [1,2,1; 2,3,0; 3,4,1; 4,1,1];
obs = [];
sbnde = [];

cd datafiles
save channel_grid.mat mv xy bound mbound outbnd grid_type outbc x y bndxy bnde obs
fprintf('generated datatfile: channel_grid.mat \n')
clear variables
load channel_grid.mat

% ------------------------------------------------------------------------%
% Build discretization
% ------------------------------------------------------------------------%
pde=3; domain=10; enclosed=0;
qmethod = 2; % Q1-Q1 Elements
[x,y,xy,xyp,mp,map] = q2q1gridx(x,y,xy,mv,bound);
[A,B,Q,G,Bx,By,f,g] = stokes_q2q1(xy,xyp,mv,mp);
save channel_stokes_nobc.mat pde domain outbnd qmethod grid_type A B Q f g xy xyp mbound bound x y
save channel_stokes_nobc.mat Bx By bndxy bnde obs -append
save channel_stokes_nobc.mat mv mp G map -append
fprintf('system matrices saved in channel_stokes_nobc.mat ...\n')
% boundary conditions
[Ast,Bst,fst,gst] = flowbc(A,B,f,g,xy,bound);
np=length(gst); nu=length(f)/2;
% ----------------------------------------------------------------------- %
% Direct Solution
% ----------------------------------------------------------------------- %
tic;
beta=0;
xst=[Ast,Bst';Bst,sparse(np,np)]\[fst;gst];
etoc=toc; fprintf('Stokes system solved in %8.3e seconds\n',etoc);
% ----------------------------------------------------------------------- %
% Error estimation
% ----------------------------------------------------------------------- %
stokespost;
% ----------------------------------------------------------------------- %
% Plot Solution
% ----------------------------------------------------------------------- %
spc=default('uniform/nonuniform streamlines 1/2 (default uniform)',1);
flowplot(qmethod,xst,By,Bx,A,xy,xyp,x,y,bound,spc,33);
fprintf('\n');
% ------------------------------------------------------------------------%
% Inf-Sup condition
% ------------------------------------------------------------------------%
[np,nu]=size(Bst);
if size(A,1) < 5000
    fprintf('inf-sup eigenvalue estimation \n')
    %leig=default('eig/eigs? 0/1 (default is eig)',0);
    leig = 0;
    S=Bst*(Ast\Bst');
    % compute discrete eigenvalues
    S=0.5*(S+S');Q=0.5*(Q+Q');
    if leig==1
        [V,D,flag] = eigs(S,Q,4,'sm');
        e=diag(D); e=sort(real(e));
        fprintf('inf-sup constant is %g \n',e(2))
    else
        [V,D]=eig(full(S),full(Q));
        e=diag(D); [e,k]=sort(real(e)); V=V(:,k); pp=length(e);
        fprintf('inf-sup constant is %g \n',e(2))
        fprintf('     upper bound is %g \n',e(pp))
    end
end

cd ..
%% AMG
%-------------------------------------------------------------------------%
% We build here an AMG Algorithm for the (1,1)-Block of our Saddle-Point
% matrix, then we try to use it in the FGMRES algorithm
%-------------------------------------------------------------------------%

% Firs we setup the grids:
i_method = 2; % i_method: interpolation method
              % 0=BHM direct, 1=Stuben direct, 2=Stuben standard (default is 1)
% Build the grids:
max_levels = 10;
tic;
grid_data = amg_grids_setup(Ast, i_method, max_levels);
fprintf("Time to build the grids: %1.2e\n",toc);
tic;
smoother_params = amg_smoother_params(grid_data, 'PGS', 3);
smoother_data = amg_smoother_setup(grid_data, smoother_params);
fprintf("Time to build the smoothers: %1.2e\n",toc);
% We can visualize the grids we have built by doing:
[last_level] = amg_coarsen_plot(grid_data);
% We now have a routine to apply a single V(1,1)-Cycle:
max_Clevels = 7;  % Maximum number of levels before coarse solve
x0 = zeros(nu,1); % Initial guess
x = amg_v_cycle(fst, grid_data, smoother_data, max_Clevels, x0);
fprintf("The residual after 1 Iteration is %1.1e\n",norm(Ast*x-fst));
fprintf("The relative residual after 1 Iteration is %1.1e\n",norm(Ast*x-fst)/norm(fst));
%% Systems with the mass matrix
% ----------------------------------------------------------------------- %
% Chebyshev iteration
% ----------------------------------------------------------------------- %
Q  = 0.5*(Q+Q');
p0 = zeros(np,1);
maxit = np;
tol = 1e-6;
lmax = eigs(Q,1,'largestabs');
lmin = eigs(Q,1,'smallestabs');
[p,it,r] = chebyshev(Q, gst, p0, maxit, tol, lmax, lmin);
% Now you can try to assemble together this routine with the other example 
% with FGMRES to build an incomplete solver.