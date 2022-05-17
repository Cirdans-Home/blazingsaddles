%% Iterative solution of Stokes with Stable Discretization
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
%% Iterative Solutions
cd ..
S = Bst*(Ast\Bst');
M = [Ast,Bst';Bst,sparse(np,np)];
b = [fst;gst];
% Block Diagonal Schur-Complement
P = [Ast,sparse(nu,np);sparse(np,nu),S];
[x,~,~,~,resvecds] = gmres(M,b,[],1e-7,100,[],P);
% Block Tridiagonal Schur-Complement
P = [Ast,Bst';sparse(np,nu),S];
[x,~,~,~,resvectds] = gmres(M,b,[],1e-7,100,[],P);
% Block Diagonal Mass-Matrix
P = [Ast,sparse(nu,np);sparse(np,nu),0.5*(Q+Q')];
[x,~,~,~,resvecdq] = gmres(M,b,[],1e-7,100,[],P);
% Block Tridiagonal Mass-Matrix
P = [Ast,Bst';sparse(np,nu),0.5*(Q+Q')];
[x,~,~,~,resvectdq] = gmres(M,b,[],1e-7,100,[],P);

%% Residuals
figure(10)
semilogy(1:length(resvecds),resvecds./resvecds(1),'o--',...
    1:length(resvectds),resvectds./resvectds(1),'x-',...
    1:length(resvecdq),resvecdq./resvecdq(1),'s--',...
    1:length(resvectdq),resvectdq./resvectdq(1),'h-','LineWidth',2);
legend('Block Diagonal Schur','Block Tridiagonal Schur',...
    'Block Diagonal Q1-Mass','Block Tridiagonal Q1-Mass');

%% Using Flexible-GMRES
% Here we use the implementation from Y. Saad: 
% https://www-users.cse.umn.edu/~saad/software/MATLAB_DEMOS.tar.gz
sol = zeros(size(x));
maxits = 100;
tolIts = 1e-7;
im = inf; % Infinite restart
PRE = [Ast,Bst';sparse(np,nu),0.5*(Q+Q')]; % Here we use one of the 
                                           % variants as preconditioner
precfun = @(PRE,x) PRE\x;                  % and for the moment give a
                                           % direct solution
[sol,res,its] = fgmres (M,PRE,precfun,b,sol,im,maxits,tolIts);