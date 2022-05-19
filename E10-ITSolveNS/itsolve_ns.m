%% Solution of the Navier-Stokes equation
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
pde=3; enclosed=0;
%% Define the geometry of the problem
% We define here the spatial domain we are going to use in the experiments!
outbnd=default('horizontal dimensions [-1,L]: L? (default L=5)',5);
fprintf('\nGrid generation for backward-facing step domain.\n')
nc=default('grid parameter: 3 for underlying 8x24 grid (default is 4)',4);
stretch=default('grid stretch factor (default is 1)',1);
n=2^(nc-2);
if n<1, error('illegal nc parameter, try again.'), end
%
% generate the inlet grid
if stretch==1
    xx=[-1:1/n:0];  yy=[0:1/n:1];
elseif stretch>1
    dy=(1-stretch)/(1-stretch^n);
    dd(1)=dy;for k=1:n-1, dd(k+1)=dd(k)+dy*stretch^k; end, dd(n)=1;
    xx=sort([0,-dd]);
    yy=[0,dd];
else, error('illegal stretch parameter, try again.')
end
[xyl,mvl,leftl,rightl,bottoml,topl,mboundl,xl,yl] = grid_xblock(xx,yy);
refleft=leftl;refright=rightl; refbottom=bottoml;reftop=topl;
bound=unique([leftl;rightl;topl;bottoml]);
%
if stretch==1;
    xx=[0:1/n:outbnd]; yy=[-1:1/n:1];
elseif stretch>1
    dd(1)=dy;for k=1:n-1, dd(k+1)=dd(k)+dy*stretch^k; end, dd(n)=1;
    ddx=dd(n)-dd(n-1);ndx=ceil(1/ddx);
    fprintf('   outlet subdivision parameter set to %d\n',ndx)
    fprintf('   ... associated uniform grid value is %d\n',n)
    xx=sort([0,dd,[1+1/ndx:1/ndx:outbnd]]);
    yy=sort([-dd,0,dd]);
end
[xyr,mvr,leftr,rightr,bottomr,topr,mboundr,xr,yr] = grid_xblock(xx,yy);
bound=unique([leftr;rightr;topr;bottomr]);
%
% merge grids
[xy,mv,left,bottom,top,right,mbound,x,y] = ...
    grid_mergeleftright(xyl,mvl,leftl,rightl,bottoml,topl,mboundl,xl,yl,...
    xyr,mvr,leftr,rightr,bottomr,topr,mboundr,xr,yr);
bound=unique([left;top;bottom]);
%macrogridplot(xy,mv,bound,mbound);
%figure(1),pause(1),set(gcf,'Visible','off'),
%figure(2), pause(1),set(gcf,'Visible','off'),
fprintf('  All done.\n\n')
%
% specify boundary information for graphics
% bndxy: (x,y)-coordinates of vertices that define the domain and
%         obstacle(s)
% bnde: boundary edges (node1 node2 1(for dirichlet)/0(for neumann))
% obs: obstacles (node1 node2 node3 node4)
% sbnde: boundary edges near which stretching is needed (edge1 edge2 ...)
% 'obs' and/or 'sbnde' can be absent if there is no obstacle in the problem
% and/or only uniform grid is needed
bndxy = [-1,-1;  0,-1;  outbnd,-1; outbnd,1; -1,1; -1,0; 0,0];
bnde = [2, 3, 1; 3, 4, 0; 4, 5, 1; 5, 6, 1; 6, 7, 1; 7, 2, 1];
obs = [1     2     7     6];
sbnde = [5 6];
cd datafiles
save step_grid.mat mv xy bound mbound outbnd stretch x y bndxy bnde obs sbnde
domain = 3; pde=3; enclosed=0;
load step_grid.mat
cd ..

%% Discretize everything
% What FEM space are we going to select? Is it stable or unstable?
q_in=default('Q1-Q1/Q1-P0/Q2-Q1/Q2-P1: 1/2/3/4? (default Q1-P0)',2);
qmethod=q_in-1;
if qmethod==2
    [x,y,xy,xyp,mp,map] = q2q1gridx(x,y,xy,mv,bound);
    [A,B,Q,G,Bx,By,f,g] = stokes_q2q1(xy,xyp,mv,mp);
elseif qmethod==3
    [x,y,xy,xyp,ee] = q2p1grid(x,y,xy,mv,bound);
    [A,B,Q,G,Bx,By,f,g] = stokes_q2p1(xy,xyp,mv);
elseif qmethod==0
    [ev,ee,ebound,xyp] = q1q1grid(x,y,xy,mv,bound,mbound);
    [A,B,Q,C,G,Bx,By,f,g] = stokes_q1q1(xy,ev);
elseif qmethod==1
    [ev,ee,ebound,xyp] = q1p0grid(x,y,xy,mv,bound,mbound);
    [A,B,Q,C,G,Bx,By,f,g] = stokes_q1p0(xy,xyp,mv,ev);
end
cd datafiles
save step_stokes_nobc.mat pde domain outbnd qmethod stretch A B Q f g xy xyp mbound bound x y
save step_stokes_nobc.mat Bx By bndxy bnde obs -append
if qmethod==1
    save step_stokes_nobc.mat C G ev ee ebound -append
elseif qmethod==0
    save step_stokes_nobc.mat C G ev ee ebound mv enclosed -append
elseif qmethod==2
    save step_stokes_nobc.mat G mv mp map -append
else
    save step_stokes_nobc.mat G ee mv  -append
end
fprintf('system matrices saved in step_stokes_nobc.mat ...\n')
cd ..

clear variables
fprintf('Incompressible flow problem on step domain ...\n')
viscosity=default('viscosity parameter (default 1/210)',1/210);
nlmethod=default('Picard/Newton/hybrid linearization 1/2/3 (default hybrid)',3);
nlmethod=nlmethod-1;
if nlmethod==0,
    maxit_p=default('number of Picard iterations (default 9)',9);
    maxit_n=0;
elseif nlmethod==1,
    maxit_p=0;
    maxit_n=default('number of Newton iterations (default 6)',6);
else
    maxit_p=default('number of Picard iterations (default 3)',3);
    maxit_n=default('number of Newton iterations (default 5)',5);
end
tol_nl=default('nonlinear tolerance (default 1.1*eps)',1.1*eps);
%
%
%% initialize for nonlinear iteration: compute Stokes solution
%% load assembled matrices
cd datafiles
load step_stokes_nobc.mat
cd ..
%
fprintf('stokes system ...\n')
%% boundary conditions
[Ast,Bst,fst,gst] = flowbc(A,B,f,g,xy,bound);
nlres0_norm = norm([fst;gst]);
%
nv=length(fst)/2; np=length(gst);
if qmethod>1,
    beta=0;
    xst=[Ast,Bst';Bst,sparse(np,np)]\[fst;gst];
elseif qmethod==1
    beta=1/4;     % default parameter
    xst=[Ast,Bst';Bst,-beta*C]\[fst;gst];
elseif qmethod==0
    fprintf('computing pressure stabilized solution...\n')
    xst=[Ast,Bst';Bst,-C]\[fst;gst]; beta=1;
end
% compute residual of Stokes solution
if qmethod>1
    N = navier_q2(xy,mv,xst);
elseif qmethod<=1,
    nubeta=beta/viscosity;
    fprintf('\n nubeta is set to %e \n',nubeta)
    N = navier_q1(xy,ev,xst);
end
Anst = viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N];
[Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
if     qmethod>1,  nlres = [Anst,Bst';Bst,sparse(np,np)]*xst-[fst;gst];
elseif qmethod<=1, nlres = [Anst,Bst';Bst,-nubeta*C]*xst-[fst;gst];
end
nlres_norm  = norm(nlres);
%%% plot solution
spc=default('uniform/nonuniform streamlines 1/2 (default uniform)',1);
contourn = default('number of contour lines (default 50)',50);
flowplot09(qmethod,xst,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,33)
%flowplotl(qmethod,xst,By,Bx,A,xy,xyp,x,y,bound,33);
pause(1)
fprintf('\n\ninitial nonlinear residual is %e ',nlres0_norm)
fprintf('\nStokes solution residual is %e\n', nlres_norm)
flowsol = xst;
%
%
pde=4;
it_p = 0;

%% Nonlinear iteration
% We collect here the different nonlinear solve for NS
% 1) Picard Iteration
% 2) Newton Iteration
% 3) Hybrid Solution: Localize with Picard kill it with Newton!

% set indicator for strategy for stopping criterion for Picard/Newton steps
exact = default('inexact or exact linear iteration 1/2 (default exact)',2);

%% Set here the parameters for the iterative solution algorithm!
if exact == 1
    ns_params.itmeth = 1;         % GMRES/Bicgstab(2)  1/2
    ns_params.tol    = 1.d-6;     % For stokes only, set below for NS
    ns_params.maxit  = 200;
    ns_params.precon = 4;         % no / BFBT / Fp / xBFBt / PCD* / LSC* 0/1/2/3/4/5
    ns_params.precon_format = 1;  % ideal / amg  preconditioning? 1/2
end
tau_safety=1e-2;

%% Picard startup step
while nlres_norm>nlres0_norm*tol_nl && it_p<maxit_p
    % Do we need to fix the tolerances for using iterative methods?
    if exact~=1
        ns_params.tol = tol_nl/10;
    else
        ns_params.tol = tau_safety*nlres_norm;
    end
    it_p = it_p+1;
    fprintf('\nPicard iteration number %g \n',it_p),
    % Compute Picard correction and update solution
    if     qmethod>1
        % We are using here a stable method:
        if exact ~= 1
            dxns = -[Anst,Bst';Bst,sparse(np,np)]\nlres;
        else
            dxns = -it_nstokes_cortona(Anst,Bst,sparse(np,np),beta,Q,G, nlres, ...
                flowsol,viscosity,ns_params, domain);
        end
    elseif qmethod<=1
        % We are using an unstable method, we put a stabilization C:
        if exact ~= 1
            dxns = -[Anst,Bst';Bst,-nubeta*C]\nlres;
        else
            dxns = -it_nstokes_cortona(Anst,Bst,C,nubeta,Q,G, nlres, ...
                flowsol,viscosity,ns_params, domain);
        end
    end
    xns = flowsol + dxns;
    % compute residual of new solution
    if     qmethod>1,  N = navier_q2(xy,mv,xns);
    elseif qmethod<=1, N = navier_q1(xy,ev,xns);
    end
    Anst = viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N];
    [Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
    if     qmethod>1,  nlres = [Anst,Bst';Bst,sparse(np,np)]*xns-[fst;gst];
    elseif qmethod<=1, nlres = [Anst,Bst';Bst,-nubeta*C]*xns-[fst;gst];
    end
    nlres_norm = norm(nlres);
    nnv=length(fst); soldiff=norm(xns(1:nnv)-flowsol(1:nnv));
    fprintf('nonlinear residual is %e',nlres_norm)
    fprintf('\n   velocity change is %e\n',soldiff)
    % plot solution
    flowplot09(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,66);drawnow
    %  flowplotz(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,66); drawnow
    pause(1)
    flowsol = xns;
    % end of Picard iteration loop
end


%% Newton iteration loop
it_nl = it_p;
it_n = 0;
while (nlres_norm > nlres0_norm*tol_nl) && (it_nl < maxit_p + maxit_n),
    % Do we need to fix the tolerances for using iterative methods?
    if exact~=1
        ns_params.tol = tol_nl/10;
    else
        ns_params.tol = tau_safety*nlres_norm;
    end
    it_n = it_n+1;
    it_nl = it_nl+1;
    fprintf('\nNewton iteration number %g \n',it_n),
    % compute Jacobian of current solution
    if     qmethod>1,  [Nxx,Nxy,Nyx,Nyy] = newton_q2(xy,mv,flowsol);
    elseif qmethod<=1, [Nxx,Nxy,Nyx,Nyy] = newton_q1(xy,ev,flowsol);
    end
    J = viscosity*A + [N + Nxx, Nxy; Nyx, N + Nyy];
    Jnst = newtonbc(J,xy,bound);
    % compute Newton correction and update solution
    if     qmethod>1
        % We are solving here with a stable method!
        if exact ~= 1
            dxns = -[Jnst,Bst';Bst,sparse(np,np)]\nlres;
        else
            dxns = -it_nstokes_cortona(Jnst,Bst,sparse(np,np),beta,Q,G, nlres,...
                flowsol,viscosity,ns_params, domain);
        end
    elseif qmethod<=1
        % We are solving here with an unstable method!
        if exact ~= 1
            dxns = -[Jnst,Bst';Bst,-nubeta*C]\nlres;
        else
            dxns = -it_nstokes_cortona(Jnst,Bst,C,nubeta,Q,G, nlres, ...
                flowsol,viscosity,ns_params, domain);
        end
    end
    xns = flowsol + dxns;
    % compute residual of new solution
    if     qmethod>1,  N = navier_q2(xy,mv,xns);
    elseif qmethod<=1, N = navier_q1(xy,ev,xns);
    end
    Anst = viscosity*A + [N, sparse(nv,nv); sparse(nv,nv), N];
    [Anst,Bst,fst,gst] = flowbc(Anst,B,f,g,xy,bound);
    if     qmethod>1,  nlres = [Anst,Bst';Bst,sparse(np,np)]*xns-[fst;gst];
    elseif qmethod<=1, nlres = [Anst,Bst';Bst,-nubeta*C]*xns-[fst;gst];
    end
    nlres_norm = norm(nlres);
    nnv=length(fst); soldiff=norm(xns(1:nnv)-flowsol(1:nnv));
    fprintf('nonlinear residual is %e',nlres_norm)
    fprintf('\n   velocity change is %e\n',soldiff)
    % plot solution
    flowplot09(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,66);drawnow
    %   flowplotz(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,66); drawnow
    pause(1)
    flowsol = xns;
end

%% Check on the convergence
if nlres_norm <= nlres0_norm * tol_nl
    fprintf('\nfinished, nonlinear convergence test satisfied\n\n');
else
    fprintf('\nfinished, stopped on iteration counts\n\n');
end

%% Estimating the errors
if qmethod==1
    [jmpx,jmpy,els] = stressjmps_q1p0(viscosity,flowsol,xy,ev,ebound);
    [error_x,error_y,fex,fey,ae] = navierpost_q1p0_p(viscosity,flowsol,jmpx,jmpy,els,xy,ev);
    [error_x,error_y] = navierpost_q1p0_bc(viscosity,ae,fex,fey,...
        error_x,error_y,xy,ev,ebound);
    error_div = q1div(xy,ev,flowsol);
    errorest=sqrt(sum(error_x.^2 + error_y.^2 + error_div.^2));
    fprintf('estimated overall error is %10.6e \n',errorest)
    ee_error=sqrt((error_x.^2 + error_y.^2 + error_div.^2));
    % plot macroelement errors
    mplotl(ee_error,ev,xy,x,y,67), title('Estimated error')
    pause(5), figure(66)
elseif qmethod==0
    error_div = q1div(xy,ev,flowsol);
elseif qmethod>1, navierpost,
end