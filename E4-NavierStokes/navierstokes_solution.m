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
% ------ Horizontal dimension [-1,outbound] ------------------------------%
outbnd = 5;
% -------Grid parameters--------------------------------------------------%
nc=4; % Number of refinement levels
stretch=1; % Grid stretching
n=2^(nc-2);
%% Grid Generation
if n<1, error('illegal nc parameter, try again.'), end
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
cd datafiles/
save step_grid.mat mv xy bound mbound outbnd stretch x y bndxy bnde obs sbnde
clear variables
%% Building Matrices

domain = 3; pde=3; enclosed=0;
load step_grid.mat
% ------ Type of discretization ----------------------------------------- %
q_in = 1; % Q1-Q1/Q1-P0/Q2-Q1/Q2-P1: 1/2/3/4
% ----------------------------------------------------------------------- %
qmethod=q_in-1;
if qmethod==2,
    [x,y,xy,xyp,mp,map] = q2q1gridx(x,y,xy,mv,bound);
    [A,B,Q,G,Bx,By,f,g] = stokes_q2q1(xy,xyp,mv,mp);
elseif qmethod==3,
    [x,y,xy,xyp,ee] = q2p1grid(x,y,xy,mv,bound);
    [A,B,Q,G,Bx,By,f,g] = stokes_q2p1(xy,xyp,mv);
elseif qmethod==0
    [ev,ee,ebound,xyp] = q1q1grid(x,y,xy,mv,bound,mbound);
    [A,B,Q,C,G,Bx,By,f,g] = stokes_q1q1(xy,ev);
elseif qmethod==1
    [ev,ee,ebound,xyp] = q1p0grid(x,y,xy,mv,bound,mbound);
    [A,B,Q,C,G,Bx,By,f,g] = stokes_q1p0(xy,xyp,mv,ev);
end
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
clear variables
%% Solve the nonlinear system
fprintf('Incompressible flow problem on step domain ...\n')
% ------------- Viscosity parameter --------------------------------------%
viscosity= 1/210; % 1/210
nlmethod = 2; % Picard/Newton/hybrid linearization 1/2/3
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
% Convergence tolerance ------------------------------------------------- %
tol_nl=1.1*eps;

%% Computes Stokes Solution as initial guess
load step_stokes_nobc.mat
fprintf('stokes system ...\n')
% boundary conditions
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
spc=1; % Uniform streamlines 2 for non-uniform
contourn = 50; % Number of contourn lines
flowplot09(qmethod,xst,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,33);
fprintf('\n\ninitial nonlinear residual is %e ',nlres0_norm)
fprintf('\nStokes solution residual is %e\n', nlres_norm)
flowsol = xst;

%% March the solution with one of the methods
pde=4;
it_p = 0;
nlrs_norm_picard = [];
%
% nonlinear iteration
%%% Picard startup step
while nlres_norm>nlres0_norm*tol_nl && it_p<maxit_p,
    it_p = it_p+1;
    fprintf('\nPicard iteration number %g \n',it_p),
    % compute Picard correction and update solution
    if     qmethod>1,  dxns = -[Anst,Bst';Bst,sparse(np,np)]\nlres;
    elseif qmethod<=1, dxns = -[Anst,Bst';Bst,-nubeta*C]\nlres;
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
    nlrs_norm_picard = [nlrs_norm_picard,nlres_norm];
    fprintf('nonlinear residual is %e',nlres_norm)
    fprintf('\n   velocity change is %e\n',soldiff)
    % plot solution
    flowplot09(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,33+it_p);drawnow
    % Save: requires export_fig
    set(gcf,'Color','None');
    export_fig(sprintf('Picard_it_%03d.png',it_p));
    %  flowplotz(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,66); drawnow
    pause(1)
    flowsol = xns;
    % end of Picard iteration loop
end

if nlmethod == 0
    try
        save residuals.mat nlrs_norm_picard -append
    catch
        save residuals.mat nlrs_norm_picard
    end
end
%%%
%
it_nl = it_p;
it_n = 0;
nlres_norm_newton = [];
%%%% Newton iteration loop
while (nlres_norm > nlres0_norm*tol_nl) && (it_nl < maxit_p + maxit_n),
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
    if     qmethod>1,  dxns = -[Jnst,Bst';Bst,sparse(np,np)]\nlres;
    elseif qmethod<=1, dxns = -[Jnst,Bst';Bst,-nubeta*C]\nlres;
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
    nlres_norm_newton = [nlres_norm_newton,nlres_norm];
    nnv=length(fst); soldiff=norm(xns(1:nnv)-flowsol(1:nnv));
    fprintf('nonlinear residual is %e',nlres_norm)
    fprintf('\n   velocity change is %e\n',soldiff)
    % plot solution
    flowplot09(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,bndxy,bnde,obs,contourn,spc,33+it_nl);drawnow
    %   flowplotz(qmethod,xns,By,Bx,A,xy,xyp,x,y,bound,66); drawnow
    set(gcf,'Color','None');
    export_fig(sprintf('Newton_it_%03d.png',it_nl));
    pause(1)
    flowsol = xns;
    % end of Newton iteration loop
end

if nlmethod == 1
    try
        save residuals.mat nlres_norm_newton -append
    catch
        save residuals.mat nlres_norm_newton
    end
end

if nlmethod == 2
    nlres_hyb = [nlrs_norm_picard(1:maxit_p),nlres_norm_newton];
    try
        save residuals.mat nlres_hyb -append
    catch
        save residuals.mat nlres_hyb 
    end
end
%
if nlres_norm <= nlres0_norm * tol_nl,
    fprintf('\nfinished, nonlinear convergence test satisfied\n\n');
else
    fprintf('\nfinished, stopped on iteration counts\n\n');
end
%
%%% estimate errors
if qmethod==1
    [jmpx,jmpy,els] = stressjmps_q1p0(viscosity,flowsol,xy,ev,ebound);
    [error_x,error_y,fex,fey,ae] = navierpost_q1p0_p(viscosity,flowsol,jmpx,jmpy,els,xy,ev);
    [error_x,error_y] = navierpost_q1p0_bc(viscosity,ae,fex,fey,...
        error_x,error_y,xy,ev,ebound);
    error_div = q1div(xy,ev,flowsol);
    errorest=sqrt(sum(error_x.^2 + error_y.^2 + error_div.^2));
    fprintf('estimated overall error is %10.6e \n',errorest)
    ee_error=sqrt((error_x.^2 + error_y.^2 + error_div.^2));
    %% plot element errors
    %eplotl(ee_error,ev,xy,x,y,67);
    %% plot macroelement errors
    mplotl(ee_error,ev,xy,x,y,67), title('Estimated error')
    pause(5), figure(66)
elseif qmethod==0
    error_div = q1div(xy,ev,flowsol);
elseif qmethod>1, navierpost,
end

