%% Stokes with Stabilized Discretization
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

%% We build Stokes test problem using IFISS
clear variables
% ----------------------------------------------------------------------- %
nc= 4; % Number of subdivisions > 2
% ----------------------------------------------------------------------- %
grid_type=1; % 2 for stretched
% ----------------------------------------------------------------------- %
n=2^nc; np=n/2; nq=n/4;
% y-direction
if grid_type==2
    hmax=nc/(2^(nc+1));
    x1=-1;x2=-2*hmax;x3=2*hmax;x4=1;nx1=2^(nc-1)-1;nx2=2;nx3=2^(nc-1)-1;
    y1=-1;y2=-2*hmax;y3=2*hmax;y4=1;ny1=2^(nc-1)-1;ny2=2;ny3=2^(nc-1)-1;
    y=subint(y1,y2,y3,y4,ny1,ny2,ny3);
    stretch=(y(3)-y(2))/(y(2)-y(1));
    x=y;
else
    yy=[1/np:1/np:1];
    ypos=[0,yy];
    yneg=-yy(length(yy):-1:1);
    y=[yneg,ypos]';
    x=y;
end
%
%% compute biquadratic element coordinates
nvtx=(n+1)*(n+1);
[X,Y]=meshgrid(x,y);
xx=reshape(X',nvtx,1);
yy=reshape(Y',nvtx,1);
xy=[xx(:),yy(:)];
%
kx = 1;
ky = 1;
mel=0;
for j=1:np
    for i=1:np
        mref=(n+1)*(ky-1)+kx;
        mel=mel+1;
        nvv(1) = mref;
        nvv(2) = mref+2;
        nvv(3) = mref+2*n+4;
        nvv(4) = mref+2*n+2;
        nvv(5) = mref+1;
        nvv(6) = mref+n+3;
        nvv(7) = mref+2*n+3;
        nvv(8)=  mref+n+1;
        nvv(9)=  mref+n+2;
        mv(mel,1:9)=nvv(1:9);
        kx = kx + 2;
    end
    ky = ky + 2;
    kx = 1;
end
%
%% compute boundary vertices and edges
% four boundary edges
k1=find( xy(:,2)==-1  );
e1=[]; for k=1:mel, if any(mv(k,5)==k1), e1=[e1,k]; end, end
ef1=ones(size(e1));
%
k2=find( xy(:,1)==1 & xy(:,2)<1   & xy(:,2) >-1);
e2=[]; for k=1:mel, if any(mv(k,6)==k2), e2=[e2,k]; end, end
ef2=2*ones(size(e2));
%
k3=find( xy(:,2)==1  );
e3=[]; for k=1:mel, if any(mv(k,7)==k3), e3=[e3,k]; end, end
ef3=3*ones(size(e3));
%
k4=find( xy(:,1)==-1 & xy(:,2)<1   & xy(:,2) >-1);
e4=[]; for k=1:mel, if any(mv(k,8)==k4), e4=[e4,k]; end, end
ef4=4*ones(size(e4));
%
bound=sort([k1;k2;k3;k4]);
mbound=[e1',ef1';e2',ef2';e3',ef3';e4',ef4'];

%% specify boundary information for graphics
% bndxy: (x,y)-coordinates of vertices that define the domain and
%         obstacle(s)
% bnde: boundary edges (node1 node2 1(for dirichlet)/0(for neumann))
% obs: obstacles (node1 node2 node3 node4)
% sbnde: boundary edges near which stretching is needed (edge1 edge2 ...)
% 'obs' and/or 'sbnde' can be absent if there is no obstacle in the problem
% and/or only uniform grid is neededbndxy = ...
bndxy = [-1,-1; 1,-1; 1,1; -1,1];
bnde = [1,2,1; 2,3,1; 3,4,1; 4,1,1];
obs = [];
sbnde = [1 2 3 4];
save cavity_grid.mat mv xy bound mbound grid_type  x y bndxy bnde obs
clear variables
load cavity_grid.mat

pde=3; domain=1; enclosed=1;
q_in=2; % Select Q1-P0 elements
qmethod=q_in-1;
if qmethod==2,
   [x,y,xy,xyp,mp,map] = q2q1gridx(x,y,xy,mv,bound);
   [A,B,Q,G,Bx,By,f,g] = stokes_q2q1(xy,xyp,mv,mp);
elseif qmethod==3,
   [x,y,xy,xyp,ee] = q2p1grid(x,y,xy,mv,bound);
   [A,B,Q,G,Bx,By,f,g] = stokes_q2p1(xy,xyp,mv);
elseif qmethod==8,
   [x,y,xy,xyp,ee] = q2p1grid(x,y,xy,mv,bound);
   [A,BB,QQ,G,Bx,By,f,gg] = stokes_q2p1(xy,xyp,mv);
   nnp=length(gg); ppk=[1:3:nnp];
   B=BB(ppk,:); Q=QQ(ppk,ppk); g=gg(ppk);
   fprintf('Reduction to Q2-P0 Stokes system done.\n')
elseif qmethod==0, 
   [ev,ee,ebound,xyp] = q1q1grid(x,y,xy,mv,bound,mbound);
   [A,B,Q,C,G,Bx,By,f,g] = stokes_q1q1(xy,ev);
elseif qmethod==1 
   [ev,ee,ebound,xyp] = q1p0grid(x,y,xy,mv,bound,mbound);
   [A,B,Q,C,G,Bx,By,f,g] = stokes_q1p0(xy,xyp,mv,ev);
end
save square_stokes_nobc.mat pde domain qmethod grid_type A B Q f g xy xyp mbound bound x y 
save square_stokes_nobc.mat Bx By bndxy bnde obs -append
if qmethod==1 
   save square_stokes_nobc.mat C G ev ee ebound -append
elseif qmethod==0 
   save square_stokes_nobc.mat C G ev ee ebound mv enclosed -append
elseif qmethod==2
   save square_stokes_nobc.mat mv mp G map -append
else
   save square_stokes_nobc.mat mv ee G -append
end
fprintf('system matrices saved in square_stokes_nobc.mat ...\n')

clear variables
load square_stokes_nobc.mat
%
fprintf('imposing (enclosed flow) boundary conditions ...\n') 
%% boundary conditions
[Ast,Bst,fst,gst] = flowbc(A,B,f,g,xy,bound);
%
np=length(gst); nu=length(f)/2;
% ----------------------------------------------------------------------- %
beta=1/4;   % Stabilization parameter
% ----------------------------------------------------------------------- %
tic;
xstz=[Ast,Bst',zeros(2*nu,1);
      Bst,-beta*C,ones(np,1)/np; ...
      zeros(1,2*nu),ones(1,np)/np,zeros(1,1)]\[fst;gst;0];
xst=xstz(1:end-1); multiplier=xstz(end);
etoc=toc; fprintf('Stokes system solved in %8.3e seconds\n',etoc) 

%% Plot the Solution
spc=1;
figure(2)
flowplot(qmethod,xst,By,Bx,A,xy,xyp,x,y,bound,spc,33);
fprintf('\n') 

%% The Saddle-Point matrix
M = [Ast,Bst',zeros(2*nu,1);Bst,-beta*C,ones(np,1)/np; ...
      zeros(1,2*nu),ones(1,np)/np,zeros(1,1)];
figure(3)
spy(M)