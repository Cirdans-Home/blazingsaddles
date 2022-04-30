%% Optimal Control of Poisson Equation with upper-bounded constraints
% This code solves the optimal control problem for the semilinear equation
% \nabla^2 y + y^3 = u, y = 0 on the unit square and an l1 term in the
% functional to optimize

clear;


%% Coordinates
% We use a finite difference mesh with n mesh points:
n = 200;
h=1/(n+1);
[x1,y1]=meshgrid(h:h:1-h,h:h:1-h);

%% Desired state
desiredstate = @(x,y) sin(2*pi*x).*sin(2*pi*y).*exp(2*x)/6;
yd=feval(desiredstate,x1,y1);
figure(1)
subplot(2,3,3);
mesh(x1,y1,yd);
title('Desired State');
yd=reshape(yd,n^2,1);

%% Optimization parameters
alpha = 1e-3;  % L1 weight 
beta = 0.008;  % Regularization parameter

%% Auxiliary matrices: Laplacian, Identity and O matrices
L = discretelaplacian(n,h);
I = speye(n^2);
O = sparse(n^2,n^2);

%% Initialization
u=sparse(n^2,1);
y=sparse(n^2,1);
p=sparse(n^2,1);
lam=sparse(n^2,1);

%% Semismooth Newton
res=1; iter=0;
while res >= 1e-10
    iter=iter+1; % Increase iter counte
    % Semismooth Newton step 
    Y=spdiags(y,0,n^2,n^2); % Auxiliary matrices for the step
    P=spdiags(p,0,n^2,n^2);
    Ximax=spdiags(spones(max(0,u+1/alpha*(lam-beta))),0,n^2,n^2);
    Ximin=spdiags(spones(min(0,u+1/alpha*(lam+beta))),0,n^2,n^2);
    Xi=Ximax+Ximin;
    A=[L+3*Y.^2 -I O O
        -I+6*Y.*P O L+3*Y.^2 O
        O alpha*I I I
        O I-Xi O -1/alpha*Xi];
    F=[ -L*y-y.^3+u
        -L*p-3*Y.^2*p+y-yd
        -p-alpha*u-lam
        -u+max(0,u+1/alpha*(lam-beta))+min(0,u+1/alpha*(lam+beta))];
    delta=A\F;
    uprev=u;
    yprev=y;
    pprev=p;
    y=y+delta(1:n^2);
    u=u+delta(n^2+1:2*n^2);
    p=p+delta(2*n^2+1:3*n^2);
    lam=lam+delta(3*n^2+1:4*n^2);
    res=norm(u-uprev)+norm(y-yprev)+norm(p-pprev);
    fprintf("\t Iter %d Improvement %e ||F|| = %e\n",iter,res,norm(F));
    figure(1)
    subplot(2,3,1);
    contour(x1,y1,reshape(y,size(x1)));
    subplot(2,3,2);
    mesh(x1,y1,reshape(y,size(x1)));
    title('State')
    subplot(2,3,4);
    contour(x1,y1,reshape(u,size(x1)));
    subplot(2,3,5);
    mesh(x1,y1,reshape(u,size(x1)));
    title('Control')
    subplot(2,3,6)
    contourf(x1,y1,abs(reshape(u,size(x1))-reshape(yd,size(x1))));
    pause()    
end

function L = discretelaplacian(n,h)
%DICRETELAPLACIAN builds the 2D discrete Laplacian on the square with
%homogeneous Dirichlet BCs.
d(n:n:n^2) = 1;
d=d';
e=sparse(n^2,1);
e(1:n:(n-1)*n+1) = 1;
b = ones(n^2,1);
a = [b,b-d,-4*b,b-e,b];
L =-1/(h^2)*spdiags(a,[-n,-1,0,1,n],n^2,n^2);
end