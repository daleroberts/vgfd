% Solve the Partial-integro differential equation for option pricing problems
% where the underlying has variance gamma distributed returns.
%
% Dale Roberts <dale.o.roberts@gmail.com>

clear all;

%% Parameters

global K sigma theta nu

T = 0.1;
K = 100;

r = 0.10;
q = 0;

sigma = 0.12;
theta = -0.14;
nu = 0.20;

xmin = log(10);
xmax = log(260);

nx = 120;
nt = 30;
M = 30;

%% Setup mesh

dx = (xmax-xmin)/(nx-1);
dt = T/(nt-1);

x = xmin:dx:xmax;
t = 0:dt:T;

%% Operator T

sigdx = 0.5*dx^3*(vgfd_k(dx) + vgfd_k(-dx));

P = (sigma^2 + sigdx)/2;
a = -P/(dx^2) - sigdx/(4*dx);
b = 1/dt + r + 2*P/(dx^2);
c = -P/(dx^2) + sigdx/(4*dx);

e = ones(nx,1);
z = zeros(nx,1);
T = spdiags([a*e b*e c*e],-1:1,nx,nx);

%% Operator B

lambdadx = quadgk(@vgfd_k, -M*dx, -dx) + quadgk(@vgfd_k, dx, M*dx);
omegadx = quadgk(@(x) (1-exp(x)).*vgfd_k(x),-M*dx,dx) ...
    + quadgk(@(x) (1-exp(x)).*vgfd_k(x),dx, M*dx);

D1 = 1/(2*dx)*spdiags([-e z e],-1:1,nx,nx);

%% Setup problem and solve

U=zeros(nx,nt);
U(:,1)=vgfd_psi(x);

for j=2:nt
    J = vgfd_J(x,U(:,j-1),M);
    B = U(:,j-1)/dt-(q-r-omegadx+sigma^2/2)*D1*U(:,j-1)-lambdadx*U(:,j-1)+J;
    U(:,j) = T\B;
    U(1,j) = U(1,j-1);
    U(end,j) = exp(x(end))*exp(-q*(j-1)*dt)-K*exp(-r*(j-1)*dt);
end

%% Plot solution

% style
set(gcf, 'Renderer', 'opengl');
shading faceted
colormap Jet

surf(exp(x),t,U','EdgeColor','None');
view(-140, 30);
xlabel('space x');
ylabel('time t');
zlabel('u');
