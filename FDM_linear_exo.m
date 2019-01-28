%% Solve the BVP with the finite difference method
%%   y'' = p(x) y' + q(x) y + r(x)
%%   y(a) = α, y(b) = β

%% Boundary value problem
a     = 0;
b     = 1;
alpha = 2;
beta  = 2;
p = @(x) -1./(x+1);
q = @(x) 2 + 0*x;
r = @(x) (1-x.^2).*exp(-x);

%% Parameters
N  = 10-1;
dx = (b-a)/(N+1);
x  = a+dx*(1:N);
% explanation
%------------
% 0      1      2                          N     N+1
% |------|------|------| ... |------|------|------|
% a     a+dx   a+2dx                              b


%% Initialization
A = diag(2+q(x)*dx^2)+ diag(-(1-p(x(1:x(1:N-1)))/2*dx),1) + diag(-(1+p(x(1:x(2:N)))/2*dx),-1);

vecB = -dx^2*r(x');
%% boundary condition
vecB(1) = vecB(1)+(1+p(x(1))/2*dx)*alpha;
vecB(N) = vecB(N)+(1-p(x(N))/2*dx)*beta;


%%-----       solve        ---%%
%%----------------------------%%
Y = A\vecB;
% add the boundary condition
y_sol = [alpha; Y; beta];


%%----- plot
plot([a x b],y_sol')
xlabel('x')
ylabel('y')
title('y''''=p y + q y'' + r, y(a)=alpha,y(b)=beta')
