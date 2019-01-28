%% Solve the BVP with the finite difference method
%%   y'' = p(x) y' + q(x) y + r(x)
%%   y(a) = α, y(b) = β

%% Boundary value problem
a     = 0;
b     = 1;
alpha = 0;
beta  = 2;
p = @(x) 0*x;
r = @(x) 4*x;
y_exact= @(x) (exp(2)/(exp(4)-1))*(exp(2*x)-exp(-2*x))+x;

%% Parameters
N  = 100-1;
dx = (b-a)/(N+1);
x  = a+dx*(1:N);
% explanation
%------------
% 0      1      2                          N     N+1
% |------|------|------| ... |------|------|------|
% a     a+dx   a+2dx                              b


%% Initialization
A = -diag(1+dx/2*p(x(2:N)),-1) + diag(2+dx^2*4*ones(1,N)) - diag(1-dx/2*p(x(1:(N-1))),1);
vecB = -dx^2*r(x');
%% boundary condition
vecB(1) = vecB(1) + (1+dx/2*p(x(1)))*alpha;
vecB(N) = vecB(N) + (1-dx/2*p(x(N)))*beta;


%%-----       solve        ---%%
%%----------------------------%%
Y = A\vecB;
% add the boundary condition
y_sol = [alpha; Y; beta];
y_sol_exact=[alpha y_exact(x) beta];
error=y_sol'-y_sol_exact;
%%----- plot
figure
hold on
plot([a x b],y_sol','r')
plot([a x b],y_sol_exact,'c')
xlabel('x')
ylabel('y')
title('y''''=p y + q y'' + r, y(a)=alpha,y(b)=beta')

figure 
loglog([a x b],error)
title('Error/Accuracy')

