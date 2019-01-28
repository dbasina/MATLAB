%% Solve the Hyperbolic PDE
%%    d^2u/dt^2 = c^2 * d^2u/dx^2    ,  0<x<l, t>0
%%       u(x,0) = sin(2*pi*x)      ,  0<x<l, t=0
%%        du(x,0)/dt = 2*pi*sin(2*pi*x),0<x<l
%%     u(0,t) = u(l,t) = 0   ,     t>0


%%-- parameters
c =  1;
l =  1;
T = 10;
u0 = @(x) sin(2*pi.*x); 
du0=@(x) 2*pi*sin(2*pi.*x);
d2u0=@(x)4*pi^2*cos(2*pi.*x);
exact=@(x,t) sin(2*pi*x)*(cos(2*pi*t)+sin(2*pi*t))
%%-- numerical paramters
%% B1
%dx = 2*10^-1;
%dt = 10^-2;
%%B2
dx = 10^-1;
dt = 2*10^-2;
%% Difference between the parameters is that The first set of parameters has a large dx which produces a non-smooth solution curve. 
%%The second set of parameters produces a fairly smooth solution curve.

lambda = c^2*dt^2/dx^2;

%% time and space discretization
%Number of Time and Space points
nTime = floor(T/dt+.5);
nX    = floor(l/dx+.5);

%Space discretization Vector
intX  = linspace(0,l,nX);

%%Initialize solution matrix points to zeros
U=zeros(nTime,nX);

%% Set boundary at t=0 to u(x,0) = sin(2*pi*x);
U(nTime,2:nX-1)=u0(intX(2:nX-1));

%%Define the solution at the first time step on the interior points of the solution matrixu1=u0+ du0*dt + 0.5*dt^2*d2u0
U(nTime-1,2:nX-1)= u0(intX(2:nX-1)) + du0(intX(2:nX-1))*dt + 0.5 * dt^2 * d2u0(intX(2:nX-1));

%% Iterate Using the Discretization formula for the rest of the interior points from bottom to top.
%%Bottom of matrix U is analogous to the base of the logically discretized grid.

for n= nTime-2:-1:2
    U_new=U; 
    for i=2:1:nX-1
        U_new(n-1,i)=2*U(n,i)-U(n+1,i)+lambda*(U(n,i-1)-2*U(n,i)+U(n,i+1))
    end
    
    % Plot is done on flipped matrix since our 0 is at the bottom of the matrix and matlab starts from the top of the matrix.
   plot(intX,flipud(U_new(n,:)))
   xlabel('x')
   axis([0 1 -75 75])
   title(['solution at time t=' num2str(n*dt,'%1.2f')])
   %% update
   U = U_new;
   pause(.0001)
end
    








