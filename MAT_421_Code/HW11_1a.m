%% Solve the parabolic PDE
%%    du/dt = alpha^2 d^2U/dt^2    ,  0<x<l, t>0
%%       u(x,0) = sin(pi*x/l)      ,  0<x<l, t=0
%%     u(0,t) = u(l,t) = 0   ,     t≥0


%%-- parameters

%%1b
%%T      = 100;
T=1
alpha2 = 0.5;
l      = 10;
u0 = @(x) sin(pi.*x/l);
%% 1C
%%u0 = @(x) 1+sin(2*pi.*x/l);

%%-- numerical paramters
dx = .5;
dt=0.005;
%%1b Parameters
%%dt = 0.25;
ld = alpha2/dx^2*dt;

%% init
nTime = floor(T/dt+.5);
nX    = floor(l/dx+.5);

%% IC
intX  = linspace(0,l,nX);
u     = u0(intX);
%%Store Mass loss
int=[];
col = linspace(1,200,200);

%%----------------------------------%%
%%-------      loop         --------%%
%%----------------------------------%%

for n=1:nTime
   %% new time step
   u_new = u;
   for i=2:(nX-1)
     u_new(i) = (1-2*ld)*u(i) + ld*(u(i-1)+u(i+1));
     
   end
   %% plot
   plot(intX,u_new)
   xlabel('x')
   axis([0 l 0 2])
   title(['solution at time t=' num2str(n*dt,'%1.2f')])
   %% update
   u = u_new;
   pause(.0001)
 end

for n=1:nTime
  %% new time step
  u = (1-2*ld)*u + ld*([u(2:end) 0] + [0 u(1:(end-1))]);
  int = [int trapz(intX,u)]; % Use this to show numerical loss of mass
 
end

%% Plot the evolution of mass over time
figure
plot(col,int,'r','linewidth',3), axis tight
title('Mass over time')





