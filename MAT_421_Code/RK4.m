function [y,t]=RK4(f,y0,Time,dt)
  %%
  %% The Runge-Kutta method of order 4
  %%    0   |
  %%   1/2  |  1/2
  %%   1/2  |      1/2
  %%    1   |           1
  %%        -------------------
  %%           1/6 1/3 1/3 1/6  
  %%   
  %%
  %% Input:
  %%        . f: vector fields
  %%        . y0: initial condition
  %%        . Time: final time of simulation
  %%        . dt: time step
  %% Output:
  %%        . t: time at which the solution is estimated
  %%        . y: estimation of the solutions. Its size is y(d,nTime+1)
  %%             with d dimension of the system and nTime number of time steps
  
  %%-----  Initialisation   -----%%
  d       = length(y0);		% dimension of y0
  t       = 0:dt:Time;
  nTime   = length(t) -1;	% number of ∆t
  y       = zeros(d,nTime+1);
  y(:,1)  = y0;


  %%---------------------------------------%%
  %%------------      Loop      -----------%%
  %%---------------------------------------%%

  for n=1:nTime
    %% init
    tn = t(n);
    yn = y(:,n);
    %% The slopes k1,k2,k3,k4
    k1 = f(tn,yn);
    k2 = f(tn+dt/2 , yn + k1*dt/2);
    k3 = f(tn+dt/2 , yn + k2*dt/2);
    k4 = f(tn+dt   , yn + k3*dt);
    %% RK4
    y(:,n+1) = yn + dt*(k1/6 + k2/3 + k3/3 + k4/6);
  end
 
  %%---------------------------------------%%
  %%---------------------------------------%%

  %%--- one last step if t(n)~Time
  if (t(end)<Time)
    %% last step
    tn = t(nTime+1);
    yn = y(:,nTime+1);
    dtLast = Time-t(end);
    %% The slopes k1,k2,k3,k4
    k1 = f(tn,yn);
    k2 = f(tn+dtLast/2 , yn + k1*dtLast/2);
    k3 = f(tn+dtLast/2 , yn + k2*dtLast/2);
    k4 = f(tn+dtLast   , yn + k3*dtLast);
    %%-- final value
    yLast = yn + dtLast*(k1/6 + k2/3 + k3/3 + k4/6);
    %% update
    t = [t Time];
    y = [y yLast];
  end

end