function ivp2(steps,t0,tf,y0)
%ivp2(steps,t0,tf,y0) solves the scalar initial value problem dy/dt = -y
%   with y(t0) = y0 using backward Euler and TR
%   Usage: ivp1(10,0,2,1)

dt = (tf-t0)/steps;
y1 = zeros(steps+1,1); % set y to be a steps+1 dimensional column vector
y2 = zeros(steps+1,1);
t_pts = linspace(t0,tf,steps+1);

y1(1) = y0; y2(1) = y0;
for n = 1:steps % timestep loop
    y1(n+1) = y1(n)/(1+dt); % backward Euler
    y2(n+1) = (1-dt/2)*y2(n)/(1+dt/2); % trapezoidal rule (TR)
end

plot(t_pts,y1,'b-',t_pts,y2,'c-',t_pts,exp(-t_pts),'k.')
xlim([t0 tf])
xlabel('t','FontSize',16)
ylabel('y','FontSize',16)
legend('BE','TR','exact','Location','NorthEast')

end