function burgers1(N,steps,CFL)
%burgers1(N,steps,CFL) solves the inviscid Burgers' equation 
%       u_t + f_x = 0, f(u) = u^2/2
%   with N+1 gridpoints for 0 <= x <= 1, h = 1/N, dt = CFL*h, and
%   t_f = steps*dt 
%   Uses conservative upwind method and through-flow BCs
%   Try: burgers1(200,120,0.9) for Riemann problem shock or rarefaction
%        burgers1(200,350,0.9) for merging shocks

h = 1/N;
x = linspace(0,1,N+1); t = 0;
u = zeros(N+1,1); % N+1 dimensional column vector

% ICs: Riemann problem shock wave
j = 1:N/2-1; u(j) = 2; 
u(N/2) = 1.5; 
j = N/2+1:N+1; u(j) = 1;

% ICs: Riemann problem rarefaction wave
% j = 1:N/2-1; u(j) = 1; 
% u(N/2) = 1.5; 
% j = N/2+1:N+1; u(j) = 2;

% ICs: merging shocks
% j = 1:N/4-1; u(j) = 3; 
% u(N/4) = 2.5;
% j = N/4+1:N/2-1; u(j) = 2;
% u(N/2) = 1.5;
% j = N/2+1:N+1; u(j) = 1;

for n = 1:steps % timestep loop
    umax = max(abs(u));
    dt = CFL*h/umax; 
    t = t + dt;
    % with thru-flow BCs
    % gridpts: 1 | 1 2 ... N N+1 | N+1
    uleft = [u(1);u(1:N)];

    % conservative upwind for u > 0: correct shock speed
    u = u - dt*(u.^2-uleft.^2)/(2*h);
    
    % nonconservative upwind: incorrect shock speed, incorrect rarefaction
    % speed, and incorrect unphysical rarefaction at CFL = 1
    % u = u - dt*u.*(u-uleft)/h;

    % plot(x,u,'r-'); axis([0 1 0 3.2]); % merging shocks
    plot(x,uexact(t),'b-',x,u,'r-'); axis([0 1 0 2.2]); % RP
    xlabel('x','FontSize',16)
    ylabel('u','FontSize',16)
    getframe;
end
tf = t

    function ue = uexact(t)
        ue = zeros(N+1,1);
        % shock wave
        s = 1.5;
        Ns = round(N/2 + s*t*N);
        Ns = min(N+1,Ns);
        i = 1:Ns; ue(i) = 2;  
        i = Ns+1:N+1; ue(i) = 1;
        % rarefaction wave
%         xleft = 0.5 + t;
%         xright  = 0.5 + 2*t;
%         for i = 1:N+1
%             if (x(i) <= xleft)
%                 ue(i) = 1;
%             elseif (x(i) >= xright)
%                 ue(i) = 2;
%             else
%                 ue(i) = (x(i)-0.5)/t;
%             end
%         end
    end

end