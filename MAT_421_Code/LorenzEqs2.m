function LorenzEqs(tf)
%LorenzEqs(tf) solves the Lorenz equations to tf 
%   with x(0) = x0, y(0) = y0, z(0) = z0 given below,
%   using fourth/fifth-order Runge-Kutta or TRBDF2.
%   Here w(t) = [x(t), y(t), z(t)] and the Lorenz equations are 
%       dw/dt = F(w)
%   The Jacobian for TRBDF2 is supplied in J
%   Try: LorenzEqs(40)

tic

tspan = [0 tf];
sigma = 10; r = 28; b = 8/3;
% equilibria
eta = sqrt(b*(r-1)); % or = -sqrt(b*(r-1))
xeq = eta; 
yeq = eta; 
zeq = r-1;

%w0 = [0; 1; 0]; % standard ICs
delta = 10^-3*ones(3,1); w0 = [xeq; yeq; zeq] + delta % equilibrium ICs
w0 = [0; eps; 0]; % near very unstable equilibrium at origin

% Runge-Kutta 4/5
options = odeset('RelTol',10^-12,'AbsTol',10^-15);
[t,w] = ode45(@F,tspan,w0,options);

% TRBDF2 with analytical Jacobian
% options = odeset('RelTol',10^-6,'AbsTol',10^-9,'Jacobian',@jacobian);
% [t,w] = ode23tb(@F,tspan,w0,options);

% TRBDF2 with MATLAB calculated "finite-difference" Jacobian:
% options = odeset('RelTol',10^-6,'AbsTol',10^-9);
% [t,w] = ode23tb(@F,tspan,w0,options);

toc

figure
plot(t,w(:,2),'b-','LineWidth',2)
xlabel('t','FontSize',24)
ylabel('y','FontSize',24)



    function wprime = F(t,w)
    wprime = [sigma*(w(2)-w(1))
         r*w(1)-w(2)-w(1)*w(3)
         w(1)*w(2)-b*w(3)];
    end

    function J = jacobian(t,w)
    J = [-sigma   sigma   0
         r-w(3)   -1      -w(1)
         w(2)     w(1)    -b   ];
    end

end