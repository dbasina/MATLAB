function RK2(w0,t0,tf)
%RK2(w0,t0,tf) solves dw/dt = f(w) from t0 to tf
%   with the column vector of ICs w(t0) = w0 
%   using second-order Runge-Kutta with a fixed timestep.
%   The function f(w) is implemented below for 2 first-order ODEs.
%   For (damped) harmonic oscillator: RK2([1; -0.05],0,8*pi) with N = 200
%   For van der Pol oscillator: RK2([1; 0],0,50) with N = 2000

tic
N = 2000; % number of timesteps
p = 0.05; omega2 = 1; % for damped harmonic oscillator

w = w0(:);
wsol = transpose(w); % store solution in rows
dt = (tf-t0)/N;
tsol = linspace(t0,tf,N+1);

for n = 1:N % timestep loop
    w = w+dt*f(w+dt*f(w)/2); % RK2
    wsol(end+1,:) = transpose(w);
end
toc
fprintf('\nNumber of steps = %d\n',N)

figure
plot(tsol,wsol(:,1),'r-','LineWidth',2)
% For damped harmonic oscillator, use:
% plot(tsol,exp(-p*tsol).*cos(sqrt((-p^2+omega2))*tsol),'r-',...
%     tsol,wsol(:,1),'b.','MarkerSize',12,'LineWidth',2)
% legend('exact','RK2','Location','NorthEast')
xlim([t0 tf])
xlabel('t','FontSize',24)
ylabel('y','FontSize',24)

figure
plot(wsol(:,1),wsol(:,2),'b-',wsol(1,1),wsol(1,2),'r.',...
    'MarkerSize',24,'LineWidth',2)
axis square
xlabel('y','FontSize',24)
ylabel('dy/dt','FontSize',24)
title('attractor','FontSize',24)

    function wprime = f(w)
    % wprime = [w(2); -w(1)]; % harmonic oscillator
    % wprime = [w(2); -omega2*w(1)-2*p*w(2)]; % damped harmonic oscillator
    wprime = [w(2); -w(1)+(4-w(1)^2)*w(2)]; % van der Pol oscillator
    end

end