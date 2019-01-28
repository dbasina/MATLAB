function newton(x0)
%newton(x0) finds a root of the nonlinear function specified by 
%   f and fprime
%   TRY: newton(1), newton(-pi/2), newton(-1.57)
N  = 100-1;

A=-diag(ones(1,N-1),-1)  + diag(2*ones(1,N)) - diag(ones(1,N-1),1);
B=diag(ones(1,N-1),-1)  + diag(-ones(1,N-1),1);
epsilon = eps; % governs precision of convergence
% sometimes can't achieve eps; can take epsilon = 10^-6 for example
x = x0;

k = 0;
error=[]
while abs(f(x)) > epsilon*abs(f(x0)) && k < 20
    k = k+1
    xprev = x;
    x = x - f(x)/fprime(x)
    change = x - xprev
    residual = f(x)
    err=abs(x-0.567143290409784);
    error=[error err];
end

    function y = f(x) % nested subfunction
    y = A*x+cos(x);
    end

    function yprime = fprime(x) % nested subfunction
    yprime = x*exp(x)+exp(x);
    end
plot(error)
end