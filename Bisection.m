function bisection(a,b)
%bisection(a,b) finds a root of the nonlinear function specified by f
%   Slow, but is guaranteed to find a root bracketed by a and b
%   TRY: bisection(0,1)

epsilon = eps; % governs precision of convergence
k = 0;
error=[]
while abs(b-a) > epsilon*(abs(a)+abs(b))/2 && k < 200
    k = k+1;
    x = (a+b)/2;
    err=abs(0.567143290409784-x)
    error=[error err];
    if sign(f(x)) == sign(f(b))
        b = x;
    else
        a = x;
    end
end
residual = f(x)

    function y = f(x)
    y = (x*exp(x)) -1;
    end
plot(error)
end