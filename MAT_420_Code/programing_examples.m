function programing_examples
clc
format compact

pb = 1;

switch pb
    case 1 % power
        format long
        set(0,'RecursionLimit',1021)
        A = [1 1;1 0];
        r = max(abs(eig(A)));
        A = A/r;
        n = 1019;
        tic, power_recursive(A,n), toc
        tic, power_divideconquer(A,n), toc
        tic, power_mixed(A,n), toc
        tic, A^n, toc % optimized in Matlab
        (eye(2)+[1,2;2,-1]/sqrt(5))/2
        tic, power_direct(A,n), toc % no recursion
        tic, power_mixed_direct(A,n), toc % no recursion
    case 2 % recursion
        set(0,'RecursionLimit',1026)
        n = 1019;
        if n<30, tic, fibonacci_recursive(n), toc, end % impractical
        tic, fibonacci_recursive2(n), toc
        tic, fibonacci_direct(n), toc
        tic, fibonacci_direct2(n), toc 
        round(((1+sqrt(5))/2)^n/sqrt(5))
    case 3 % polynomial evaluation
        format long
        a = [0,1./(1:100000).^2];
        x = [-1,-.5,.5,1];
        tic, peval(a,x), toc
        tic, horner(a,x), toc
        tic, polyval(a(end:-1:1),x), toc
        disp(['-pi**2/12 = ',num2str(-pi^2/12,'%17.15f'),'  ', ...
            'pi**2/6 = ',num2str(pi^2/6,'%17.15f')])
    case 4 % bisection on f(x)=0
        format long
        f = @(x)x^2-2;
        ab = [1,2]; % assumes fa*fb<0
        tol = 1e-15;
        tic, bisect(f,ab,tol,0.5), toc
        tic, bisect(f,ab,tol,0.25), toc
        tic, bisect(f,ab,tol,0.75), toc
        tic, bisect(f,ab,tol,Inf), toc
        disp(['sqrt(2) = ',num2str(sqrt(2),'%17.15f')])
    case 5 % sorting algorithms
        x = round(10*rand(1,10))
        insertion(x)
        insertion_recursive(x)
        mergesort(x)
        quicksort(x)
        sort(x)
        % timing
        x = rand(1,20000);
        tic, insertion(x); toc
        tic, mergesort(x); toc
        tic, quicksort(x); toc
        tic, sort(x); toc
    case 6 % koch fractal
        ex = 5;
        maxlev = 7;
        switch ex
            case 1
                p = [0,1;0,0];
                t = [0,1/3,1/2,2/3,1;0,0,sqrt(3)/6,0,0];
            case 2
                p = [0,1;0,0];
                t = [0,1/3,1/3,2/3,2/3,1;0,0,1/3,1/3,0,0];
            case 3
                p = [0,1,.5,0;0,0,sqrt(3)/2,0];
                t = [0,1/3,1/2,2/3,1;0,0,sqrt(3)/6,0,0];
            case 4
                p = [0,.5,1,0;0,sqrt(3)/2,0,0]; % same as case 3 reversed
                t = [0,1/3,1/2,2/3,1;0,0,sqrt(3)/6,0,0];
            case 5
                npts = 9;
                a = linspace(0,2*pi,npts);
                p = [cos(a);sin(a)];                
                t = [0,1/3,1/3,2/3,2/3,1;0,0,1/3,1/3,0,0];
            case 6
                npts = 7;
                a = linspace(0,2*pi,npts);
                p = [cos(a);sin(a)];                
                t = [0,1/3,1/2,2/3,1;0,.1,.6,.1,0];
            case 7
                npts = 6;
                a = linspace(0,2*pi,npts);
                p = [cos(a);sin(a)];                
                x = linspace(0,1,7);
                t = [x;4*x.*(1-x)];
        end
        figure(1)
        set(gcf,'color','white')
        pp = plot(p(1,:),p(2,:),'b');
        axis equal
        hold on
        %%plot(.2*t(1,:),.2*t(2,:)+.7,'r');
        hold off
        axis off%%
        %%pause(5)
        for lev = 1:maxlev
            q = koch(lev,p,t);
            set(pp,'xdata',q(1,:),'ydata',q(2,:))
            %%pause(2)
        end
    case 7 % generic partitioning
        p = 8;
        n = 2^p;
        m = 100;
        F = linspace(0,1,m+1);
        T = zeros(m+1,3);
        for beta = 0:2
            for k = 0:m
                t = tic;
                partition(n,k/m,beta);
                T(k+1,beta+1) = toc(t);
            end
        end
        figure(1)
        set(gcf,'DefaultAxesfontSize',18)
        semilogy(F,T,'.-','markersize',20)
        ylim([1e-3,1e0])
        grid on
        legend({'$\beta = 0$','$\beta = 1$','$\beta = 2$'},'location','Best','interpreter','latex')
        xlabel('$\alpha$','interpreter','latex')
        ylabel('$W_n$ with $n=2^8$','interpreter','latex')
        
        n = 2.^(4:12);
        T = zeros(length(n),3);
        for beta = 0:2
            for k =1:length(n)
                nk = n(k);
                t = tic;
                partition(nk,.5,beta);
                T(k,beta+1) = toc(t);
            end
        end
        figure(2)
        set(gcf,'DefaultAxesfontSize',18)
        loglog(n,T,'.-','markersize',30)
        grid on
        legend({'$\beta = 0$','$\beta = 1$','$\beta = 2$'},'location','Best','interpreter','latex')
        xlabel('$n$','interpreter','latex')
        ylabel('$W_n$','interpreter','latex')
    case 8 % tridiagonal linear system
        n = 7;
        e = ones(n,1);
        A = spdiags([-e,2*e,-e],[-1,0,1],n,n);
        b = e;
        [solve_direct(A,b), solve_recursive(A,b), solveRB_recursive(A,b), A\b]
        % timing
        n = 2^10-1;
        e = ones(n,1);
        A = spdiags([-e,2*e,-e],[-1,0,1],n,n);
        b = e;
        set(0,'RecursionLimit',1025)
        tic, solve_direct(A,b); toc
        tic, solve_recursive(A,b); toc
        tic, solveRB_recursive(A,b); toc
        tic, A\b; toc
    case 9 % FFT
        format long
        x = (0:7)';
        [dft_direct(x), dft_recursive(x), fft(x)]
        % timing
        n = 2^10;
        x = rand(n,1);
        tic, dft_direct(x); toc
        tic, dft_recursive(x); toc
        tic, fft(x); toc
end

function a = power_recursive(a,n)
if n==1, return; end
a = a*power_recursive(a,n-1);

function a = power_divideconquer(a,n)
if n==1, return; end
if mod(n,2)
    a = power_divideconquer(a,(n-1)/2)*power_divideconquer(a,(n+1)/2); % n odd
else
    a = power_divideconquer(a,n/2)^2; % n even
end

function a = power_mixed(a,n)
if n==1, return; end
if mod(n,2)
    a = a*power_mixed(a,n-1); % n odd
else
    a = power_mixed(a,n/2)^2; % n even
end

function an = power_direct(a,n)
an = a;
for k = 1:n-1
    an = a*an;
end

function an = power_mixed_direct(a,n)
steps = [];
while n>1
    if mod(n,2)
        n = n-1; steps = [steps, 1]; % n odd
    else
        n = n/2; steps = [steps, 0]; % n even
    end
end
an = a;
for k = length(steps):-1:1
    if steps(k)
        an = a*an; % odd step
    else
        an = an^2; % even step
    end
end

function F = fibonacci_recursive(n)
if n==0
    F = 0;
elseif n==1
    F = 1;
else
    F = fibonacci_recursive(n-1)+fibonacci_recursive(n-2);
end

function F = fibonacci_recursive2(n)
if n==0
    F = 0;
elseif n==1 || n==2
    F = 1;
else
    if mod(n,2)
        F = fibonacci_recursive2((n-1)/2)^2+fibonacci_recursive2((n+1)/2)^2;
    else
        F = fibonacci_recursive2(n/2+1)^2-fibonacci_recursive2(n/2-1)^2;
    end
        
end

function F = fibonacci_direct(n)
F = zeros(n+1,1);
F(2) = 1;
for k = 3:n+1
    F(k) = F(k-1)+F(k-2);
end
F = F(n+1); % storing all previous values not practical for large n

function F = fibonacci_direct2(n)
if n==0
    F = 0;
else
    Fn1 = 0; F = 1;
    for k = 2:n
        Fn2 = Fn1; Fn1 = F;
        F = Fn1+Fn2;
    end
end

function p = peval(a,x)
% naive evaluation of polynomial p = a(1)+a(2)*x+...+a(n+1)*x^n
n = length(a)-1;
p = 0*x+a(1);
for k = 1:n
    p = p + a(k+1)*x.^k;
end

function p = horner(a,x)
% horner's rule to evaluate p = a(1)+a(2)*x+...+a(n+1)*x^n
n = length(a)-1;
p = 0*x+a(n+1);
for k = n:-1:1
    p = p.*x+a(k);
end

function ab = bisect(f,ab,tol,alfa)
a = ab(1); b = ab(2);
if b-a<tol, return; end
fa = f(a); % assumes f(a)*f(b)<0
v = alfa;
if isinf(alfa), v = (1+2*rand)/4; end % random number in [1/4,3/4]
x = a+v*(b-a); fx = f(x);
if fx*fa < 0
    ab = [a,x];
else
    ab = [x,b];
end
ab = bisect(f,ab,tol,alfa);

function x = insertion(x)
n = length(x);
for i = 2:n
    j = i-1;
    xi = x(i);
    while j>=1 && xi<x(j), j = j-1; end
    x = [x(1:j),xi,x(j+1:i-1),x(i+1:end)];
end

function x = insertion_recursive(x)
if isempty(x), return; end
xi = x(end);
x = [insertion_recursive(x(1:end-1)), xi];
j = length(x)-1;
while j>=1 && xi<x(j), j = j-1; end
x = [x(1:j),xi,x(j+1:end-1)];

function x = mergesort(x)
n = length(x);
if n>1
    m = fix(n/2);
    x1 = mergesort(x(1:m));
    x2 = mergesort(x(m+1:end));
    x = merge(x1,x2);
end

function x = merge(x1,x2)
n1 = length(x1);
n2 = length(x2);
i1 = 0;
i2 = 0;
x = [x1,x2];
while i1<n1 || i2<n2
    if i1<n1 && i2<n2
        if x1(i1+1)<=x2(i2+1)
            x(i1+i2+1) = x1(i1+1);
            i1 = i1+1;
        else
            x(i1+i2+1) = x2(i2+1);
            i2 = i2+1;
        end
    elseif i1<n1
        x(i1+i2+1) = x1(i1+1);
        i1 = i1+1;
    elseif i2<n2
        x(i1+i2+1) = x2(i2+1);
        i2 = i2+1;
    end
end

function x = quicksort(x)
n = length(x);
if n>1
    m = floor(n/2);
    i1 = 1;
    i2 = 1;
    x1 = [];
    x2 = [];
    xm = x(m);
    for i = 1:n
        if i~=m
            if x(i)<xm
                x1(i1) = x(i);
                i1 = i1+1;
            else
                x2(i2) = x(i);
                i2 = i2+1;
            end
        end
    end
    x1 = quicksort(x1);
    x2 = quicksort(x2);
    x = [x1,xm,x2];
end

function p = koch(lev,p,t)
if lev==0, return ; end
p = koch(lev-1,p,t); % get pattern from lower level
n = size(p,2); m = size(t,2);
pp = repmat(p,m,1); % initialize coordinates of inserted points
for k = 2:n
    u = diff(p(:,[k-1,k]),1,2); % vector between consecutive points
    a = angle(complex(u(1),u(2))); % angle between vector and horizontal line
    c = cos(a); s = sin(a);
    R = [c,-s;s,c]; % rotation matrix
    q = R*(t*norm(u)); % scale and rotate line segment
    pp(:,k-1) = pp(:,k-1)+q(:); % shift line segment
end
p = [p(:,1),reshape(pp(3:2*m,1:n-1),2,(m-1)*(n-1))]; % delete duplicate points

function x = partition(n,alfa,beta)
if n>1
    n1 = round(alfa*(n-1)+1-alfa);
    partition(n1,alfa,beta);
    partition(n-n1,alfa,beta);
end
overhead(n,beta);

function overhead(n,beta)
s = 0;
for k = 1:n^beta
    s = s+(1/k)^2;
end
   
function b = solve_direct(A,b)
n = length(b);
for k = 1:n-1
    i = k; j = k+1:n;
    D = A(i,i); % D is 1x1!
    A(j,j) = A(j,j)-A(j,i)*(D\A(i,j));
    b(j) = b(j)-A(j,i)*(D\b(i));
end
x(n) = b(n)/A(n,n);
for k = n:-1:1
    i = k; j = k+1:n;
    D = A(i,i);
    b(i) = D\(b(i)-A(i,j)*b(j));
end

function x = solve_recursive(A,b)
n = length(b);
if n==1, x = b/A; return; end
i = 1; j = 2:n;
D = A(i,i); % D is 1x1!
B = A(j,j)-A(j,i)*(D\A(i,j));
d = b(j)-A(j,i)*(D\b(i));
y = solve_recursive(B,d);
x = [D\(b(i)-A(i,j)*y); y];

function x = solveRB_recursive(A,b)
n = length(b);
if n==1, x = b/A; return; end
i = 1:2:n; j = 2:2:n; % red/black ordering
D = A(i,i); % D is diagonal when A is tridiagonal
B = A(j,j)-A(j,i)*(D\A(i,j));
d = b(j)-A(j,i)*(D\b(i));
y = solveRB_recursive(B,d);
x = [D\(b(i)-A(i,j)*y); y]; % solution of Ax=b in red/black ordering
x([i,j]) = x; % re-order unknowns back

function x = dft_direct(x)
n = length(x);
z = exp(complex(0,-1)*2*pi*(0:n-1)/n);
F = vander(z);
x = F*x(end:-1:1);

function x = dft_recursive(x)
n = length(x);
if n==1, return; end
y1 = dft_recursive(x(1:2:n)); % F_{n/2}x_{odd}
y2 = exp(complex(0,-1)*pi*(0:2:n-2)/n).'.*dft_recursive(x(2:2:n)); % D_{n/2}F_{n/2}x_{even}
x = [y1+y2;y1-y2];
