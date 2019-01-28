function fpm=fixed(f,x,tol,N) 
i=1 
y=feval(f,x) 
if y==x
    fprintf('fixed point is %f', y)
end
while abs(x-y)>tol && i+1<=N 
    i=i+1
    x=y;
    y=feval(f,x)
end
end