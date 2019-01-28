function x = tridisolve(a,b,c,d)
x = d;
n = length(x);
% forward elimination
for j = 1:n-1                   
    mu = a(j)/b(j);             %(n-1) [divisions]
    b(j+1) = b(j+1) - mu*c(j);  % (n-1)  [ substractions and multiplications ]
    x(j+1) = x(j+1) - mu*x(j);  % (n-1) [ substractions and multiplications ]
end
% back solve
x(n) = x(n)/b(n);               %n [divisions]                       
for j = n-1:-1:1                
    x(j) = (x(j)-c(j)*x(j+1))/b(j); % (n-1) operations
end
end                                 %Total 4(n-1)+n =5n-1 ~ 5n
