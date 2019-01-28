%% BVP y''= y^3-y*y'
tic
N=10
h = 1/N; % h = dx
x = linspace(1,2,N+1);
y = zeros(N-1,1);

e = ones(N-1,1);
D2 = spdiags([e -2*e e],[-1 0 1],N-1,N-1)/h^2;
D1 = spdiags([-e e],[-1 1],N-1,N-1)/(2*h);

normresidual = Inf;
iter = 0;
while iter<5
    F = -D2*y +D1*y.^2 -diag(ones(1,N-1))*y.^3;
    J = -D2 + 2*D1*y - spdiags(3*y.^2,0,N-1,N-1);
    y = y - J\F;
    iter = iter + 1;
end
y = [1/2; y; 1/3];
toc
figure
plot(x,y,'b-','LineWidth',2) % Use in turn 'b-','c-','r-' for 3 plots
xlabel('x','FontSize',24)
ylabel('y','FontSize',24)


