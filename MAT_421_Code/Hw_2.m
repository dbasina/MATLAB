%% BVP y''= y^3-y*y'
%% We need to zoom in a lot to be able to see all the graphs
tic
colors=['b-','c-','r-','-g'];
iter=4
for i= 1:iter
N=10^(i)
h = 1/N; % h = dx
x = linspace(0,1,N+1);
y = zeros(N-1,1);

e = ones(N-1,1);
D2 = spdiags([e -2*e e],[-1 0 1],N-1,N-1)/h^2;

normresidual = Inf;
iter = 0;
while iter<5
    F = D2*y - diag(ones(1,N-1))*cos(y);
    J = D2  - spdiags(-sin(y),0,N-1,N-1);
    y = y - J\F;
    iter = iter + 1;
end
y = [0; y; 0]
toc
hold on
plot(x,y,colors(i),'LineWidth',2) % Use in turn 'b-','c-','r-' for 3 plots
xlabel('x','FontSize',24)
ylabel('y','FontSize',24)
end

