function bvp(N)
%bvp(N) solves the linear boundary value problem 
%       y'' + 2y' + y = 0, y(0) = 1, y(1) = 0 
%   using central differences with N dx.
%   The exact solution is y(x) = (1 - x) exp(-x).
%   Try bvp(40).  For displaying full matrices, try bvp(4).

h = 1/N; % h = dx
y = zeros(N-1,1); % set y = N-1 dimensional column vector (interior pts)
xpts = linspace(0,1,N+1);

e = ones(N-1,1);
D2 = spdiags([e -2*e e],[-1 0 1],N-1,N-1)/h^2;
% full_h2_D2 = full(h^2*D2)
D1 = spdiags([-e e],[-1 1],N-1,N-1)/(2*h);
% full_2h_D1 = full(2*h*D1)
A = D2 + 2*D1 + speye(N-1);
% full_A = full(A)
b = [-1/h^2 + 1/h; zeros(N-3,1); 0];
y = A\b;

y = [1; y; 0]; % add in BCs to y
figure
plot(xpts,y,'r.',xpts,(1-xpts).*exp(-xpts),'b-','MarkerSize',24,...
    'LineWidth',2)
xlim([0 1])
xlabel('x','FontSize',24)
ylabel('y','FontSize',24)
legend('comp','exact','Location','NorthEast')