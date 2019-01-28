%% Solve the BVP with the finite difference method
%%   y'' = f(y)
%%   y(a) = α, y(b) = β

%% Boundary value problem
a = 1;
b = 2;
alpha = 1/2;
beta  = 1/3;
N  = 100-1;
dx = 1/(N+1);
x  = a+dx*(1:N);
F  = @(y) y.^3-y.*y';
DF = @(y) diag(3*y.^2-y.*F(y)+F(y));

%% Parameters

% explanation
%------------
% 0      1      2                          N     N+1
% |------|------|------| ... |------|------|------|
% a     a+dx   a+2dx                              b

%% numerical parameter
nbIterations = 5;
stockY = zeros(nbIterations+1,N);

%% Initialization
Y = 1+x';			% choice first iteration
L = 2*diag(ones(1,N)) - diag(ones(1,N-1),1) - diag(ones(1,N-1),-1);
vecBC = [alpha; zeros(N-2,1); beta];
stockY(1,:) = Y';


%%-----       loop        ---%%
%%---------------------------%%
C=F(Y)
for k=1:nbIterations
  %% Newton method
  DJ = L + dx^2*DF(Y);
  Y  = Y - DJ\ (L*Y + dx^2*C- vecBC);
  %% save  
  stockY(k+1,:) = Y';
end

%% plot
figure(1)
plot([a x b],[alpha Y' beta])
xlabel('x')
legend('y5','location','northwest')
title('y''=y^2, y(0)=1,y(1)=2')

figure(2)
plot([a x b],[alpha stockY(1,:) beta],...
     [a x b],[alpha stockY(2,:) beta],...
     [a x b],[alpha stockY(3,:) beta],...
     [a x b],[alpha stockY(4,:) beta])
xlabel('x')
title('y''''=y^2, y(0)=1,y(1)=2')
legend('Y0','Y1','Y2','Y3','location','northwest')
