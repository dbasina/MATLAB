%% Solve the BVP using Finite-Element Method
%%   -(py?)? + q(x) y = f(x)
%%   y(0)=y(1)=0
%%---- The Boundary Value Probem (BVP)
p = @(x) x.^2
q = @(x) 2+(x*0)
f = @(x) -4*x.^2
%%--- numerical parameters: vector space V
N = 10
dx = 1/(N+1);
%%--- estimation matrix
A = zeros(N,N);
b = zeros(N,1);
%%-- construction of the matrix ?A? and vector ?b? to estimate the integrals
intX = linspace(0,1,10000); 
for i=1:N
    % a) phi_i and phi_i?
    [phi_i,dphi_i] = Basis_Phi(intX,i,N); 
    % b) b_i = L(phi_i)
    b(i) = trapz(intX,f(intX).*phi_i);
    % c) estimation of A(phi_i,phi_j) 
    for j=1:N
        % phi_j and phi_j?
       [phi_j,dphi_j] = Basis_Phi(intX,j,N); 
       % a_ij = A(phi_i,phi_j)
       A(i,j) = trapz(intX,p(intX).*dphi_i.*dphi_j +q(intX).*phi_i.*phi_j);
    end
     
    
end
%%-------------------------------%%
%%--- Resolution system ---%%
%%-------------------------------%% 
c = A\b;
%%-- construction solution
%% y = c(1)*phi_1 + ... + c(N)*phi_N 
ySol = zeros(1,length(intX));
for i=1:N
[phi_i,dphi_i]= Basis_Phi(intX,i,N);
ySol= ySol + c(i)*phi_i;
end
%%-- plot
yExact = @(x) x.*(x-1);
plot(intX,ySol,intX,yExact(intX))
legend('Numerical Solutio','Exact','location','southwest')
function [phi,dphi]=Basis_Phi(x, i,N)
%% The piece-wise linear basis function %%
%%                  |  (x-x_im)/dx  if x_im < x <= x_i
%%    phi(x)  =     |  (x_ip-x)/dx  if x_i <x< x_ip 
%%                  |       0       otherwise
%%
%%                  |      1/dx     if x_im < x <= x_i
%%    phi?(x) =     |     -1/dx     if x_i <x< x_ip 
%%                  |       0       otherwise
%% init
  dx = 1/(N+1);
  x_im = (i-1)*dx;
  x_i  =     i*dx;
  x_ip = (i+1)*dx;
%% value
phi = (x-x_im)/dx.*((x-x_im)>0).*((x-x_i)<=0) ...
     +(x_ip-x)/dx.*((x-x_i)>0).*((x-x_ip)<0);
dphi = 1/dx.*((x-x_im)>0).*((x-x_i)<=0) ...
      -1/dx.*((x-x_i)>0).*((x-x_ip)<0); 
end