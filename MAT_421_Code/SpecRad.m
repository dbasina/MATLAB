%% For N=4,h=1/5 the spectral rad computation is =cos(pi/5) as shown below.

%setup the diagonal matrix
h=1/5
ones=[-1,-1,-1,-1];
twos=[2,2,2,2,2];
A=diag(ones,-1)+diag(ones,1)+diag(twos)

%LDU Factorization of A, extract D from U
[L,U]=lu(A);       
D=diag(diag(U)); 

% Set Mj(Jacobi Matrix) 
Mj=D;

%Compute B, compute eigenvalues and check that max of eigenvlaues is equal
%to cos(pi/5)
B=inv(Mj)*(Mj-A);
[Z,E]=eig(B);
SpectralRadius=E(1)   %Max Eigen val of B is the Jacobi spectral radius.Cos Pi is given under for comparision
cosine=cos(pi/5)