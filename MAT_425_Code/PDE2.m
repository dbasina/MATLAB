%Discretized Laplacian Solution
% Problem Statement: 
% d^2v/dr^2 + 1/r dv/dr - v/r^2 + d^2v/dz^2 =0
% BC:
% V(0,Z)=V(1,Z)=0   0<=Z<=h
% V(r,h)=0          0<=r<=1
% V(r,0)=r          0<=r<=1

%Variables
Nr=40; %number of grid points in r
Nz=30; % number of grid points in z
dr=1/(Nr-1); %discretization of r and z
dz=1/(Nz-1);

%AU= d^2v/dr^2 + 1/r dv/dr - v/r^2
%UB= d^2v/dz^2 +F
% AU+UB=F

A=zeros(Nr);  
B=zeros(Nz);

%Construct A-Matrix using equations obtained for tridiagonal form
%A middle Diagonal
for i=1:Nr
    A(i,i)=-2-(1/i^2);
end
%A Upper Diagonal
for i= 2:Nr
    A(i-1,i)=1+(1/(2*(i-1)));           % no dr in there???
end
%A lower Diagonal
for i=1:Nr-1
    A(i+1,i)=1-(1/(2*(i+1)));           % same question
end
A=A/dr^2;

%Construct B-Matrix using constants obtained for tridiagonal form
%B middle Diagonal 
for j=1:Nz
    B(j,j)=-2;
end
%B Upper Diagonal
for j=2:Nz
    B(j-1,j)=-1;
end
%B Lower Diagonal
for j=1:Nz-1
    B(j+1,j)=-1;
end
B=B/(dz^2);

% AU+UB=F where F represents the boundary conditions in the first column.
F=zeros(Nr,Nz);
for k=1:Nr
    F(k,1) = k*dr;
end 
% Decompose A: A=PDP^1 thus
%   DV + VB = P^-1 F
% with V=P^-1 U
[Z,E]=eig(A);
H = transpose( inv(Z)*F );
%B=transpose(B);
 
% find V_transpose
V_t = zeros(Nz,Nr);
for i=1:Nr
    %V(:,i) = linsolve(B+(E(i,i)*eye(Nz)),H(:,i));
    V_t(:,i) = (B+(E(i,i)*eye(Nz)))\H(:,i)
end
% deduce U (U=PV)
U = Z*transpose(V_t);
% plot
intR = (1:Nr)*dr;
intZ = (1:Nz)*dz;
surf(intR,intZ,U')
xlabel('r')
ylabel('z')
