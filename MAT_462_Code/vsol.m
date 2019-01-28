%%Grid initialization. 
Nr=10; %number of grid points in r
Nz=8; % number of grid points in z
dr=1/Nr; %discretization of r and z
dz=1/Nz;

%% Function for d(v)/dt 
v_t=@(v,psi,i,j) (((psi(i,j+1)*v(i+1,j))-(psi(i,j+1)*v(i-1,j))-(psi(i,j-1)*v(i+1,j))+(psi(i,j-1)*v(i-1,j)))/2*dz^2*dr)... %% Term 1
                  - (((psi(i+1,j)*v(i,j+1))-(psi(i+1,j)*v(i,j-1))-(psi(i-1,j)*v(i,j-1))+(psi(i-1,j)*v(i,j-1)))/2*dz*dr^2)... %% Term 1
                  +(v(i,j)*(psi(i,j+1)-psi(i,j-1))/2*dz*(i*dr)^2)... %% Term 2
                  + (((v(i,j+1)-2*v(i,j)+v(i,j-1))/dz^2)+((v(i+1,j)-2*v(i,j)+v(i-1,j))/dr^2)+((v(i+1,j)-v(i-1,j))/2*i*dr^2)-(v(i,j)/(i*dr)^2))/Re  %% Term 3
%% Function for d(ita)/dt
ita_t=@(v,psi,ita,i,j) ((psi(i,j+1)-psi(i,j-1))/2*dz)* (((psi(i+1,j)-psi(i-1,j))/2*dz)-(ita(i,j)/(i*dr)^2))... %% Term 1
                     - (((psi(i+1,j)-psi(i-1,j))/2*dr)*((ita(i,j+1)-ita(i,j-1))/2*i*dr*dz))...%% Term 1
                    + (2*v(i,j)*(v(i,j+1)-v(i,j-1))/2*i*dr*dz)...% Term 2
                    +(((ita(i,j-1)-2*ita(i,j)+ita(i,j+1))/dz^2)+((ita(i-1,j)-2*ita(i,j)+ita(i+1,j))/dr^2)+((ita(i+1,j)-ita(i-1,j))/2*dr^2*i)-(ita(i,j)/(i*dr)^2))/Re %% Term 4

%% Set up psi function 

%Apsi= d^2psi/dr^2 - 1/r d(psi)/dr 
%psiB= d^2v/dz^2
%F=-r*ita
%Apsi + psiB= -r*ita
% AU+UB=-r*ita

A=zeros(Nr-1);  
B=zeros(Nz-1);

%Construct A-Matrix using equations obtained for tridiagonal form
%A middle Diagonal
for i=1:(Nr-1)
    A(i,i)=-2;
end
%A Upper Diagonal
for i=1:(Nr-2)
    A(i,i+1)=1+(1/(2*i));           % no dr in there???
end
%A lower Diagonal
for i=2:(Nr-1)
    A(i,i-1)=1-(1/(2*i));           % same question
end
A=A/(dr^2);

%Construct B-Matrix using constants obtained for tridiagonal form
%B middle Diagonal 
for j=1:(Nz-1)
    B(j,j)=-2;
end
%B Upper Diagonal
for j=2:(Nz-1)
    B(j-1,j)=1;
end
%B Lower Diagonal
for j=1:(Nz-2)
    B(j+1,j)=1;
end
B=B/(dz^2);

% AU+UB=F where F represents the boundary conditions in the first column.
F=zeros(Nr-1,Nz-1);
for k=1:(Nr-1)
    F(k,1) = k*dr.*ita(k,:);
end 
F=-F/(dz^2);
% Decompose A: A=PDP^1 thus
%   DV + VB = P^-1 F
% with V=P^-1 U
[Z,E]=eig(A);
H = transpose( inv(Z)*F );
%B=transpose(B);
 
% find V_transpose
V_t = zeros(Nz-1,Nr-1);
for i=1:(Nr-1)
    %V(:,i) = linsolve(B+(E(i,i)*eye(Nz)),H(:,i));
    V_t(:,i) = (B+(E(i,i)*eye(Nz-1)))\H(:,i);
end
% deduce U (U=PV)
U = Z*transpose(V_t);
        

