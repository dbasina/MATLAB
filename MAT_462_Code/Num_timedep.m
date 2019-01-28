clear,clc;
clear all;
h = 1;
a = 1;
nr = 100;
nz = 80;
% use the length scale a to diagonlize lengths this should give (h/a)
dr = 1/(nr); % small values of r between points
dz = (h/a)/(nz); % small values of z between points
n = nr-1; % interior in the r direction
m = nz-1; % interior in the z direction

% diagonal part of the matrix A for the r part of the PDE
    for i = 1:n
    A(i,i) = -2 -1/i.^2;
    end
% lower diagonal of the matrix A for the r part of the PDE
    for i = 2:n
    A(i,i-1) = 1-1/(2*i);
    end

% upper diagonal of the matrix A for the r part of the PDE
    for i  = 1:n-1
    A(i,i+1) = 1+1/(2*i);
    end

    A = A/(dr^2); % this is the complete matrix for the r part of the PDE

% diagonal matrix of B
    for j = 1:m
    B(j,j) = -2;
    end

% lower diagonal of matrix B
    for j = 2:m
    B(j,j-1) = 1;
    end

% upper diagonal of matrix B
    for j = 1:m-1
    B(j,j+1) = 1;
    end

    B = B/(dz^2); % this is the complete matrix for B


    F = zeros(n,m); % the forcing term

% first column is due to the boundary condition.
    for i = 1:n
     F(i,1) = dr*i;
    end
    F = -F/(dz^2);
    
 [Z,E] = eig(A);
 Zinv = inv(Z)';
 Ft = F';
H = Ft*Zinv; 
for i = 1:n
    v(:,i) = (B+E(i,i)*eye(m))\H(:,i);
end 

for i = 1:n
    bc(i) = i*dr;
end
u = Z*v';
nt = 100;
tspan = linspace(0,0.0025,nt);
dt = tspan(2)-tspan(1);
U = zeros(n,m,nt-1);
for i = 2:nt
    U(:,:,i) = U(:,:,i-1)+ dt*(A*U(:,:,i-1)+U(:,:,i-1)*B+F);
    Unew = U(:,:,i-1)+(dt/2)*( (A*U(:,:,i)+U(:,:,i)*B+F));
end
% h = uicontrol(gcf,'style','slider','units','pix','position',[100 5 300 20]);
% set(h,'min',min(tspan),'max',max(tspan) );
% set(h,'callback','i=find(t==nearest(get(h,"value"))')
contour(U(:,:,40)')
surfc(-1*U(:,:,100)')
