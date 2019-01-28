
nr = 80; % number of grids in the r direction
nz = 100; % nuber of grids in the z direction
h = 1; % height of the container 
a = 1; % radius of the container the 
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

 %Solution for Au+uB=F
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
U = Z*v';
bc = bc;


function dvdt=f(v)
nr = 100; % number of grids in the r direction
nz = 100; % nuber of grids in the z direction
h = 1; % height of the container 
a = 1; % radius of the container the 
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

F = zeros(n,1); % the forcing term

% first column is due to the boundary condition.
 for i = 1:n
     F(i,1) = dr*i;
 end
 F = -F/(dz^2);
 dvdt=A*v+B*v; 
 
end



    





    
