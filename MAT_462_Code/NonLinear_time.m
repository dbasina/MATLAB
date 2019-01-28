%% Main Script
clear,clc;
clear all;
tic
nr = 60;
nz = 60;
ar = 2.5;
dr = 1/(nr+1);
dz = ar/(nz+1);
Re = 1000;
dt = 5*10^(-3);
nt = 30000; % number of iterations
% Making the discretize matrices
% matrix A represents the r part of the pde while matrix
% B represents the z part of the pde

% diag
for i = 1:nr
    A(i,i) = -2;
end
% lower diag
for i = 2:nr
    A(i,i-1) = 1+1/(2*i);
end
% upper diag
for i = 1:nr-1
    A(i,i+1) = 1-1/(2*i);
end

A = A/(dr^2);


% diag
for j = 1:nz
    B(j,j) = -2;
end

% lower diag
for j = 2:nz
    B(j,j-1) = 1;
end

% upper diag
for j = 1:nz-1
    B(j,j+1) = 1;
end
B = B/(dz^2);
[Z,E] = eig(A);
Z = -Z;
Zinv = inv(Z);
for i = 1:nr
[L(:,:,i),U(:,:,i)] = lu(B + E(i,i)*eye(nz));
end
% initialized eta, v, and psi with size nr+2 by nz+2
eta = zeros(nr+2,nz+2);
v = zeros(nr+2, nz+2);
for i = 1:nr+2
    v(i,1) =(i-1)*dr;
end
psi = zeros(nr+2,nz+2);
psi_pre = zeros(nr+2,nz+2);

nor = [];
for k = 1:nt    
    % Euler Method (The predictor)
    v_pre = v+dt*vrhs(psi,v,nr,nz,Re,dr,dz);
    eta_pre = eta + dt*etarhs(psi,v,eta,nr,nz,Re,dr,dz);
    psi_pre(2:nr+1,2:nz+1) = psisolve(eta_pre(2:nr+1,2:nz+1),L,U,Z,Zinv,nr,dr);
    
    % boundary conditions
    eta_pre = etabc(psi_pre,eta_pre,nr,nz,dr,dz);
     vtemp = v;
%     etatemp = eta;
%     psitemp = psi;
    
    % Heun's Method (The corrector)
    v = v + (dt/2)*(vrhs(psi,v,nr,nz,Re,dr,dz) +vrhs(psi_pre,v_pre,nr,nz,Re,dr,dz));
    eta = eta +(dt/2) * (etarhs(psi,vtemp,eta,nr,nz,Re,dr,dz) + etarhs(psi_pre,v_pre,eta_pre,nr,nz,Re,dr,dz));
    psi(2:nr+1,2:nz+1) = psisolve(eta(2:nr+1,2:nz+1),L,U,Z,Zinv,nr,dr);
    eta = etabc(psi,eta,nr,nz,dr,dz);
    
    if mod(k,30) == 0
        figure(1)
        contour(linspace(0,1,nr+2), linspace(0,ar,nz+2),eta',300)
        figure(2)
        contour(linspace(0,1,nr+2), linspace(0,ar,nz+2),psi',60)
        figure(3)
        contour(linspace(0,1,nr+2), linspace(0,ar,nz+2),v',60)
        drawnow;
        nor = [nor norm(v')];
    end
    
end
% 
figure(4)
plot(nor)
toc






