Nr=60; %number of grid points in r
Nz=80; % number of grid points in z
dr=1/Nr-1; %discretization of r and z
dz=1/Nz-1;

%AU= d^2v/dr^2 + 1/r dv/dr - v/r^2
%UB= d^2v/dz^2 +F
% AU+UB=F

A=zeros(Nr-1);  
B=zeros(Nz-1);

%Construct A-Matrix using equations obtained for tridiagonal form
%A middle Diagonal
for i=1:(Nr-1)
    A(i,i)=-2-(1/i^2);
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
    F(k,1) = k*dr;
end
F=-F/(dz^2);
 %time discretization
 t0=0;
 dt=0.01;
 tf=500;
 Tspan=[t0:dt:tf];
 nt=length(Tspan);
 U=zeros(Nr-1,Nz-1,nt);
 tic
for i= [1:nt]
 Uint= U(:,:,i)+(dt.*(A*U(:,:,i)+U(:,:,i)*B+F));
 U(:,:,i+1)= U(:,:,i)+((dt*0.5).*((A*U(:,:,i)+U(:,:,i)*B+F)+(A*Uint+Uint*B+F)));
end
h = uicontrol(gcf,'style','slider','units','pix','position',[100 5 300 20]);
set(h,'min',min(Tspan),'max',max(Tspan) );
set(h,'callback','i=find(t==nearest(get(h,"value"))')
 for i = 1:nt
     figure(1)
     contour(U(:,:,i)')
     pause(0);
     set(h,'value',Tspan(i));
 end
