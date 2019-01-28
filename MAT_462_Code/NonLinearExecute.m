clear;
clc;
%%Initialize the Matrices V,ita,psi
nr=60;
nz=60;
Re=1000;
ar=1;
dr=1/(nr+1);
dz=ar/(nz+1);
dt=5*(10^-3);
Tn=1000;
v=zeros(nr+2,nz+2,Tn);
ita=zeros(nr+2,nz+2,Tn);
psi= zeros(nr+2,nz+2,Tn);
% f1=figure('Name','Velocity');
% f2=figure('Name','Vorticity');
% f3=figure('Name','Stream Lines');
% f4=figure('Name','Velocity Surface');
% f5=figure('Name','Vorticity Surface');
 f6=figure('Name','Stream Surface');

%% Add boundary conditions to v.
for k=1:Tn
    for i=1:nr+2
        
        v(i,1,k)= (i-1)*dr;
    end
end
%% Set up the PSI derivative function using tri diagonals on the interior points
A=zeros(nr);
%A Middle Diagonal
for i=1:(nr)
    A(i,i)=-2;
end
%A Upper Diagonal
for i=1:(nr-1)
    A(i,i+1)=(1-(1/(2*i)));           
end
%A lower Diagonal
for i=2:(nr)
    A(i,i-1)=(1+(1/(2*i)));
end
A=A/(dr^2);

B=zeros(nz);
%B middle Diagonal 
for j=1:(nz)
    B(j,j)=-2;
end
%B Upper Diagonal
for j=1:(nz-1)
    B(j,j+1)=1;
end
%B Lower Diagonal
for j=2:(nz)
    B(j,j-1)=1;
end
B=B/(dz^2);

%% Heunn's Method
tic
for k=1:Tn-1
   
    %Euler Method
    v_temp= v(:,:,k)+(dt*dv_dt(v(:,:,k),psi(:,:,k),Re,nr,nz,ar));
    ita_temp = ita(:,:,k)+(dt*dita_dt(v(:,:,k),ita(:,:,k),psi(:,:,k),Re,nr,nz,ar));
    ita_interior= ita_temp(2:nr+1,2:nz+1);
    for i=1:nr
        r=dr*(i+1);
        F_temp(i,:)=-r.*ita_interior(i,:);
    end
    psi_temp_interior= sylvester(A,B,F_temp);
    psi_temp=padarray(psi_temp_interior,[1 1]);
    ita_temp=ita_bc(ita_temp,psi_temp,nr,nz,dr,dz);
    
    %%Correct    
    v(:,:,k+1)=v(:,:,k)+(dt/2)*(dv_dt(v(:,:,k),psi(:,:,k),Re,nr,nz,ar)+ dv_dt(v_temp,psi_temp,Re,nr,nz,ar));
    ita(:,:,k+1)=ita(:,:,k)+ (dt/2)*(dita_dt(v(:,:,k),ita(:,:,k),psi(:,:,k),Re,nr,nz,ar)+dita_dt(v_temp,ita_temp,psi_temp,Re,nr,nz,ar));
    ita_interior_C=ita(2:nr+1,2:nz+1,k+1);
    for i=1:nr
        r=dr*(i+1);
        F(i,:)=-r.*ita_interior_C(i,:);
    end
    psi_interior=sylvester(A,B,F);
    psi(:,:,k+1)=padarray(psi_interior,[1 1]);
    ita(:,:,k+1)=ita_bc(ita(:,:,k+1),psi(:,:,k+1),nr,nz,dr,dz);
    
     drawnow
%      figure(f1);
%      contour(transpose(v(:,:,k+1)));
%      figure(f2);
%      contour(transpose(ita(:,:,k)),60);
%      figure(f3);
%      contour(transpose(psi(:,:,k)));
%      figure (f4)
%      surfc(transpose(v(:,:,k)));
%      figure (f5)
%      surfc(transpose(ita(:,:,k)));
      figure (f6)
      surfc(transpose(psi(:,:,k)));
      
end
toc

