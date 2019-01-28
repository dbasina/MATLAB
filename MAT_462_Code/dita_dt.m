function ditadt=dita_dt(v,ita,psi,Re,nr,nz,ar)
ditadt=zeros(nr+2,nz+2);
dr=1/(nr+1);
dz=ar/(nz+1);
for i=2:nr+1
    r=(i-1)*dr;
    for j= 2:nz+1
    jacobian= (D1z(psi,i,j,dz)*((D1r(ita,i,j,dr)*(1/r)) -(ita(i,j)*(1/r^2)))) - (D1r(psi,i,j,dr)*((1/r).*D1z(ita,i,j,dz)));
    second= (2*(v(i,j)/r))*D1z(v,i,j,dz);
    reynolds= (D2z(ita,i,j,dz)+D2r(ita,i,j,dr)+(D1r(ita,i,j,dr)/r) - (ita(i,j)/r^2)).*(1/Re);
    ditadt(i,j)= jacobian+second+reynolds;
    end
end
    
