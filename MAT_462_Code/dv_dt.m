function dvdt=dv_dt(v,psi,Re,nr,nz,ar)
dvdt=zeros(nr+2,nz+2); %% should dvdt be the full matrix or just interior?
dr=1/(nr+1);
dz=ar/(nz+1);
for j= 2:nz+1      
      for i=2:nr+1
        r=(i-1)*dr;
        jacobian =((D1z(psi,i,j,dz)*D1r(v,i,j,dr))- (D1r(psi,i,j,dr)*D1z(v,i,j,dz)))./r;
        second= (v(i,j)/(r^2)).*(D1z(psi,i,j,dz));
        reynolds= (D2z(v,i,j,dz) + D2r(v,i,j,dr)+(D1r(v,i,j,dr)./r) - (v(i,j)./r^2)).*(1/Re);
        dvdt(i,j)=jacobian+second+reynolds;
    end
end




