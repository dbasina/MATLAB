function itabc=ita_bc(ita,psi,nr,nz,dr,dz)
    for j=1:nz+2
        ita(nr+2,j)=-2*psi(nr+1,j)/(dr^2);
    end
    
    for i=1:nr+2
        r=(i)*dr;
        ita(i,1)=-2*psi(i,2)/(r*dz^2);
    end
    
    for i=1:nr+2
        r=(i)*dr;
        ita(i,nz+2)=-2*psi(i,nz+1)/(r*dz^2);
    end
    itabc=ita;
end