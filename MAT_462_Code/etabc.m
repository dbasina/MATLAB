function eta_bc= etabc(psi,eta,nr,nz,dr,dz)
%    for i = 1:nr+2
%        r = (i-1)*dr;
%    end
   for j = 1:nz+2
       eta(nr+2,j) =  -2*psi(nr+1,j)/(dr^2); % at wall r = 1
   end
   for i = 1:nr+2
        r = (i-1)*dr;
        eta(i,1) = -2*psi(i,2)/(r*dz^2); % at bottom z = 0
   end 
   for i = 1:nr+2
        r = (i-1)*dr;
        eta(i,nz+2) = -2*psi(i,nz+1)/(r*dz^2); % at the top z = 1
   end
    
   
    eta_bc = eta;

end
