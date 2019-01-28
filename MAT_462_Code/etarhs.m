function deta_dt = etarhs(psi,v,eta,nr,nz,Re,dr,dz)
deta_dt = zeros(nr+2,nz+2); % initialize deta_dt no change at t = 0 full matrix

% for i = 2:nr+1
%     r = (i-1)*dr;
% end

for i = 2:nr+1
        r = (i-1)*dr;
    for j = 2:nz+1
        deta_dt(i,j) = (1/r)*((psi(i,j+1)-psi(i,j-1))/(2*dz))*((eta(i+1,j)-eta(i-1,j))/(2*dr)) -...
                       (eta(i,j)/r^2)*((psi(i,j+1)-psi(i,j-1))/(2*dz)) - (1/r)*((psi(i+1,j)-psi(i-1,j))/(2*dr))*((eta(i,j+1)-eta(i,j-1))/(2*dz)) +...
                       (2*v(i,j)/r)*((v(i,j+1)-v(i,j-1))/(2*dz)) +...
                       (1/Re)*( (eta(i,j+1)-2*eta(i,j)+eta(i,j-1))/(dz^2) + (eta(i+1,j)-2*eta(i,j) + eta(i-1,j))/(dr^2) +...
                       (1/r)*((eta(i+1,j)-eta(i-1,j))/(2*dr)) - (eta(i,j)/r^2) );
    end
end
            
end