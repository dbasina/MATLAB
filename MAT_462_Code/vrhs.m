function dv_dt = vrhs(psi,v,nr,nz,Re,dr,dz)
% This function file computes the velocity change of the fluid in the
% cylindrical chamber using the finite difference method.


dv_dt = zeros(nr+2,nz+2); % initialize dv_dt no change at t = 0 full matrix
% for i = 1:nr+2
%     r = (i-1)*dr;
% end
for i = 2: nr+1
    r = (i-1)*dr;
    for j = 2:nz+1
        % This is the pde of dv/dt
        dv_dt(i,j) = (1/r)*((psi(i,j+1)-psi(i,j-1))*(v(i+1,j)-v(i-1,j))/(4*dz*dr)) - (1/r)*((psi(i+1,j)-psi(i-1,j))*(v(i,j+1)-v(i,j-1))/(4*dz*dr)) +...
            (v(i,j)/r^2)*((psi(i,j+1)-psi(i,j-1))/(2*dz)) + (1/Re)*((v(i,j+1)-2*v(i,j)+v(i,j-1))/(dz^2) + (v(i+1,j)-2*v(i,j) + v(i-1,j))/(dr^2) +...
            (1/r)*(v(i+1,j)-v(i-1,j))/(2*dr) - (v(i,j)/r^2));
    end
end

end




