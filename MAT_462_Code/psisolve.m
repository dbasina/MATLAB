function psi_in = psisolve(eta,L,U,Z,Zinv,nr,dr)

     for i = 1:nr
         F(i,:) = -(i*dr)*eta(i,:);
     end
        H = transpose(F)*transpose(Zinv);
     for i = 1:nr
      X(i,:) = U(:,:,i)\(L(:,:,i)\(H(:,i)));
      
     end
     psi_in = Z*X;
end
 