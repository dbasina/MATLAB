function psisolve = psi_sol(ita,L,U,Z,Zinv,dr,nr)
    for i=1:nr
       r=i*dr;
       F(i,:)=-r*ita(i,:);
    end
    H=transpose(F)*transpose(Zinv);
    for i=1:nr
        X(i,:)=U(:,:,i)\L(:,:,i)\H(:,i);
    end    
   psisolve=Z*X;       
end

    