function sdz=D2z(matrix,i,j,h)
    sdz= (matrix(i,j+1)-2*matrix(i,j)+matrix(i,j-1))/h^2;
end