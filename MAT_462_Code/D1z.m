function fdz=D1z(matrix,i,j,h)    
    fdz= (matrix(i,j+1)- matrix(i,j-1))/ 2*h;
end