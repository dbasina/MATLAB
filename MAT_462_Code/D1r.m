function fdr=D1r(matrix,i,j,h)   
    fdr= (matrix(i+1,j)-matrix(i-1,j))/2*h;
end