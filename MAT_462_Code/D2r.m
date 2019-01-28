function sdr=D2r(matrix,i,j,h)
    sdr= (matrix(i+1,j)-2*matrix(i,j)+matrix(i-1,j))/ h^2;
end