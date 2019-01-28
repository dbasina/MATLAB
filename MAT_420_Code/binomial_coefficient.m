function C = binomial_coefficient(n,k)
    if k==0
        C=1;
    elseif k==n
        C= 1;
    else
        C= binomial_coefficient(n-1,k)+binomial_coefficient(n-1,k-1);
end

        