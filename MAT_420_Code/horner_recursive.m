%% Input 1 for tracker. tracker goes from 1,2,3....,length(a).
function P= horner_recursive(a,tracker,x)
    n=tracker;
    if n == length(a)
        P=a(n);    
    else
        P=a(n)+horner_recursive(a,n+1,x).*x; 
    end
        
    