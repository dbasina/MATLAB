function v = polyinterp(x,y,u) 
%polyinterp(x,y,u) calculates the full interpolating polynomial by Lagrange
%   interpolation
% Example:
% x = 1:6; 
% y = [16 18 21 17 15 12]; 
% disp([x; y]) 
% u = 0.75:0.05:6.25; % dense set of x points
% v = polyinterp(x,y,u); 
% plot(x,y,'o',u,v,'r-')

% code processes all components of u simultaneously using MATLAB array
% operations

% v = P(u) = sum_k product_{j != k} (u-x(j))/(x(k)-x(j))*y(k)
% != is "not equal to"

n = length(x); % polynomial is of degree n-1
v = zeros(size(u)); % allocate array of dense set of y points
for k = 1:n
    w = ones(size(u)); % allocate array for terms in product
    for j = [1:k-1 k+1:n]
        w = (u-x(j))./(x(k)-x(j)).*w; 
    end
    v = v + w*y(k); % terms in sum
end

end