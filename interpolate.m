%interpolate finds the polynomial, piecewise linear, spline, and pchip fits

% pchip() and spline() are available as built-in functions
% piecelin() and polyinterp() are in Moler's NCM and can be downloaded
% polyinterp.m is on the course web page
% Here we use interp1()

y = [80; 125; 170; 215; 350];
x = [1; 2; 3; 4; 5];
p = polyfit(x,y,5)
xpts = 1:0.01:5;
ypoly = polyval(p,xpts);
ylin = interp1(x,y,xpts,'linear');
yspline = interp1(x,y,xpts,'spline');
ypchip = interp1(x,y,xpts,'pchip');


plot(x,y,'o',xpts,ypoly,'r-',xpts,ylin,'g-',xpts,yspline,'c-',xpts,ypchip,'b-')
xlabel('x','FontSize',16)
ylabel('y','FontSize',16)
legend('data','poly','piecelin','spline','pchip','Location','South')

ypchip0 = interp1(x,y,-0.3,'pchip') % to find pchip value at x = -0.3