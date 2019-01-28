function ifs(mfile,max_iterations)

%
% ifs     Iterated function system fractal drawer
%
%    Calling sequence    ifs('mfilename',max_iterations)
%                  or    ifs mfilename max_iterations
%
%    where mfilename.m is an m-file containing the parameters of an
%    affine iterated function system, and max_iterations is an integer.
%    The second argument is optional, defaults to 10 000, and is
%    capped at 1 million. The default value is chosen for speed of 
%    execution, rather than resolution; 50 000 appears to be a good
%    compromise between the two.
%
%    The input m-file must define an (n x 6)-matrix called "parameters",
%    with one row per function of the IFS. The format for each row is
%
%       xscale    yscale    shear    theta    xshift    yshift
%
%    where
%         xscale, yscale     are scale factors in the x and y directions,
%                              and must lie between -1 and 1;
%         shear              is the amount of shear parallel to the x-axis
%                              (i.e. multiplication by [ 1 shear ; 0 1 ])
%         theta              is an angle of rotation, anti-clockwise and
%                              measured in degrees;
%         [xshift ; yshift ] is the translation vector of the affine map.
%
%    Each function of the IFS is constructed by applying the transformations 
%    in the order given, i.e. scale, then shear, then rotate, then translate.
%    Each function is chosen with probability approximately proportional
%    to its determinant, i.e. to |(xscale)(yscale)|.
%
%    A sample parameter file is
%
%               parameters=[
%               0.95	0.95	0	 10	1	-1/2
%               1/4	-1/8	1	-40	0	4.8
%               ];
%
%    ifs called with no arguments defaults to 10 000 iterations of this
%    system. 
%
%    The input m-file may optionally define an additional positive integer 
%    parameter "order", specifying an order of rotational symmetry. n copies
%    of the fractal will be plotted, rotated through angles 360k/n degrees
%    about the origin, k=0:n-1.

order=1;
if (nargin==0)
parameters =[
0.95	0.95	0	 10	1	-1/2
1/4	-1/8	1	-40	0	4.8
];
else
eval(mfile)
end

if (nargin <2)
     max_iterations = 10000;
else
     if ischar(max_iterations)
            max_iterations=str2num(max_iterations);
     end
     max_iterations = min(max_iterations,1000000);
end

[n,m]  = size(parameters);
if m<6
error('Parameter matrix has fewer than 6 columns')
elseif m>6
warning('Parameter matrix has more than 6 columns--additional columns ignored')
end

xscale = parameters(:,1);
yscale = parameters(:,2);
shear  = parameters(:,3);
theta  = parameters(:,4);
shift  = [ parameters(:,5)'; parameters(:,6)' ];

if max(abs([xscale ; yscale])) >= 1
       error('All scale factors must lie between -1 and 1')
end

Map = zeros(2,2,n);

for i=1:n

rotation = [ cosd(theta(i)) -sind(theta(i)) ; sind(theta(i)) cosd(theta(i)) ];
Shear    = [ 1 shear(i) ; 0 1 ];
scale    = [ xscale(i) 0 ; 0 yscale(i) ];
Map(:,:,i) = rotation*Shear*scale;

end


d1=sprintf('\n%d iterations of an iterated function system with %d functions',max_iterations,n);
disp(d1)

detvec = abs(xscale.*yscale);
tol = 0.25*min(detvec(detvec>0));
detvec = detvec+tol;
weights = detvec/sum(detvec)
cweights = cumsum(weights);


x = zeros(2,max_iterations);
y = [0;0];

for i = 1:100
map = sample(cweights);
     y = Map(:,:,map)*y + shift(:,map);
end

x(:,1) = y;

for i = 2:max_iterations
map = sample(cweights);
     x(:,i) = Map(:,:,map)*x(:,i-1) + shift(:,map);
end

plot(x(1,:),x(2,:),'b.','MarkerSize',0.1)
axis equal
%axis off

%rotational symmetry
order = floor(order);
if order > 1
    angle = 360/order;
    rot = [ cosd(angle) -sind(angle) ; sind(angle) cosd(angle) ];
    hold on
    for i = 1:order
        x = rot*x;
        plot(x(1,:),x(2,:),'b.','MarkerSize',0.1)
    end
    hold off
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sampler
function choice = sample(cumweight)
     r = rand;
     choice = min(find(cumweight>r));
