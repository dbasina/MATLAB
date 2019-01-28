% IFS parameter file fern.m
tilt = 5;
stalk = 0.05;
left = 0.3;
right = 0.3;
body =0.85;
height = 1;


xscale = [ 0 -left right 0.95*body ];
yscale = [ 3.8*stalk left right body ];
angles = [ 0 100 -85 tilt ];
shear  = [ 0 0 0 0];
translates = [ [0;0] [ 0; 1 ] [0; 1.3 ]  [ 0 ; 1.8 ]  ];

parameters = [ xscale ;yscale ; shear ;angles ;translates ]';

