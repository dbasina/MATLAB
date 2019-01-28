% IFS parameter file snowflake.m
xscale = [ 1/3 1/3 0.35 0.35 1/3 ];
yscale = xscale; %[ 0.25 0.25 0.4 0.4 0.25 0.25 ];
angles = [ 0 0 60 -60 0  ];
shear  = [ 0 0 0 0 0 ];
translates = [ [0;0] [ 0; 1 ] [0; 2]  [ 0 ; 2] [0; 2]];

parameters = [ xscale ; yscale ; shear ; angles ; translates ]';
order = 6;


