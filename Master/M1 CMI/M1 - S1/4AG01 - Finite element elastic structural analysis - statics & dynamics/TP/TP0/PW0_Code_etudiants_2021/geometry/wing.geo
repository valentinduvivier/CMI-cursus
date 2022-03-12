// global dimensions
L1 = 48.0;
L2 = 44.0;
L3 = 60.0;

// mesh density
lc1 = 5.;

// points coordinates
Point(1) = {0.0,0.0,0.0,lc1};
Point(2) = {L1,L2,0.0,lc1};
Point(3) = {L1,L3,0.0,lc1};
Point(4) = {0.0,L2,0.0,lc1};

// lines and line loops
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(5) = {1,2,3,4};

// surface
Plane Surface(1) = {5};

// physical entities
Physical Line("right") = {2}; 
Physical Line("left") = {4}; 
Physical Surface("wing") = {1}; 
