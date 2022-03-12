// global dimensions
h = 15.;

// mesh density
lc1 = 10^-2;
lc2 = 1;

// points coordinates
Point(1) = {0, 0, 0, lc1};
Point(2) = {h, 0, 0, lc2};
Point(3) = {0, h, 0, lc2};

// lines and line loops
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line Loop(4) = {1, 2, 3};

// surface
Plane Surface(1) = {4};

// physical entities
Physical Line("bottom") = {1}; 
Physical Line("left")   = {3}; 
Physical Surface("dam") = {1}; 
