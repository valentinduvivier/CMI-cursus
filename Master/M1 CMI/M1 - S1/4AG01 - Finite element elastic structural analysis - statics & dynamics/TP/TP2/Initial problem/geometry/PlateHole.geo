// global dimensions

R    = 0.025; 
L    = 0.4; 

// mesh density

l     = 1.;

// points coordinates

Point(1) = {-L/2, -L/2, 0, l};
Point(2) = { L/2, -L/2, 0, l};
Point(3) = { L/2,  L/2, 0, l};
Point(4) = {-L/2,  L/2, 0, l};

// lines and line loops

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};


// Points coordinates arround the circle

Point(10) = { 0, 0, 0, l};
Point(11) = { R, 0, 0, l};
Point(12) = { 0, R, 0, l};
Point(13) = {-R, 0, 0, l};
Point(14) = { 0,-R, 0, l};

// Circle

Circle(15) = { 11, 10, 12};
Circle(16) = { 12, 10, 13};
Circle(17) = { 13, 10, 14};
Circle(18) = { 14, 10, 11};

Line Loop (6) = {15, 16, 17, 18};

// physical entities

// lines

Physical Line("right")   = {2}; 
Physical Line("left")    = {4}; 
Physical Line("top")     = {3}; 
Physical Line("bottom")  = {1}; 

// Surfaces 

Plane Surface(1) = {5, 6};
Physical Surface("plate") = {1}; 
