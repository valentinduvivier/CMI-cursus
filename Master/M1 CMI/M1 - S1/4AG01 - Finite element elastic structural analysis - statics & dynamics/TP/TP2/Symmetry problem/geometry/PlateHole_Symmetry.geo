// global dimensions
R   = 0.025; 
L2  = 0.4;
L1  = L2; 

// mesh density
lc1     = 10^-2;
lc2     = 10^-3;

// points coordinates
Point(1) = {   0,    0, 0, lc1};
Point(2) = {L1/2,    0, 0, lc1};
Point(3) = {L1/2, L2/2, 0, lc1};
Point(4) = {   0, L2/2, 0, lc1};

// Points coordinates arround the circle
Point(5) = {    L1/4, 0, 0, lc2};
Point(6) = {L1/4 - R, 0, 0, lc2};
Point(7) = {L1/4 + R, 0, 0, lc2};

// lines square
Line(1) = {7, 2}; // top left
Line(2) = {2, 3}; // top
Line(3) = {3, 4}; // left
Line(4) = {4, 1}; // bottom left
Line(5) = {1, 6}; // bottom right

// Circle
Circle(6) = {7, 5, 6};
//Ellipse(6) = {6, 2, 5, 5};

// Connection lines and circle
Line Loop (7) = {1, 2, 3, 4, 5, -6};

// Surfaces 
Plane Surface(1) = {7};

// physical entities
Physical Line("right")   = {2}; 
Physical Line("left")    = {4}; 
Physical Line("top")     = {3}; 
Physical Line("bottom")  = {1}; 
Physical Line("hole")    = {5};

Physical Surface("plate_ellipse") = {1}; 

