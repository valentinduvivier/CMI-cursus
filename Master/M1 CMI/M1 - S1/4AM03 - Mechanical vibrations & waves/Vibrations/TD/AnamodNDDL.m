

function [X, f] = AnamodNDDL(M, K);

disp('Matrix of inertia')
M
disp('Matrix of stiffness')
K

disp('Heart : A = M-1*K')
A = (K/M)'

[X, w2] = eig(A)

% -------------------------------------------------

disp('Eigenvalues : ')
vp = abs(diag(w2)')

disp('Eigenvectors : ')
f = sqrt(vp)/(2*pi)

disp('Eigenmodes : ')
X

% -------------------------------------------------

disp('orthogonality of the modes : ');
disp('Modal masses');
XtMx = X'*M*X
disp('Modal stiffness')
XtKx = X'*K*X
