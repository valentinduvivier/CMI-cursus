% Exo 4
m = 1;
k = 1;

m1 = m;
m2 = m;
m3 = m;

k1 = k;
k2 = k;

K = [k1 -k1   0;
    -k1 k1-k2 -k2;
    0   -k2   k2];

M = [m1 0 0; 0 m2 0; 0 0 m3];

[V,f] = AnamodNDDL(M, K);

%AnimNDDL(V,f)