t = linspace(0,2*pi,100);

y_sin = sin(t);

y_LU_n10 = 0.283845 - 0.018279*t - 0.946138*(t.^2) + 0.920037 * (t.^3) - 0.335589*(t.^4) + 0.053360*(t.^5) - 0.003107*(t.^6);
y_Chol_n10 = 0.436856 - 1.137458*t + 0.85094*(t.^2) - 0.2155*(t.^3) - 0.000773*(t.^4) + 0.007078*(t.^5) - 0.000684*(t.^6);

y_LU_n1 = 0.108879 + 0.070487*t + 1.650026*(t.^2) - 1.303064*(t.^3) + 0.362983*(t.^4) - 0.044481*(t.^5) + 0.002072*(t.^6);
y_Chol_n1  = 0.108879 + 0.070487*t + 1.650026*(t.^2) - 1.303064*(t.^3) + 0.362983*(t.^4) - 0.044481*(t.^5) + 0.002072*(t.^6);

%subplot(1,2,1)
%plot(t, y_sin, 'r-', t, y_Chol_n10, 'b-*', t, y_LU_n10, 'g-*');
%title("Courbe Sinus(10*t) (r) et approximation Cholesky (b) et LU (g) n=10");

%subplot(1,2,2)
plot(t, y_sin, 'r-', t, y_LU_n1, 'b-', t, y_Chol_n1, 'g-');
title("Courbe Sinus(t) (r) et approximation Cholesky (b) == LU (g) n=1");
