% TO : 
% to = 1 ? tr = 6,994s
% to = 10 ? tr = 30.970s
% to = 2 ? tr = 4.090s
% 
% 
% pulsation de coupure :
% w = 2*pi*f
% to = 1 ? 1*2pi = 2pi
% to = 10 ? 0,1*2pi = 0,628
% to = 2 ? 0,5*2pi = pi

to = {6.994, 4.090, 30.970};
w = {2*pi, pi, 0.628};

x = linspace(0,20,50);

subplot(1,2,1)
plot(4.090, 0.725, 'r*')
hold on
plot(6.994, 0.43, 'r*')
plot(30.97, 0.10, 'r*')
hold off
xlabel('tr 5%'); ylabel('pulsation de coupure - w');
title('Courbe w = f(tr)');
axis([0 32 0 1])

subplot(1,2,2)
plot(1.363, 0.725, 'r*')
hold on
plot(2.331, 0.43, 'r*')
plot(10.32, 0.10, 'r*')
plot(x, 1./x, 'b-');
hold off
xlabel('tau'); ylabel('pulsation de coupure - w');
title('Concordance y = f(1/x) et w = f(tau)');
axis([0 32 0 1])