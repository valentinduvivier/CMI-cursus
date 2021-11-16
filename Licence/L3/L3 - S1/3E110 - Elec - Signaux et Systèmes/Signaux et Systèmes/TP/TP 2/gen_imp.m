A=2;
a=10;
t1=0.1;
t2=0.4;

Fe = 512;
N = 256;
Te = 1/Fe;
t = 0:Te:(N-1)*Te;
tt = 0:Te:(2*N-2)*Te;

x = zeros(1, length(t)); 
ind = find(and((t>t1),(t<t2)));
x(ind) = A;

f = exp(-a*t);

figure(1);
% subplot(2,2,1)
% plot(t,x)
% title(['Fonction f(t) avec A=',num2str(A), ', t1=',num2str(t1), ', t2=', num2str(t2)]);
% subplot(2,2,2)
% plot(t,f,'r-')
% title(['Fonction g(t) pour a=', num2str(a)]);

freq = 0:(Fe/N):(N-1)*Fe/N;

X = abs(fft(x))/N;
subplot(2,1,1)
plot(freq, X);
title('TF pour la fonction f(t)');
F = abs(fft(f))/N;

subplot(2,1,2)
plot(freq, F, 'r-');
title('TF pour la fonction g(t)');

prod_conv = F.*X;
%prod_tf = 
figure(2);
%plot(tt,abs(fft(conv(x,f))));
%title(['TF de la convolution de f(t) et g(t)']);
plot(freq, prod_conv)
title(['Produit des TF de f(t) et g(t)']);
