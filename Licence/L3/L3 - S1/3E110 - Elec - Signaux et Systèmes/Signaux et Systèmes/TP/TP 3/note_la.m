% Note  LA 440Hz

Fe=44100;   % Frequence d'echantillonnage 44.1kHz (standard audio)
Fee = Fe/2;
Te=1/Fe;

f0= 293.7;    % Fréquence du signal de la note LA 440
f1 = 440;
fp = 3000;      

filtrage = 400;

t=0:Te:3;   % Intervalle de temps de 0 à 3 secondes
N = length(t);
f=0:Fe/N:Fe/N*(N-1);
fc = -Fee/ N*(N-1) : Fe/N : Fee/N*(N-1);

tt=-3:Te:3;   % Intervalle de temps de 0 à 3 secondes


note=sin(2*pi*f0*t);   % Signal correspondant à la note
note2 = sin(2*pi*f1*t);

p = 2*sin(2*pi*fp*t);

%wavwrite(note/max(note),Fe,'note.wav');

audiowrite('note.wav',note/max(note),Fe);

%sound(note,Fe); 
%sound(p,Fe); 

H = abs(fft(note))/N;
G = abs(fft(p))/N;
M = (note + note2) .* p;

figure (1)
plot(f,H);
axis([0 5000 0 2])

figure (2)
plot(f,G);
axis([0 5000 0 2])

sig = note.*(p/2);
P = sig.*(p/2);

x = zeros(1, length(fc)); 
ind = find(and((fc>-filtrage),(fc<filtrage)));
x(ind) = 1;

yr = note - note.*cos(4*pi*fp*t);

figure (3)
plot(f,abs(fft(M))/N);
axis([0 5000 0 2])
%sound(M,Fe);

figure (4)
plot(f,abs(fft(P))/N);
axis([0 5000 0 2])

Z = abs(fftshift(fft(sig))); 

figure (5)
plot(fc, Z .* x, '*r');
axis([-7000 7000 0 2])

figure (6)
plot(fc, x);
axis([-7000 7000 0 2])

xlabel('Frequence');
ylabel('Amplitude');