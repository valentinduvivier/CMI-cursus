
clc
clear
close all

set(0, 'DefaultLineLineWidth', 2)

L = 1;
H = 0.5;
m = 100;

k = 5e3;

% ---------------------------------------
% ---------------------------------------
% Matrices

% Matrix of inertia
M_mat = 1/6*m*[2 1 ; ...
               1 2];
      
% ---------------------------------------
           
% Stiffness matrix
K_mat = k*[1 0 ; ...
           0 1];
       
% ---------------------------------------
% ---------------------------------------
% Eigenvalues / eigenvectors

[V, D] = eig(K_mat, M_mat);

Freq = diag(sqrt(D)/2/pi); % w/2pi = f

% Normalisation
V = V./(ones(2,1)*max(abs(V)));

X = 0:2;

% Carter
z10 = H;
z20 = H;

figure(1)

for m = 1:2
    
    subplot(1,2,m)
    
    % Bati
    U0(m) = line(1.2*[-L L], [0 0], 'Linewidth', 20, 'color', 'k');
    hold on
    U1(m) = line([-L L], [z10 z10], 'Linewidth', 15); % body
    U2(m) = line([-L -L]*0.95, [0 z10], 'Linewidth', 5, 'color', 'k'); % spring 1
    U3(m) = line([+L +L]*0.95, [0 z20], 'Linewidth', 5, 'color', 'k'); % spring 2

    axis equal
    xlim([-L L]*1.2);
    ylim([-0.1 1]);
    axis off
    hold off
end

T = 1;
amp = 0.05;

pause

for t = 0:0.001:T
    for m = 1:2
        
        subplot(1,2,m)
        
        z1 = z10 + 3*amp * V(1, m)*sin(2*pi*Freq(m)*t);
        z2 = z20 + 3*amp * V(2, m)*sin(2*pi*Freq(m)*t);
        TH = (z2 - z1)/(2*L);
        
        set(U1(m), 'Xdata', [-L L], 'Ydata', [z1 z2]);
        set(U2(m), 'Ydata', [0 z1]);
        set(U3(m), 'Ydata', [0 z2]);
       
    end
    drawnow
end
  
hold off

%%

F = [1; 0];
c = 1;
C = 1;
C_mat = c*K_mat/k + C*M_mat/m; % Rayleigh damping

% Modal matrices
Km = V'*K_mat*V;
Mm = V'*M_mat*V;
Cm = V'*C_mat*V;
Fm = V'*F;

xi1 = 0.05;
xi2 = 0.05;
xi  = diag([xi1, xi2]); % Modal damping

freq = 0:0.0001:max(Freq)*1.1;
HRayleigh   = zeros(2, size(freq, 2));
HAmodal     = zeros(2, size(freq, 2));

ii = 0;

for omega = freq*2*pi
    ii = ii + 1;
    HRayleigh(:,ii) = (Km - omega^2*Mm + 1i*omega*Cm)\(Fm);
    HAmodal(:,ii)   = (Km - omega^2*Mm + 1i*omega*xi)\(Fm);
end
    
%% 

figure(2)
subplot(2,1,1)

H1 = abs(V(1,:)*HRayleigh); H1 = H1./max(H1);
H2 = abs(V(2,:)*HRayleigh); H2 = H2./max(H2);

plot(freq, 20*log(H1), freq, 20*log(H2))

xlim([0 max(Freq)*1.1])
ylim([-150 10])

title('Transfer functions - Rayleigh damping or proportionnal damping')
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')

grid('on'); grid minor
legend('z__1/F', 'z__2/F')
  
% ---------------------------------------

subplot(2,1,2)

plot(freq, 180/pi*angle(V(1,:)*HRayleigh), ...
        freq, 180/pi*angle(V(2,:)*HRayleigh))
    
xlim([0 max(Freq)*1.1])

title('Transfer functions - Modal damping')
xlabel('Frequency [Hz]')
ylabel('Phase [deg°]')

grid('on'); grid minor
    
% ---------------------------------------
% ---------------------------------------

%%
figure(3)
subplot(2,1,1)

H1 = abs(V(1,:)*HAmodal); H1 = H1./max(H1);
H2 = abs(V(2,:)*HAmodal); H2 = H2./max(H2);

plot(freq, 20*log(H1), freq, 20*log(H2))

xlim([0 max(Freq)*1.1])
ylim([-150 10])

title('Transfer functions - Modal damping')
xlabel('Frequency [Hz]')
ylabel('|H| [dB]')

grid('on'); grid minor
legend('z__1/F', 'z__2/F')

% ---------------------------------------

subplot(2,1,2)

plot(freq, 180/pi*angle(V(1,:)*HAmodal), ...
        freq, 180/pi*angle(V(2,:)*HAmodal))
    
xlim([0 max(Freq)*1.1])

title('Transfer functions - Modal damping')
xlabel('Frequency [Hz]')
ylabel('Phase |H| [deg°]')

grid('on'); grid minor

% ---------------------------------------
% ---------------------------------------

%%

h   = 1;
g   = 1;

t   = 0:0.01:10;
v   = 1;
t0  = 2*L/v;

HSA         = 1*t./t;
HSB         = HSA;
HSB(t<t0)   = 0;

m1 = Mm(1,1);
m2 = Mm(2,2);

o1  = 2*pi*Freq(1);
o2  = 2*pi*Freq(2);
o1d = o1*sqrt(1-xi1^2);
o2d = o2*sqrt(1-xi2^2);

p1 = g*h/m1*HSA.*exp(-xi1*o1*t).*sin(o1d*t)+ ...
        g*h/m1*HSB.*exp(-xi1*o1*(t-t0)).*sin(o1d*(t-t0));
    
p2 = g*h/m2*HSA.*exp(-xi2*o2*t).*sin(o2d*t)+ ...
        g*h/m2*HSB.*exp(-xi2*o2*(t-t0)).*sin(o2d*(t-t0));
    
zA = p1 + p2;
zB = p1 - p2;

% ---------------------------------------
% ---------------------------------------


figure(4)
subplot(2,1,1)
plot(t, p1, t, p2)

xlabel('temps [s]')
ylabel('Mode amplitude')

xlim([0 max(t)])
grid('on'); grid minor
legend('p1', 'p2')

% ---------------------------------------

subplot(2,1,2)
plot(t, zA, t, zB)

xlabel('temps [s]')
ylabel('Amplitude [m]')

xlim([0 max(t)])
grid('on'); grid minor
legend('zA', 'zB')


