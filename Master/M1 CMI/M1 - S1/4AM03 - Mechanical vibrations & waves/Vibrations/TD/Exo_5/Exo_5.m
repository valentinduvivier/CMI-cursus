% Modelisation vibration ship

clc,
clear all
close all

srcsz = get(0, 'ScreenSize');
set(0, 'defaultlinelinewidth', 1.5)

k = 10e3 * (75*60); % rho * S

M = 10e6;

kf = k*...
    [1 0 0 0;
     0 1 0 0;
     0 0 1 0;
     0 0 0 1];
 
 L = 75;
 K = 10e11; % Nm/rad
 
 ks = K/L^2 * ...
    [ 2 -2  0  0;
     -2  4 -2  0;
      0 -2  4 -2;
      0  0 -2  2];
  
  m = M*eye(4);
  
  k = ks + kf;
  
  [v, d] = eig(k,m);
  
  fr = diag(sqrt(d)/2/pi);

  V = v./repmat(max(abs(v)), [4,1]);
  
  
  % ------------------------------------------- %
  % Calculation of frequential response with Rayleigh's damping
  
  LSwell = 45; % m
  CSwell = 1;  % m.s-1
  FSwell = CSwell/LSwell;
  
  x = L*(0:3)';
  f = exp(1i*2*pi*FSwell*x);
  c1 = 0.1;
  c2 = 0.1;
  c = c1*k/k + c2*m/m; % Rayleigh's damping
  
  km = v'*k*v;
  mm = v'*m*v;
  cm = v'*c*v;
  fm = v'*f;