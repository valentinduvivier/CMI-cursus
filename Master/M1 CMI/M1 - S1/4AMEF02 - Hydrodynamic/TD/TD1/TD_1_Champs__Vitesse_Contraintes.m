% Code to display the speed field as well as the contraint field in the
% case of a sheared flow.

clear all; clc; close all;

% Constant variables

    % Tube caracteristics
        h = 1;  % Height of the pipe
        
    % Flow caractesristics
        nu_w    = 10^-6;    % [Pa.s]    % at T° = 20°C
        rho_w   = 1000;     % [kg.m-3]    
        

    % Problem caracteristics
        % Time of study
            T   = 100;           % Time of the study

        % Precision Fourier series
            N   = 500;           % Nb n considered for the precision with Fourier sum
        
        % Vertical position
            N_y     = 100;      % Nb y position for y => [0,h]
            Dy      = h/N_y;    % Space step
        
            y_dim   = Dy:Dy:h;
            
        % Transitory time
            tau     = (h^2/nu_w)/(3*10^3);             % Transitory time
    
        % Upper side speed
            U = 10;     % [m.s-1]
        
% Flow_field & Constraint_field

    Sum_u = 0;    % Term for Fourier series
    Sum_sig = 0;

    u = zeros(N_y, T);    % Initial flow field
    sig = zeros(N_y, T);  % Initial constraint field
    
    for t = 1:T
       
        for y = 1:N_y
            
            for n = 1:N     % For each t, we estimate that we need N terms on the sum to get a precise enough result
                Sum_u = Sum_u + 2/(n*pi)*((-1)^n)*sin(n*pi*(Dy*y)/h)*exp(-((n*pi)^2)*t/tau);
                Sum_sig = Sum_sig + ((-1)^n)*cos(n*pi*(Dy*y)/h)*exp(-((n*pi)^2)*t/tau);
            end
            
            u(y,t)      = U*((Dy*y)/h + Sum_u);
            sig(y,t)    = (U*(nu_w*rho_w)/h)*(1 + 2*Sum_sig);
            
            Sum_u = 0;
            Sum_sig = 0;
        end
    end
    
    
%% Displaying results
    
    figure(1)
        for t = 1:T
            txt = ['',num2str(t/tau)];
            plot(u(:,t), y_dim', 'color', t*[1/(1*T), 1/(2*T), 1/(10*T)], 'DisplayName', txt);
            xlabel('Speed'); ylabel('Height [0,h]'); title('Speed field');
            hold on
            grid on
        end
        hold off
        legend show

        
    figure(2)
        for t = 1:T
            txt = ['',num2str(t/tau)];
            plot(sig(:,t), y_dim', 'color', t*[1/(1*T), 1/(2*T), 1/(10*T)], 'DisplayName', txt);
            xlabel('Constraint'); ylabel('Height [0,h]'); title('Constraint field');
            hold on
            grid on
        end        
        hold off
        legend show
        
        