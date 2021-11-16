% Homework Assignement 2 - Hyperbolic Equations
    % Throughout this code we will put in action the theory seen during the
    % course and try to have a better understanding on how to deal with
    % resolution of non-linear equation.
    
% Objectives : 
%   - Have more indsight on the effect of inhomogeneous term in a numerical
%       scheme
%   - Deal with different initial boundary for hyperbolic systems

clear all; clc; close all;

%% 2 Problem's variables

    % We consider the following fixed variables : 
        L = 10;     % [m]       % Length of the mesh
        H = 1;      % [m]       % Initial wave height, that stays constant throughout the entire motion
        g = 9.61;   % [m/s^2]   % Gravity's acceleration
        w = 0.4;    % [m]       % "width" of the hill
        
%% 2.1 - Numerical Solution
    
    % During this part we consider :
        epsilon = 0.1;
      
    % Method used
        % We use Lax-Friedrichs method : u = (h,hv)^T :
        %   uj^n+1 - uj^n + (Delta_t/Delta_x)*(F(u_j+1,u_j)^(F*L_x) -
        %   F(u_j,u_j-1)^(F*L_x)) = 0

        % With : 
        %   - F(u_j+1,u_j)^(F*L_x) := 1/2[f(u_j+1) + f(u_j) - alpha(u_j+1 -
        %   u_j)]
        %   - alpha = max|f'(u)| on u, for j = 0, +/-1, +/-2, ... .
    
% ------------------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------------------ %

% Q 2.1 - 1 : resolution system

    % We here consider : 
        %   - alpha = Delta_x/Delta_t
    
            N = 100;        % Number of cells for the square mesh - has impact on the length
                            % to cross and the precision for the curve
                                    
            N_t = 800;      % Number of iteration - define the size of the loop that gives U

            T = 5;          % Time of the analysis (sec)
            
            Dx = L/N;       % Space-step
            Dt = 0.3*Dx;    % T/N_t;     % Time-step such that we ensure at the maximum stability - Delta_t/Delta_x < 1
            
            x = linspace(0, L, N);  % Vector x to be used for initial wave's height

            alpha = Dx/Dt;    % Parameter for discrete equation
            
        % Vector containing u = (h,rho)^T   % rho = hv
            u = zeros(2,N+2); % We consider h the upper vector and hv the lower one
        
% ------------------------------------------------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------------------------------------------------ %
        
    % Conditions : initial and boundary
        
        % Speed
            % Initial condition
            for k = 1:N+2   
                u(2,k) = 0;     % Initial moment m = h*u = 0
            end
     
        % Height
            % Initial condition
            % We keep h(1) and h(n+2) at zero and the rest considering the
            % initial height
            for k = 1:N
               u(1,k+1) = H + epsilon*exp(-((x(k) - L/2)^2)/w^2); 
            end
  
            figure(1)     % Graphic of the initial height 
            plot(x, u(1,2:N+1), 'r'); % We represent only the point of h for the considered initial value (size(x,2) points)
            xlabel('Position x (m)'); ylabel('Initial height (m)'); title('Graphic representing the function h = f(x) at t=0s');
            
            
% ------------------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------------------ %

        
    % Flux : 
        %   F(u_j+1,u_j)^(F*L_x) - F(u_j,u_j-1)^(F*L_x)
        %
        %   --> 1/2[f(u_j+1) + f(u_j-1) - alpha(u_j+1 - 2*u_j + u_j-1)]
        %
        
        
        % Vector to keep in memory the iteration of the height and speed
            U = zeros(2*(N_t+1),N+2);
            U = [u]; 

        for t = 1:N_t
            
            % Reflecting boundary for the height + boundary for the moment
                u(1,1) = u(1,2);    u(1,N+2) = u(1,N+1);    % Same height
                
                u(2,1) = -u(2,2);   u(2,N+2) = -u(2,N+1);   % Opposed momentum (and thus opposed speed)       
                
            % Variables from general equation's form 
                % Function f to work on F (Lax-Friedrichs flux) :
                    f = [u(2,:) ; (u(2,:).^2)./u(1,:) + (1/2)*g*(u(1,:).^2)];
                
            % We calulate the speed and height
                % What is at u(:,1) and u(:,n+2) is already defined : 
                %   - v is always zero at u(:,1) and u(:,n+2)
                %   - h(1) and h(n+2) as well as m(1) and m(n+2) (m = h*v) are defined above at the begining of the loop
                
                % Lax-Friedrichs flux
                for k = 1:N+1
                    F1(:,k) = (1/2)*(f(:,k+1) + f(:,k) - alpha*(u(:,k+1) - u(:,k)));
                end
                
                % Resolution of the system to get h and h*v
                for k = 2:N+1
                    u(:,k) = u(:,k) - (Dt/Dx)*(F1(:,k) - F1(:,k-1));
                end

            % We stock the value for the different iteration (time
            % depending function U = (h,h*v))
                 U = [U; u];
        end
      
        
% ------------------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------------------ %

%% Display of the wave

        % U(2*k,2:end-1)./U(2*(k-1) + 1,2:end-1) to get v
        % U(2*(k-1)+1,2:end-1) to get h
        
        figure(2)
        for k = 1:N_t+1
            plot(x, U(2*(k-1)+1,2:end-1), 'r');
            axis([0 10 0.95 1.1]);
            xlabel('Position on the mesh (m)'); ylabel('Wave''s height (m)'); title(['Wave at t = ',num2str((k-1)*Dt,2),' s']);
            drawnow
        end
        
        % Graphs of the height at different timeline
        t1 = 1; t2 = 50; t3 = 100; t4 = 101;    % Timelines
        
        figure(3)
            subplot(2,2,1)
            plot(x, U(2*t1-1,2:end-1), 'r');
            axis([0 10 1 1.1]);
            xlabel('Position on the mesh (m)'); ylabel('Wave height (m)'); title(['Wave at t = ',num2str(t1*Dt,2),' s']);

            subplot(2,2,2)
            plot(x, U(2*t2-1,2:end-1), 'r');
            axis([0 10 1 1.1]);
            xlabel('Position on the mesh (m)'); ylabel('Wave height (m)'); title(['Wave at t = ',num2str(t2*Dt,2),' s']);

            subplot(2,2,3)
            plot(x, U(2*t3-1,2:end-1), 'r');
            axis([0 10 1 1.1]);
            xlabel('Position on the mesh (m)'); ylabel('Wave height (m)'); title(['Wave at t = ',num2str(t3*Dt,2),' s']);

            subplot(2,2,4)
            plot(x, U(2*t4-1,2:end-1), 'r');
            axis([0 10 1 1.1]);
            xlabel('Position on the mesh (m)'); ylabel('Wave height (m)'); title(['Wave at t = ',num2str(t4*Dt,2),' s']);

