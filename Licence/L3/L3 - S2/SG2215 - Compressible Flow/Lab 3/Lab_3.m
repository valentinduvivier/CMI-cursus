%% Lab 3 : Pressure distribution on a biconvex airfoil in supersonic flow

    clear all; close all; clc;

    % Theory
    
        % DATA 
            % Flow data
                M_inf   = 2.0;  % Given freestream Mach number

                AOA     = 1.5;  % Given angle of attack converted in radians
                AOA     = AOA*pi/180;
                
                g = 1.4;    % Heat capacity ratio
                
            % Experiment data
                x = [0, 7, 17, 27, 37, 47, 57, 67, 70]; % Pressure taps positions
                
                r = 210;    % Radius for bi-convex airfoil
                
                c = 70;     % Chord length
    
        %% Linearized theory
        
            % We calculate the angles phi conresponding to each pressure
            % taps on the lower side (as we have an infinity of shock waves
            % we can have as much angles phi as we have of considered points x
            % i.e. size(phi) = size(x)
                phi = asin((c/2 - x)/r);                    
            
            % As well we deduce the angles Theta
                Theta_up = phi - AOA; 
                Theta_low = phi + AOA; 
                
                Theta = [Theta_low; Theta_up];
                
            % Cp profile
            
                Cp_linear = 2*Theta/sqrt(M_inf^2 - 1);  % Cp_linear_low = Cp_linear(1) & Cp_linear_up = Cp_linear(2)
                
               figure(1)
                plot(x/c, Cp_linear(1,:), 'r--', x/c, Cp_linear(2,:), 'k--', 'linewidth', 2);
                hold on
                legend({'Linear theory lower surface', 'Linear theory upper surface'}, 'location', 'northwest')
                grid on

                
        %% Shock-Expansion theory
        
            % Data at LE
                M_1 = M_inf;
                
                Theta_1 = Theta(:,1);     % Initial angle Theta for each side
                
                % p_inf = p1 % Useful relation once we calculate the pressure ratios
            
        %% Lower & Upper surfaces

            % Numerical resolution
                x0 = 0.7; % Value close from the intended one
                
            % We loop to get Beta for both lower and upper surfaces
                for k = 1:2
                    Beta(k) = fsolve(@(x) 2*cot(x).*(M_1^2*sin(x).^2-1)./(M_1^2*(g + cos(2*x)) + 2) - tan(Theta_1(k,1)), x0);
                end
                
                Beta = Beta';

                % Beta_low = Beta(1) & Beta_up = Beta(2); and so for the
                % following variables, i.e. (1) = low & (2) = up
                
                
                % Normal Mach number M_n1
                    M_n1 = M_1.*sin(Beta);
                
                % Pressure ratio from normal shock relation
                    p_ratio_2_1 = 1 + ((2*g)/(g+1)).*(M_n1.^2 - 1);
                    
                % Normal Mach number M_n2 from normal shock relation
                    M_n2 = sqrt((2 + (g-1)*M_n1.^2)./(2*g.*M_n1.^2 - (g-1)));
                
                % M_2
                    M_2 = M_n2./sin(Beta - Theta_1);
                
                % v_2 - Prandtl Meyer angle (rad) from M_2
                    v_2 = sqrt((g+1)/(g-1))*atan(sqrt((g-1)/(g+1).*(M_2.^2 - 1))) - atan(sqrt(M_2.^2 - 1));
                
                % v_k - Prandtl Meyer angle for every other pressure tap location
                    v = abs(Theta - Theta_1) + v_2;
                    
            %% Function to get the associated Mach number M_k
                
                x0 = 2;
             
            % M_k from Prandtle Meyer function
                for i = 1:2
                    for k = 1:size(x,2)
                        M(i,k) = fsolve(@(x) - v(i,k) + sqrt((g+1)/(g-1))*atan(sqrt(((g-1)/(g+1))*(x^2 - 1))) - atan(sqrt(x^2 - 1)), x0);
                    end
                end

%% Pressure ratios
                % Pressure ratio p_k/p_2
                    p_ratio_2_k = ((2 + (g-1).*M.^2)./(2 + (g-1).*M_2.^2)).^(g/(g-1));
                 
                % Pressure ratio p_k/p_inf
                    p_ratio_2_inf = p_ratio_2_1;    % as p_inf = p_1

                    p_ratio_k_inf = p_ratio_2_inf./p_ratio_2_k;

                % Pressure profile
                    Cp_expansion = ((2./(g*M_inf^2)).*(p_ratio_k_inf - 1));
                
                % Displaying
                   figure(2)
                    plot(x/c, Cp_expansion(1,:), 'b*-', x/c, Cp_expansion(2,:), 'k-');
                    grid on
                    xlabel('Cp'); ylabel('x/c');
                    legend({'Shock-Expansion theory : Lower surface', 'Shock-Expansion theory : Upper surface'}, 'location', 'northwest');
                    hold on
                

%% Displaying results - Curves comparison

                   figure(3)

                    plot(x/c, Cp_linear(1,:), 'r--', x/c, Cp_linear(2,:), 'k--', x/c, Cp_expansion(1,:), 'r-', x/c, Cp_expansion(2,:), 'k-');
                    grid on
                    xlabel('x/c'); ylabel('Cp');
                    legend({'Linear theory lower surface', 'Linear theory upper surface', 'Shock-Expansion theory : Lower surface', 'Shock-Expansion theory : Upper surface'}, 'location', 'northeast');
                    hold off
    
    