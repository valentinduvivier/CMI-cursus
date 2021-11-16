            function F = fonction(X)
                     
            % Heat capacity ratio
                g_air = 1.4;    % [-]
                g_hel = 1.667;  % [-]

            % Gas specific constant
                r_air = 287.058;    % [J/(kg*K)]
                r_hel = 2076.90;    % [J/(kg*K)]
            % Temperature at driven and driver sections
                T_1 = 300;  % [K]
                T_4 = 300;  % [K]

            % Accélerations
                a_1 = (g_air*r_air*T_1)^(1/2);
                a_4 = (g_hel*r_hel*T_4)^(1/2);
                
                M_2 = 1.04;
                
                % theoretical_opti_M_s
                    F = 0.31*X^3 - 2.6*X^2 + 8.1*X - 5.5  -  3.2;
                
                % theoretical Mach number M_s_theory_2
%                     F = (2*(X.^2 - 1))./(((2*g_air*X.^2 - (g_air - 1)).^(1/2)) * (((g_air - 1)*X.^2 + 2).^(1/2)))  -  M_2;            
            end
            