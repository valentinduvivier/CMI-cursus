
            fun = @fonction;

            % Initial point close to zero to facilitate/ensure convergence
            % toward a physical result
                x0 = 1;

            % M_s (see the function to change the points (x2 - x1 considered)
                M_s_theory_opti_2 = fsolve(fun, x0);  
                
                M_2_theory_opti_2 = (2*(M_s_theory_opti_2^2 - 1))/(((2*g_air*M_s_theory_opti_2^2 - (g_air - 1))^(1/2)) * (((g_air - 1)*M_s_theory_opti_2^2 + 2)^(1/2)));            


                
%% Flow in shock tube

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                  %
%                 General Variables                %
%                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    % Mach number
        M_2 = [0.46, 1.04, 1.34];
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                  %
%%                    Experiment 2                 %
%                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % We use here M_2 = M_2(2)

    % 1
        % a - to get M_s out of M_2 through numeric method
        
            % Function giving M_s if provided M_2
                fun = @fonction;

            % Initial point close to the final result to facilitate/ensure convergence
            % toward a physical result
                x0 = 1;

            % M_s
                M_s_theory_2 = fsolve(fun, x0);  
                
        % b - theoretical ratio p4/p1 from M_s
            ratio_P4_P1_theory_2  = ((2*g_air*M_s_theory_2^2 - (g_air - 1))/(g_air + 1))*...
                ((1 - (((g_hel-1)*a_1)/((g_air+1)*a_4))*(M_s_theory_2 - M_s_theory_2^-1))^(-2*g_hel/(g_hel - 1)));
            log_ratio_P4_P1_theory_2 = log(ratio_P4_P1_theory_2);  
            
        % c - "experimental" ratio p4/p1 from M_s
            M_s_corr_2 = 0.31*M_s_theory_2^3 - 2.6*M_s_theory_2^2 + 8.1*M_s_theory_2 - 5.5;
            ratio_P4_P1_exp_2     = exp(M_s_corr_2);     % From measurement in figure 6
            
    % 2
        % a - Experimental M_s from "Oscilloscope_trace.pdf"
            P_1_2 = 10*10^3;    % [kPa]
            P_2_2 = 0.52*10^5;  % [kPa]
                
            M_s_exp_2 = (((g_air + 1)/(2*g_air))*(P_2_2/P_1_2 - 1) + 1)^(1/2);   % [-]
            
        % b - Experimental M_2 deduced from M_s_exp_2
            M_2_exp_2 = (2*(M_s_exp_2^2 - 1))/(((2*g_air*M_s_exp_2^2 - (g_air - 1))^(1/2)) * (((g_air - 1)*M_s_exp_2^2 + 2)^(1/2)));            
        
            absolute_error_M_2_2 = abs(M_2_exp_2 - M_2(2))/M_2(2);
                
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                  %
%%                    Experiment 3                 %
%                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % M_2 = M_2(3);

    % 1
        % a - to get M_s out of M_2 through numeric method
        
            % Function giving M_s if provided M_2
                fun = @fonction;

            % Initial point close to zero to facilitate/ensure convergence
            % toward a physical result
                x0 = 1;

            % M_s (see the function to change the points (x2 - x1 considered)
                M_s_theory_3 = fsolve(fun, x0);  
                
        % b - theoretical ratio p4/p1 from M_s
            ratio_P4_P1_theory_3  = ((2*g_air*M_s_theory_3^2 - (g_air - 1))/(g_air + 1))...
                *((1 - (((g_hel-1)*a_1)/((g_air+1)*a_4))*(M_s_theory_3 - M_s_theory_3^-1))^(-2*g_hel/(g_hel - 1)));
            log_ratio_P4_P1_theory_3 = log(ratio_P4_P1_theory_3);  
            
        % c - "experimental" ratio p4/p1 from M_s
            M_s_corr_3 = 0.31*M_s_theory_3^3 - 2.6*M_s_theory_3^2 + 8.1*M_s_theory_3 - 5.5;
            ratio_P4_P1_exp_3     = exp(M_s_corr_3);     % From measurement in figure 6
             
    % 2
        % a - Experimental M_s from "Oscilloscope_trace.pdf"
            P_1_3 = 10*10^3;    % [kPa]
            P_2_3 = 1*10^5;     % [kPa]
                
            M_s_exp_3 = (((g_air + 1)/(2*g_air))*(P_2_3/P_1_3 - 1) + 1)^(1/2);   % [-]
            
        % b - Experimental M_2 deduced from M_s_exp_3
            M_2_exp_3 = (2*(M_s_exp_3^2 - 1))/(((2*g_air*M_s_exp_3^2 - (g_air - 1))^(1/2)) * (((g_air - 1)*M_s_exp_3^2 + 2)^(1/2)));            
        
            absolute_error_M_2_3 = abs(M_2_exp_3 - M_2(3))/M_2(3);               
                
                
        % c - Experimental M_2 obtained from shadow graph
            % Variables
                % TOP
                    Beta_top    = 65.8*pi/180;     % [ °]
                    Theta_top   = 8*pi/180;        % [ °]

                % BOTTOM
                    Beta_bot    = 55*pi/180;       % [ °]
                    Theta_bot   = 4*pi/180;        % [ °]
                
            % Mach numbers
                % Function giving M_2 if provided Theta and beta
                    fun = @fonction;

                % Initial point close to zero to facilitate/ensure convergence
                % toward a physical result
                    x0 = 1;

                % M_s (see the function to change the points (x2 - x1 considered)
%                     M_2_exp_3_sg_top = fsolve(fun, x0);         
%                     M_2_exp_3_sg_bot = fsolve(fun, x0);         
                
                    
            absolute_error_M_2_3_sg_top = abs(M_2_exp_3_sg_top - M_2(3))/M_2(3);               
            absolute_error_M_2_3_sg_bot = abs(M_2_exp_3_sg_bot - M_2(3))/M_2(3);               
