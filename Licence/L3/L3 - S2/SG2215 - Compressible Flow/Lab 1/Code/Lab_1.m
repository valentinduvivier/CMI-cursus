
clear all; close all; clc;

%% Compressible pipe flow

    % Variables
        % Tube
            cf = 0.005; 
            pipe_length             = 2670;     % [mm]
            maximal_mass_flow_rate  = 25;       % [g.s-1] 

            d_1     = 19.5; % [mm]  % Inlet tube's diameter
            d_2     = 10;   % [mm]  % Outlet tube's diameter
            D       = 6.4;  % [mm]  % Diameter of the pipe tube
            
        % Gas
            R   = 8.134;    % [J/(mol.K)]   % Universal perfect gases' constant
            r   = 287;      % [J/(kg.K)]    % Specific perfect gases' constant
        
        % General
            pi  = 3.14;     % []
            g   = 1.4;      % []    % g = gamma
            T0  = 20;       % [°C]   % Stagnation inlet temperature
            
        % We create the vector containing the position of the taps along the
        % tube
        	x = [0, 76, 760, 1520, 2280, 2400, 2530, 2600, 2630, 2670];  % [mm]
        
        % We load the three files from experiment
            load('Data.mat');

            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                           SI Units                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Tube
            pipe_length     = pipe_length*10^-3;        % [m]
            
            d_1     = d_1*10^-3;    % [m]
            d_2     = d_2*10^-3;    % [m]
            D       = D*10^-3;      % [m]

        % General
            T0  = T0 + 273;     % [K]   % Stagnation inlet temperature
            
        	x   = x*10^-3;      % [m]
            
        
%% 4 - Home assignement

    % a
    
        M_2__HA = 1;
    
        % We get M_1 from a numerical approximation of our nonlinear equation
        
            % Link to the function called "fonction"
                fun = @fonction;

                % Initial point close to zero to facilitate/ensure convergence
                % toward a physical result
                x0 = 0.1;

                % M_1 (see the function to change the points (x2 - x1 considered)
                    M_HA = fsolve(fun, x0);  
        
            % We get that :
                M_1__HA = 0.2532;   % M_1__HA = Mach number at the pipe inlet
                disp(['The inlet Mach number to get to M2 = 1 is M1 = ', num2str(M_1__HA),'']);
                               
    % b
    
        % If p2 = 100kPa & M2 = 1 (we use non-isentropic relation) :
            
            % Given
                p2__HA = 100;  % [kPa]

            % Pressure ratio for an non-isentropic tube
                p1__HA = (M_2__HA/M_1__HA)*(((2 + (g-1)*M_2__HA^2)/(2 + (g-1)*M_1__HA^2))^(1/2))*p2__HA;

            disp(['The inlet pressure to get to M2 = 1 is P0 = ', num2str(p1__HA),'kPa']);
        
    % c 
    
        % We are willing to determine Mi to then get the associated pressures
            
            % Mach repartition over the tube
                M_1__HA = [0.2532, 0.2562, 0.2892, 0.3494, 0.4970, 0.5517, 0.6544, 0.7780, 0.9949, 1.0];        
            
        % Work on the ratio pi/p0inlet
            
            % Stagnation state relation (for Venturi tube)
                T_1__HA = T0/(1+((g-1)*M_1__HA(1)^2)/2);
            
            % Stagnation state relation (for Venturi tube)
                p_0_inlet = p1__HA*(((1 + ((g-1)*M_1__HA(1)^2)/2))^(g/(g-1)));
            
            % Pressure ratio (for the pipe tube)
                pressure_ratio = (M_2__HA./M_1__HA(2:end-1)).*(((2 + (g-1)*M_2__HA^2)./(2 + (g-1).*M_1__HA(2:end-1).^2)).^(1/2))*(p2__HA/p_0_inlet);
            
            
%% 6 - Evaluation of Data
   
    % 1 
    
    % We get T0_2 form the perfect gas law for an adiabatic flow (existence of gamma) : P(V^gamma) = cst 
    
    % The pressure measured at the Venturi section are obtained from the
    % file Exp1.txt
    
    % We create the vectors of pressure from the file        
        p0 = Exp1(1,3:end)*10^3;        % [Pa]      % Stagnation pressure at the inlet
        
        pn = Exp1(2:3,3:end)*10^3;      % [Pa]      % pn1 = pn(1,:) & pn2 = pn(2,:)         % Inlet and outlet pressure at the Venturi section
        pp = Exp1(4:end-2,3:end)*10^3;  % [Pa]      % pp1 = pp(1,:) & pp2 = pp(2,:) & etc   % Pressure along the pipe
        
        pB = Exp1(end-1,3:end)*10^3;    % [Pa]      % Stagnation outlet pressure
        ppitot = Exp1(end,3:end)*10^3;  % [Pa]      % Pressure form the Pitot tube
        
                
        % a - mass flow rate
            
        % We look at the pressure pn2  to compare the resultats with pB (?)
        
            % Calculation of the Mach number at the inlet & outlet
                M_2 = ((2*((p0./pn(2,:)).^((g-1)/g) - 1))/(g-1)).^(1/2);

            % We deduce T_V form M_V :
                T_2 = (T0./(1 + ((g-1).*M_2.^2)/2));

            % We the get rho_V form : P/rT = rho
                rho_2 =  pn(2,:)./(r.*T_2);

            % Necessary variables to get the mass flow rate :
                A_1 = pi*(d_1/2)^2;
                A_2 = pi*(d_2/2)^2;
                A = [A_1; A_2];
                
                u_2 = M_2.*((g*r.*T_2).^(1/2));
                
            mass_flow_Venturi = rho_2.*A_2.*u_2;

            
        % Plot
        % We plot the mass flow rate over pB (outlet stagnation pressure)
            figure(1)
            plot(pB*10^-3, mass_flow_Venturi*10^3, 'c*-', 'linewidth', 2);   % We worked above with SI units but we display the results in kPa & g.s-1
            xlabel('Back pressure pB [kPa]'); ylabel('Mass flow rate [kg.s-1]'); title('Mass flow rate over pB');
            
        % b
      
            % Different valve pressure = different 
            figure(2)
            for k = 1:size(pp,2)
                l2(k) = plot(x(2:end-1)*10^-3, pp(:,k)*10^-3, '*-', 'color', k*[1/17 1/20 1/80]);
                hold on
            end
            xlabel('pB [kPa]'); ylabel('Pressure pp [kPa]'); title('Pressure along the pipe over pB');
            legend([l2],{['pB = ', num2str(round(pB(1)*10^-3)),' kPa'], ['pB = ', num2str(round(pB(2)*10^-3)),' kPa'], ['pB = ', num2str(round(pB(3)*10^-3)),' kPa'], ...
            ['pB = ', num2str(round(pB(4)*10^-3)),' kPa'], ['pB = ', num2str(round(pB(5)*10^-3)),' kPa'], ['pB = ', num2str(round(pB(6)*10^-3)),' kPa'], ['pB = ', num2str(round(pB(7)*10^-3)),' kPa'], ...
            ['pB = ', num2str(round(pB(8)*10^-3)),' kPa'], ['pB = ', num2str(round(pB(9)*10^-3)),' kPa'], ['pB = ', num2str(round(pB(10)*10^-3)),' kPa'], ['pB = ', num2str(round(pB(11)*10^-3)),' kPa'], ...
            ['pB = ', num2str(round(pB(12)*10^-3)),' kPa'], ['pB = ', num2str(round(pB(13)*10^-3)),' kPa'], ['pB = ', num2str(round(pB(14)*10^-3)),' kPa'], ['pB = ', num2str(round(pB(15)*10^-3)),' kPa']}, 'location', 'southwest');
          
            
%% ----------------------------------------------------------------------------------------------           
   
    % 2
    
    % Non-dimensional value for position
    % To use those will allow us to consider the pressure regarding to a
    % normalised length (so considering the diameter at this position (?))
    
        x_D = x(2:end-1)./D;     % inner_diameter = 6.4
        
        p0 = Exp2(1,3);         % Stagnation pressure at the inlet
        
        pn = Exp2(2:3,3);       % pn1 = pn(1,:) & pn2 = pn(2,:)         % Inlet and outlet pressure at the Venturi section
        pp = Exp2(4:end-2,3);   % pp1 = pp(1,:) & pp2 = pp(2,:) & etc   % Pressure along the pipe
        
        pB = Exp2(end-1,3);     % Stagnation outlet pressure
        ppitot = Exp2(end,3);   % Pressure form the Pitot tube
        
        % a - Plot Pp over x/D
            figure(3)
            plot(x_D, (pp/p0)', '*-', 'linewidth', 2);
            xlabel('Normalized position x/D [ ]'); ylabel('Pressure ratio pp/p0 [ ]'); title('Pressure ratio along the pipe for chocked conditions');
    
        % b - Calculation Mach number & Plot Mach over x/D
             
            % 1 - Calculation Mach number
                M_2 = 1;
                
                M = fsolve(fun, x0); % We get the beloww mach number using numerical approximation 
                
                % By opening the function "fonction.m", you have a
                % dedicated line for this part and you can modify the
                % considered pressure pp_i you want
                
                M_1 = [0.2589, 0.2886, 0.3464, 0.4889, 0.5254, 0.5865, 0.6332, 0.7518];
            
            % 2 - Plot Mach over position difference (xfinal - x1,x2,...,x8)
            
                figure(4)
                plot(x_D, M_1, 'r*-', 'linewidth', 2);
                xlabel('Normalized position x/D [ ]'); ylabel('Mach number [ ]'); title('Mach number at each tap against normalized distance from inlet');

            
        % c - Comparison with home assignement

            figure(5)
            l5_1 = plot(x_D, pp/p0, 'b*-', 'linewidth', 2);
            hold on
            l5_2 = plot(x_D, pressure_ratio, 'g*-', 'linewidth', 2);
            hold off
            xlabel('Normalized position x/D [ ]'); ylabel('Pressure ratio [ ]'); title('Comparison pressure ratio for theory and experiment');
            legend([l5_1, l5_2], {'Experiment', 'Theory'}, 'location', 'southwest')
            
            figure(6)
           	l6_1 = plot(x_D, M_1, 'r*-', 'linewidth', 2);
            hold on
            l6_2 = plot(x_D, M_1__HA(2:end-1), 'm*-', 'linewidth', 2);
            hold off
            xlabel('Normalized position x/D [ ]'); ylabel('Mach number [ ]'); title('Comparison Mach number for theory and experiment');
            legend([l6_1, l6_2], {'Experiment', 'Theory'}, 'location', 'northwest')        

            
%% ----------------------------------------------------------------------------------------------           

    % 3

        % Defining pressures from Exp 3 
            p0_3 = Exp3(1,:);         % Stagnation pressure at the inlet

            pn_3 = Exp3(2:3,:);       % pn1 = pn(1,:) & pn2 = pn(2,:)         % Inlet and outlet pressure at the Venturi section
            pp_3 = Exp3(4:end-2,:);   % pp1 = pp(1,:) & pp2 = pp(2,:) & etc   % Pressure along the pipe

            pB_3 = Exp3(end-1,:);     % Stagnation outlet pressure
            ppitot_3 = Exp3(end,:);   % Pressure form the Pitot tube

        % a    
            IR = ppitot_3./pB_3;      % Isentropic pressure ratio

            % Mach number found with the isentropic ratio
                M1_IR = (2/(g-1)*((IR).^((g-1)/g) - 1)).^(1/2);

            % Mach number found with the Rayleigh Pitot tube formula
                M1_RP = [0,1.0196,1.121,1.1472,1.15316,1.13535,1.07172,1.00754,0.8418,0,0];

                
            x = linspace(0,1,11);     % Radial position array

            figure(7)
            plot(x,real(M1_IR),'b+-','LineWidth',1)
            hold on
            plot(x,IR1,'ro--','LineWidth',1)
            grid on
            xlabel('Radial position [-]')
            ylabel('Mach number [-]')
            legend('Isentropic ratio','Rayleigh Pitot tube formula','Location','South')
            hold off