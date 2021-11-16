%% 2.2 - Tests with non-horizontal bottom - Type 4 - entropy term

    % B = B(x) --> s(h,m,x) != 0
    
    % General variables
        
        L = 10;         % Length of the mesh
        T = 20;         % Time defining the study (we look at the function for T second(s))
        
        N_x = 200;      % Number of point defining the mesh
        N_t = 3000;     % Number of iteration
        
        H = 1;          % Initial height of the wave
        w = 0.1*L;      % Width of the Gaussian pulse
        a = H/5;        % Coeff to ensure a stable scheme (low right hand side with respect to H)
        
        g = 9.81;       % Gravity's acceleration
        
    % Bottom
        B0  = H/10;
        r   = L/6;
        
    % Entropy
        delta = 0.1;
        
    % Mesh
        x = linspace(0, L, N_x);    % Position on the mesh [0;10] m
        
        Delta_x = L/(N_x-1);        % Space step (with respect to the one implicitly used on the linspace)
        Delta_t = T/N_t;            % Time step
        
    % Checking stability (CFL condition)
        CFL = (4*Delta_t)/Delta_x;
    
        if CFL > 1
            disp('error CFL conition not fulfilled');
        end 
       
% ----------------------------------------------------------------------------------------------- %
% ----------------------------------------------------------------------------------------------- %
        
    % Vector solution
        u = zeros(2, N_x+2);    % 2 lines (h & m) and N_x columns, one for each points of the mesh
        B = zeros(1, N_x);      % Function of the bottom
        S = zeros(1, N_x);      % Function s(x,m,h) -> consideration of a non-horizontal bottom
        
    % IC :
        % Bottom
            for k = 1:N_x
                if abs(x(k) - L/2) < r
                    B(k) = B0*cos((pi*(x(k) - L/2))/(2*r))^2; 
                else
                    B(k) = 0;
                end
            end
         
        % Height
            for k = 1:N_x
                u(1,k+1) = 1.5;
            end
            
        % Moment
            for k = 2:N_x+1
                u(2,k) = 5.5;
            end
            
        % Dirichlet condition for inflow
            Inflow_Condition = 5.5;

        % Variables to stock the iteration of u
        % Even lines = m ; Odd lines = h
            U_2 = zeros(2*(N_t+1), N_x+2);
            U_2 = [u]; 

        % Variables to stock the last iteration of U
        % Even lines = m ; Odd lines = h
        % Purpose it to compare graphs for several delta
            U_delta = zeros(size(delta,2), size(U_2,2));
            U_delta = [U_2];

        for t = 1:N_t            
                
            % BC
                % Height
                    u(1,1)      = 1.5;    % Dirichlet BC for ghost cell (j =  -1) 
                    u(1,N_x+2)  = 1.5;    % Dirichlet BC for ghost cell (j = N+2)

                % Moment 
                    % Inflow by Dirichlet BC at the ghost cell (j = -1)
                        u(2,1) = Inflow_Condition;  

                    % Outflow by extrapolation (j = N+2)
                        u(2,N_x+2) = 2*u(2,N_x+1) - u(2,N_x);

            % Function f to work on F (Roe numerical flux) :
                f = [u(2,:) ; (u(2,:).^2)./u(1,:) + (1/2)*g*(u(1,:).^2)];
                
            % Eigen-values & Eigen-vectors
                for k = 1:N_x+1
                    h_tilde(k) = (1/2)*(u(1,k+1) + u(1,k));
                    u_tilde(k) = ((u(1,k+1)^(1/2))*u(2,k+1)/u(1,k+1) + (u(1,k)^(1/2))*u(2,k)/u(1,k))/(u(1,k+1)^(1/2) + u(1,k)^(1/2));
                    c_tilde(k) = (g*h_tilde(k)).^(1/2);

                    l_1(t,k) = u_tilde(k) - c_tilde(k);
                    l_2(t,k) = u_tilde(k) + c_tilde(k);
                    
                    a_1(k) = ((u_tilde(k) + c_tilde(k)).*(u(1,k+1) - u(1,k)) - (u(2,k+1) - u(2,k)))./(2*c_tilde(k));
                    a_2(k) = ((u(2,k+1) - u(2,k)) - (u_tilde(k) - c_tilde(k)).*(u(1,k+1) - u(1,k)))./(2*c_tilde(k));
                end
                    
            % Simplifications
                One = linspace(1,1,size(l_1,2));    % Recall : size(l_1,2) = size(l_2,2)
                
                r_1 = [One; l_1(t,:)];
                r_2 = [One; l_2(t,:)];
                
                W_1 = r_1.*a_1;
                W_2 = r_2.*a_2;
                
            % Harten's entropy fix
                for k = 1:N_x+1
                    % lambda 1 --> Phi(1)
                        if abs(l_1(t,k)) >= delta
                            Phi(1,k) = abs(l_1(t,k));
                        else
                            Phi(1,k) = (l_1(t,k)^2 + delta^2)/(2*delta);
                        end
                    % lambda 2 --> Phi(2)
                        if abs(l_2(t,k)) >= delta
                            Phi(2,k) = abs(l_2(t,k));
                        else
                            Phi(2,k) = (l_2(t,k)^2 + delta^2)/(2*delta);
                        end
                end
                
                % Viscous terms
                    Visc_entropy_minus = [(1/2)*(l_1(t,:) - Phi(1,:)); (1/2)*(l_2(t,:) - Phi(2,:))];
                    Visc_entropy_plus  = [(1/2)*(l_1(t,:) + Phi(1,:)); (1/2)*(l_2(t,:) + Phi(2,:))];

                    Visc_entropy = Visc_entropy_plus - Visc_entropy_minus;

                    Visc_1 = Visc_entropy(1,:).*W_1;
                    Visc_2 = Visc_entropy(2,:).*W_2;

                % Roe flux
                    for k = 1:N_x+1                    
                        F(:,k) = (1/2)*(f(:,k+1) + f(:,k)) - (1/2)*(Visc_1(:,k) + Visc_2(:,k));
                    end
                
                % Source term
                    for k = 1:N_x
                        if abs(x(k) - L/2) < r
                            S(2,k) = (g*u(1,k)*B0*pi/(2*r))*sin((pi*(x(k) - L/2))/r);
                        else
                            S(2,k) = 0;
                        end
                    end
                
                % Roe scheme
                    for k = 2:N_x+1
                        u(:,k) = u(:,k) - (Delta_t/Delta_x)*(F(:,k) - F(:,k-1)) + Delta_t*S(:,k-1);
                    end

            % We stock the value for the different iteration (time
            % depending function U = (h,h*v))
                U_2 = [U_2; u];
        end
        
        U_delta = [U_delta; U_2(end-1:end,:)];    
    
%% Display of the wave - drawnow

        % U_2(2*k + 2,2:end-1)./U_2(2*(k-1) + 1,2:end-1)    to get u
        % U_2(2*k + 1,2:end-1)                              to get h
        
        % Variables axis
            % Height 
                x_lim_inf_h = x(1);     x_lim_sup_h = x(end);
                y_lim_inf_h = 0;        y_lim_sup_h = 3;  
                
            % Speed 
                x_lim_inf_s = x(1);     x_lim_sup_s = x(end);
                y_lim_inf_s = 0;        y_lim_sup_s = 2;
        
        % Height
%             figure(110)
%             for k = 1:N_t/2
%                 plot(x, U_2(2*k + 1,2:end-1) + B, 'r');
%                 axis([x_lim_inf_h x_lim_sup_h y_lim_inf_h y_lim_sup_h]);
%                 xlabel('Position on the mesh [m]'); ylabel('Wave''s height + B [m]'); title(['t = ',num2str(round((T/N_t)*k,2)),' s']);
%                 drawnow
%             end

        % Speed
%             figure(111)
%             for k = 1:N_t/2
%                 plot(x, U_2(2*k + 2,2:end-1)./U_2(2*k + 1,2:end-1), 'r');
%                 axis([0 10 0 2]);
%                 xlabel('Position on the mesh [m]'); ylabel('Wave''s speed [m.s^{-1}]'); title(['t = ',num2str(round((T/N_t)*k,2)),' s']);
%                 drawnow
%             end
           
        % Eigen-values
%             figure(112)
%             for k = 1:N_t/2
%                 plot(x, l_2(k,1:end-1), 'r');
%                 axis([0 10 3 4]);
%                 xlabel('Position on the mesh [m]'); ylabel('Wave''s low eigenvalue [m.s^{-1}]'); title(['t = ',num2str(round((T/N_t)*k,2)),' s']);
%                 drawnow
%             end

%% Presentation of results
    
    figure(113)
        k = N_t;
        p_113(1:3) = plot(x, B, 'r--', x, U_2(2*k + 1, 2:end-1) + B, 'b--', x, U_2(2*k + 2,2:end-1)./U_2(2*k + 1,2:end-1), 'g-');
        hold on;
        p_113(4) = plot([0,10], [3.83,3.83], 'r-');
        axis([0 10 -0.2 5.5]);
        xlabel('Position on the mesh [m]'); title(['t = ',num2str(round(k*T/N_t,2)),' s']);
        grid on;
        legend([p_113], {'B', 'h+B', 'u', 'Sonic speed'}, 'location', 'northwest');
        hold off;