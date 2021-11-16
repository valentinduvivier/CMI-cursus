%% 2.1 - Tests with flat bottom

    % B = 0 --> s(h,m,x) = 0 as flat bottom
    
    % General variables
        
        L = 10;         % Length of the mesh
        T = 20;         % Time defining the study (we look at the function for T second(s))
        
        N_x = 320;      % Number of point defining the mesh
        N_t = 3000;     % Number of iteration
        
        H = 1;          % Initial height of the wave
        w = 0.1*L;      % Width of the Gaussian pulse
        a = H/5;        % Coeff to ensure a stable scheme (low right hand side with respect to H)
        
        g = 9.81;       % Gravity's acceleration
        
    % Mesh
        x = linspace(0, L, N_x);
        
        Delta_x = L/(N_x-1);    % Space step
        Delta_t = T/N_t;        % Time step
        
    % Checking stability (CFL condition)
        CFL = (4*Delta_t)/Delta_x;

        if CFL > 1
            disp('error CFL condition not fulfilled');
        end
        
% ----------------------------------------------------------------------------------------------- %
% ----------------------------------------------------------------------------------------------- %

    % Vector solution
        u = zeros(2, N_x+2);  % 2 lines (h & m) and N_x columns, one for each points of the mesh
        
    % IC :
        % Height
            for k = 1:N_x
                u(1,k+1) = H + a*exp(-(x(k)-L/2)^2)/(w^2);
            end
        
        % Moment
            for k = 2:N_x+1
                u(2,k) = (u(1,k) - H)*sqrt(g*u(1,k));
            end
            
    % Display initial data to make sure the initialization went well
%         figure(1)
%         plot(x, u(1,2:end-1),'r-', 'linewidth', 2);
%         xlabel('Position on the mesh [m]'); ylabel('Wave''s height'); title('Initial height of the wave on the mesh');
%             
        % Variables to stock the iteration of m & h
        % Even lines = m; Odd lines = h
            U_2 = zeros(2*(N_t+1), N_x+2);
            U_2 = [u]; 

        for t = 1:N_t
            
            % "Wall" (bouncing) BC :
                u(1,1) = u(1,2);    u(1,N_x+2) = u(1,N_x+1);
                
                u(2,1) = -u(2,2);   u(2,N_x+2) = -u(2,N_x+1);            
                
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
                
                Visc_1 = abs(l_1(t,:)).*W_1;
                Visc_2 = abs(l_2(t,:)).*W_2;

            % Roe flux
                for k = 1:N_x+1                    
                    F(:,k) = (1/2)*(f(:,k+1) + f(:,k)) - (1/2)*(Visc_1(:,k) + Visc_2(:,k));
                end
                
            % Roe scheme
                for k = 2:N_x+1
                    u(:,k) = u(:,k) - (Delta_t/Delta_x)*(F(:,k) - F(:,k-1));
                end

            % We stock the value for the different iteration (time
            % depending function U = (h, m=h*v))
                U_2 = [U_2; u];
        end
                
%% Display of the wave - drawnow

        % U(2*k + 2,2:end-1)./U(2*(k-1) + 1,2:end-1)    to get u
        % U(2*k + 1,2:end-1)                            to get h
        
        % Variables axis
            % Height
                x_lim_inf_h = x(1);   x_lim_sup_h = x(end);
                y_lim_inf_h = 0;      y_lim_sup_h = 2;
                
            % Speed
                x_lim_inf_s = x(1);   x_lim_sup_s = x(end);
                y_lim_inf_s = -2;      y_lim_sup_s = 2;     
                
        % Height
%             figure(2)
%             for k = 1:N_t/2
%                 plot(x, U_2(2*k + 1,2:end-1), 'r');
%                 axis([x_lim_inf_h x_lim_sup_h y_lim_inf_h y_lim_sup_h]);
%                 xlabel('Position on the mesh [m]'); ylabel('Wave''s height [m]'); title(['t = ',num2str(round((T/N_t)*k,2)),' s']);
%                 drawnow
%             end        

        % Speed
%             figure(3)
%             for k = 1:N_t/2
%                 plot(x, U_2(2*k + 2,2:end-1)./U_2(2*k + 1,2:end-1), 'r');
%                 axis([x_lim_inf_s x_lim_sup_s y_lim_inf_s y_lim_sup_s]);
%                 xlabel('Position on the mesh [m]'); ylabel('Wave''s speed [m.s^{-1}]'); title(['t = ',num2str(round((T/N_t)*k,2)),' s']);
%                 drawnow
%             end
 
%% Presentation results

        % Variables axis
            x_lim_inf_h = x(1);     x_lim_sup_h = x(end);
            y_lim_inf_h = 0.9;      y_lim_sup_h = 1.3;
            
        figure(4)
            subplot(1,3,1)
            k = round(N_t/20);
            plot(x, U_2(2*k + 1, 2:end-1), 'r');
            axis([x_lim_inf_h x_lim_sup_h y_lim_inf_h y_lim_sup_h]);
            ylabel('Wave''s height [m]'); title(['t = ',num2str(round(k*T/N_t,2)),' s']);

            subplot(1,3,2)
            k = round(N_t/15);
            plot(x, U_2(2*k + 1, 2:end-1), 'r');
            axis([x_lim_inf_h x_lim_sup_h y_lim_inf_h y_lim_sup_h]);
            xlabel('Position on the mesh [m]'); title(['t = ',num2str(round(k*T/N_t,2)),' s']);

            subplot(1,3,3)
            k = round(N_t/7);
            plot(x, U_2(2*k + 1, 2:end-1), 'r');
            axis([x_lim_inf_h x_lim_sup_h y_lim_inf_h y_lim_sup_h]);
            title(['t = ',num2str(round(k*T/N_t,2)),' s']);
