%% 2.2.3 - Lax-Friedrich scheme

    % Problem's variables
        L = 10;
        T = 20;
        
        H = 1;
        g = 9.81;
        w = 0.1*L;
        epsilon = H/5;
        
    % Bottom
        B0  = H/10;
        r   = L/6;
        
    % Mesh
        n = 80;        
        iter_t = 3000;           
        
        x = linspace(0, L, n);
            
        Delta_x = L/((n-1));
        Delta_t = T/iter_t;
            
        alpha = Delta_x/Delta_t;
            
    % Vector containing u = (h,hv)^T
        u = zeros(2, n+2);
        
% ------------------------------------------------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------------------------------------------------ %
        
    % Conditions : initial and boundary
        
        % Bottom
            for k = 1:n
                if abs(x(k) - L/2) < r
                    B(k) = B0*cos((pi*(x(k) - L/2))/(2*r))^2; 
                else
                    B(k) = 0;
                end
            end    
        
        % Type 3
            ii = 1;
            jj = 2.1;
            
        % Type 4
%             ii = 1.5;
%             jj = 5.5;
        
        % Height
            for k = 1:n
               u(1,k+1) = ii; 
            end
        
        % Moment
            for k = 1:n
                u(2,k+1) = jj;
            end
            
            U = zeros(2*(iter_t+1),n+2);
            U = [u];

            Inflow_condition = jj;
            
        for t = 1:iter_t
            % BC
                % Height
                    u(1,1)   = 2*u(1,2) - u(1,3);   % Extrapolation for ghost cell (j = -1)  
                    u(1,n+2) = ii;                 % Dirichlet BC for ghost cell  (j = N+2)

                % Moment
                    u(2,1) = Inflow_condition;   
                    u(2,n+2) = 2*u(2,n+1) - u(2,n);            

            % Function f to work on F (Lax-Friedrichs flux) :
                f = [u(2,:) ; (u(2,:).^2)./u(1,:) + (1/2)*g*(u(1,:).^2)];
                
            % Source term
                for k = 1:n
                    if abs(x(k) - L/2) < r
                        S(2,k) = (g*u(1,k)*B0*pi/(2*r))*sin((pi*(x(k) - L/2))/r);
                    else
                        S(2,k) = 0;
                    end
                end
              
            % Lax-Friedrichs flux
                for k = 1:n+1
                    F(:,k) = (1/2)*(f(:,k+1) + f(:,k) - alpha*(u(:,k+1) - u(:,k)));
                end
                
            % Resolution of the system to get h and h*v
                for k = 2:n+1
                    u(:,k) = u(:,k) - (Delta_t/Delta_x)*(F(:,k) - F(:,k-1)) + Delta_t*S(:,k-1);
                end
                
            U = [U; u];
        end
        
%% STABILITY

        % Variables axis
            x_lim_inf_h = x(1);   x_lim_sup_h = x(end);
            y_lim_inf_h = 1;      y_lim_sup_h = 1.5; 

%         figure(120)
%             subplot(2,2,1)
%             k = 1;
%             plot(x, U(2*k + 1, 2:end-1), 'r');
%             axis([x_lim_inf_h x_lim_sup_h y_lim_inf_h y_lim_sup_h]);
%             xlabel('position on the mesh [m]'); ylabel('Wave ''s height (m)'); title(['t = ',num2str(round(k*T/iter_t,2)),' s']);
% 
%             subplot(2,2,2)
%             k = 50;
%             plot(x, U(2*k + 1, 2:end-1), 'r');
%             axis([x_lim_inf_h x_lim_sup_h y_lim_inf_h y_lim_sup_h]);
%             xlabel('position on the mesh [m]'); ylabel('Wave ''s height (m)'); title(['t = ',num2str(round(k*T/iter_t,2)),' s']);
% 
%             subplot(2,2,3)
%             k = iter_t;
%             plot(x, U(2*k + 1, 2:end-1), 'r');
%             axis([x_lim_inf_h x_lim_sup_h y_lim_inf_h y_lim_sup_h]);
%             xlabel('position on the mesh [m]'); ylabel('Wave ''s height (m)'); title(['t = ',num2str(round(k*T/iter_t,2)),' s']);
% 
%             subplot(2,2,4)
%             k = 1;
%             plot(x, U(2*k + 1, 2:end-1), 'r');
%             axis([x_lim_inf_h x_lim_sup_h y_lim_inf_h y_lim_sup_h]);
%             xlabel('position on the mesh [m]'); ylabel('Wave ''s height (m)'); title('t = 150*Delta_t s = 3 bouncings');

%% Presentation results

    figure(121)
        k = iter_t;
        p_121 = plot(x, B, 'r--', x, U(2*k + 1, 2:end-1), 'b--', x, U(2*k + 2,2:end-1)./U(2*k + 1,2:end-1), 'g-');
        axis([0 10 -0.2 3]);
        xlabel('Position on the mesh [m]'); title(['t = ',num2str(round(k*T/iter_t,2)),' s']);
        grid on;
        legend([p_121], {'B', 'h + B', 'u'}, 'location', 'northwest');
  
