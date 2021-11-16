%% 2.1.3 - Lax-Friedrich scheme

    % Problem's variables
        L = 10;
        T = 100;
        
        H = 1;
        g = 9.81;
        w = 0.1*L;
        epsilon = H/5;

    % Mesh
        n = 80;
        iter_t = 5000; % 620; 1240; 2480          
        
        x = linspace(0, L, n);
            
        Delta_x = L/((n-1));
        Delta_t = T/iter_t;
            
        alpha = Delta_x/Delta_t;
            
    % Vector containing u = (h,hv)^T
        u = zeros(2,n+2);
        
% ------------------------------------------------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------------------------------------------------ %
        
    % Conditions : initial and boundary
        
        % Height
            for k = 1:n
               u(1,k+1) = H + epsilon*exp(-((x(k) - L/2)^2)/w^2); 
            end
        
        % Moment
            for k = 1:n
                u(2,k+1) = (u(1,k+1) - H)*sqrt(g*u(1,k+1));     
            end
                 
            U = zeros(2*(iter_t+1),n+2);
            U = [u];

        for t = 1:iter_t
            
            % Reflecting boundary for height + moment
                u(1,1) = u(1,2);    u(1,n+2) = u(1,n+1);
                
                u(2,1) = -u(2,2);   u(2,n+2) = -u(2,n+1);            
                
                % Function f to work on F (Lax-Friedrichs flux) :
                    f = [u(2,:) ; (u(2,:).^2)./u(1,:) + (1/2)*g*(u(1,:).^2)];
                
                % Lax-Friedrichs flux
                    for k = 1:n+1
                        F(:,k) = (1/2)*(f(:,k+1) + f(:,k) - alpha*(u(:,k+1) - u(:,k)));
                    end
                
                % Lax-Friedrichs scheme
                    for k = 2:n+1
                        u(:,k) = u(:,k) - (Delta_t/Delta_x)*(F(:,k) - F(:,k-1));
                    end
                
            U = [U; u];
        end
        
%% Presentation results

    % Height for bounce
        figure(11)
            subplot(2,2,1)
            k = 1;
            plot(x, U(2*k + 1, 2:end-1), 'r');
            axis([0 10 1 1.5]);
            xlabel('position on the mesh [m]'); ylabel('Wave ''s height (m)'); title(['t = ',num2str(round(k*T/iter_t,2)),' s']);

            subplot(2,2,2)
            k = 50;
            plot(x, U(2*k + 1, 2:end-1), 'r');
            axis([0 10 1 1.5]);
            xlabel('position on the mesh [m]'); ylabel('Wave ''s height (m)'); title(['t = ',num2str(round(k*T/iter_t,2)),' s']);

            subplot(2,2,3)
            k = 500;
            plot(x, U(2*k + 1, 2:end-1), 'r');
            axis([0 10 1 1.5]);
            xlabel('position on the mesh [m]'); ylabel('Wave ''s height (m)'); title(['t = ',num2str(round(k*T/iter_t,2)),' s']);

            subplot(2,2,4)
            k = iter_t;
            plot(x, U(2*k + 1, 2:end-1), 'r');
            axis([0 10 1 1.5]);
            xlabel('position on the mesh [m]'); ylabel('Wave ''s height (m)'); title(['t = ',num2str(round(k*T/iter_t,2)),' s']);

