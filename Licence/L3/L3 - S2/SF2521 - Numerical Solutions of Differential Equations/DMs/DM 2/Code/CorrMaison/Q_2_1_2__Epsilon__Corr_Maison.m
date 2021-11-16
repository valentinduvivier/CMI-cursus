%% 2.1 - Numerical Solution - continuation

% Q 2.1 - 2 : Influence of epslion on height and speed
    % During this part we will consider consider :        
        for k = 1:10    % number of epsilon you want to test
            epsilon_2(k) = k*0.01;     % 0.4 = epsilon step
        end
            
        % U still stocks the variable h and m
            U_2 = zeros(size(u,1)*size(N_t,2),size(u,2));
            
        % U_epsilon stocks U for different epsilon
            U_epsilon = zeros(size(U_2,1)*size(epsilon_2,2),size(U_2,2));
        
for i = 1:size(epsilon_2,2)
    
    % Vector containing u = (h,hv)^T
        u = zeros(2, N+2); % We consider h the upper vector and hv the lower ones
                
        
    % Conditions : initial and boundary
        
        % Speed
            % Initial condition
            for k = 1:N+2   
                u(2,k) = 0;     % Initial speed (or initial moment m = h*u = 0)
            end
          
        % Height
            for k = 1:N
               u(1,k+1) = H + epsilon_2(i)*exp(-((x(k) - L/2)^2)/w^2); 
            end
        
        % Vector to keep in memory the iteration of the height and speed
            U_2 = [u]; 
            
        for t = 1:N_t     % We iterate iter_t times
            
            % Reflecting boundary for the height
                u(1,1) = u(1,2);    u(1,N+2) = u(1,N+1);
                
                u(2,1) = -u(2,2);   u(2,N+2) = -u(2,N+1); 

            % Variables from general equation's form 
                % Functions to work on F (Lax-Friedrichs flux) :
                    f = [u(2,:) ; (u(2,:).^2)./u(1,:) + (1/2)*g*(u(1,:).^2)];
                
            % We calulate the speed and height
                % What is at u(:,1) and u(:,n+2) is already defined : 
                %   - v is always zero at u(:,1) and u(:,n+2)
                %   - h(1) and h(n+2) are defined above at the begining of the loop
                
                for k = 1:N+1
                    F1(:,k) = (1/2)*(f(:,k+1) + f(:,k) - alpha*(u(:,k+1) - u(:,k)));
                end
                
                for k = 2:N+1
                    u(:,k) = u(:,k) - (Dt/Dx)*(F1(:,k) - F1(:,k-1));
                end

            % We stock the value for the different iteration (time
            % depending function u)
                 
                U_2 = [U_2; u];
        end
        
        % We now stock u for different epsilon
            U_epsilon = [U_epsilon; U_2];
end      
 
        U_epsilon = U_epsilon(size(u,1)*size(epsilon_2,2) + 1:end,:);

% ------------------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------------------ %
           
        % Display of the initial height 
            
            t1 = 451;       % Up to 601 (only odd number as it's height)
            t2 = 452;       % Up to 602 (only even number as it's height)
            t3 = t2 - 1;    % Time before t2 which represents the moment (to get the speed)
            
            figure(6)
            for i = 1:size(epsilon_2,2)
                %Height%
                p(i) = plot(x, U_epsilon(2*(N_t+1)*(i-1) + t1, 2:N+1), 'color', i*[0.1,0.1,0]); % We represent only the point of h for the considered initial value (size(x,2) points)
                hold on
                grid on
                xlabel('Position x (m)'); ylabel('Initial height (m)'); title(['Graphic representing the function h = f(x) at t = ',num2str(((t1-1)/2)*Dt),' s']); % 401 == 200 iterations so 2/3 * T = 2s.
            end
            hold off
            legend([p], {'e = 0.4', 'e = 0.8', 'e = 1.2', 'e = 1.6', 'e = 2.0', 'e = 2.4', 'e = 2.8', 'e = 3.2', 'e = 3.6', 'e = 4.0'}, 'location', 'northwest');

            figure(7)
            for i = 1:size(epsilon_2,2)
                %Speed%
                q(i) = plot(x, U_epsilon(2*(N_t+1)*(i-1) + t2, 2:N+1)./U_epsilon(2*(N_t+1)*(i-1) + t3, 2:N+1),'color', i*[0.1,0.1,0]); % We represent only the point of h for the considered initial value (size(x,2) points)
                hold on
                grid on
                xlabel('Position x (m)'); ylabel('Speed (m)'); title(['Graphic representing the function h = f(x) at t = ',num2str((t2/2)*Dt),' s']);
            end
            hold off
            legend([q], {'e = 0.4', 'e = 0.8', 'e = 1.2', 'e = 1.6', 'e = 2.0', 'e = 2.4', 'e = 2.8', 'e = 3.2', 'e = 3.6', 'e = 4.0'}, 'location', 'northwest');

            
% Q 2.1 - 2 : Influence of epsilon on wave collisions
%%
    figure(8)
    
        % Time to consider for the waves' collision
            t1 = 431; t2 = 601; t3 = 801; t4 = 1001;
    
        subplot(2,2,1)
            for i = 1:size(epsilon_2,2)
                %Height%
                p(i) = plot(x, U_epsilon(2*(N_t+1)*(i-1) + t1, 2:N+1), 'color', i*[0.1,0.1,0]); % We represent only the point of h for the considered initial value (size(x,2) points)
                hold on
                grid on
                xlabel('Position x (m)'); ylabel('Initial height (m)'); title(['Graphic representing the function h = f(x) at t = ',num2str(((t1-1)/2)*Dt),' s']); % 401 == 200 iterations so 2/3 * T = 2s.
            end
            hold off
            legend([p], {'e = 0.4', 'e = 0.8', 'e = 1.2', 'e = 1.6', 'e = 2.0', 'e = 2.4', 'e = 2.8', 'e = 3.2', 'e = 3.6', 'e = 4.0'}, 'location', 'northwest');

        subplot(2,2,2)
            for i = 1:size(epsilon_2,2)
                %Height%
                p(i) = plot(x, U_epsilon(2*(N_t+1)*(i-1) + t2, 2:N+1), 'color', i*[0.1,0.1,0]); % We represent only the point of h for the considered initial value (size(x,2) points)
                hold on
                grid on
                xlabel('Position x (m)'); ylabel('Initial height (m)'); title(['Graphic representing the function h = f(x) at t = ',num2str(((t2-1)/2)*Dt),' s']); % 401 == 200 iterations so 2/3 * T = 2s.
            end
            hold off
            legend([p], {'e = 0.4', 'e = 0.8', 'e = 1.2', 'e = 1.6', 'e = 2.0', 'e = 2.4', 'e = 2.8', 'e = 3.2', 'e = 3.6', 'e = 4.0'}, 'location', 'northwest');

        subplot(2,2,3)
            for i = 1:size(epsilon_2,2)
                %Height%
                p(i) = plot(x, U_epsilon(2*(N_t+1)*(i-1) + t3, 2:N+1), 'color', i*[0.1,0.1,0]); % We represent only the point of h for the considered initial value (size(x,2) points)
                hold on
                grid on
                xlabel('Position x (m)'); ylabel('Initial height (m)'); title(['Graphic representing the function h = f(x) at t = ',num2str(((t3-1)/2)*Dt),' s']); % 401 == 200 iterations so 2/3 * T = 2s.
            end
            hold off
            legend([p], {'e = 0.4', 'e = 0.8', 'e = 1.2', 'e = 1.6', 'e = 2.0', 'e = 2.4', 'e = 2.8', 'e = 3.2', 'e = 3.6', 'e = 4.0'}, 'location', 'northwest');

        subplot(2,2,4)
            for i = 1:size(epsilon_2,2)
                %Height%
                p(i) = plot(x, U_epsilon(2*(N_t+1)*(i-1) + t4, 2:N+1), 'color', i*[0.1,0.1,0]); % We represent only the point of h for the considered initial value (size(x,2) points)
                hold on
                grid on
                xlabel('Position x (m)'); ylabel('Initial height (m)'); title(['Graphic representing the function h = f(x) at t = ',num2str(((t4-1)/2)*Dt),' s']); % 401 == 200 iterations so 2/3 * T = 2s.
            end
            hold off
            legend([p], {'e = 0.4', 'e = 0.8', 'e = 1.2', 'e = 1.6', 'e = 2.0', 'e = 2.4', 'e = 2.8', 'e = 3.2', 'e = 3.6', 'e = 4.0'}, 'location', 'northwest');

            
            % ---------------------------------------------------------------------- %
            % As a conclusion, we can easily see that amplitude is widely
            % influenced by epsilon. Thereby, bigger epsilon is, bigger
            % the apmlitude gets. 
            % Moreover, the epsilon influences the width of the initial
            % wall which, together with the defined mesh, makes that the
            % wave bounds and fall sometimes as a block. 
            % To finish, we see that different epsilons bring different
            % inertias and behaviors for the waves, and so they collide at 
            % different moment and with different force.
            % FFor the speed, as mentioned the inertia isn't the same when
            % bouncing on the walls or at the middle, etc, so not the same
            % speed.
            % ---------------------------------------------------------------------- %

            