%% 2.1 - Numerical Solution - continuation

% Q 2.1 - 3 : Influence of alpha
    % During this part we consider :
        
        epsilon_3 = 0.1;
        
        C = linspace(0.5, 1.2, 10);
        alpha_3 = (Dx/Dt)*C;    
        
        x = linspace(0, L, N);  % Vector x to be used for initial wave's height

        % U_3 stocks the variable h and m
            U_3 = zeros(size(u,1)*size(N_t,2),size(u,2));
        
        % U_alpha stocks U_3 for different alpha
            U_alpha = zeros(size(U_3,1)*size(alpha_3,2),size(U_3,2));
            
% Q 2.1 - 2

for i = 1:size(alpha_3,2)
    
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
               u(1,k+1) = H + epsilon_3*exp(-((x(k) - L/2)^2)/w^2); 
            end
                    
        % Vector to keep in memory the iteration of the height and speed
            U_3 = [u]; 
            
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
                    F1(:,k) = (1/2)*(f(:,k+1) + f(:,k) - alpha_3(i)*(u(:,k+1) - u(:,k)));
                end
                
                for k = 2:N+1
                    u(:,k) = u(:,k) - (Dt/Dx)*(F1(:,k) - F1(:,k-1));
                end

            % We stock the value for the different iteration (time
            % depending function u)
                 
                 U_3 = [U_3; u];
        end
        
        % We now stock u for different epsilon
        U_alpha = [U_alpha; U_3];
end  

        U_alpha = U_alpha(size(u,1)*size(alpha_3,2) + 1:end,:);

%% ------------------------------------------------------------------------------------------ %
           
        % Display of the initial height 
            
            t1 = 251;       % Up to 601 (only odd number as it's height)
            t2 = 452;       % Up to 602 (only even number as it's height)
            t3 = t2 - 1;    % Time before t2 which represents the moment (to get the speed)
            
            figure(9)
            for i = 1:size(alpha_3,2)
                %Height%
                p(i) = plot(x, U_alpha(2*(N_t+1)*(i-1) + t1, 2:N+1), 'LineWidth', 1.2, 'color', i*[0.09,0.07,0.1]); % We represent only the point of h for the considered initial value (size(x,2) points)
                hold on
                grid on
                xlabel('Position x (m)'); ylabel('Initial height (m)'); title(['Graphic representing the function h = f(x) at t = ',num2str(((t1-1)/2)*Dt),' s']); % ex : 401 == 200 iterations so 2/3 * T = 2s.
            end
            
            p(i+1) = plot(x, U(t1,2:end-1), 'r');
            axis([0 10 1 1.5]);
            xlabel('Position on the mesh (m)'); ylabel('Wave height (m)'); title(['Wave at t = ',num2str(t1*Dt,2),' s']);
            hold off
            legend([p], {['\alpha = ',num2str(alpha_3(1)),''], ['\alpha = ',num2str(alpha_3(2)),''], ['\alpha = ',num2str(alpha_3(3)),''], ['\alpha = ',num2str(alpha_3(4)),''], ['\alpha = ',num2str(alpha_3(5)),''], ...
                ['\alpha = ',num2str(alpha_3(6)),''], ['\alpha = ',num2str(alpha_3(7)),''], ['\alpha = ',num2str(alpha_3(8)),''], ['\alpha = ',num2str(alpha_3(9)),''], ['\alpha = ',num2str(alpha_3(10)),''], ['\alpha = ',num2str(Dx/Dt,2),' = reference']}, 'location', 'north');
            
            % ---------------------------------------------------------------------- %
            % As a conclusion, for too big alpha, the Lax-Friedrich flux gets unstable
            % and so we get no solution for too big alpha. On the other
            % hand, for the one converging, we see that bigger the alpha,
            % higher the amplitude of u.
            % These unstabilities confirm the idea that to solve PDE
            % necessitates to look at conditions of stability. Here we have
            % one bounding the variables, alpha, Dx and Dt.
            % ---------------------------------------------------------------------- %

            