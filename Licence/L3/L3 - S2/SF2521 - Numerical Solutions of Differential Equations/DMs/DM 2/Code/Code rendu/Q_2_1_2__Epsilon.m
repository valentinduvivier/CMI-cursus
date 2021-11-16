%% 2.1 - Numerical Solution - continuation

% Q 2.1 - 2 : Influence of epslion
    % During this part we consider :
        E = 10;     % number of epsilon you want to test
        
        for k = 1:E    
            epsilon_2(k) = k*0.4;     % 0.4 = epsilon step
        end
            
        % U still stocks the variable h and m
        U_2 = zeros(size(u,1)*size(t_n,2),size(u,2));
        % U_epsilon stocks U for different epsilon
        U_epsilon = zeros(size(U_2,1)*size(epsilon_2,2),size(U_2,2));
        
% Q 2.1 - 2

for i = 1:size(epsilon_2,2)
    
    % Vector containing u = (h,hv)^T
        u = zeros(2, n+2); % We consider h the upper vector and hv the lower ones
                
    % Conditions : initial and boundary
        
        % Speed
            % Initial condition
            for k = 1:n+2   
                u(2,k) = 0;     % Initial speed (or initial moment m = h*u = 0)
            end
          
        % Height
            for k = 1:n
               u(1,k+1) = H + epsilon_2(i)*exp(-((x(k) - L/2)^2)/w^2); 
            end
        
        % Vector to keep in memory the iteration of the height and speed
            U_2 = [u]; 
            
        for t = t_n     % We iterate iter_t times
            
            % Reflecting boundary for the height
                u(1,1) = u(1,2);    u(1,n+2) = u(1,n+1);
                
                u(2,1) = -u(2,2);   u(2,n+2) = -u(2,n+1); 

            
            % Variables from general equation's form 
                % Functions to work on F (Lax-Friedrichs flux) :
                    f = [u(2,:) ; (u(2,:).^2)./u(1,:) + (1/2)*g*(u(1,:).^2)];
                
            % We calulate the speed and height
                % What is at u(:,1) and u(:,n+2) is already defined : 
                %   - v is always zero at u(:,1) and u(:,n+2)
                %   - h(1) and h(n+2) are defined above at the begining of the loop
                
                for k = 1:n+1
                    F1(:,k) = (1/2)*(f(:,k+1) + f(:,k) - alpha*(u(:,k+1) - u(:,k)));
                end
                
                for k = 2:n+1
                    u(:,k) = u(:,k) - (Delta_t/Delta_x)*(F1(:,k) - F1(:,k-1));
                end

            % We stock the value for the different iteration (time
            % depending function u)
                 
                U_2 = [U_2; u];
        end
        
        % We now stock u for different epsilon
        U_epsilon = [U_epsilon; U_2];
end      
 
    % Re-sizing of U_epsilon
        % Has we couldn't as usual declare the first value of U_epsilon, we have
        % to remove the extra data (that are basically equal to zero)

            U_epsilon = U_epsilon((size(U_2,1)-2)*size(epsilon_2,2) + 1:end,:);
         
% ------------------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------------------ %
 
        % Initial position depending on epsilon 
          
        % Display of the initial height 
        figure(6)
        for i = 1:size(epsilon_2,2)
            p(i) = plot(x, U_epsilon(19 + 2*iter_t*(i-1),2:n+1)); % We represent only the point of h for the considered initial value (size(x,2) points)
            hold on
            xlabel('Position x (m)'); ylabel('Initial height (m)'); title('Graphic representing the function h = f(x)');
        end
        hold off
        legend([p], {'e = 0.4', 'e = 0.8', 'e = 1.2', 'e = 1.6', 'e = 2.0', 'e = 2.4', 'e = 2.8', 'e = 3.2', 'e = 3.6', 'e = 4.0'}, 'location', 'northwest');
           
     
%% Influence on the wave's shape with respect to epsilon
        % We calculate the maximum taken by the wave, and thus for every
        % point(which depends on Delta_t, etc) and all the epsilon
        
        % The purpose is to compare the position in x of the maximums and
        % to see if they are contant over the epsilon : 
        %   - if we have the same position in x (not the amplitude), then epsilon as no impact on
        %   the shape ;
        %   - if we don't : wve's shape is influenced by epsilon
        
        
        % Variable gathering the maximums for all the epsilon (each line = a
        % different epsilon)
            maximum = zeros(size(epsilon_2,2),iter_t+1);
        
        % Corresponding position (note : 2 maximums as the functions is
        % symmetrical)
            x_epsilon = zeros(2*size(epsilon_2,2),iter_t+1);
            x_max_2 = zeros(2,iter_t+1);
            
        % We calculate these maximum and their associated position
        for i = 1:size(epsilon_2,2)
           for k = 1:iter_t+1
              maximum(i,k) = max(U_epsilon(2*k-1 + size(U_2,1)*(i-1),2:end-1)); 
              x_max_2(:,k) = x(U_epsilon(2*k-1 + size(U_2,1)*(i-1),2:end-1) == max(U_epsilon(2*k-1 + size(U_2,1)*(i-1),2:end-1)));
           end
           x_epsilon = [x_epsilon; x_max_2];
        end
        
% Re-sizing of x_epsilon
        x_epsilon = x_epsilon(size(x_epsilon)/2+1:end,:);

%% Approach 1 - Display of the position in one graph
        figure(7)
        for k = 1:size(maximum,1)
           plot(x_epsilon(2*(k-1) + 1,:), maximum(k,:), '*', x_epsilon(2*k,:), maximum(k,:), '*', 'color', k*[0.05 0.1 0.04]);
           hold on
           xlabel('Position x [0,L]'); ylabel('maximum of the function'); title('Graphic of the maximum in function of the position x');
        end
        hold off

% Conclusion --> the way of displaying isn't consistent. We decide an other
% approach

%% Approach 2 - Display of the shape in funtion of the crossed distance
    
    i = 0; % Variable to know when we did one loop
    
    x_length2 = zeros(size(epsilon_2,2),size(x_epsilon,2)+1);     % Distance crossed
    x_length2(:,1) = 0;
    
    for k = 1:size(x_epsilon,2)
        for e = 1:size(epsilon_2,2)
            if x_epsilon(2,k) == 10
                i = i + 1;      % We are at x = 10*i, we continue from this value
            end
            x_length2(e,k+1) = x_length2(e,k) + x_epsilon(2*e,k)+10*i - 5;
        end
    end
    
    % NOTE : we only worked on x > 5 as we are on a symmetrical mesh.
    % Moreover, the "-5" term comes from the fact that we offset our value
    % so we begin at 0.
    
        x_length2 = x_length2(:,2:end);   % Re-sizing
    
    % We now display the new curves that are not restricted by the size of
    % the mesh anymore
    figure(8)
        rr = plot(x_length2', maximum', '*');
        axis([0 10*10^5 1 2]);
        xlabel('Distance theoretically crossed by te function (m)'); ylabel('Maximal Height of the wave (m)'); title('Behavior of the wave depending on epsilon');
        legend([rr], {'e = 0.4', 'e = 0.8', 'e = 1.2', 'e = 1.6', 'e = 2.0', 'e = 2.4', 'e = 2.8', 'e = 3.2', 'e = 3.6', 'e = 4.0'}, 'location', 'northeast');

    
%% Display of the graphs showing the mean speed over one loop

% We here pay attention to the amplitude taken by the speed's graphs rather
% than their shape (as the shape are ~~ identical)
    figure(9)
        for i = 1:size(epsilon_2,2)
            subplot(4,3,i)
            for k = 1:iter_t+1
                plot(x, U_epsilon(2*k + size(U_2,1)*(i-1),2:end-1)./U_epsilon(2*(k-1) ...
                    + 1 + size(U_2,1)*(i-1),2:end-1), 'color', k*[1/305 1/305 1/600]);
                hold on
                xlabel('position on the mesh [0;L](m)'); ylabel('Wave speed (m.s-1)');
            end
            title(['Wave at epsilon = ',num2str(i*0.4),'']);
        end
        hold off

    