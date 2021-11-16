%% 2.1 - Numerical Solution - continuation

% Q 2.1 - 3 : Influence of alpha
    % During this part we consider :
        
        epsilon_3 = 0.1;
        
        C = linspace(0.5,1.5,10);
        alpha_3 = (Delta_x/Delta_t)*C;    
        
        x = linspace(0, L, n);  % Vector x to be used for initial wave's height

        % U_3 stocks the variable h and m
            U_3 = zeros(size(u,1)*size(t_n,2),size(u,2));
        % U_alpha stocks U_3 for different alpha
            U_alpha = zeros(size(U_3,1)*size(alpha_3,2),size(U_3,2));
            
% Q 2.1 - 2

for i = 1:size(alpha_3,2)
    
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
               u(1,k+1) = H + epsilon_3*exp(-((x(k) - L/2)^2)/w^2); 
            end
                    
        % Vector to keep in memory the iteration of the height and speed
            U_3 = [u]; 
            
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
                    F1(:,k) = (1/2)*(f(:,k+1) + f(:,k) - alpha_3(i)*(u(:,k+1) - u(:,k)));
                end
                
                for k = 2:n+1
                    u(:,k) = u(:,k) - (Delta_t/Delta_x)*(F1(:,k) - F1(:,k-1));
                end

            % We stock the value for the different iteration (time
            % depending function u)
                 
                 U_3 = [U_3; u];
        end
        
        % We now stock u for different epsilon
        U_alpha = [U_alpha; U_3];
end  

    % Re-sizing of U_alpha
        % Has we couldn't as usual declare the first value of EE, we have
        % to remove the extra data (that are basically equal to zero)
        
            U_alpha = U_alpha((size(U_3,1)-2)*size(alpha_3,2) + 1:end,:);
            
% ------------------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------------------ %
         
        
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
            maximum = zeros(size(alpha_3,2),iter_t+1);
        
        % Corresponding position (note : 2 maximums as the functions is
        % symmetrical)
            x_alpha = zeros(2*size(alpha_3,2),iter_t+1);
            x_max_3 = zeros(2,iter_t+1);

            U_alpha(isnan(U_alpha)) = 0;
            U_alpha(isinf(U_alpha)) = 0;
            
%         % We calculate these maximum and their associated position
        for i = 1:size(alpha_3,2)
           for k = 1:iter_t+1
              maximum(i,k) = max(U_alpha(2*k-1 + size(U_3,1)*(i-1),2:end-1)); 
              if U_alpha(2*k-1 + size(U_3,1)*(i-1),2:end-1) ~= 0
                  x_max_3(:,k) = x(U_alpha(2*k-1 + size(U_3,1)*(i-1),2:end-1) == max(U_alpha(2*k-1 + size(U_3,1)*(i-1),2:end-1)));
              else
                  x_max_3(:,k) = [0,0];
              end
           end
           x_alpha = [x_alpha; x_max_3];
        end

        % Re-sizing of x_alpha
            x_alpha = x_alpha(size(x_alpha)/2+1:end,:);

            
%% Display of the shape in function of the crossed distance
    
    i = 0; % Variable to know when we completed one loop (= crossing 10 meters)
    
    x_length_3 = zeros(size(alpha_3,2),size(x_alpha,2)+1);     % Distance crossed
    x_length_3(:,1) = 0;
    
    for k = 1:size(x_alpha,2)
        for a = 1:size(alpha_3,2)
            if x_alpha(2,k) == 10
                i = i + 1;      % We are at x = 10*i, we continue from this value
            end
            x_length_3(a,k+1) = x_length_3(a,k) + x_alpha(2*a,k)+10*i - 5;
        end
    end
    
    % NOTE : we only worked on x > 5 as we are on a symmetrical mesh.
    % Moreover, the "-5" term comes from the fact that we offset our value
    % so we begin at 0.
    
        x_length_3 = x_length_3(:,2:end);   % Re-sizing
    
    % We now display the new curves that are not restricted by the size of
    % the mesh anymore
    figure(12)
    for k = 1:size(maximum,1)
        rr(k) = plot(x_length_3(k,:), maximum(k,:), '-');
        axis([0 10*10^5 1 1.04]);
        hold on
        xlabel('Distance theoretically crossed by the function (m)'); ylabel('Maximal Height of the wave (m)'); title('Behavior of the wave depending on alpha');
    end
    legend([rr], {'a = 5.0', 'a = 6.1', 'a = 7.2', 'a = 8.3', ' = 9.4', 'a = 10.5', 'a = 11.6', 'a = 12.7', 'a = 13.8', 'a = 15.0'}, 'location', 'northeast');

%%
    % Display of an instable solution for alpha = 10.5 > Delta_x/Delta_t
        figure(13)
        subplot(2,3,1)
        plot(x, U_alpha(2*(30-1)+1 + size(U_3,1)*(6-1),2:end-1), 'r');
        axis([0 10 1 1.15]);
        xlabel('position on the mesh [0;L](m)'); ylabel('Wave height (m)'); title('Wave at t = 30*Delta_t');
        
        subplot(2,3,2)
        plot(x, U_alpha(2*(200-1)+1 + size(U_3,1)*(6-1),2:end-1), 'r');
        axis([0 10 1 1.15]);
        xlabel('position on the mesh [0;L](m)'); ylabel('Wave height (m)'); title('Wave at t = 200*Delta_t');
        
        subplot(2,3,3)
        plot(x, U_alpha(2*(360-1)+1 + size(U_3,1)*(6-1),2:end-1), 'r');
        axis([0 10 1 1.15]);
        xlabel('position on the mesh [0;L](m)'); ylabel('Wave height (m)'); title('Wave at t = 360*Delta_t');
        
        subplot(2,3,4)
        plot(x, U_alpha(2*(375-1)+1 + size(U_3,1)*(6-1),2:end-1), 'r');
        axis([0 10 1 1.15]);
        xlabel('position on the mesh [0;L](m)'); ylabel('Wave height (m)'); title('Wave at t = 375*Delta_t');
        
        subplot(2,3,5)
        plot(x, U_alpha(2*(383-1)+1 + size(U_3,1)*(6-1),2:end-1), 'r');
        axis([0 10 1 1.15]);
        xlabel('position on the mesh [0;L](m)'); ylabel('Wave height (m)'); title('Wave at t = 383*Delta_t');
    
        subplot(2,3,6)
        plot(x, U_alpha(2*(388-1)+1 + size(U_3,1)*(6-1),2:end-1), 'r');
        axis([0 10 1 1.15]);
        xlabel('position on the mesh [0;L](m)'); ylabel('Wave height (m)'); title('Wave at t = 388*Delta_t');
        
        
% Display of the graphs showing the mean speed over one loop

% We mainly pay attention to the fact that the shape are no more the same
% once we passed the range of stable alpha
%%
%     figure(14)
%     for i = 1:size(alpha_3,2)
%         subplot(4,3,i)
%         for k = 1:iter_t+1
%             plot(x, U_alpha(2*k + size(U_3,1)*(i-1),2:end-1)/U_alpha(2*(k-1) + 1 + size(U_3,1)*(i-1),2:end-1), 'color', k*[1/305 1/305 1/600]);
%             hold on
%         end
% %         xlabel('position on the mesh [0;L](m)'); ylabel('Wave speed (m.s-1)'); title(['Wave at alpha = ',num2str(l),'']);
%         hold off
%     end
