%% We measure in this code the consistency of the considered point in the Q_2_2_1 over Delta_x

        % Vector containing u = (h,hv)^T
            u_5 = zeros(2,n+2); % We consider h the upper vector and hv the lower ones

    % Conditions : initial and boundary
            Delta_t_5 = Delta_t;
            Delta_x_5 = linspace(0.11,1.1,10);
            
            alpha_5 = alpha;
            
        % Speed
            % Initial condition
            for k = 1:n+2   
                u_5(2,k) = 0;     % Initial speed (or initial moment m = h*u = 0)
            end
     
        % Height
            % Initial condition
            % We keep h(1) and h(n+2) at zero and the rest considering the
            % initial height
            for k = 1:n
               u_5(1,k+1) = H + epsilon*exp(-((x(k) - L/2)^2)/w^2); % 
            end
  
%             figure(18)     % Graphic of the initial height 
%             plot(x, u_4(1,2:n+1), 'r'); % We represent only the point of h for the considered initial value (size(x,2) points)
%             xlabel('Position x (m)'); ylabel('Initial height (m)'); title('Graphic representing the function h = f(x) at t=0s');
%             
            
%%            
        
        % Vector to keep in memory the iteration of the height and speed
        % for the different Delta_x_5
            U_5 = zeros(2* (size(Delta_t_5,2)*(iter_t+1) +1),n+2);
            U_5 = [u_5]; 
            
            f = zeros(2,n+2);
            FF_5 = zeros(2*( size(Delta_t_5,2)*(iter_t+1) +1),n+2);
            FF_5 = [f];
            
            tau_5 = (zeros(size(Delta_t_5,2),size(t_n,2)));
            
for i = 1:size(Delta_x_5,2)    
        for t = t_n
            
            % Reflecting boundary for the height + boundary for the moment
                u_5(1,1) = u_5(1,2);    u_5(1,n+2) = u_5(1,n+1);
                
                u_5(2,1) = u_5(2,2);    u_5(2,n+2) = u_5(2,n+1);            
                
            % Variables from general equation's form 
                % Function f to work on F (Lax-Friedrichs flux) :
                    f = [u_5(2,:) ; (u_5(2,:).^2)./u_5(1,:) + (1/2)*g*(u_5(1,:).^2)];
            % We calulate the speed and height
                % What is at u(:,1) and u(:,n+2) is already defined : 
                %   - v is always zero at u(:,1) and u(:,n+2)
                %   - h(1) and h(n+2) as well as m(1) and m(n+2) (m = h*v) are defined above at the begining of the loop
                
                % Lax-Friedrichs flux
                for k = 1:n+1
                    F1(:,k) = (1/2)*(f(:,k+1) + f(:,k) - alpha*(u_5(:,k+1) - u_5(:,k)));
                end
                
                % Resolution of the system to get h and h*v
                for k = 2:n+1
                    u_5(:,k) = u_5(:,k) - (Delta_t_5/Delta_x_5(i))*(F1(:,k) - F1(:,k-1));
                end
                
            % We stock the value for the different iteration (time
            % depending function U = (h,h*v))
                 U_5 = [U_5; u_5];
                 FF_5 = [FF_5; f];
        end
end


        
%%        
        % Calculation tau for x = n/2
        for t = 1:size(Delta_x_5,2)
            for i = 1:iter_t-1
                tau_5(t,i) = (U_5(2*(i-1)+3 + 2*iter_t*(t-1), n/2) - U_5(2*(i-1)+1 + 2*iter_t*(t-1), n/2))/Delta_t_5 + ...
                    (FF_5(2*(i-1)+1 + 2*iter_t*(t-1), n/2+1) - FF_5(2*(i-1)+1 + 2*iter_t*(t-1), n/2-1) - alpha_5*( U_5(2*(i-1)+1 + 2*iter_t*(t-1), n/2+1) ...
                    - 2*U_5(2*(i-1)+1 + 2*iter_t*(t-1), n/2) + U_5(2*(i-1)+1 + 2*iter_t*(t-1), n/2-1) ))/(Delta_x_5(t)*2);
            end
        end
        
        % We here apply the L_1 norm to tau
        for t = 1:size(Delta_x_5,2)
            Sum = 0;
            for j = 1:iter_t-1
                Sum = Sum + abs(tau_5(t,j));
            end
            tau_norm_space(1,t) = Sum*Delta_x_5(t);
        end
        
        %%
        
        % We display the value of norm L_1 of tau depending on the value of
        % Delta_x_5
        % We thus get the type of dependence and a coeff for it
        figure(19)
        for t = 1:size(Delta_x_5,2)-2
            plot(Delta_x_5, (tau_norm_space), 'r*');
            xlabel('Delta_x'); ylabel('Norm L_1 of the Truncation error'); title('Graph representing tau = a*Delta_x + b');
        end
        
        
        
        