%% We consider BC such taht the flux can cross the left and right boundaries
   
    % Conditions : initial and boundary
        % Speed
            % Initial condition
            for k = 1:N+2   
                u_no_boundary(2, k) = 0;     % Initial speed (or initial moment m = h*u = 0)
            end
     
        % Height
            % Initial condition
            % We keep h(1) and h(n+2) at zero and the rest considering the
            % initial height
            for k = 1:N
               u_no_boundary(1,k+1) = H + epsilon*exp(-((x(k) - L/2)^2)/w^2); % 
            end
               
            
% ------------------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------------------ %
        
        % Vector to keep in memory the iteration of the height and speed
            U_no_boundary = [u_no_boundary]; 

        for t = 1:N_t
            
            % Reflecting boundary for the height + boundary for the moment
                u_no_boundary(1,1) = u_no_boundary(1,2);    u_no_boundary(1,N+2) = u_no_boundary(1,N+1);
                
                u_no_boundary(2,1) = u_no_boundary(2,2);    u_no_boundary(2,N+2) = u_no_boundary(2,N+1);            
                
            % Variables from general equation's form 
                % Function f to work on F (Lax-Friedrichs flux) :
                    f = [u_no_boundary(2,:) ; (u_no_boundary(2,:).^2)./u_no_boundary(1,:) + (1/2)*g*(u_no_boundary(1,:).^2)];
                
            % We calulate the speed and height
                % What is at u(:,1) and u(:,n+2) is already defined : 
                %   - v is always zero at u(:,1) and u(:,n+2)
                %   - h(1) and h(n+2) as well as m(1) and m(n+2) (m = h*v) are defined above at the begining of the loop
                
                % Lax-Friedrichs flux
                for k = 1:N+1
                    F1(:,k) = (1/2)*(f(:,k+1) + f(:,k) - alpha*(u_no_boundary(:,k+1) - u_no_boundary(:,k)));
                end
                
                % Resolution of the system to get h and h*v
                for k = 2:N+1
                    u_no_boundary(:,k) = u_no_boundary(:,k) - (Dt/Dx)*(F1(:,k) - F1(:,k-1));
                end

            % We stock the value for the different iteration (time
            % depending function U = (h,h*v))
                 U_no_boundary = [U_no_boundary; u_no_boundary];
        end
      
        
        figure(11)
        for k = 1:N_t+1
            plot(x, U_no_boundary(2*(k-1)+1,2:end-1), 'r');
            axis([0 10 0.95 1.1]);
            xlabel('Position on the mesh (m)'); ylabel('Wave''s height (m)'); title(['Wave at t = ',num2str((k-1)*Dt,2),' s']);
            drawnow
        end
        
% ------------------------------------------------------------------------------------------ %
% ------------------------------------------------------------------------------------------ %

%% Quantification of reflection
    
    % Working on the case of non reflecting boundary
        Sum1 = zeros(N_t,1);

        for t = 1:N_t

            for i = 2:N+1
               Sum1(t,1) = Sum1(t,1) + U_no_boundary(2*(t-1) + 1, i); 
            end
        end
        
    % Working on the case of reflecting boundary
        Sum2 = zeros(N_t,1);

        for t = 1:N_t

            for i = 2:N+1
               Sum2(t,1) = Sum2(t,1) + U(2*(t-1) + 1, i); 
            end
        end
      
    % To mesure how well the reflection er non reflection went :
    
        % Reflection case :
        
            Reflection_rate = abs(Sum2(end,1) - Sum2(1,1))/Sum2(1,1) * 100;
        
            % --> 10^-14 % of reflection ==> NOTHING
            
        % Non-reflecting case
            
            Non_reflection_rate = abs(Sum1(end,1) - 100)/100 * 100;    % 100 as we have a normal height of water of 1m + wave. At normal, Volume = 100 ==> Sum = 100 as reference
        
            % 7.019 * 10^-4 % of non reflection. A bigger time T might ameliorate
            % this result even though it implies longer calculations. To
            % get a better N_t (more precise time step) might help for the precision of the
            % results as well, but the negative aspect is the same.
            
    % Displaying height sum variation
    figure(12)
    
        for k = 1:N_t/10
           plot(Dt*k, Sum1(k,1), 'b*-');
           hold on
           grid on
        end
        xlabel('time (s)'); ylabel('Sum of height'); title('Sum of height C = 0.3');

        % ------------------------------------------------------------------- %
        % As a conclusion, we do have the wave leaving the mesh as we see
        % that the overall sum of height diminishes just before 1.5 s.
        %
        % ------------------------------------------------------------------- %
        