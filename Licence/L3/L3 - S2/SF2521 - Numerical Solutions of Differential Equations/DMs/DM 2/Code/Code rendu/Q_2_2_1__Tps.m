
    % Vector containing u = (h,hv)^T
            u_4 = zeros(2,n+2); % We consider h the upper vector and hv the lower ones

    % Conditions : initial and boundary

            % Variables for dependence quantification
            h1 = Delta_t; h2 = h1/2; h3 = h1/4;
            
            h_1 = [h1, h2, h3];
            
            h_2 = 0.3*h_1; h_3 = 0.4*h_1; h_4 = 0.7*h_1;
            h_5 = 0.47*h_1; h_6 = 0.12*h_1; h_7 = 0.83*h_1;
            h_8 = 0.98*h_1;
            
            h = [h_1, h_2, h_3, h_4, h_5, h_6, h_7, h_8];
            Delta_t_4 = h;
            
%             Delta_t_4 = linspace(0.001,0.01,10);
            
%             Delta_x_4 = linspace(0.011,0.13,10);
            Delta_x_4 = Delta_x;
            
            alpha_4 = alpha;
            
        % Speed
            % Initial condition
            for k = 1:n+2   
                u_4(2,k) = 0;     % Initial speed (or initial moment m = h*u = 0)
            end
     
        % Height
            % Initial condition
            % We keep h(1) and h(n+2) at zero and the rest considering the
            % initial height
            for k = 1:n
               u_4(1,k+1) = H + epsilon*exp(-((x(k) - L/2)^2)/w^2); % 
            end
  
%             figure(15)     % Graphic of the initial height 
%             plot(x, u_4(1,2:n+1), 'r'); % We represent only the point of h for the considered initial value (size(x,2) points)
%             xlabel('Position x (m)'); ylabel('Initial height (m)'); title('Graphic representing the function h = f(x) at t=0s');
%             
            
%%            
        
        % Vector to keep in memory the iteration of the height and speed
            U_4 = zeros(2* (size(Delta_t_4,2)*(iter_t+1) +1),n+2);
            U_4 = [u_4]; 
            
            f = zeros(2,n+2);
            FF_4 = zeros(2*( size(Delta_t_4,2)*(iter_t+1) +1),n+2);
            FF_4 = [f];
            
            tau_4 = (zeros(size(Delta_t_4,2),size(t_n,2)));
            
for i = 1:size(Delta_t_4,2)    
        for t = t_n
            
            % Reflecting boundary for the height + boundary for the moment
                u_4(1,1) = u_4(1,2);    u_4(1,n+2) = u_4(1,n+1);
                
                u_4(2,1) = u_4(2,2);    u_4(2,n+2) = u_4(2,n+1);            
                
            % Variables from general equation's form 
                % Function f to work on F (Lax-Friedrichs flux) :
                    f = [u_4(2,:) ; (u_4(2,:).^2)./u_4(1,:) + (1/2)*g*(u_4(1,:).^2)];
            % We calulate the speed and height
                % What is at u(:,1) and u(:,n+2) is already defined : 
                %   - v is always zero at u(:,1) and u(:,n+2)
                %   - h(1) and h(n+2) as well as m(1) and m(n+2) (m = h*v) are defined above at the begining of the loop
                
                % Lax-Friedrichs flux
                for k = 1:n+1
                    F1(:,k) = (1/2)*(f(:,k+1) + f(:,k) - alpha*(u_4(:,k+1) - u_4(:,k)));
                end
                
                % Resolution of the system to get h and h*v
                for k = 2:n+1
                    u_4(:,k) = u_4(:,k) - (Delta_t_4(i)/Delta_x_4)*(F1(:,k) - F1(:,k-1));
                end
                
            % We stock the value for the different iteration (time
            % depending function U = (h,h*v))
                 U_4 = [U_4; u_4];
                 FF_4 = [FF_4; f];
        end
end


        
%%        
        % The formula used for tau is given on the report for more clarity
        % We calculate tau which is the truncation error at the points [h_0,v_0]
        % We consider the point x = n/2, middle of the mesh.
        for t = 1:size(Delta_t_4,2)
            for i = 1:iter_t-1
                tau_4(t,i) = (U_4(2*(i-1)+3 + 2*iter_t*(t-1), n/2) - U_4(2*(i-1)+1 + 2*iter_t*(t-1), n/2))/Delta_t_4(t) + ...
                    (FF_4(2*(i-1)+1 + 2*iter_t*(t-1), n/2+1) - FF_4(2*(i-1)+1 + 2*iter_t*(t-1), n/2-1) - alpha_4*( U_4(2*(i-1)+1 + 2*iter_t*(t-1), n/2+1) ...
                    - 2*U_4(2*(i-1)+1 + 2*iter_t*(t-1), n/2) + U_4(2*(i-1)+1 + 2*iter_t*(t-1), n/2-1) ))/(Delta_x_4*2);
            end
        end
        
        % We here apply the L_1 norm to tau
        for t = 1:size(Delta_t_4,2)
            Sum = 0;
            for j = 1:iter_t-1
                Sum = Sum + abs(tau_4(t,j));
            end
            tau_norm_tps(1,t) = Sum*Delta_x_4;
        end
        
        %%
        
        P1 = (tau_norm_tps(1) - tau_norm_tps(2))/(tau_norm_tps(2) - tau_norm_tps(3));
        P2 = (tau_norm_tps(4) - tau_norm_tps(5))/(tau_norm_tps(5) - tau_norm_tps(6));
        P3 = (tau_norm_tps(7) - tau_norm_tps(8))/(tau_norm_tps(8) - tau_norm_tps(9));
        P4 = (tau_norm_tps(10) - tau_norm_tps(11))/(tau_norm_tps(11) - tau_norm_tps(12));
        
        P5 = (tau_norm_tps(13) - tau_norm_tps(14))/(tau_norm_tps(14) - tau_norm_tps(15));
        P6 = (tau_norm_tps(16) - tau_norm_tps(17))/(tau_norm_tps(17) - tau_norm_tps(18));
        P7 = (tau_norm_tps(19) - tau_norm_tps(20))/(tau_norm_tps(20) - tau_norm_tps(21));
        P8 = (tau_norm_tps(22) - tau_norm_tps(23))/(tau_norm_tps(23) - tau_norm_tps(24));

        figure(16)
        plot(h1, log(P1)/log(2), '*', 'Color', [0 0.5 0.5]);
        hold on
        plot(h1*0.3, log(P2)/log(2), '*', 'Color', [0 0.5 0.5]);
        plot(h1*0.4, log(P3)/log(2), '*', 'Color', [0.7 0 0.3]);
        plot(h1*0.7, log(P4)/log(2), '*', 'Color', [0.7 0.3 0]);
        plot(h1*0.47, log(P5)/log(2), '*', 'Color', [0 0.7 0.3]);
        plot(h1*0.12, log(P6)/log(2), '*', 'Color', [0.5 0.5 0]);
        plot(h1*0.83, log(P7)/log(2), '*', 'Color', [0.5 0 0.5]);
        plot(h1*0.98, log(P8)/log(2), '*', 'Color', [0.5 0.5 0.5]);
        
        %%
        
        % We display the value of norm L_1 of tau depending on the value of
        % Delta_t_4
        % We thus get the type of dependence and a coeff for it
        figure(17)
        for t = 1:size(Delta_t_4,2)
            plot(Delta_t_4, (tau_norm_tps), '*'); 
            axis([0 0.012 -0.02 0.02]);
            xlabel('Delta_t'); ylabel('(Norm L_1 of the Truncation error)'); title('Graph representing tau = a*Delta_t + b');
        end
        
        
        
        