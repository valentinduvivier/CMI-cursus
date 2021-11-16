%% 2.1 - Tests with flat bottom - a = 2*H

    % B = 0 --> s(h,m,x) = 0 as flat bottom
    
    % General variables
        
        L = 10;         % Length of the mesh
        T = 20;         % Time defining the study (we look at the function for T second(s))
        
        N_x = 80;       % Number of point defining the mesh
        iter_t = 1321;  % Number of iteration
        
        H = 1;          % Initial height of the wave
        w = 0.1*L;      % Width of the Gaussian pulse
        a = 2*H;        % Coeff to ensure a stable scheme (low right hand side with respect to H)
        
        g = 9.81;       % Gravity's acceleration
        
    % Mesh
        x = linspace(0, L, N_x);
        
        Delta_x = L/(N_x-1);    % Space step
        Delta_t = T/iter_t;     % Time step
        
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
            for k = 1:size(x,2)
                u(1,k+1) = H + a*exp(-(x(k)-L/2)^2)/(w^2);
            end
        
        % Speed
            % To be determined to have a single pulse
            for k = 1:size(u,2)-2
                u(2,k+1) = (u(1,k) - H)*sqrt(g*u(1,k));
            end

        % Variables to stock the iteration of u
        % Even lines = m; Odd lines = height
            U_3 = zeros(2*(iter_t+1), N_x+2);
            U_3 = [u]; 

        for t = 1:iter_t
            
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
            % depending function U = (h,h*v))
                U_3 = [U_3; u];
        end
                
%% Display of the wave - drawnow

        % U_3(2*k + 2, 2:end-1)./U(2*k + 1, 2:end-1)  to get v
        % U_3(2*k + 1, 2:end-1)                       to get h
        
%         figure(40)
%         for k = 1:N_t
%             plot(x, U_3(2*k + 1, 2:end-1), 'r');
%             axis([0 10 0 3.5]);
%             xlabel('Position on the mesh [m]'); ylabel('Wave''s height (m)'); title(['Wave at t = ',num2str(round((T/N_t)*k,2)),' s']);
%             drawnow
%         end
%         
%         figure(41)
%         for k = 1:T+1
%             plot(x, U_3(2*(k+1),2:end-1), 'r');
% %             axis([0 10 0.95 1.1]);
%             xlabel('Position on the mesh [m]'); ylabel('Wave''s height (m)'); title(['Wave at t = ',num2str(k),'*Delta_t ~~ 1 bouncing']);
%             drawnow
%         end
                
%% STABILITY

%         figure(42)
%             subplot(2,2,1)
%             k = round(iter_t/iter_t);
%             plot(x, U_3(2*k + 1, 2:end-1), 'r');
%             axis([0 10 0.5 3.2]);
%             xlabel('Position on the mesh [m]'); ylabel('Wave''s height [m]'); title(['t = ',num2str(round(k*T/iter_t,2)),' s']);
% 
%             subplot(2,2,2)
%             k = round(iter_t/15);
%             plot(x, U_3(2*k + 1, 2:end-1), 'r');
%             axis([0 10 0.5 3.2]);
%             xlabel('Position on the mesh [m]'); ylabel('Wave''s height [m]'); title(['t = ',num2str(round(k*T/iter_t,2)),' s']);
% 
%             subplot(2,2,3)
%             k = round(iter_t/7);
%             plot(x, U_3(2*k + 1, 2:end-1), 'r');
%             axis([0 10 0.5 3.2]);
%             xlabel('Position on the mesh [m]'); ylabel('Wave''s height [m]'); title(['t = ',num2str(round(k*T/iter_t,2)),' s']);
% 
%             subplot(2,2,4)
%             k = round(iter_t/1);
%             plot(x, U_3(2*k + 1, 2:end-1), 'r');
%             axis([0 10 0.5 3.2]);
%             xlabel('Position on the mesh [m]'); ylabel('Wave''s height [m]'); title(['t = ',num2str(round(k*T/iter_t,2)),' s']);

 %% Presentation results

        figure(43)
            subplot(1,3,1)
            k = round(iter_t/50);
            plot(x, U_3(2*k + 1, 2:end-1), 'r');
            axis([0 10 0.5 3.2]);
            % xlabel('Position on the mesh [m]'); 
            ylabel('Wave''s height [m]'); title(['t = ',num2str(round(k*T/iter_t,2)),' s']);

            subplot(1,3,2)
            k = round(iter_t/20);
            plot(x, U_3(2*k + 1, 2:end-1), 'r');
            axis([0 10 0.5 3.2]);
            xlabel('Position on the mesh [m]'); % ylabel('Wave''s height [m]'); 
            title(['t = ',num2str(round(k*T/iter_t,2)),' s']);

            subplot(1,3,3)
            k = round(iter_t/10);
            plot(x, U_3(2*k + 1, 2:end-1), 'r');
            axis([0 10 0.5 3.2]);
            % xlabel('Position on the mesh [m]'); ylabel('Wave''s height [m]'); 
            title(['t = ',num2str(round(k*T/iter_t,2)),' s']);
