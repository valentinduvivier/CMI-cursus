%% 2.1.4 - Order of Accuracy - Delta x

    % B = 0 --> s(h,m,x) = 0 as flat bottom
    
    % General variables
        L = 10;         % Length of the mesh
        T = 20;         % Time defining the study (we look at the function for T second(s))
        
        N_x = 80;       % Number of point defining the mesh
        N_t = 3000;     % Number of iteration
        
        H = 1;          % Initial height of the wave
        w = 0.1*L;      % Width of the Gaussian pulse
        a = H/5;        % Coeff to ensure a stable scheme (low right hand side with respect to H)
        
        g = 9.81;       % Gravity's acceleration
        
    % Mesh
        x = linspace(0, L, N_x);
        
        Delta_x = L/N_x;    % Space step
        Delta_t = T/N_t;        % Time step
        
    % Checking stability (CFL condition)
        CFL = (4*Delta_t)/Delta_x;
    
        if CFL > 1
            disp('error CFL conition not fulfilled');
        end 

% ----------------------------------------------------------------------------------------------- %
% ----------------------------------------------------------------------------------------------- %
        
        % Variables for accuracy
            h1 = 4; h2 = h1/2; h3 = h2/2; h4 = h3/2; h5 = h4/2;
            h6 = h5/2; h7 = h6/2;

            h = [h1, h2, h3, h4, h5, h6, h7];

        U_j = zeros(1,size(h,2));
        
    for j = 1:size(h,2)
        
        N_x = round(L / h(j));
        x   = linspace(0, L, N_x);
        u   = zeros(2, N_x+2);  
        
        % IC
            % Height
                for k = 1:N_x
                    u(1,k+1) = H + a*exp(-(x(k)-L/2)^2)/(w^2);
                end

            % Moment
                for k = 2:N_x+1
                    u(2,k) = (u(1,k) - H)*sqrt(g*u(1,k));
                end
        
        U_2 = zeros(2*(N_t+1), N_x+2);
        U_2 = [u];
        
        for t = 1:N_t
            
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
                    u(:,k) = u(:,k) - (Delta_t/h(j))*(F(:,k) - F(:,k-1));
                end

                U_2 = [U_2; u];
        end
        
        indice = find(x<L/2);
        indice
        h(j)
        U_j(j) = U_2(end-1,5);
    end
    
%% Presentation results
 
    % (u_h - u_h/2)/(u_h/2 - u_h/4)
        P1 = abs((U_j(1) - U_j(2))/(U_j(2) - U_j(3)));
        P2 = abs((U_j(2) - U_j(3))/(U_j(3) - U_j(4)));
        P3 = abs((U_j(3) - U_j(4))/(U_j(4) - U_j(5)));
        P4 = abs((U_j(4) - U_j(5))/(U_j(5) - U_j(6)));
        P5 = abs((U_j(5) - U_j(6))/(U_j(6) - U_j(7)));

    % log_2((u_h - u_h/2)/(u_h/2 - u_h/4)) --> relation sought
        M1 = log2(P1);
        M2 = log2(P2);
        M3 = log2(P3);
        M4 = log2(P4);
        M5 = log2(P5);
       
        
% A bit fixed but still not very good