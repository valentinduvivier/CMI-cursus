%% 2.1.4 - Order of Accuracy - Delta t

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
        
        Delta_x = L/(N_x-1);    % Space step
        Delta_t = T/N_t;        % Time step
        
    % Checking stability (CFL condition)
        CFL = (4*Delta_t)/Delta_x;
    
        if CFL > 1
            disp('error CFL conition not fulfilled');
        end 
        
% ----------------------------------------------------------------------------------------------- %
% ----------------------------------------------------------------------------------------------- %

        % Variables for accuracy
            h1 = 0.001; h2 = (1/2)*h1; h3 = (1/2)*h2; h4 = (1/2)*h3; h5 = (1/2)*h4;
            h6 = (1/2)*h5; h7 = (1/2)*h6; h8 = (1/2)*h7; h9 = (1/2)*h8;

            h10 = (1/2)*h9; h11 = (1/2)*h10; h12 = (1/2)*h11; h13 = (1/2)*h12; h14 = (1/2)*h13;
            h15 = (1/2)*h14; h16 = (1/2)*h15; h17 = (1/2)*h16; h18 = (1/2)*h17;

            h19 = (1/2)*h18; h20 = (1/2)*h19; h21 = (1/2)*h20; h22 = (1/2)*h21; h23 = (1/2)*h22;
            h24 = (1/2)*h23; h25 = (1/2)*h24; h26 = (1/2)*h25; h27 = (1/2)*h26;        

            h = [h1, h2, h3, h4, h5, h6, h7, h8, h9, ...
                h10, h11, h12, h13, h14, h15, h16, h17, h18, ...
                h19, h20, h21, h22, h23, h24, h25, h26, h27];  
    
        U_j = zeros(1,size(h,2));
        
    for j = 1:size(h,2)
        
        u = zeros(2, N_x+2);  
        
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
                    u(:,k) = u(:,k) - (h(j)/Delta_x)*(F(:,k) - F(:,k-1));
                end

                U_2 = [U_2; u];
        end
        
        U_j(j) = U_2(end-1,5);
    end       
                
%% Presentation results
 
    % u_h - u_h/2
        P1 = (U_j(1) - U_j(2))/(U_j(2) - U_j(3));
        P2 = (U_j(2) - U_j(3))/(U_j(3) - U_j(4));
        P3 = (U_j(3) - U_j(4))/(U_j(4) - U_j(5));
        P4 = (U_j(4) - U_j(5))/(U_j(5) - U_j(6));
        P5 = (U_j(5) - U_j(6))/(U_j(6) - U_j(7));
        P6 = (U_j(6) - U_j(7))/(U_j(7) - U_j(8));
        P7 = (U_j(7) - U_j(8))/(U_j(8) - U_j(9));
        P8 = (U_j(8) - U_j(9))/(U_j(9) - U_j(10));
        P9 = (U_j(9) - U_j(10))/(U_j(10) - U_j(11));
        P10 = (U_j(10) - U_j(11))/(U_j(11) - U_j(12));
        P11 = (U_j(11) - U_j(12))/(U_j(12) - U_j(13));
        P12 = (U_j(12) - U_j(13))/(U_j(13) - U_j(14));
        P13 = (U_j(13) - U_j(14))/(U_j(14) - U_j(15));
        P14 = (U_j(14) - U_j(15))/(U_j(15) - U_j(16));
        P15 = (U_j(15) - U_j(16))/(U_j(16) - U_j(17));
        P16 = (U_j(16) - U_j(17))/(U_j(17) - U_j(18));
        P17 = (U_j(17) - U_j(18))/(U_j(18) - U_j(19));
        P18 = (U_j(18) - U_j(19))/(U_j(19) - U_j(20));
        P19 = (U_j(19) - U_j(20))/(U_j(20) - U_j(21));
        P20 = (U_j(20) - U_j(21))/(U_j(21) - U_j(22));
        P21 = (U_j(21) - U_j(22))/(U_j(22) - U_j(23));
        P22 = (U_j(22) - U_j(23))/(U_j(23) - U_j(24));
        P23 = (U_j(23) - U_j(24))/(U_j(24) - U_j(25));
        P24 = (U_j(24) - U_j(25))/(U_j(25) - U_j(26));
        P25 = (U_j(25) - U_j(26))/(U_j(26) - U_j(27));

    % log_2((u_h - u_h/2)/(u_h/2 - u_h/4)) --> relation sought
        M1 = log(P1)/log(2);
        M2 = log(P2)/log(2);
        M3 = log(P3)/log(2);
        M4 = log(P4)/log(2);
        M5 = log(P5)/log(2);
        M6 = log(P6)/log(2);
        M7 = log(P7)/log(2);
        M8 = log(P8)/log(2);
        M9 = log(P9)/log(2);
        M10 = log(P10)/log(2);
        M11 = log(P11)/log(2);
        M12 = log(P12)/log(2);
        M13 = log(P13)/log(2);
        M14 = log(P14)/log(2);
        M15 = log(P15)/log(2);
        M16 = log(P16)/log(2);
        M17 = log(P17)/log(2);
        M18 = log(P18)/log(2);
        M19 = log(P19)/log(2);
        M20 = log(P20)/log(2);
        M21 = log(P21)/log(2);
        M22 = log(P22)/log(2);
        M23 = log(P23)/log(2);
        M24 = log(P24)/log(2);
        M25 = log(P25)/log(2);
      
        