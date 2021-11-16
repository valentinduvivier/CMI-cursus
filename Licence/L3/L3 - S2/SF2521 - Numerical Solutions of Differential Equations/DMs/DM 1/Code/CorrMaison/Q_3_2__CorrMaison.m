
%% Q_3.2 
    % I. Calculation of the errors
    
%  % ---------------------------------------------------------------------------------------------------------------------------- %        
    
    % Initialisation
    
    % Space variables
        Lx = 1; Ly = 1;         % Dimensions of the area over which we descretize
        
    % We consider a square mesh for now on
        N = 10;    
        M = N;

        Dx = Lx/M; Dy = Ly/N;   % Space step

    % Time variables
        T = 10;         % The code will be running for 10s
        N_t = 100;      % Time scale - Nb of iter°

        Dt = T/N_t;     % Time step

    % Diverse variables
        w = 1/4;
        
    % We create point to define the mesh for the dicretization
        x = linspace(0, Lx, M);
        y = linspace(0, Ly, N);

        Q = zeros(M,N);     % As i = 1 --> M and j = 1 --> N


%  % ---------------------------------------------------------------------------------------------------------------------------- %        

      
    % Ax
        Ax = zeros(M);

        % Diagonal
            % (See BC)
                Ax(1,1) = -1;
                Ax(M,M) = -1;

            for k = 2:M-1
                Ax(k,k) = -2;
            end

        % Diagonal superior & inferior
            for k = 1:M-1
                Ax(k,k+1) = 1;

                Ax(k+1,k) = 1;
            end
        
    % Ay
        Ay = zeros(N);

        % Diagonal
            % (See BC)
                Ay(1,1) = -1;
                Ay(N,N) = -1;

            for k = 2:N-1
                Ay(k,k) = -2;
            end

        % Diagonal superior & inferior
            for k = 1:N-1
                Ay(k,k+1) = 1;

                Ay(k+1,k) = 1;
            end

%  % ---------------------------------------------------------------------------------------------------------------------------- %        
        
  
%     % PART 1 - Influence of space step h
% 
%     [X3,Y3] = meshgrid(x,y);
%     
%     for k = 1:18
%         h(k) = (1/(2^(k-1))) * Dx;
%     end
%
%     Q_p = zeros(1,size(h,2));
%     
%     % We loop for the space-step defined upper
%     for p = 1:size(h,2)
%         xs = 1/2; ys = 1/2;
%         w = 1/4;
%         
%         S1 = exp(-((X3 - xs).^2 + (Y3 - ys).^2)/(w^2));
% 
%         r = sqrt(((X3 - xs).^2 + (Y3 - ys).^2));
%         eps = sqrt(h(p));
% 
%         S2 = zeros(N);  % M = N
% 
%         for i = 1:M
%             for j = 1:N
%                 if r(j,i) < eps
%                     S2(j,i) = 2*(pi/(eps^2*((pi^2) - 4))) * (1 + cos((pi/eps)*r(j,i)));
%                 end
%             end
%         end
%         
% %  % ---------------------------------------------------------------------------------------------------------------------------- %
%         
%         
%         % We calculate Tx and Ty for a different space-step at each loop
%             Tx = kron(Ax,eye(N))/h(p)^2;
%             Ty = kron(eye(M),Ay)/h(p)^2;
% 
%             I = eye(M*N);
% 
%         % Constant matrices 
%             A = (I - Dt*(Tx+Ty));
%             [L,U] = lu(A);
% 
%             Q = zeros(M);     % M = N
%             
%         % Reorganisation of matrix Q as a column vector
%             QQ = sparse(reshape(Q,[],1));
%         
%         % Matrix gathering QQ over time. This matrix will allow us to deal
%         % with the heat function over time
%             FF = sparse([QQ]);
% 
%  % ---------------------------------------------------------------------------------------------------------------------------- %        
% 
% 
% % S1 - Uncomment all if you're willing to test S1   
% 
%         % Source vector
%         SS = sparse(reshape(S1,[],1));
% 
%         for t = 1:N_t     % We loop for as many itera° as we have 
%             Y = L^-1*(QQ + Dt*SS);
%             QQ = U^-1*Y;  % Resolution of heat equation in 2D with source
%             
%             % We ensure that the temperature inside the mesh dosn't go
%             % above the source heat
%             for k = 1:size(QQ,1)
%                if QQ(k,1) > max(SS)
%                   QQ(k,1) = max(SS);
%                end
%             end
%             
%         FF = [FF QQ];
%             
%         end
%         
%         
% % S2 - Uncomment all if you're willing to test S2 
% 
% %         SS = reshape(S2,[],1);
% % 
% %         for t = 1:N_t
% % 
% %             while Dt*t < 0.25
% %                 Y = L^-1*(QQ + Dt*SS);
% %                 QQ = U^-1*Y;                % Resolution of heat equation in 2D with source
% %                 
% %                 % We ensure that the temperature inside the mesh dosn't go
% %                 % above the source heat
% %                 for k = 1:size(QQ,1)
% %                    if QQ(k,1) > max(SS)
% %                       QQ(k,1) = max(SS);
% %                    end
% %                 end
% %                 
% %                 FF = [FF QQ];
% % 
% %                 t = t + 1;
% %             end
% %             
% %             Y = L^-1*QQ;
% %             QQ = U^-1*Y;            
% %              
% %             % We ensure that the temperature inside the mesh dosn't go
% %             % above the source heat
% %             for k = 1:size(QQ,1)
% %                 if QQ(k,1) > max(SS)
% %                   QQ(k,1) = max(SS);
% %                 end
% %             end
% % 
% %         FF = [FF QQ];
% % 
% %         end
%         
% 
%         FF = FF(:,2:end);
% 
% %  % ---------------------------------------------------------------------------------------------------------------------------- %        
% 
%         Q_p(p) = FF(4,5);   % Random point to check for space and time dependance
%         % Some columns of FF might not work as we define that given a
%         % certain value (equal to the max of the source heat), we have the
%         % heat that is stagnating (100°C + 100°C = 100°C). Thereby, the
%         % difference between our value is then compromised . We thus
%         % suggest to pay attention to these column and to restrict the
%         % nalysis to the "good" columns, i.e. the ones where the heat isn't to
%         % its maximum.
%     end
%       
%  
% % Presentation results
%  
%     % (u_h - u_h/2) / (u_h/2 - u_h/4)
%         P1 = (Q_p(1) - Q_p(2))/(Q_p(2) - Q_p(3));
%         P2 = (Q_p(2) - Q_p(3))/(Q_p(3) - Q_p(4));
%         P3 = (Q_p(3) - Q_p(4))/(Q_p(4) - Q_p(5));
%         P4 = (Q_p(4) - Q_p(5))/(Q_p(5) - Q_p(6));
%         P5 = (Q_p(5) - Q_p(6))/(Q_p(6) - Q_p(7));
%         P6 = (Q_p(6) - Q_p(7))/(Q_p(7) - Q_p(8));
%         P7 = (Q_p(7) - Q_p(8))/(Q_p(8) - Q_p(9));
%         P8 = (Q_p(8) - Q_p(9))/(Q_p(9) - Q_p(10));
%         P9 = (Q_p(9) - Q_p(10))/(Q_p(10) - Q_p(11));
%         P10 = (Q_p(10) - Q_p(11))/(Q_p(11) - Q_p(12));   
%         
%         P11 = (Q_p(11) - Q_p(12))/(Q_p(12) - Q_p(13));
%         P12 = (Q_p(12) - Q_p(13))/(Q_p(13) - Q_p(14));
%         P13 = (Q_p(13) - Q_p(14))/(Q_p(14) - Q_p(15));
%         P14 = (Q_p(14) - Q_p(15))/(Q_p(15) - Q_p(16));
%         P15 = (Q_p(15) - Q_p(16))/(Q_p(16) - Q_p(17));
%         P16 = (Q_p(16) - Q_p(17))/(Q_p(17) - Q_p(18));
% %         P17 = (Q_p(17) - Q_p(18))/(Q_p(18) - Q_p(19));
% %         P18 = (Q_p(18) - Q_p(19))/(Q_p(19) - Q_p(20));
% %         P19 = (Q_p(19) - Q_p(20))/(Q_p(20) - Q_p(21));
% %         P20 = (Q_p(20) - Q_p(21))/(Q_p(21) - Q_p(22));
% %         
% %         P21 = (Q_p(21) - Q_p(22))/(Q_p(22) - Q_p(23));
% %         P22 = (Q_p(22) - Q_p(23))/(Q_p(23) - Q_p(24));
% %         
%     % log_2((u_h - u_h/2)/(u_h/2 - u_h/4)) --> relation sought
%         M1 = log(P1)/log(2);
%         M2 = log(P2)/log(2);
%         M3 = log(P3)/log(2);
%         M4 = log(P4)/log(2);
%         M5 = log(P5)/log(2);
%         M6 = log(P6)/log(2);
%         M7 = log(P7)/log(2);
%         M8 = log(P8)/log(2);
%         M9 = log(P9)/log(2);
%         M10 = log(P10)/log(2);
% 
%         M11 = log(P11)/log(2);
%         M12 = log(P12)/log(2);
%         M13 = log(P13)/log(2);
%         M14 = log(P14)/log(2);
%         M15 = log(P15)/log(2);
%         M16 = log(P16)/log(2);
% %         M17 = log(P17)/log(2);
% %         M18 = log(P18)/log(2);
% %         M19 = log(P19)/log(2);
% %         M20 = log(P20)/log(2);
%         
%         % => The tendency of M gives the power of influence of h
% 
%         figure(7)
%         
%         xx = linspace(10^-4,1,100);
%         yy = xx.^2;
%         
%         for k = 1:size(h,2)-1
%             loglog(h(k), abs(Q_p(k+1) - Q_p(k)), '*r', xx, yy, 'b');
%             xlabel('h'); ylabel('|Q_{p+1} - Q_{p}|'); title('|Q_{p+1} - Q_{p}| in function of h');
%             legend('Q_{p+1} - Q_{p}', 'y = x^2', 'location', 'northwest');
%             grid on
%             hold on
%         end
   

% ------------------------------------------------------------------------------------------- %
% ------------------------------------------------------------------------------------------- %


%     % PART 2 - Influence of time step dt
%     
%     % We consider a square mesh for now on
%         M = N;
%         h = Dx;     % Dx = Dy = h;
%     
%     % We define a few space-step that are gonna be used to get the factor of
%     % influence of the space-step over the heat equation (Q_3.2)
%     
%     
%     % ---------------------------------------------------------------------------- %
%     %    We define our k such that following k are multiple of 2. We seek
%     %    the factor of influence of h. The Q_p chosen is arbitrary, and has
%     %    the only purpose to get the coefficient p (see ConvRate paper + subject DM1).
%     % ---------------------------------------------------------------------------- %
% 
%     for k = 1:18
%        dt(k) = (1/(2^(k-1))) * (T/N_t); 
%     end
% 
%     
% %  % ---------------------------------------------------------------------------------------------------------------------------- %        
%     
% 
%     [X4,Y4] = meshgrid(x,y);
%     
%     Q_p = zeros(1,size(dt,2));
%         
%     % We loop for the space-step defined upper
%     for p = 1:size(dt,2)
%         xs = 1/2; ys = 1/2;
%         w = 1/4;
%         
%         S1 = exp(-((X4 - xs).^2 + (Y4 - ys).^2)/(w^2));
% 
%         r = sqrt(((X4 - xs).^2 + (Y4 - ys).^2));
%         eps = sqrt(h);
% 
%         S2 = zeros(N,M);
% 
%         for i = 1:M
%             for j = 1:N
%                 if r(j,i) < eps
%                     S2(j,i) = 2*(pi/(eps^2*((pi^2) - 4))) * (1 + cos((pi/eps)*r(j,i)));
%                 end
%             end
%         end
%         
% %  % ---------------------------------------------------------------------------------------------------------------------------- %
%         
%         
%         % We calculate Tx and Ty for a different space-step at each loop
%             Tx = kron(Ax,eye(N))/h^2;
%             Ty = kron(eye(M),Ay)/h^2;
% 
%             I = eye(M*N);
% 
%         % Constant matrices 
%             A = (I - dt(p)*(Tx+Ty));
%             [L,U] = lu(A);
% 
%             Q = zeros(M);     % M = N
%             
%         % Reorganisation of matrix Q as a column vector
%             QQ = sparse(reshape(Q,[],1));
%         
% %         % Matrix gathering QQ over time. This matrix will allow us to deal
% %         % with the heat function over time
%             FF = sparse([QQ]);
% 
% %  % ---------------------------------------------------------------------------------------------------------------------------- %        
% 
% 
% % S1 - Uncomment all if you're willing to test S1   
% 
%         % Source vector
%         SS = sparse(reshape(S1,[],1));
% 
%         for t = 1:N_t     % We loop for as many itera° as we have 
%             Y = L^-1*(QQ + dt(p)*SS);
%             QQ = U^-1*Y;  % Resolution of heat equation in 2D with source
%             
%             % We ensure that the temperature inside the mesh dosn't go
%             % above the source heat
%             for k = 1:size(QQ,1)
%                if QQ(k,1) > max(SS)
%                   QQ(k,1) = max(SS);
%                end
%             end
%             
%         FF = [FF QQ];
%             
%         end
%         
%         
% % S2 - Uncomment all if you're willing to test S2 
% 
% %         SS = reshape(S2,[],1);
% % 
% %         for t = 1:N_t
% 
% %             while dt(p)*t < 0.25
% %                 Y = L^-1*(QQ + dt(p)*SS);
% %                 QQ = U^-1*Y;                % Resolution of heat equation in 2D with source
% %                 
% %                 % We ensure that the temperature inside the mesh dosn't go
% %                 % above the source heat
% %                 for k = 1:size(QQ,1)
% %                    if QQ(k,1) > max(SS)
% %                       QQ(k,1) = max(SS);
% %                    end
% %                 end
% %                 
% %                 FF = [FF QQ];
% % 
% %                 t = t + 1;
% %             end
% %             
% %             Y = L^-1*QQ;
% %             QQ = U^-1*Y;            
% %              
% %             % We ensure that the temperature inside the mesh dosn't go
% %             % above the source heat
% %             for k = 1:size(QQ,1)
% %                 if QQ(k,1) > max(SS)
% %                   QQ(k,1) = max(SS);
% %                 end
% %             end
% % 
% %         FF = [FF QQ];
% % 
% %         end
%         
% 
%         FF = FF(:,2:end);
% 
% %  % ---------------------------------------------------------------------------------------------------------------------------- %        
% 
%         Q_p(p) = FF(8,8);   % Random point to check for space and time dependance
%                 
%     end
%       
%  
%     % Presentation results
%  
%     % (u_h - u_h/2) / (u_h/2 - u_h/4)
%         P1 = (Q_p(1) - Q_p(2))/(Q_p(2) - Q_p(3));
%         P2 = (Q_p(2) - Q_p(3))/(Q_p(3) - Q_p(4));
%         P3 = (Q_p(3) - Q_p(4))/(Q_p(4) - Q_p(5));
%         P4 = (Q_p(4) - Q_p(5))/(Q_p(5) - Q_p(6));
%         P5 = (Q_p(5) - Q_p(6))/(Q_p(6) - Q_p(7));
%         P6 = (Q_p(6) - Q_p(7))/(Q_p(7) - Q_p(8));
%         P7 = (Q_p(7) - Q_p(8))/(Q_p(8) - Q_p(9));
%         P8 = (Q_p(8) - Q_p(9))/(Q_p(9) - Q_p(10));
%         P9 = (Q_p(9) - Q_p(10))/(Q_p(10) - Q_p(11));
%         P10 = (Q_p(10) - Q_p(11))/(Q_p(11) - Q_p(12));   
%         
%         P11 = (Q_p(11) - Q_p(12))/(Q_p(12) - Q_p(13));
%         P12 = (Q_p(12) - Q_p(13))/(Q_p(13) - Q_p(14));
%         P13 = (Q_p(13) - Q_p(14))/(Q_p(14) - Q_p(15));
%         P14 = (Q_p(14) - Q_p(15))/(Q_p(15) - Q_p(16));
%         P15 = (Q_p(15) - Q_p(16))/(Q_p(16) - Q_p(17));
%         P16 = (Q_p(16) - Q_p(17))/(Q_p(17) - Q_p(18));
% %         P17 = (Q_p(17) - Q_p(18))/(Q_p(18) - Q_p(19));
% %         P18 = (Q_p(18) - Q_p(19))/(Q_p(19) - Q_p(20));
% %         P19 = (Q_p(19) - Q_p(20))/(Q_p(20) - Q_p(21));
% %         P20 = (Q_p(20) - Q_p(21))/(Q_p(21) - Q_p(22));
% %         
% %         P21 = (Q_p(21) - Q_p(22))/(Q_p(22) - Q_p(23));
% %         P22 = (Q_p(22) - Q_p(23))/(Q_p(23) - Q_p(24));
% %         
%     % log_2((u_h - u_h/2)/(u_h/2 - u_h/4)) --> relation sought
%         M1 = log(P1)/log(2);
%         M2 = log(P2)/log(2);
%         M3 = log(P3)/log(2);
%         M4 = log(P4)/log(2);
%         M5 = log(P5)/log(2);
%         M6 = log(P6)/log(2);
%         M7 = log(P7)/log(2);
%         M8 = log(P8)/log(2);
%         M9 = log(P9)/log(2);
%         M10 = log(P10)/log(2);
% 
%         M11 = log(P11)/log(2);
%         M12 = log(P12)/log(2);
%         M13 = log(P13)/log(2);
%         M14 = log(P14)/log(2);
%         M15 = log(P15)/log(2);
%         M16 = log(P16)/log(2);
% %         M17 = log(P17)/log(2);
% %         M18 = log(P18)/log(2);
% %         M19 = log(P19)/log(2);
% %         M20 = log(P20)/log(2);
%         
%         % => The tendency of M gives the power of influence of h
% 
%         figure(8)
%         
%         xx = linspace(10^-4,1,100);
%         yy = xx.^2;
%         
%         for k = 1:size(dt,2)-1
%             loglog(dt(k), abs(Q_p(k+1) - Q_p(k)), 'r*-', xx, yy, 'b');
%             xlabel('\Delta_t'); ylabel('|Q_{p+1} - Q_{p}|'); title('|Q_{p+1} - Q_{p}| in function of \Delta_t');
%             legend('Q_{p+1} - Q_{p}', 'y = x^2', 'location', 'northwest');
%             grid on
%             hold on
%         end
%         
 % ---------------------------------------------------------------------------------------------------------------------------- %        
 % ---------------------------------------------------------------------------------------------------------------------------- %        


        % II. Influence of epsilon on source S2

        h = Dx;

        % We test our epsilon for different alpha
        alpha = [1 0.1 0.01 0.001 0.0001 10e-10];

        [X5,Y5] = meshgrid(x,y);

        Q_p = zeros(1,size(h,2));
    
        % We loop for the space-step defined upper
        for p = 1:size(alpha,2)
            xs = 1/2; ys = 1/2;
            w = 1/4;

            r = sqrt(((X5 - xs).^2 + (Y5 - ys).^2));
            eps = alpha(p)*sqrt(h);

            S2 = zeros(N);  % M = N

            for i = 1:M
                for j = 1:N
                    if r(j,i) < eps
                        S2(j,i) = 2*(pi/(eps^2*((pi^2) - 4))) * (1 + cos((pi/eps)*r(j,i)));
                    end
                end
            end

    %  % ---------------------------------------------------------------------------------------------------------------------------- %


            % We calculate Tx and Ty for a different space-step at each loop
                Tx = kron(Ax,eye(N))/h^2;
                Ty = kron(eye(M),Ay)/h^2;

                I = eye(M*N);

            % Constant matrices 
                A = (I - Dt*(Tx+Ty));
                [L,U] = lu(A);

                Q = zeros(M);     % M = N

            % Reorganisation of matrix Q as a column vector
                QQ = sparse(reshape(Q,[],1));

            % Matrix gathering QQ over time. This matrix will allow us to deal
            % with the heat function over time
                FF = sparse([QQ]);

     % ---------------------------------------------------------------------------------------------------------------------------- %        
        
        
            % S2 - Uncomment all if you're willing to test S2 

            SS = reshape(S2,[],1);

            for t = 1:N_t

                while Dt*t < 0.25
                    Y = L^-1*(QQ + Dt*SS);
                    QQ = U^-1*Y;                % Resolution of heat equation in 2D with source

                    % We ensure that the temperature inside the mesh dosn't go
                    % above the source heat
                    for k = 1:size(QQ,1)
                       if QQ(k,1) > max(SS)
                          QQ(k,1) = max(SS);
                       end
                    end

                    FF = [FF QQ];

                    t = t + 1;
                    
                    if t >= N_t
                        break
                    end
                end

                if t ~= N_t
                    Y = L^-1*QQ;
                    QQ = U^-1*Y;            

                    % We ensure that the temperature inside the mesh dosn't go
                    % above the source heat
                    for k = 1:size(QQ,1)
                        if QQ(k,1) > max(SS)
                          QQ(k,1) = max(SS);
                        end
                    end

                    FF = [FF QQ];
                else
                    break
                end
            end

            FF = FF(:,2:end);
  
            Q_p(p) = FF(5,5);
        end

% Presentation results
 
    % (u_h - u_h/2) / (u_h/2 - u_h/4)
        P1 = (Q_p(1) - Q_p(2))/(Q_p(2) - Q_p(3));
        P2 = (Q_p(2) - Q_p(3))/(Q_p(3) - Q_p(4));
        P3 = (Q_p(3) - Q_p(4))/(Q_p(4) - Q_p(5));
        P4 = (Q_p(4) - Q_p(5))/(Q_p(5) - Q_p(6));
        
    % log_2((u_h - u_h/2)/(u_h/2 - u_h/4)) --> relation sought
        M1 = log(P1)/log(2);
        M2 = log(P2)/log(2);
        M3 = log(P3)/log(2);
        M4 = log(P4)/log(2);
        
        % => The tendency of M gives the power of influence of h

        figure(9)
        
        xx = linspace(10^-4,1,100);
        yy = xx.^2;
        
        for k = 1:size(alpha,2)-1
            loglog(alpha(k), abs(Q_p(k+1) - Q_p(k)), 'r*-', xx, yy, 'b');
            xlabel('\alpha'); ylabel('|Q_{p+1} - Q_{p}|'); title('|Q_{p+1} - Q_{p}| in function of \alpha');
            legend('Q_{p+1} - Q_{p}', 'y = x^2', 'location', 'northwest');
            grid on
            hold on
        end
        
        % ------------------------------------------------------------------------- %
        % As a conclusion we see that to decrease epsilon makes the dirac
        % closer from theory (unit infinity pic), but it as well makes that
        % the condition r < eps isn't fulfilled anymore, and so that the 
        % equation gets no result in the end. 
        % ------------------------------------------------------------------------- %
