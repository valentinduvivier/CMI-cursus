        
%% Q_3.3 Numerical conservation
    
    % Integral
            
    % We consider a square mesh for now on
        Lx = 1; Ly = 1;
    
        N = 10;    
        M = N;

        Dx = Lx/M; Dy = Ly/N;   % Space step
        h = Dx;
        
    % Time variables
        T = 2;         % The code will be running for 10s
        N_t = 100;      % Time scale - Nb of iter°

        Dt = T/N_t;     % Time step
        
%     We create point to define the mesh for the dicretization
        x = linspace(0, Lx, M);
        y = linspace(0, Ly, N);

        Q = zeros(M,N);     % As i = 1 --> M and j = 1 --> N


 % ---------------------------------------------------------------------------------------------------------------------------- %        

      
    % Ax
        Ax = zeros(M);

%         Diagonal
%             (See BC)
                Ax(1,1) = -1;
                Ax(M,M) = -1;

            for k = 2:M-1
                Ax(k,k) = -2;
            end

%         Diagonal superior & inferior
            for k = 1:M-1
                Ax(k,k+1) = 1;

                Ax(k+1,k) = 1;
            end
        
%     Ay
        Ay = zeros(N);

%         Diagonal
%             (See BC)
                Ay(1,1) = -1;
                Ay(N,N) = -1;

            for k = 2:N-1
                Ay(k,k) = -2;
            end

%         Diagonal superior & inferior
            for k = 1:N-1
                Ay(k,k+1) = 1;

                Ay(k+1,k) = 1;
            end

 % ---------------------------------------------------------------------------------------------------------------------------- %        
        
  
    % PART 1 - Influence of space step h

    [X3,Y3] = meshgrid(x,y);
    
    Q_p = zeros(1,size(h,2));
    
    % We loop for the space-step defined upper
        xs = 1/2; ys = 1/2;
        w = 1/4;
        
        S1 = exp(-((X3 - xs).^2 + (Y3 - ys).^2)/(w^2));

        r = sqrt(((X3 - xs).^2 + (Y3 - ys).^2));
        eps = sqrt(h);

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


% S1 - Uncomment all if you're willing to test S1   

        % Source vector
        SS = sparse(reshape(S1,[],1));

        Sum = 0;
        SUM = [Sum];
        
        for t = 1:N_t     % We loop for as many itera° as we have 
            Y = L^-1*(QQ + Dt*SS);
            QQ = U^-1*Y;  % Resolution of heat equation in 2D with source
            
            % We ensure that the temperature inside the mesh dosn't go
            % above the source heat
            for k = 1:size(QQ,1)
               if QQ(k,1) > max(SS)
                  QQ(k,1) = max(SS);
               end
            end
            
            % Sum
            Sum = 0;

            for i = 1:M*N
                Sum = Sum + QQ(i,1);
            end
            
        SUM = [SUM, Sum];
            
        FF = [FF QQ];
            
        end
        
        
        SUM = SUM * h^2;
        SUM = SUM(2:end);
        
% S2 - Uncomment all if you're willing to test S2 

%         SS = reshape(S2,[],1);
% 
%         for t = 1:N_t
% 
%             while Dt*t < 0.25
%                 Y = L^-1*(QQ + Dt*SS);
%                 QQ = U^-1*Y;                % Resolution of heat equation in 2D with source
%                 
%                 % We ensure that the temperature inside the mesh dosn't go
%                 % above the source heat
%                 for k = 1:size(QQ,1)
%                    if QQ(k,1) > max(SS)
%                       QQ(k,1) = max(SS);
%                    end
%                 end
%                 
%                 FF = [FF QQ];
% 
%                 t = t + 1;
%             end
%             
%             Y = L^-1*QQ;
%             QQ = U^-1*Y;            
%              
%             % We ensure that the temperature inside the mesh dosn't go
%             % above the source heat
%             for k = 1:size(QQ,1)
%                 if QQ(k,1) > max(SS)
%                   QQ(k,1) = max(SS);
%                 end
%             end
% 
%         FF = [FF QQ];
% 
%         end
        

        FF = FF(:,2:end);

        
        for k = 1:size(FF,2)
           INT(k) = (FF(end,k) + FF(1,k))/2; 
        end
            
%  % ---------------------------------------------------------------------------------------------------------------------------- %        


         % Not sure about the veracity of the results plotted but it looks
         % like the numerical conservation is true up to an error of 10e-3
         ERROR = abs(SUM - INT);
         
         figure(10)
         plot(1:N_t, ERROR, '*-b');     % We note that when we are at N_t (100) we have T = 2s ans that at 1 we have T > 0s
         xlabel('Iteration'); ylabel('Error = |Sum - Int|'); title('Quantification / Verification of numerical conservation'); 
         
         