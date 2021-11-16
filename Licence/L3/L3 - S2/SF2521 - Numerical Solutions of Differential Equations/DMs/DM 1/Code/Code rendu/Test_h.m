%%
Lx = 1; Ly = 1;     % Dimensions of the area over which we descretize
Nx = 30; Ny = 30;   % Number of point considered for the discretization

Dx = 1/Nx; Dy = 1/Ny; % Step for the position

Dt = 0.2;          % Step for the time (any as we are in implicit)
% As the time's step isn't a condition of the convergence, we advice
% to consider TPS as an integer and better change Dt if you feel a
% need to check for other range of time.

TPS = 40;            

% We create point to define the mesh for the dicretization
x = linspace(0,1,Nx);
y = linspace(0,1,Ny);

Q = zeros(Nx,Ny);     % As i = 1 --> M and j = 1 --> N


%  % ---------------------------------------------------------------------------------------------------------------------------- %        


% Boundary condition : Nabla Q * n = 0 :

    % The boundar condition (no flux oat the boundaries) is sumed up in the
    % matrix A as the no flux as an impact on the second derivative
    % difference operator's matrix

    
%  % ---------------------------------------------------------------------------------------------------------------------------- %         
    

% Drawing Sources

    % Source S1
    
    % Origin point for the source (xs,ys) :
    xs_1 = 1/2; ys_1 = 1/2;

    [X1,Y1] = meshgrid(x,y);
    S1 = exp(-((X1 - xs_1).^2 + (Y1 - ys_1).^2)/(0.25^2));

    figure(1)
    subplot(1,2,1);
    surf(X1,Y1,S1);
    shading interp
    title('S1 for a mesh of 30-by-30 points, at t < 1/4');
    xlabel('x'); ylabel('y');
    

    % Source S2

    % Origin point for the source (xs,ys) :
    xs_2 = 1/2; ys_2 = 1/2;
         
    [X2,Y2] = meshgrid(x,y);
    
    r = sqrt(((X2 - xs_2).^2 + (Y2 - ys_2).^2));
    eps = sqrt(max(Dx,Dy));
    
    S2 = zeros(Ny,Nx);
    
    for j = 1:Nx
        for i = 1:Ny
            if r(i,j) < eps
            	S2(i,j) = 2*(pi/(eps^2*((pi^2) - 4))) * (1 + cos((pi/eps).*r(i,j)));
            end
        end
    end
    
    figure(1)
    subplot(1,2,2);
    surf(X2,Y2,S2);
    shading interp
    title('S2 for a mesh of 30-by-30 points, at t < 1/4');
    xlabel('x'); ylabel('y');


%  % ---------------------------------------------------------------------------------------------------------------------------- %        


% Creation matrix Tx/Ty

%     % Matrix A (leads to Tx or Ty)
%     % We consider two different A :
%     %   - one of the size of M
%     %   - one of the size of N
%     
%     % Thus, we will be able to consider every cases :
%     %   - the ones where there is a different number of column and line
%     %   - the particular case where we will have as much lines as columns
%     
%     

%  % ---------------------------------------------------------------------------------------------------------------------------- %        
%      
       
         Ax = zeros(Nx,Nx);

        % Diagonal
         Ax(1,1) = -1;
         Ax(Nx,Nx) = -1;

        for k = 2:Nx-1
            Ax(k,k) = -2;
        end

        % Diagonal superior and inferior
        for k = 1:Nx-1
            Ax(k,k+1) = 1;
      
            Ax(k+1,k) = 1;
        end
        
 
        Ay = zeros(Ny,Ny);

        % Diagonal
         Ay(1,1) = -1;
         Ay(Ny,Ny) = -1;

        for k = 2:Ny-1
            Ay(k,k) = -2;
        end

        % Diagonal superior
        for k = 1:Ny-1
            Ay(k,k+1) = 1;
        
            Ay(k+1,k) = 1;
        end

        
%  % ---------------------------------------------------------------------------------------------------------------------------- %        
        
        
% Heat function
        % We apply the Kronecker product to get matrices of the same size
        Tx = kron(Ax,eye(Ny,Ny))/Dx^2;
        Ty = kron(eye(Nx,Nx),Ay)/Dy^2;

        % We introduce a matrix identity of the size of the one jus above
        I = eye(Nx*Ny,Nx*Ny);
        
        % Constant matrices
        A = (I - Dt*(Tx+Ty));
        A_inverse = A^-1;
         
        % Reorganisation of matrix Q as a column vector
        QQ = reshape(Q,[],1);
        
        % Matrix gathering QQ over time. This matrix will allow us to deal
        % with the heat function over time
        FF = zeros(size(QQ,1),TPS);
        FF = [QQ];

% S1 - Uncomment all if you're willing to test S1   

%         % Source vector
%         SS = reshape(S1,[],1);
%         
%         for tps = Dt:Dt:TPS     % We loop for different time 
%             for i = 1:(M*N)
%                 if QQ(i,1) < max(SS)    % If the heat function is not at the maximum
%                    QQ(i,1) = A_inverse(i,:)*QQ(:,1) + Dt*A_inverse(i,:)*SS(:,1);  % Resolution of heat equation in 2D with source
%                 else        % otherwise we have the heat function at her maximum, and thus we stop implementing the heat. 
%                             % We only have the heat that will spread to its
%                             % neighbors.
%                     QQ(i,1) = max(SS);
%                 end
%             end
%             FF = [FF QQ];       % We concatenate the value of QQ over time in FF
%         end
        
        
% S2 - Uncomment all if you're willing to test S2 
% 
        SS = reshape(S2,[],1);

        for tps = Dt:Dt:TPS
            if tps < 0.25
               QQ = A_inverse*QQ + Dt*A_inverse*SS;  % Resolution of heat equation in 2D with source
            else
                QQ = A_inverse*QQ;
            end
            FF = [FF QQ];       % We concatenate the value of QQ over time in FF
        end
        

%  % ---------------------------------------------------------------------------------------------------------------------------- %        


% Display of the results
        % We reshape once again to have the results under matrix form : 
        S = reshape(SS,Nx,Ny);
        Q = reshape(QQ,Nx,Ny);
%%
        Ny = Nx;

    % Point considered : 
        X0 = 1/4; Y0 = 1/4;
    
    % We define a few space-step that are gonna be used to get the factor of
    % influence of the space-step over the heat equation (Q_3.2)
        h1 = 1; h2 = (1/2)*h1; h3 = (1/2)*h2; h4 = (1/2)*h3; h5 = (1/2)*h4;
        h6 = (1/2)*h5; h7 = (1/2)*h6; h8 = (1/2)*h7; h9 = (1/2)*h8;
        
        h10 = (1/2)*h9; h11 = (1/2)*h10; h12 = (1/2)*h11; h13 = (1/2)*h12; h14 = (1/2)*h13;
        h15 = (1/2)*h14; h16 = (1/2)*h15; h17 = (1/2)*h16; h18 = (1/2)*h17;
       
        h19 = (1/2)*h18; h20 = (1/2)*h19; h21 = (1/2)*h20; h22 = (1/2)*h21; h23 = (1/2)*h22;
        h24 = (1/2)*h23; h25 = (1/2)*h24; h26 = (1/2)*h25; h27 = (1/2)*h26;        

        h = [h1, h2, h3, h4, h5, h6, h7, h8, h9, ...
            h10, h11, h12, h13, h14, h15, h16, h17, h18, ...
            h19, h20, h21, h22, h23, h24, h25, h26, h27];
    
%  % ---------------------------------------------------------------------------------------------------------------------------- %        

    [X,Y] = meshgrid(x,y);
    
    Ax = sparse(zeros(Nx,Nx));
    Ay = sparse(zeros(Ny,Ny));
    
    % Ax (--> Tx)
    Ax(1,1)   = -1;
    Ax(Nx,Nx) = -1;
  
    for k = 1:Nx-1
       Ax(k,k+1) =  1;
       
       Ax(k+1,k) =  1;       
    end
    
    for k = 2:Nx-1
        Ax(k,k)  = -2;
    end
    
    % Ay (--> Ty)
    Ay(1,1)   = -1;
    Ay(Ny,Ny) = -1;
  
    for k = 1:Ny-1
       Ay(k,k+1) =  1;
       
       Ay(k+1,k) =  1;       
    end
    
    for k = 2:Ny-1
        Ay(k,k)   = -2;
    end
    
    % We loop for the space-step defined upper
    for k = 1:size(h,2)
        xs = 1/2; ys = 1/2;

        S1 = exp(-((X - xs).^2 + (Y - ys).^2)/(0.25^2));

        r = sqrt(((X - xs).^2 + (Y - ys).^2));
        eps = sqrt(h(k));

        S2 = zeros(Nx,Ny);

        for j = 1:Nx
            for i = 1:Ny
                if r(i,j) < eps
                    S2(i,j) = 2*(pi/(eps^2*((pi^2) - 4))) * (1 + cos((pi/eps).*r(i,j)));
                end
            end
        end
        
%  % ---------------------------------------------------------------------------------------------------------------------------- %
        
        % We calculate Tx and Ty for a different space-step at each loop
            Tx = kron(Ax,eye(Ny,Nx))/h(k)^2;
            Ty = kron(eye(Nx,Ny),Ax)/h(k)^2;

            I = eye(Nx*Ny,Nx*Ny);

        % Constant matrices 
            A = (I - Dt*(Tx+Ty));
            A_inverse = A^-1;

            Q_1 = zeros(Nx,Ny);
            
        % Reorganisation of matrix Q as a column vector
            QQ = reshape(Q_1,[],1);

        % Matrix gathering QQ over time
            FF = zeros(size(QQ,1),TPS);
        
%  % ---------------------------------------------------------------------------------------------------------------------------- %        

    % S1 - Uncomment all if you're willing to test S1   

        % Source vector
            SS = reshape(S1,[],1);
        
        % Same loop as seen before
            for tps = Dt:Dt:TPS
                for i = 1:Nx*Ny
                    if QQ(i,1) < max(SS)
                        QQ(i,1) = A_inverse(i,:)*QQ(:,1) + Dt*A_inverse(i,:)*SS(:,1);  % Resolution of heat equation in 2D with source
                    else
                        QQ(i,1) = A_inverse(i,:)*QQ(:,1);
                    end
                end
                FF = [FF QQ];
            end

            S = reshape(SS,Nx,Ny);
            Q_1 = reshape(QQ,Nx,Ny);

            
    % S2 - Uncomment all if you're willing to test S2            
            
          % Source vector
%         SS = reshape(S2,[],1);
% 
%         for tps = Dt:Dt:TPS
%             if tps < 0.25
%                QQ = A_inverse*QQ + Dt*A_inverse*SS;  % Resolution of heat equation in 2D with source
%             else
%                 QQ = A_inverse*QQ;
%             end
%                 
%             FF = [FF QQ];
%         end
% 
%         S = reshape(SS,Nx,Ny);
%         Q_1 = reshape(QQ,Nx,Ny);


%  % ---------------------------------------------------------------------------------------------------------------------------- %        


% Q 3.2 b)
            if (mod(Ny,2) == 1)
                Uh(k) = (Q_1((Ny-1)/4,(Ny-1)/4) + Q_1((Ny+3)/4,(Ny+3)/4))/2;
            else if (mod(Ny,4) == 0)
                Uh(k) = Q_1(Ny/4,Ny/4);
                else
                Uh(k) = (Q_1((Ny-2)/4,(Ny-2)/4) + Q_1((Ny+2)/4,(Ny+2)/4))/2;
                end
            end
            
            Uh_2(k) = Q_1(7,7);
    end
      
            
%  We plot the error for the different h considered
%%
        % Those calculation are base on the method of the link given on Cansas
        % We calculate here the log_2.
        P1 = ((Uh(1) - Uh(2))/(Uh(2) - Uh(3))); P2 = ((Uh(4) - Uh(5))/(Uh(5) - Uh(6)));
        P3 = ((Uh(7) - Uh(8))/(Uh(8) - Uh(9))); P4 = ((Uh(10) - Uh(11))/(Uh(11) - Uh(12)));
        P5 = ((Uh(13) - Uh(14))/(Uh(14) - Uh(15))); P6 = ((Uh(16) - Uh(17))/(Uh(17) - Uh(18)));
        P7 = ((Uh(19) - Uh(20))/(Uh(20) - Uh(21))); P8 = ((Uh(22) - Uh(23))/(Uh(23) - Uh(24)));
        P9 = ((Uh(25) - Uh(26))/(Uh(26) - Uh(27)));
        % We plot the different point to see if there is a possible coefficient
        % defining the error bounded to h
        

%%
        
        M1 = log10(P1)/log10(2);
        M2 = log10(P2)/log10(2);
        M3 = log10(P3)/log10(2);
        M4 = log10(P4)/log10(2);
        M5 = log10(P5)/log10(2);
        M6 = log10(P6)/log10(2);
        M7 = log10(P7)/log10(2);
        M8 = log10(P8)/log10(2);
        M9 = log10(P9)/log10(2);

%%
        x = -5:1;
        y = 1.2*x-3;

        figure(2)
        plot(h(1), log10(P1)/log10(2), 'r*');
        hold on
        plot(h(4), log10(P2)/log10(2), 'g*');
        plot((h(7)), log10(P3)/log10(2), 'c*');
        plot((h(10)), log10(P4)/log10(2), 'm*');
        plot((h(13)), log10(P5)/log10(2), 'b*');
        plot((h(16)), log10(P6)/log10(2), 'y*');
        plot((h(19)), log10(P7)/log10(2), 'r*');
        plot((h(22)), log10(P8)/log10(2), 'g*');
        plot((h(25)), log10(P9)/log10(2), 'c*');

        hold off

        figure(3)
        plot(log10(h(1)), log10(abs(Uh(1) - Uh(2))), 'r*');
        hold on
        plot(h(4), log10(abs(Uh(4) - Uh(5))), 'g*');
        plot(log10(h(7)), log10(abs(Uh(7) - Uh(8))), 'c*');
        plot(log10(h(10)), log10(abs(Uh(10) - Uh(11))), 'm*');
        plot(log10(h(13)), log10(abs(Uh(13) - Uh(14))), 'b*');
        plot(log10(h(16)), log10(abs(Uh(16) - Uh(17))), 'r*')
        plot(log10(h(19)), log10(abs(Uh(19) - Uh(20))), 'b*')
        plot(log10(h(22)), log10(abs(Uh(22) - Uh(23))), 'y*')
        plot(log10(h(25)), log10(abs(Uh(25) - Uh(26))), 'r*')

        plot(x, y, 'g-');
        hold off
        
%  % ---------------------------------------------------------------------------------------------------------------------------- %        
%  % ---------------------------------------------------------------------------------------------------------------------------- %        
