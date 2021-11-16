
% Elementary information : 
%     - each major question using Matlab has its own section ;
%     - you should still run them in the order of appearance
%     to make shure every variable are computed at least once ;
% 

%% Q_2.5

% Initialisation data

clear all

Lx = 1; Ly = 1;     % Dimensions of the area over which we descretize
M = 30; N = 20;   % Number of point considered for the discretization

Dx = 1/M; Dy = 1/N; % Step for the position

Dt = 0.02;          % Step for the time (any as we are in implicit)
% As the time's step isn't a condition of the convergence, we advice
% to consider TPS as an integer and better change Dt if you feel a
% need to check for other range of time.

TPS = 4;            

% We create point to define the mesh for the dicretization
x = linspace(0,1,M);
y = linspace(0,1,N);

Q = zeros(M,N);     % As i = 1 --> M and j = 1 --> N


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

    % To uncomment if you want to display the source S1
    subplot(1,2,1);
    surf(X1,Y1,S1);
    title('S1 for a mesh of 30-by-30 points, at t < 1/4');
    xlabel('x'); ylabel('y');
    

    % Source S2

    % Origin point for the source (xs,ys) :
    xs_2 = 1/2; ys_2 = 1/2;
         
    [X2,Y2] = meshgrid(x,y);
    
    r = sqrt(((X2 - xs_2).^2 + (Y2 - ys_2).^2));
    eps = sqrt(max(Dx,Dy));
    
    S2 = zeros(N,M);
    
    for j = 1:M
        for i = 1:N
            if r(i,j) < eps
            	S2(i,j) = 2*(pi/(eps^2*((pi^2) - 4))) * (1 + cos((pi/eps).*r(i,j)));
            end
        end
    end
    
	% To uncomment if you want to display the source S2
    subplot(1,2,2);
    surf(X2,Y2,S2);
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
       
         Ax = zeros(M,M);

        % Diagonal
         Ax(1,1) = -1;
         Ax(M,M) = -1;

        for k = 2:M-1
            Ax(k,k) = -2;
        end

        % Diagonal superior and inferior
        for k = 1:M-1
            Ax(k,k+1) = 1;
      
            Ax(k+1,k) = 1;
        end
        
 
        Ay = zeros(N,N);

        % Diagonal
         Ay(1,1) = -1;
         Ay(N,N) = -1;

        for k = 2:N-1
            Ay(k,k) = -2;
        end

        % Diagonal superior
        for k = 1:N-1
            Ay(k,k+1) = 1;
        
            Ay(k+1,k) = 1;
        end

        
%  % ---------------------------------------------------------------------------------------------------------------------------- %        
        
        
% Heat function
        % We apply the Kronecker product to get matrices of the same size
        Tx = kron(Ax,eye(N,N))/Dx^2;
        Ty = kron(eye(M,M),Ay)/Dy^2;

        % We introduce a matrix identity of the size of the one jus above
        I = eye(M*N,M*N);
        
        % Constant matrices
        A = (I - Dt*(Tx+Ty));
        A_inverse = A^-1;
        
        % Reorganisation of matrix Q as a column vector
        QQ = reshape(Q,[],1);
        
        % Matrix gathering QQ over time. This matrix will allow us to deal
        % with the heat function over time
        FF = zeros(size(QQ,1),TPS);
        

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
        S = reshape(SS,M,N);
        Q = reshape(QQ,M,N);

% We are now able to get the heat over the mesh using Q and the heat over
% the mesh and the time using reshape(FF).
   

%% Q_3.1 

% To decoment if you want to display the result of the question Q_3.1

        % S1
%         subplot(1,4,1)
%         HM1 = imagesc(reshape(FF(:,5),M,N)); 
%         title('Heat repartition, t = 0.10s for S1');
%         xlabel('x'); ylabel('y');
%         
%         subplot(1,4,2)
%         HM1 = imagesc(reshape(FF(:,13),M,N)); 
%         title('Heat repartition, t = 0.26s for S1');
%         xlabel('x'); ylabel('y');
%         
%         subplot(1,4,3)
%         HM1 = imagesc(reshape(FF(:,137),M,N)); 
%         title('Heat repartition, t = 3.48s for S1');
%         xlabel('x'); ylabel('y');
%         
%         subplot(1,4,4)
%         HM1 = imagesc(reshape(FF(:,145),M,N)); 
%         title('Heat repartition, t = 3.5s for S1');
%         xlabel('x'); ylabel('y');
%         
        % S2  
        subplot(1,3,1)
        HM1 = imagesc(reshape(FF(:,5),M,N)); 
        title('Heat repartition, t = 0.10s for S2');
        xlabel('x'); ylabel('y');
        
        subplot(1,3,2)
        HM1 = imagesc(reshape(FF(:,13),M,N)); 
        title('Heat repartition, t = 0.26s for S2');
        xlabel('x'); ylabel('y');
        
        subplot(1,3,3)
        HM1 = imagesc(reshape(FF(:,66),M,N)); 
        title('Heat repartition, t = 1.32s for S2');
        xlabel('x'); ylabel('y');
        


%% Q_3.2

% PART 1 : Calculation of the errors

    % Point considered : 
    X0 = 1/4; Y0 = 1/4;
    
    % We consider a square mesh for now on
    M = N;
    % Dx = Dy = h;
    
    % We define a few space-step that are gonna be used to get the factor of
    % influence of the space-step over the heat equation (Q_3.2)
    h1 = 1/M; h2 = (1/2)*h1; h3 = (1/4)*h1;
    h = [h1, h2, h3];
    h_2 = h/2; h_3 = h/4; h_4 = h/7;

    h = [h, h_2, h_3, h_4];

    
%  % ---------------------------------------------------------------------------------------------------------------------------- %        
    

    [X,Y] = meshgrid(x,y);
    
    % We loop for the space-step defined upper
    for k = 1:size(h,2)
        xs = 1/2; ys = 1/2;

        S1 = exp(-((X - xs).^2 + (Y - ys).^2)/(0.25^2));


        r = sqrt(((X - xs).^2 + (Y - ys).^2));
        eps = sqrt(h(k));

        S2 = zeros(M,N);

        for j = 1:M
            for i = 1:N
                if r(i,j) < eps
                    S2(i,j) = 2*(pi/(eps^2*((pi^2) - 4))) * (1 + cos((pi/eps).*r(i,j)));
                end
            end
        end
        
        
%  % ---------------------------------------------------------------------------------------------------------------------------- %
        
        
        % We calculate Tx and Ty for a different space-step at each loop
        Tx = kron(Ax,eye(N,M))/h(k)^2;
        Ty = kron(eye(M,N),Ax)/h(k)^2;

        I = eye(M*N,M*N);

        % Constant matrices 
        A = (I - Dt*(Tx+Ty));
        A_inverse = A^-1;

        Q_1 = zeros(M,N);
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
            for i = 1:M*N
                if QQ(i,1) < max(SS)
                    QQ(i,1) = A_inverse(i,:)*QQ(:,1) + Dt*A_inverse(i,:)*SS(:,1);  % Resolution of heat equation in 2D with source
                else
                    QQ(i,1) = A_inverse(i,:)*QQ(:,1);
                end
            end
            FF = [FF QQ];
        end
        
        S = reshape(SS,M,N);
        Q_1 = reshape(QQ,M,N);

            
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
%         S = reshape(SS,M,N);
%         Q_1 = reshape(QQ,M,N);


%  % ---------------------------------------------------------------------------------------------------------------------------- %        


% Q 3.2 b)
            % Value at (x,y) = (1/4,1/4); to get the coordinates (i,j)
            % associated and thus the value of the heat Q. 
            
            % This part is the product of hand calculation to then 
            % schematize/generalize the position of any heat value needed
            if (mod(N,2) == 1)
                Uh(k) = (Q_1((N-1)/4,(N-1)/4) + Q_1((N+3)/4,(N+3)/4))/2;
            else if (mod(N,4) == 0)
                Uh(k) = Q_1(N/4,N/4);
                else
                Uh(k) = (Q_1((N-2)/4,(N-2)/4) + Q_1((N+2)/4,(N+2)/4))/2;
                end
            end
        
    end
      
            
%  We plot the error for the different h considered

        % Those calculation are base on the method of the link given on Cansas
        % We calculate here the log_2.
        P1 = sqrt((Uh(1) - Uh(2))/(Uh(2) - Uh(3))); P2 = sqrt((Uh(4) - Uh(5))/(Uh(5) - Uh(6)));
        P3 = sqrt((Uh(7) - Uh(8))/(Uh(8) - Uh(9))); P4 = sqrt((Uh(10) - Uh(11))/(Uh(11) - Uh(12)));

        % We plot the different point to see if there is a possible coefficient
        % defining the error bounded to h
 %%
 t = linspace(0,0.02,100);
        y = t*200;
        figure(2)
        plot(h1, P2, 'r*');
        hold on
        plot(h1/2, P1, 'r*');
        plot(h1/4, P3, 'r*');
        plot(h1/7, P4, 'r*');
        plot(t,y,'g-'); title('y(t) = 200*t'); xlabel('h'); ylabel('log_2');
        hold off


%  % ---------------------------------------------------------------------------------------------------------------------------- %        
%  % ---------------------------------------------------------------------------------------------------------------------------- %        


%% PART 2 : Influence of epsilon on source S2

        h = 1/30;
        
        % We test our epsilon for different alpha
        alpha = [1 0.1 0.01 0.001 0.0001 10e-10];
        epsilon = alpha*h;
        r = sqrt(((X - xs).^2 + (Y - ys).^2));

        S2 = zeros(M,N);
        SS = reshape(S2,[],1);
        
        %R = reshape(r,[],1);
        
        %Once again, we creat a matrix gathering the value of an other
        %matrix through time, to then be able do deal with thos values.
        FFF = zeros(size(SS,1),size(epsilon,2));
        
        for t = 1:size(epsilon,2)
            S2 = reshape(SS,M,N);
                for i = 1:M
                    for j = 1:N
                        if r(i,j) < epsilon(t)
                            S2(i,j) = 2*(pi/(epsilon(t)^2*((pi^2) - 4))) * (1 + cos((pi/epsilon(t))*r(i,j)));
                        end
                    end
                end
            SS = reshape(S2,[],1);
            FFF = [FFF SS];
        end
  
        
%%        
% To decomment if you want to display the result of Q_3.2       
     
%      subplot(2,3,1);
%      surf(X,Y,reshape(FFF(:,7),M,N));
%      title(['eps alpha = 1']);
%      xlabel('x'); ylabel('y');
%      
%      subplot(2,3,2);
%      surf(X,Y,reshape(FFF(:,8),M,N));
%      title(['eps alpha = 0.1']);
%      xlabel('x'); ylabel('y');
%      
%      subplot(2,3,3);
%      surf(X,Y,reshape(FFF(:,9),M,N));
%      title(['eps alpha = 0.01']);
%      xlabel('x'); ylabel('y');
%      
%      subplot(2,3,4);
%      surf(X,Y,reshape(FFF(:,10),M,N));
%      title(['eps alpha = 0.001']);
%      xlabel('x'); ylabel('y');
%      
%      subplot(2,3,5);
%      surf(X,Y,reshape(FFF(:,11),M,N));
%      title(['eps alpha = 0.0001']);
%      xlabel('x'); ylabel('y');
%      
%      subplot(2,3,6);
%      surf(X,Y,reshape(FFF(:,12),M,N));
%      title(['eps alpha = 10e-10']);
%      xlabel('x'); ylabel('y');
     
        
%% Q_3.3 Numerical conservation
    
    % Sum
    Sum = 0;

    for i = 1:M
        for j = 1:N
            Sum = Sum + Q(i,j);
        end
    end

    Sum = (h^2)*Sum;
    
    % Integral
    QQ = reshape(Q,[],1);
    
    Integral = 0;
    
%     for k = 1:size(QQ,1)-1
%         Integral = Integral + (h/2)*(QQ(k+1,1) + QQ(k,1));
%     end
    
Integral = (1/2)*(Q(1,1) + Q(M,N));
    