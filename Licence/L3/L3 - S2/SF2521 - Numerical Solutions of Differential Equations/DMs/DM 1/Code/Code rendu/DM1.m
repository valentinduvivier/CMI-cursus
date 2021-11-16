%% Q 2.5

clear all

Lx = 1; Ly = 1;     % Dimensions of the area over which we descretize
M = 40; N = 40;       % Number of point considered for the discretization

Dx = 1/M; Dy = 1/N; % Pas pour la position
Dt = 0.01;             % Pas pour le temps (qlconque)

% Number of iteration for the time. As the time's step isn't a condition of
% the convergence, we advice to consider TPS as an integer and better
% change Dt if you feel a need to check for other range of time.


TPS = 2;            

x = linspace(0,1,M);
y = linspace(0,1,N);

Q = zeros(M,N);     % As i = 1 --> M and j = 1 --> N

% Boundary condition : Nabla Q * n = 0 :

    % The boundar condition (no flux oat the boundaries) is sumed up in the heat function part
 
    
% Drawing of the source S1

    % Point d'origine de la source (xs,ys) :
    xs_1 = 1/2; ys_1 = 1/2;

    [X1,Y1] = meshgrid(x,y);
    S1 = exp(-((X1 - xs_1).^2 + (Y1 - ys_1).^2)/(0.2^2));

    % To uncomment if you want to display the source S1
%     figure(1)
%     surf(X1,Y1,S1);
%     title('S1');
    
% Drawing of the source S2

    % Point d'origine de la source (xs,ys) :
    xs_2 = 1/2; ys_2 = 1/2;
         
    [X2,Y2] = meshgrid(x,y);
    
    r = sqrt(((X2 - xs_2).^2 + (Y2 - ys_2).^2));
    eps = sqrt(max(Dx,Dy));
    
    S2 = zeros(M,N)';
    
    for j = 1:M
        for i = 1:N
            if r(i,j) < eps
            	S2(i,j) = 2*(pi/(eps^2*((pi^2) - 4))) * (1 + cos((pi/eps).*r(i,j)));
            end
        end
    end
    
	% To uncomment if you want to display the source S1
%     figure(2)
%     surf(X2,Y2,S2);
%     title('S2');
    

% Heat function
    % Matrix D (leads to Tx or Ty)
    % We consider two different D :
    %   - one of the size of M
    %   - one of the size of N
    
    % Thus, we will be able to consider every cases :
    %   - the one were there is a different number of column and line
    %   - the particular case where we will have as much lines as columns
    
        
 % ---------------------------------------------------------------------------------------------------------------------------- %        
     
        Ax = zeros(M,M);

        % Diagonal
        Ax(1,1) = -1;
        Ax(M,M) = -1;

        for k = 1:M
            Ax(k,k) = -2;
        end

        % Diagonal superior and inferior
        for k = 1:M-1
            Ax(k,k+1) = 1;
      
            Ax(k+1,k) = 1;
        end
        
 % ---------------------------------------------------------------------------------------------------------------------------- %
 
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


        Tx = kron(Ax,eye(N,N))/Dx^2;
        Ty = kron(eye(M,M),Ay)/Dy^2;

        I = eye(M*N,M*N);
        
        A = (I - Dt*(Tx+Ty));
        [L,U] = lu(A);
        A_inverse = A^-1;
        
        QQ = reshape(Q,[],1);
        SS = reshape(S1,[],1);
        
        FF = zeros(size(QQ,1),TPS);

        for tps = 0:Dt:TPS
            QQ = A_inverse*QQ + Dt*A_inverse*SS;
            FF = [FF QQ];
        end
         
        FF(:,2);
        
        % We reshape once again to have the results under matrix form : 
        S = reshape(SS,M,N);
        Q = reshape(QQ,M,N);
        
%         for k = 1:size(FF,2)
%             HM1 = heatmap(reshape(FF(:,k),M,N),'GridVisible','off')
%         end

        for k = 1:size(FF,2)
            HM1 = imagesc(reshape(FF(:,k),M,N));
        end
            
        
        

