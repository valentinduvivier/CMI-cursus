clear all

%% Q_2.5

% Mesh
    % lengths
    Lx = 1; Ly = 1;
    
    % Number of point
    Nx = 30;    Ny = 30;
    
    % Space-step
    Dx = Lx/Nx;   Dy = Ly/Ny;
    
    % Time-step
    Dt = 0.2;
 
    % Coordinates and variables fo sources
    x_s = Lx/2;  y_s = Ly/2;
    w = 0.2;
        
% Functions    
    % Heat function
    Q = sparse(zeros(Nx,Ny));
    
    % Sources
    S_1 = zeros(Nx,Ny);
    S_2 = sparse(zeros(Nx,Ny));
    
    x = Dx:Dx:Lx;     y = Dy:Dy:Ly;
    
    [X,Y] = meshgrid(x,y);
    
    % Source S1
        
        S_1 = exp(-((X - x_s).^2 + (Y - y_s).^2)/w^2);

        % Displaying of the source S1 in 3D
        figure(1);
        subplot(1,2,1);
        surf(X,Y,S_1);
        title('Source S1'); xlabel('x'); ylabel('y'); zlabel('Intensity of S1');
    
    % Source S2
        g = 2; % g = 2 for t < 1/4
        
        eps = sqrt(max(Dx,Dy));
        
        [X2,Y2] = meshgrid(x,y);       
        
        % Dirac function
        r = sqrt((X2 - x_s).^2 + (Y2 - y_s).^2);
        dirac_e = sparse(zeros(Ny,Nx));

        K = (pi/(eps^2*(pi^2 - 4))); % Constant appearing in S2
        
        for k = 1:Nx
            for l = 1:Ny
                if r(l,k) < eps
                    dirac_e(l,k) = K * (1 + cos((pi/eps) * r(l,k)));
                end
            end
        end
        
        S_2 = sparse(g*dirac_e);

        % Display of the source S2 in 3D
        subplot(1,2,2);
        surf(X2,Y2,S_2);
        title('Source S2'); xlabel('x'); ylabel('y'); zlabel('Intensity S2');
        
%% Second order difference operator's associated matrix

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
       

%% System resolution

    Tx = sparse(kron(Ax,eye(Ny,Ny)))/Dx^2;
    Ty = sparse(kron(eye(Nx,Nx),Ay))/Dy^2;

    I = speye(Ny*Nx,Nx*Ny);
    
    A = sparse(I - Dt*(Tx + Ty));
    
    T = 10; % number of iteration
    
    % Reshape
    QQ = reshape(Q,[],1);
    
    SS_1 = reshape(S_1,[],1);
    SS_2 = reshape(S_2,[],1);
    
    FF = zeros(size(QQ,1),T);
    
    B_1 = Dt*(A\SS_1);
    B_2 = Dt*(A\SS_2);
   
    
    
    % S_1
%     for tps = Dt:Dt:T
%         QQ = A\QQ + B_1;  % Resolution of heat equation in 2D with source
%         for i = 1:Nx*Ny
%              if QQ(i,1) > max(max(SS_1))
%             	QQ(i,1) = max(max(SS_1));
%              end
%         end
%         FF = [FF QQ];       % We concatenate the value of QQ over time in FF
%     end
    
    
    % S_2
    for tps = Dt:Dt:T
        if tps < 0.25
            QQ = A\QQ + B_2;  % Resolution of heat equation in 2D with source
        else
            QQ = A\QQ;
        end
        FF = sparse([FF QQ]);       % We concatenate the value of QQ over time in FF
    end
    
    % Used for Q_3.3
    Q = reshape(QQ,Nx,Ny);
    
%% 2D not that relevant

%         figure(2);
%         
%         subplot(2,2,1);
%         imagesc(reshape(FF(:,T+1),Nx,Ny)); 
%         title('Heat repartition, t = 0.10s for S2');
%         xlabel('x'); ylabel('y');    
%     
%         subplot(2,2,2);
%         imagesc(reshape(FF(:,(T+15)),Nx,Ny)); 
%         title('Heat repartition, t = 0.10s for S2');
%         xlabel('x'); ylabel('y');   
%         
%         subplot(2,2,3);
%         imagesc(reshape(FF(:,T+19),Nx,Ny)); 
%         title('Heat repartition, t = 0.10s for S2');
%         xlabel('x'); ylabel('y');   
%         
%         subplot(2,2,4);
%         imagesc(reshape(FF(:,T+500),Nx,Ny));
%         title('Heat repartition, t = 0.10s for S2');
%         xlabel('x'); ylabel('y');   


%% 3D representation
    Nb_graphique = 5;
    
    figure(3);
    for k = 1:Nb_graphique
        subplot(2,3,k);
        surf(X,Y,reshape(FF(:,(T+2*k)),Nx,Ny));
%         subplot(2,3,2);
%         surf(X,Y,reshape(FF(:,(T+13)),Nx,Ny));
%         subplot(2,3,3);
%         surf(X,Y,reshape(FF(:,(T+15)),Nx,Ny));
%         subplot(2,3,4);
%         surf(X,Y,reshape(FF(:,(T+19)),Nx,Ny));
%         subplot(2,3,5);
%         surf(X,Y,reshape(FF(:,(T+80)),Nx,Ny));
    end
    
    
%% Q_3.3 Numerical conservation
    
    % Sum
    Sum = 0;

    for i = 1:Nx
        for j = 1:Ny
            Sum = Sum + Q(i,j);
        end
    end

    Sum = (Dx*Dy)*Sum;
    
    % Integral  
    Integral = 0;
    
%     for k = 1:size(QQ,1)-1
%         Integral = Integral + (h/2)*(QQ(k+1,1) + QQ(k,1));
%     end
    
Integral = (1/2)*(Q(1,1) + Q(Nx,Ny));         
    