%% Q_3.1

clear all

% PART 1 : Calculation of the errors

    % Point considered : 
    X0 = 1/4; Y0 = 1/4;
    
    % We consider a square mesh for now on
    N = 50;     % Nx = Ny = N
    
    L = 1;       % Lx = Ly
    D = L/N;    % Dx = Dy
    Dt = 0.2;
    
    x = D:D:L;
    y = D:D:L;
    
    % We define a few space-step that are gonna be used to get the factor of
    % influence of the space-step over the heat equation (Q_3.2)
    h1 = 1/N; h2 = 1/(2*N); h3 = 1/(4*N);    
    h_1 = [h1, h2, h3];
    h_2 = h_1/2; h_3 = h_1/4; h_4 = h_1/7; h_5 = h_1/15; h_6 = h_1/23; h_7 = h_1/30;

    h = [h_1, h_2, h_3, h_4, h_5, h_6, h_7];

%  % ---------------------------------------------------------------------------------------------------------------------------- %        
    
    % Constant variables
        % Sources
        [X3,Y3] = meshgrid(x,y);
        x_s = 1/2; y_s = 1/2; w = 0.2;

        S_1 = sparse(zeros(N));
        S_2 = sparse(zeros(N));
        
        dirac_3_1 = sparse(zeros(N));
        r = sparse(zeros(N));
        
        S_1 = sparse(exp(-((X3 - x_s).^2 + (Y3 - y_s).^2)/(w^2)));
        r = sqrt(((X3 - x_s).^2 + (Y3 - y_s).^2));

    T = 5;
    
    Uh = zeros(1,size(h,2));
        
    % We calculate Tx and Ty for a different space-step at each loop
    
        % Ax (--> Tx & Ty)
        Ax(1,1) = -1;
        Ax(N,N) = -1;

        for k = 1:N-1
            Ax(k,k+1) =  1;

            Ax(k+1,k) =  1;
        end

        for k = 2:N-1
            Ax(k,k)  = -2;
        end


        % Ay (--> Ty)
        Ay(1,1)   = -1;
        Ay(N,N) = -1;

        for k = 1:N-1
            Ay(k,k+1) =  1;

            Ay(k+1,k) =  1;
        end

        for k = 2:N-1
            Ay(k,k) = -2;
        end

        
%  % ---------------------------------------------------------------------------------------------------------------------------- %

        
    % We loop for the space-step defined upper
    for k = 1:size(h,2)
        
        eps = sqrt(h(k));
        g = 2;
        
        for j = 1:N
            for i = 1:N
                if r(j,i) < eps
                    dirac_3_1(j,i) = (pi/(eps^2*((pi^2) - 4))) * (1 + cos((pi/eps)*r(j,i)));
                end
            end
        end
        
        S_2 = sparse(g*dirac_3_1');
        
%  % ---------------------------------------------------------------------------------------------------------------------------- %
        
        Tx = sparse(kron(Ax,eye(N)))/h(k)^2;
        Ty = sparse(kron(eye(N),Ay))/h(k)^2;

        I = eye(N*N);

        % Constant matrices 
        A = sparse(I - Dt*(Tx+Ty));               
                        
        SS_1 = reshape(S_1,[],1);
        SS_2 = sparse(reshape(S_2,[],1));
        
        B_1 = Dt*(A\SS_1);
        B_2 = Dt*(A\SS_2);
        
%  % ---------------------------------------------------------------------------------------------------------------------------- %        
    
    % Heat equation  
        Q = zeros(N);
        QQ = reshape(Q,[],1);
        FF = sparse(zeros(size(QQ,1),T));
        
    % S_1
    for tps = Dt:Dt:T
        QQ = A\QQ + B_1;  % Resolution of heat equation in 2D with source
        for i = 1:N*N
             if QQ(i,1) > max(max(SS_1))
            	QQ(i,1) = max(max(SS_1));
             end
        end
        FF = [FF QQ];       % We concatenate the value of QQ over time in FF
    end
    
    
    % S_2
%     for tps = Dt:Dt:T
%         if tps < 0.25
%             QQ = A\QQ + B_2;  % Resolution of heat equation in 2D with source
%         else
%             QQ = A\QQ;
%         end
%         FF = [FF QQ];       % We concatenate the value of QQ over time in FF
%     end
%             


%  % ---------------------------------------------------------------------------------------------------------------------------- %        

        Q = reshape(QQ,N,N);
        
% Q 3.2 b)
            % Value at (x,y) = (1/4,1/4); to get the coordinates (i,j)
            % associated and thus the value of the heat Q. 
            
            % This part is the product of hand calculation to then 
            % schematize/generalize the position of any heat value needed
            if (mod(N,2) == 1)
                Uh(k) = (Q((N-1)/4,(N-1)/4) + Q((N+3)/4,(N+3)/4))/2;
            else if (mod(N,4) == 0)
                Uh(k) = Q(N/4,N/4);
                else
                Uh(k) = (Q((N-2)/4,(N-2)/4) + Q((N+2)/4,(N+2)/4))/2;
                end
            end
    end
      
            
%  We plot the error for the different h considered
%%
        % Those calculation are base on the method of the link given on Cansas
        % We calculate here the log_2.
        P1 = (Uh(1) - Uh(2))/(Uh(2) - Uh(3)); P2 = (Uh(4) - Uh(5))/(Uh(5) - Uh(6));
        P3 = (Uh(7) - Uh(8))/(Uh(8) - Uh(9)); P4 = (Uh(10) - Uh(11))/(Uh(11) - Uh(12));
        P5 = (Uh(13) - Uh(14))/(Uh(14) - Uh(15)); P6 = (Uh(16) - Uh(17))/(Uh(17) - Uh(18)); P7 = (Uh(19) - Uh(20))/(Uh(20) - Uh(21));
        
        % We plot the different point to see if there is a possible coefficient
        % defining the error bounded to h
        
 %%
 t = linspace(-10,2,100);
        y = t*0.025 + 0.75;
        figure(2)
        plot((h1), (P1), 'r*');
        hold on
        plot(log(h1/2), log(P2)/2, 'r*');
        plot(log(h1/4), log(P3)/2, 'r*');
        plot(log(h1/7), log(P4)/2, 'r*');

        plot(log(h1/15), log(P5)/2, 'r*');
        plot(log(h1/23), log(P6)/2, 'r*');
        plot(log(h1/30), log(P7)/2, 'r*');
        
        plot(t,y,'g-*'); title('y(t) = 0.025*t'); xlabel('h'); ylabel('log_2');
        axis([-8 -4 0 1]);
        hold off


%  % ---------------------------------------------------------------------------------------------------------------------------- %        
%  % ---------------------------------------------------------------------------------------------------------------------------- %        


%% PART 2 : Influence of epsilon on source S2

        h = 1/30;
        
        % We test our epsilon for different alpha
        alpha = [1 0.1 0.01 0.001 0.0001 10e-10];
        epsilon = alpha*h;
        r = sqrt(((X3 - x_s).^2 + (Y3 - y_s).^2));

        S2 = zeros(N);
        SS = reshape(S2,[],1);
        
        %R = reshape(r,[],1);
        
        %Once again, we creat a matrix gathering the value of an other
        %matrix through time, to then be able do deal with thos values.
        FFF = zeros(size(SS,1),size(epsilon,2));
        
        for t = 1:size(epsilon,2)
            S2 = reshape(SS,N,N);
                for i = 1:N
                    for j = 1:N
                        if r(j,i) < epsilon(t)
                            S2(j,i) = 2*(pi/(epsilon(t)^2*((pi^2) - 4))) * (1 + cos((pi/epsilon(t))*r(j,i)));
                        end
                    end
                end
            SS = reshape(S2',[],1);
            FFF = [FFF SS];
        end
  
        
%%        
% To decomment if you want to display the result of Q_3.2       
     
     subplot(2,3,1);
     surf(X3,Y3,reshape(FFF(:,7),N,N));
     title('eps alpha = 1');
     xlabel('x'); ylabel('y');
     
     subplot(2,3,2);
     surf(X3,Y3,reshape(FFF(:,8),N,N));
     title(['eps alpha = 0.1']);
     xlabel('x'); ylabel('y');
     
     subplot(2,3,3);
     surf(X3,Y3,reshape(FFF(:,9),N,N));
     title(['eps alpha = 0.01']);
     xlabel('x'); ylabel('y');
     
     subplot(2,3,4);
     surf(X3,Y3,reshape(FFF(:,10),N,N));
     title(['eps alpha = 0.001']);
     xlabel('x'); ylabel('y');
     
     subplot(2,3,5);
     surf(X3,Y3,reshape(FFF(:,11),N,N));
     title(['eps alpha = 0.0001']);
     xlabel('x'); ylabel('y');
     
     subplot(2,3,6);
     surf(X3,Y3,reshape(FFF(:,12),N,N));
     title(['eps alpha = 10e-10']);
     xlabel('x'); ylabel('y');
