  
    ALPHA = -0.08;
    BETA = 0;
        
    N = 500;
    T = 50;
    
    DELTA_X = 0.1;   
    DELTA_T = 0.01;   

    U = zeros(1,N+2);
    
    X = linspace(0,N,N);
    
    for k = 1:N
        U(1,k+1) = sin(X(k));
    end
    
    UU = [U];
    
    for t = 1:T
    
%         u(1,1) = u(1,2);    u(1,1) = u(T,1);
%                 
%         u(2,1) = -u(2,2);   u(2,n+2) = -u(2,n+1);
        
        for k = 2:N+1
        U(1,k) = U(1,k)*(1 - 2*ALPHA*DELTA_T/DELTA_X^2 + 6*BETA*DELTA_T/DELTA_X^2) ...
            + (U(1,k) + U(1,k-1))*(ALPHA*DELTA_T/DELTA_X^2 - 4*BETA*DELTA_T/DELTA_X^4);% + ...
%             (U(1,k+2) + U(1,k-2))*(BETA*DELTA_T/DELTA_X^4);
        end
        
        UU = [UU; U];
    end
       
    figure(2)
    for k = 1:size(UU,1)
        plot(X, UU(k,2:end-1), 'r');
        axis([0 50 -1.2 1.2]);
        drawnow
    end    

    subplot(2,2,1)
    plot(X, UU(5,2:end-1), 'r');
    xlabel('position [m]'); ylabel('Height [m]'); title('Scheme for the well posed case at t = 1s');
        axis([0 50 -8 8]);
    subplot(2,2,2)
    plot(X, UU(10,2:end-1), 'r');
    xlabel('position [m]'); ylabel('Height [m]'); title('Scheme for the well posed case at t = 2s');
        axis([0 50 -8 8]);
    subplot(2,2,3)
    plot(X, UU(25,2:end-1), 'r');
    xlabel('position [m]'); ylabel('Height [m]'); title('Scheme for the well posed case at t = 5s');
        axis([0 50 -8 8]);
    subplot(2,2,4)
    plot(X, UU(40,2:end-1), 'r');
    xlabel('position [m]'); ylabel('Height [m]'); title('Scheme for the well posed case at t = 8s');
        axis([0 50 -8 8]);
