% We are willing to approximate our non-linear function using a linearized
% function, that would stand as an approximation. We are thus here gonna
% measure and show how well the linearized equation fit the non-linear one

% Initial constant state considered :
    % U1 = h
        h_0 = 1;
        
    % U2 = h*v
        v_0 = 0; % h_0*v_0 = 0 as we choose v_0 = 0
    
    % Vector containig the varible u_tilde appearing in the formule Q_2_2_2
        U_tilde = zeros(size(U,1),size(U,2));
    
    % We create the vector U_tilde for which we apply the formula 
    % U_tilde = (U - v(x))/epsilon, v(x) being the constant vector we put
    % [h_0, v_0] for us
    
    for k = 1:iter_t+1
        U_tilde(2*(k-1) + 1, :) = ( U(2*(k-1) + 1, :) - h_0 )/epsilon;
        U_tilde(2*k,:)          = ( U(2*k,:)          - v_0 )/epsilon;
    end

    
%%
        % Drawing of the curve representing the wave's height for the
        % u_tilde function (part of the linearized function)
        figure(20)
        for k = 1:iter_t+1
            plot(x, U_tilde(2*(k-1)+1, 2:end-1), 'r');
            axis([0 10 0 1.0]);
            xlabel('position on the mesh [0;L](m)'); ylabel('Wave height (m)'); title(['Wave at t = ',num2str(k),'*Delta_t ~~ 1 bouncing']);
            drawnow
        end


%%
        
        % Display of some defined points showing the wave's state
        figure(21)
        subplot(2,2,1)
        plot(x, U_tilde(2*(3-1)+1, 2:end-1), 'r');
        axis([0 10 0 1.0]);
        xlabel('position on the mesh [0;L](m)'); ylabel('Wave height (m)'); title(['Wave at t = ',num2str(3),'*Delta_t (1 bouncing = 150)']);
              
        subplot(2,2,2)
        plot(x, U_tilde(2*(75-1)+1, 2:end-1), 'r');
        axis([0 10 0 1.0]);
        xlabel('position on the mesh [0;L](m)'); ylabel('Wave height (m)'); title(['Wave at t = ',num2str(75),'*Delta_t']);
        
        subplot(2,2,3)
        plot(x, U_tilde(2*(180-1)+1, 2:end-1), 'r');
        axis([0 10 0 1.0]);
        xlabel('position on the mesh [0;L](m)'); ylabel('Wave height (m)'); title(['Wave at t = ',num2str(180),'*Delta_t --> 1 bouncing ']);
        
        subplot(2,2,4)
        plot(x, U_tilde(2*(240-1)+1, 2:end-1), 'r');
        axis([0 10 0 1.0]);
        xlabel('position on the mesh [0;L](m)'); ylabel('Wave height (m)'); title(['Wave at t = ',num2str(240),'*Delta_t ~~ 2 bouncings(= 301)']);
        
