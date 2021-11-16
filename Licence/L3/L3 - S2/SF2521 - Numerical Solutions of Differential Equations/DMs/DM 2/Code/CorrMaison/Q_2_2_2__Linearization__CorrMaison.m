% We are willing to approximate our non-linear function using a linearized
% function, that would stand as an approximation. We are thus here gonna
% measure and show how well the linearized equation fit the non-linear one


% ATTENTION : run code Q_2_1_1__Wave__CorrMaison.m with Dt = 0.3*Dx in order to run this part of the code in an
% efficient way

        H = 1;

    % Function linearized
        u_linear = zeros(2,N);
    
    % v_c
        v_c = [H; 0];
    
    % Function of variation
        u_tilde = zeros(2,N);
    
    time = 1.7;
        
    for k = 1:N
        u_tilde(:,k) = [exp(-(((x(k) + sqrt(g*H)*time - L/2)/w)^2)) + exp(-((x(k) - sqrt(g*H)*time - L/2)/w)^2); sqrt(g*H)*( -exp(-((x(k) + sqrt(g*H)*time - L/2)/w)^2) + exp(-((x(k) - sqrt(g*H)*time - L/2)/w)^2) )];
        u_linear(:,k) = v_c + epsilon * u_tilde(:,k);
    end
    
    figure(10)
        subplot(1,2,1)
            plot(x, U(2*(round(time/Dt) - 1) + 1,2:end-1), 'r');  % 1/Dt as we have Nt*Dt = 1 (we consider t = 1s)
            axis([0 10 0.8 1.4]);
            grid on
            hold on
            plot(x, u_linear(1,:), 'b-');
            xlabel('x position'); ylabel('Wave''s height'); title('Height t = 1s & C = 0.3 (\Delta_t = 0.3*\Delta_x)');
            hold off
            
        subplot(1,2,2)
            plot(x, U(2*(round(time/Dt) - 1) + 2, 2:end-1)./U(2*(round(time/Dt) - 1) + 1, 2:end-1), 'r');  % 1/Dt as we have Nt*Dt = 1 (we consider t = 1s)
            axis([0 10 -0.2 0.2]);
            grid on
            hold on
            plot(x, u_linear(2,:)./u_linear(1,:), 'b');
            xlabel('x position'); ylabel('Wave''s speed'); title('Speed t = 1s & C = 0.3 (\Delta_t = 0.3*\Delta_x)');
            hold off
            
            
            % --------------------------------------------------------------- %
            % As a conslusion, we have that the CFL conditions implies to
            % define Dt such that we are at the limit of study. We can
            % thereby conclude on both solutions by looking at them before
            % the "critical time". The particularity of the solutions are
            % just volontarily enhanced.
            % In the end, we see that the linear case follows well theory,
            % except for what is up to the sharpness due to physical
            % application (bigger wave --> bigger speed and sharper look).
            % The linear theory reproduces as well the solution for the
            % wave equation, except tat it is physically less reliable.
            % --------------------------------------------------------------- %

            