%% Influence of the angle : DELTA P --> FREQ = 10 & B0

% Example for the blade b0 with : Freq = 10, angle = [10 20 30]

% We take from the above matrix LOG_Delta_P the values corresponding to the
% data we selected (f = 10 and angle = [10 20 30])
    for k = 1:15
        for i = 1:3
            % This matrix will be 15x3 : 15 ligns for each sensors on the blade and 3 columns for each data where f = 20 
            log_delta_p_b0(k,i) = LOG_Delta_P(30+k,1+3*(i-1)); 
        end
    end

% We create a matrix of the form of 'log_delta_p_b0' that contain the
% values [10 20 30]
% This is to then quantify the influence of the angle over the Aero loading
    for k = 1:15
       for i = 1:3
          % Matrix of 15x3
          angulus(k,i) = 10*i;
       end
    end

% We are using the function log(Delta_P/P) = a*log(Angle) + b,
% With b the ordinate at the origin and a the slope.
% As we keep the freq cst, it will eventually give the following relation :
%                   Delta_P/P = Cst*(Angle^a);

% We thus are searching for which sensor this slope is the same to quantify
% definitely and precisely the influence of the angle on the loading

% We calculate the logarithm of the angles
    log_angle = log10(angulus);

% We display the function log(Delta_P/P) = a*log(Angle) + b in portions
% to visually simplify the observations
    figure(22)
    subplot(1,3,1)
    for k = 1:5
        plot(log_angle(k,(1:3)), log_delta_p_b0(k,(1:3)), 'g*-', 'color', [k/7 k/15 k/7]);
        hold on 
    end
    hold off
    title('positions 1:5');
    xlabel('log(angle)'); ylabel('log(\DeltaP/P)');

    subplot(1,3,2)
    for k = 6:8
        plot(log_angle(k,(1:3)), log_delta_p_b0(k,(1:3)), 'g*-', 'color', [k/17 k/15 k/17]);
        hold on 
    end
    hold off
    title('positions 6:8');
    xlabel('log(angle)'); ylabel('log(\DeltaP/P)');

    subplot(1,3,3)
    for k = 9:15
        plot(log_angle(k,(1:3)), log_delta_p_b0(k,(1:3)), 'g*-', 'color', [k/17 k/15 k/17]);
        hold on 
    end
    hold off
    title('positions 9:15');
    xlabel('log(angle)'); ylabel('log(\DeltaP/P)');

    
    figure(23)
    for k = 1:15
        plot(log_angle(k,(1:3)), log_delta_p_b0(k,(1:3)), 'g*-', 'color', [k/17 0 k/35]);
        hold on 
    end
    hold off
    title('Function log(\DeltaP/P) = a*log(angle) + b, frequency = 10');
    xlabel('log(angle)'); ylabel('log(\DeltaP/P)');


%% Influence of the angle : P --> FREQ = 10 & B0

% Concatenation into one vector to simplify the code
    P = [A_a10_f10 A_a10_f15 A_a10_f20 A_a20_f10 A_a20_f15 A_a20_f20 A_a30_f10 A_a30_f15 A_a30_f20];

% Creation of the logarithm of the above ratio
    LOG_P = log10(P);

% Example for the blade b0 with : Freq = 10, angle = [10 20 30]

% We take from the above matrix LOG_Delta_P the values corresponding to the
% data we selected (f = 10 and angle = [10 20 30])
    for k = 1:15
        for i = 1:3
            % This matrix will be 15x3 : 15 ligns for each sensors on the blade and 3 columns for each data where f = 20 
            log_P_b0(k,i) = LOG_P(30+k,1+3*(i-1)); 
        end
    end

% We create a matrix of the form of 'log_p_b0' that contain the
% values [10 20 30]
% This is to then quantify the influence of the angle over the Aero loading
    for k = 1:15
       for i = 1:3
          % Matrix of 15x3
          angulus(k,i) = 10*i;
       end
    end

% We are using the function log(P) = a*log(Angle) + b,
% With b the ordinate at the origin and a the slope.
% As we keep the freq cst, it will eventually give the following relation :
%                   P = Cst*(Angle^a);

% We thus are searching for which sensor this slope is the same to quantify
% definitely and precisely the influence of the angle on the loading

% We calculate the logarithm of the angles
    log_angle = log10(angulus);

% We display the function log(P) = a*log(Angle) + b in portions
% to visually simplify the observations

    A = 0; % Variable used to get the mean slope
    i = 0;

    figure(24)
    subplot(1,3,1)
    for k = 1:5
        plot(log_angle(k,(1:3)), log_P_b0(k,(1:3)), 'g*-', 'color', [k/7 k/15 k/7]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_angle(k,3) - log_angle(k,1));
        i = i + 1;
    end
    hold off
    title('positions 1:5');
    xlabel('log(angle)'); ylabel('log(P)');

    slope_1 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation 
    
    subplot(1,3,2)
    for k = 6:8
        plot(log_angle(k,(1:3)), log_P_b0(k,(1:3)), 'g*-', 'color', [k/17 k/15 k/17]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_angle(k,3) - log_angle(k,1));
        i = i + 1;
    end
    hold off
    title('positions 6:8');
    xlabel('log(angle)'); ylabel('log(P)');

    slope_2 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation 
    
    subplot(1,3,3)
    for k = 9:15
        plot(log_angle(k,(1:3)), log_P_b0(k,(1:3)), 'g*-', 'color', [k/17 k/15 k/17]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_angle(k,3) - log_angle(k,1));
        i = i + 1;
    end
    hold off
    title('positions 9:15');
    xlabel('log(angle)'); ylabel('log(P)');

    slope_3 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation 
    
    figure(25)
    for k = 1:15
        plot(log_angle(k,(1:3)), log_P_b0(k,(1:3)), 'g*-', 'color', [k/17 0 k/35]);
        hold on 
    end
    hold off
    title('Function log(P) = a*log(angle) + b, frequency = 10');
    xlabel('log(angle)'); ylabel('log(P)');


%% Influence of the angle : FREQ = 15 & B0

% Concatenation into one vector to simplify the code
    P = [A_a10_f10 A_a10_f15 A_a10_f20 A_a20_f10 A_a20_f15 A_a20_f20 A_a30_f10 A_a30_f15 A_a30_f20];

% Creation of the logarithm of the above ratio
    LOG_P = log10(P);

% Example for the blade b0 with : Freq = 15, angle = [10 20 30]

% We take from the above matrix LOG_Delta_P the values corresponding to the
% data we selected (f = 15 and angle = [10 20 30])
    for k = 1:15
        for i = 1:3
            % This matrix will be 15x3 : 15 ligns for each sensors on the blade and 3 columns for each data where f = 20 
            log_P_b0(k,i) = LOG_P(30+k,2+3*(i-1)); 
        end
    end

% We create a matrix of the form of 'log_p_b0' that contain the
% values [10 20 30]
% This is to then quantify the influence of the angle over the Aero loading
    for k = 1:15
       for i = 1:3
          % Matrix of 15x3
          angulus(k,i) = 10*i;
       end
    end

% We are using the function log(P) = a*log(Angle) + b,
% With b the ordinate at the origin and a the slope.
% As we keep the freq cst, it will eventually give the following relation :
%                   P = Cst*(Angle^a);

% We thus are searching for which sensor this slope is the same to quantify
% definitely and precisely the influence of the angle on the loading

% We calculate the logarithm of the angles
    log_angle = log10(angulus);

% We display the function log(P) = a*log(Angle) + b in portions
% to visually simplify the observations

    A = 0; i = 0;

    figure(26)
    subplot(1,3,1)
    for k = 1:5
        plot(log_angle(k,(1:3)), log_P_b0(k,(1:3)), 'g*-', 'color', [k/7 k/35 k/7]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_angle(k,3) - log_angle(k,1));
        i = i + 1;
    end
    hold off
    title('positions 1:5');
    xlabel('log(angle)'); ylabel('log(P)');

    slope_4 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation 
    
    subplot(1,3,2)
    for k = 6:9
        plot(log_angle(k,(1:3)), log_P_b0(k,(1:3)), 'g*-', 'color', [k/17 k/35 k/17]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_angle(k,3) - log_angle(k,1));
        i = i + 1;
    end
    hold off
    title('positions 6:9');
    xlabel('log(angle)'); ylabel('log(P)');

    slope_5 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation     
    
    subplot(1,3,3)
    for k = 10:15
        plot(log_angle(k,(1:3)), log_P_b0(k,(1:3)), 'g*-', 'color', [k/17 k/35 k/17]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_angle(k,3) - log_angle(k,1));
        i = i + 1;
    end
    hold off
    title('positions 10:15');
    xlabel('log(angle)'); ylabel('log(P)');

    slope_6 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation 
    
    figure(27)
    for k = 1:15
        plot(log_angle(k,(1:3)), log_P_b0(k,(1:3)), 'g*-', 'color', [k/17 0 k/35]);
        hold on 
    end
    hold off
    title('Function log(P) = a*log(angle) + b, frequency = 15');
    xlabel('log(angle)'); ylabel('log(P)');
 

%% Influence of the angle : FREQ = 20 & B0

% Concatenation into one vector to simplify the code
    P = [A_a10_f10 A_a10_f15 A_a10_f20 A_a20_f10 A_a20_f15 A_a20_f20 A_a30_f10 A_a30_f15 A_a30_f20];

% Creation of the logarithm of the above ratio
    LOG_P = log10(P);

% Example for the blade b0 with : Freq = 20, angle = [10 20 30]

% We take from the above matrix LOG_Delta_P the values corresponding to the
% data we selected (f = 20 and angle = [10 20 30])
    for k = 1:15
        for i = 1:3
            % This matrix will be 15x3 : 15 ligns for each sensors on the blade and 3 columns for each data where f = 20 
            log_P_b0(k,i) = LOG_P(30+k,3+3*(i-1)); 
        end
    end

% We create a matrix of the form of 'log_p_b0' that contain the
% values [10 20 30]
% This is to then quantify the influence of the angle over the Aero loading
    for k = 1:15
       for i = 1:3
          % Matrix of 15x3
          angulus(k,i) = 10*i;
       end
    end

% We are using the function log(P) = a*log(Angle) + b,
% With b the ordinate at the origin and a the slope.
% As we keep the freq cst, it will eventually give the following relation :
%                   P = Cst*(Angle^a);

% We thus are searching for which sensor this slope is the same to quantify
% definitely and precisely the influence of the angle on the loading

% We calculate the logarithm of the angles
    log_angle = log10(angulus);

% We display the function log(P) = a*log(Angle) + b in portions
% to visually simplify the observations

    A = 0; % Variable used to get the mean slope
    i = 0;
    
    figure(28)
    subplot(1,3,1)
    for k = 1:4
        plot(log_angle(k,(1:3)), log_P_b0(k,(1:3)), 'g*-', 'color', [k/7 k/15 k/7]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_angle(k,3) - log_angle(k,1));
        i = i + 1;
    end
    hold off
    title('positions 1:4');
    xlabel('log(angle)'); ylabel('log(P)');
    
    slope_7 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation 
    
    subplot(1,3,2)
    for k = 5:9
        plot(log_angle(k,(1:3)), log_P_b0(k,(1:3)), 'g*-' , 'color', [k/17 k/15 k/17]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_angle(k,3) - log_angle(k,1));
        i = i + 1;
    end
    hold off
    title('positions 5:9');
    xlabel('log(angle)'); ylabel('log(P)');

    slope_8 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation 
    
    subplot(1,3,3)
    for k = 10:15
        plot(log_angle(k,(1:3)), log_P_b0(k,(1:3)), 'g*-' , 'color', [k/17 k/15 k/17]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_angle(k,3) - log_angle(k,1));
        i = i + 1;
    end
    hold off
    title('positions 10:15');
    xlabel('log(angle)'); ylabel('log(P)');

    slope_9 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    
    figure(29)
    for k = 1:15
        plot(log_angle(k,(1:3)), log_P_b0(k,(1:3)), 'g*-', 'color', [k/17 0 k/35]);
        hold on 
    end
    hold off
    title('Function log(P) = a*log(angle) + b, frequency = 20');
    xlabel('log(angle)'); ylabel('log(P)');
     
