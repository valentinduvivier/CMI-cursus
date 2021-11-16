%% Influence of the frequency : ANGLE = 10 & B0

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
            log_P_b0(k,i) = LOG_P(30+k, i); 
        end
    end

% We create a matrix of the form of 'log_p_b0' that contain the
% values [10 20 30]
% This is to then quantify the influence of the angle over the Aero loading
    for k = 1:15
       for i = 1:3
          % Matrix of 15x3
          frequen(k,i) = 10+5*(i-1);
       end
    end
    
% We are using the function log(P) = a*log(Frequency) + b,
% With b the ordinate at the origin and a the slope.
% As we keep the freq cst, it will eventually give the following relation :
%                   P = Cst*(Frequency^a);

% We thus are searching for which sensor this slope is the same to quantify
% definitely and precisely the influence of the angle on the loading

% We calculate the logarithm of the angles
    log_frequence = log10(frequen);

% We display the function log(P) = a*log(Angle) + b in portions
% to visually simplify the observations

    A = 0; % Variable used to get the mean slope
    i = 0;
    
    figure(16)
    subplot(1,3,1)
    for k = 1:5
        plot(log_frequence(k,(1:3)), log_P_b0(k,(1:3)), 'g-*', 'color', [k/7 k/15 k/7]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_frequence(k,3) - log_frequence(k,1));
        i = i + 1;
    end
    hold off
    title('positions 1:5');
    xlabel('log(frequency)'); ylabel('log(P)');
    
    slope_10 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation 
    
    subplot(1,3,2)
    for k = 6:7
        plot(log_frequence(k,(1:3)), log_P_b0(k,(1:3)), 'g-*', 'color', [k/17 k/15 k/17]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_frequence(k,3) - log_frequence(k,1));
        i = i + 1;
    end
    hold off
    title('positions 6:7');
    xlabel('log(frequency)'); ylabel('log(P)');

    slope_11 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation 
    
    subplot(1,3,3)
    for k = 8:15
        plot(log_frequence(k,(1:3)), log_P_b0(k,(1:3)), 'g-*', 'color', [k/17 k/15 k/17]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_frequence(k,3) - log_frequence(k,1));
        i = i + 1;
    end
    hold off
    title('positions 8:15');
    xlabel('log(frequency)'); ylabel('log(P)');

    slope_12 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;
    
    figure(17)
    for k = 1:15
        plot(log_frequence(k,(1:3)), log_P_b0(k,(1:3)), 'g-*', 'color', [k/17 0 k/35]);
        hold on 
    end
    hold off
    title('Function log(P) = a*log(frequency) + b, angle = 10');
    xlabel('log(frequency)'); ylabel('log(P)');
     
    
%% Influence of the frequency : ANGLE = 20 & B0

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
            log_P_b0(k,i) = LOG_P(30+k, 3+i); 
        end
    end

% We create a matrix of the form of 'log_p_b0' that contain the
% values [10 20 30]
% This is to then quantify the influence of the angle over the Aero loading
    for k = 1:15
       for i = 1:3
          % Matrix of 15x3
          frequen(k,i) = 10+5*(i-1);
       end
    end
    
% We are using the function log(P) = a*log(Frequency) + b,
% With b the ordinate at the origin and a the slope.
% As we keep the freq cst, it will eventually give the following relation :
%                   P = Cst*(Frequency^a);

% We thus are searching for which sensor this slope is the same to quantify
% definitely and precisely the influence of the angle on the loading

% We calculate the logarithm of the angles
    log_frequence = log10(frequen);

% We display the function log(P) = a*log(Angle) + b in portions
% to visually simplify the observations

    A = 0; % Variable used to get the mean slope
    i = 0;
    
    figure(18)
    subplot(1,3,1)
    for k = 1:6
        plot(log_frequence(k,(1:3)), log_P_b0(k,(1:3)), 'g-*', 'color', [k/7 k/15 k/7]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_frequence(k,3) - log_frequence(k,1));
        i = i + 1;
    end
    hold off
    title('positions 1:6');
    xlabel('log(frequency)'); ylabel('log(P)');
    
    slope_13 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation 
    
    subplot(1,3,2)
    for k = 7
        plot(log_frequence(k,(1:3)), log_P_b0(k,(1:3)), 'g-*', 'color', [k/17 k/15 k/17]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_frequence(k,3) - log_frequence(k,1));
        i = i + 1;
    end
    hold off
    title('position 7');
    xlabel('log(frequency)'); ylabel('log(P)');

    slope_14 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation 
    
    subplot(1,3,3)
    for k = 8:15
        plot(log_frequence(k,(1:3)), log_P_b0(k,(1:3)), 'g-*', 'color', [k/17 k/15 k/17]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_frequence(k,3) - log_frequence(k,1));
        i = i + 1;
    end
    hold off
    title('positions 8:15');
    xlabel('log(frequency)'); ylabel('log(P)');

    slope_15 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;
    
    figure(19)
    for k = 1:15
        plot(log_frequence(k,(1:3)), log_P_b0(k,(1:3)), 'g-*', 'color', [k/17 0 k/35]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_frequence(k,3) - log_frequence(k,1));
        i = i + 1;
    end
    hold off
    title('Function log(P) = a*log(frequency) + b, angle = 20');
    xlabel('log(frequency)'); ylabel('log(P)');
    
    slope_25 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;
    
%% Influence of the frequency : ANGLE = 30 & B0

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
            log_P_b0(k,i) = LOG_P(30+k, 6+i); 
        end
    end

% We create a matrix of the form of 'log_p_b0' that contain the
% values [10 20 30]
% This is to then quantify the influence of the angle over the Aero loading
    for k = 1:15
       for i = 1:3
          % Matrix of 15x3
          frequen(k,i) = 10+5*(i-1);
       end
    end
    
% We are using the function log(P) = a*log(Frequency) + b,
% With b the ordinate at the origin and a the slope.
% As we keep the freq cst, it will eventually give the following relation :
%                   P = Cst*(Frequency^a);

% We thus are searching for which sensor this slope is the same to quantify
% definitely and precisely the influence of the angle on the loading

% We calculate the logarithm of the angles
    log_frequence = log10(frequen);

% We display the function log(P) = a*log(Angle) + b in portions
% to visually simplify the observations

    A = 0; % Variable used to get the mean slope
    i = 0;
    
    figure(20)
    subplot(1,3,1)
    for k = 1:6
        plot(log_frequence(k,(1:3)), log_P_b0(k,(1:3)), 'g-*', 'color', [k/7 k/15 k/7]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_frequence(k,3) - log_frequence(k,1));
        i = i + 1;
    end
    hold off
    title('positions 1:6');
    xlabel('log(frequency)'); ylabel('log(P)');
    
    slope_16 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation 
    
    subplot(1,3,2)
    for k = 7
        plot(log_frequence(k,(1:3)), log_P_b0(k,(1:3)), 'g-*', 'color', [k/17 k/15 k/17]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_frequence(k,3) - log_frequence(k,1));
        i = i + 1;
    end
    hold off
    title('position 7');
    xlabel('log(frequency)'); ylabel('log(P)');

    slope_17 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;   % We reload the variables Angle and number of iteation 
    
    subplot(1,3,3)
    for k = 8:15
        plot(log_frequence(k,(1:3)), log_P_b0(k,(1:3)), 'g-*', 'color', [k/17 k/15 k/17]);
        hold on 
        A = A + (log_P_b0(k,3) - log_P_b0(k,1))/(log_frequence(k,3) - log_frequence(k,1));
        i = i + 1;
    end
    hold off
    title('positions 8:15');
    xlabel('log(frequency)'); ylabel('log(P)');

    slope_18 = A/i;  % Slope for the four first points : cumulated angle / Nb of iteration 
    A = 0; i = 0;
    
    figure(21)
    for k = 1:15
        plot(log_frequence(k,(1:3)), log_P_b0(k,(1:3)), 'g-*', 'color', [k/17 0 k/35]);
        hold on 
    end
    hold off
    title('Function log(P) = a*log(frequency) + b, angle = 30');
    xlabel('log(frequency)'); ylabel('log(P)');
     
