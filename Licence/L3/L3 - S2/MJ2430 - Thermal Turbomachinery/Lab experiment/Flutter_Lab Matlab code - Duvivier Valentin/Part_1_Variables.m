%% Steady Aerodynamic Loading

% For each angle, there is 3 associated freq
    % Angles considered during the experiment
    angle = [10, 20, 30];

    % Associated frequenceies
    freq = [10, 15, 20];

    % Cordinnates of the sensor on a blade (located throughout the boundary
        coordinates = [
        -0.52       % Point 1
        -0.4        % Point 2
        -0.3        % Point 3
        -0.22       % Point 4
        -0.18       % Point 5
        -0.1        % Point 6
        -0.05       % Point 7
        -0.01       % Point 8
        0           % Point 9
        0.02        % Point 10
        0.08        % Point 11
        0.14        % Point 12
        0.2         % Point 13
        0.3         % Point 14
        0.41];      % Point 15

    
% We take the values of loading from the files received by mail
    for k = 3:77
        % Example : we have here a vector for the loading corresponding to
        % an angle of 30° and a frequence of 20
        A_a30_f20(k-2,1) = Angle30Frequency20SpringXX0108(k,1)*10^-6;

        A_a30_f15(k-2,1) = Angle30Frequency15SpringXX0108(k,1)*10^-6;

        A_a30_f10(k-2,1) = Angle30Frequency10SpringXX0108(k,1)*10^-6;

        A_a20_f20(k-2,1) = Angle20Frequency20SpringXX0108(k,1)*10^-6;

        A_a20_f15(k-2,1) = Angle20Frequency15SpringXX0108(k,1)*10^-6;

        A_a20_f10(k-2,1) = Angle20Frequency10SpringXX0108(k,1)*10^-6;

        A_a10_f20(k-2,1) = Angle10Frequency20SpringXX0108(k,1)*10^-6;

        A_a10_f15(k-2,1) = Angle10Frequency15SpringXX0108(k,1)*10^-6;

        A_a10_f10(k-2,1) = Angle10Frequency10SpringXX0108(k,1)*10^-6;
        
        % Variable for unsteady study
        Flutter_B0(k-2,1) = FlutterAngle10Frequency32Spring530109(k,1)*10^-6;
    end

% ---------------------------------------------------------------------------------------------------------------------- %

% CASE where we study Delta_P/P

% We then see the influence of each term on the aero loading by looking to
% the ratio between the diffrence (loading - P0) with P0

% Each Delta_P = abs(P0 - P[kPa])/P0, and this for each couple (angle,freq)
    for k = 1:75
            Delta_P1(k,1) = abs(Angle10Frequency10SpringXX0108(1,6)*10^-6 - A_a10_f10(k,1))/(Angle10Frequency10SpringXX0108(1,6)*10^-6);
            
            Delta_P2(k,1) = abs(Angle10Frequency15SpringXX0108(1,6)*10^-6 - A_a10_f15(k,1))/(Angle10Frequency15SpringXX0108(1,6)*10^-6);

            Delta_P3(k,1) = abs(Angle10Frequency20SpringXX0108(1,6)*10^-6 - A_a10_f20(k,1))/(Angle10Frequency20SpringXX0108(1,6)*10^-6);

            Delta_P4(k,1) = abs(Angle20Frequency10SpringXX0108(1,6)*10^-6 - A_a20_f10(k,1))/(Angle20Frequency10SpringXX0108(1,6)*10^-6);

            Delta_P5(k,1) = abs(Angle20Frequency15SpringXX0108(1,6)*10^-6 - A_a20_f15(k,1))/(Angle20Frequency15SpringXX0108(1,6)*10^-6);

            Delta_P6(k,1) = abs(Angle20Frequency20SpringXX0108(1,6)*10^-6 - A_a20_f20(k,1))/(Angle20Frequency20SpringXX0108(1,6)*10^-6);

            Delta_P7(k,1) = abs(Angle30Frequency10SpringXX0108(1,6)*10^-6 - A_a30_f10(k,1))/(Angle30Frequency10SpringXX0108(1,6)*10^-6);

            Delta_P8(k,1) = abs((Angle30Frequency15SpringXX0108(1,6)*10^-6) - A_a30_f15(k,1))/(Angle30Frequency15SpringXX0108(1,6)*10^-6);

            Delta_P9(k,1) = abs((Angle30Frequency20SpringXX0108(1,6)*10^-6) - A_a30_f20(k,1))/(Angle30Frequency20SpringXX0108(1,6)*10^-6);
    end

% Concatenation into one vector to simplify the code
    Delta_P = [Delta_P1 Delta_P2 Delta_P3 Delta_P4 Delta_P5 Delta_P6 Delta_P7 Delta_P8 Delta_P9];

% Creation of the logarithm of the above ratio
    LOG_Delta_P = log10(Delta_P);
    

