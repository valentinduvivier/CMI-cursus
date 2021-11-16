% Geometry parameters
    % Origin: middle point of leading edge of inner wing
    % x axis: from origin to tail
    % y axis: from origin to right
    % z axis: from origin to up 

clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Inner Wing         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % basic given parameters:
        TR_w1 = 0.84;

    % input variables:
        c_w1r = 0.85;           % c_wing_right = c_w1r
        b_w1  = 0.73;           % b_wing_1 = b_w1

    % output geometry parameters:
        b_w1half       = b_w1/2;                                           % b_wing_1_half = b_w1_half
        c_w1t          = c_w1r*TR_w1;
        S_w1           = (c_w1r + c_w1t)*b_w1half;
        swept_w1       = atand((0.25*c_w1r - 0.25*c_w1t)/b_w1half);
        AR_wing_1      = (b_w1^2)/S_w1;                                    % Aspect Ratio

        Aero_mean_chord = 2/3*c_w1r*(1 + TR_w1 + TR_w1^2)/(1 + TR_w1);     % Aerodyanic mean chord
        
        y_w1ac         = b_w1half * (1-(Aero_mean_chord - c_w1t)/(c_w1r - c_w1t));             % y position of aerodynamic center of inner right half wing
        x_w1ac         = y_w1ac * tand(swept_w1) + Aero_mean_chord/4;

        X_w1 = [0, 0, b_w1half, b_w1half, 0, 0, -b_w1half, -b_w1half, 0];
        Y_w1 = [0, -c_w1r, -b_w1half*tand(swept_w1) - c_w1t, -b_w1half*tand(swept_w1)...
            , 0, -c_w1r, -b_w1half*tand(swept_w1) - c_w1t, -b_w1half*tand(swept_w1), 0];
        Z_w1 = zeros(9);


        % We begin the plot for the structure
        % This drawing is to be used for geometric and dynamic valuess
        figure(1)
        grid on
        plot3(X_w1,Y_w1,Z_w1,'*')
        hold on
        plot3(X_w1,Y_w1,Z_w1);
        axis equal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Outer Wing         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % basic given parameters:
        swept_w2 = -5;       % swept angle of leading edge
        TR_w2    = 0.35;
        CR       = 0.7;      % ratio of root chord of outer wing to tip chord of inner wing 

    % input variables:
        b_w = 1.7;

    % output geometry parameters:
        c_w2r       = CR*c_w1t;
        c_w2t       = c_w2r*TR_w2;
        b_w2half    = (b_w - b_w1)/2;

        S_w2        = (c_w2r + c_w2t)*b_w2half;    % total area of double outer wings
        MAC_w2      = 2/3*c_w2r*(1 + TR_w2 + TR_w2^2)/(1 + TR_w2);

        y_w2ac      = b_w2half*(1-(MAC_w2 - c_w2t)/(c_w2r - c_w2t)) + b_w1half;            % y position of aerodynamic center of inner right half wing
        x_w2ac      = (y_w2ac - b_w1half)*tand(swept_w2) + MAC_w2/4 + b_w1half*tand(swept_w1);

        X_w2r       = [b_w1half, b_w1half, b_w1half+b_w2half, b_w1half+b_w2half, b_w1half];
        Y_w2r       = [-b_w1half*tand(swept_w1), -b_w1half*tand(swept_w1) - c_w2r, ...
                     -b_w1half*tand(swept_w1) - b_w2half*tand(swept_w2) - c_w2t, ...
                     -b_w1half*tand(swept_w1) - b_w2half*tand(swept_w2), -b_w1half*tand(swept_w1)];
        Z_w2r       = zeros(5);

        
    % Plot for the outer wing
        plot3(X_w2r, Y_w2r, Z_w2r, '*');
        plot3(X_w2r, Y_w2r, Z_w2r);

        X_w2l = [-b_w1half, -b_w1half, -b_w1half-b_w2half, -b_w1half-b_w2half, -b_w1half];
        Y_w2l = [-b_w1half*tand(swept_w1), -b_w1half*tand(swept_w1) - c_w2r ...
                , -b_w1half*tand(swept_w1) - b_w2half*tand(swept_w2) - c_w2t ...
                , -b_w1half*tand(swept_w1) - b_w2half*tand(swept_w2), -b_w1half*tand(swept_w1)];
        Z_w2l = zeros(5);

        plot3(X_w2l, Y_w2l, Z_w2l,'*');
        plot3(X_w2l, Y_w2l, Z_w2l)



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         Whole Wing         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        S_w   = S_w1 + S_w2;
        MAC_w = (Aero_mean_chord*S_w1 + MAC_w2*S_w2)/(S_w1 + S_w2);     % Mach speed of the wing
        x_wac = (x_w1ac*S_w1 + x_w2ac*S_w2)/(S_w1 + S_w2);
        AR_w  = (b_w^2)/S_w;                                            % Aspect ratio for the whole wing

        x_cg  = 0.25*c_w1r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Vertical Tail       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % basic given parameters:
        TR_v    = 0.75;
        V_v     = 0.1;
        l_v     = 1.45*c_w1r;
        b_v     = 0.42*c_w1r;
        swept_v = 44;

    % input variables:
    
    % output geometry parameters:
        S_v     = 1.4 * V_v * b_w*S_w/l_v;                  % area of double vertical tail in total
        c_vr    = S_v/((1 + TR_v)*b_v);
        c_vt    = TR_v*c_vr;
        
        MAC_v   = 1.1 * 2/3 * c_vr * (1 + TR_v + TR_v^2)/(1 + TR_v);
        AR_v    = ((b_v^2)*2)/S_v;

        x_vac   = x_cg + l_v;
        z_vac   = b_v*(1-(MAC_v - c_vt)/(c_vr - c_vt));
        
        x_vt    = x_vac - MAC_v/4 + (b_v - z_vac)*tand(swept_v);                               % x coordinate of tip chord leading edge
        x_vr    = x_vac - MAC_v/4 - z_vac*tand(swept_v);
        
    % Position of the tail in space
        X_vr=[b_w1half, b_w1half, b_w1half, b_w1half, b_w1half,];
        Y_vr=[-x_vr, -x_vr-c_vr, -x_vt-c_vt, -x_vt, -x_vr];
        Z_vr=[0, 0, b_v, b_v, 0];

    % Plot of the vertical Tail
        plot3(X_vr, Y_vr, Z_vr, '*');
        plot3(X_vr, Y_vr, Z_vr);

        
    % Position of the symetrical Tail    
        X_vl = -X_vr;
        Y_vl = Y_vr;
        Z_vl = Z_vr;

    % Plot of the symetric vertical Tail
        plot3(X_vl, Y_vl, Z_vl, '*');
        plot3(X_vl, Y_vl, Z_vl);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Horizon Tail        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % basic given parameters:
        b_h = b_w1;
        c_h = c_vt;
        l_h = sqrt((x_vt + 0.25*c_vt - x_cg)^2 + b_v^2);

    % input variables:
    
    % output geometry parameters:
        S_h     = b_h*c_h;
        b_hhalf = b_h/2;
        V_h     = S_h * l_h / S_w / MAC_w;
        MAC_h   = c_h;
        y_hac   = 0;
        x_hac   = l_h + x_cg;
        AR_h    = b_h^2/S_h;            % Aspect ratio

    % Position of the fin (horizontal tail)
        X_hr = [-b_w1half, -b_w1half, b_w1half, b_w1half, -b_w1half];
        Y_hr = [-x_vt, -x_vt-c_h, -x_vt-c_h, -x_vt, -x_vt];
        Z_hr = [b_v, b_v, b_v, b_v, b_v];

    % Plot of the fin (horizontal tail)
        plot3(X_hr, Y_hr, Z_hr, '*');
        plot3(X_hr, Y_hr, Z_hr);
        hold off;                       % end of the plot

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Fuselage          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    l_f = 1.4*l_h;
    r_f = 0.05;
    S_f = 2*l_f*3.14*r_f*2;       % total skin area of double fuselage


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Whole plane         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S_skin = 2.2*(S_w1 + S_w2 + S_h + S_v) + S_f;         % total skin area

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         Displays           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    S_skin
    swept_w1
    S_h
    S_w1
    S_w2
    S_h
    S_v
    S_f