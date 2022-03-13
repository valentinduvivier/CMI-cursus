% ================================================================================== %
%                                                                                    %
%        Time-marching solution of the SUPERsonic Panel Flutter problem              %
%                             FINITE ELEMENT METHOD                                  %
%                                                                                    %
%                          ---------------------------                               %
%                                                                                    %
%   structure     : Flat plate 1D, circular bending, "2-point" bending element       %
%   aerodynamics  : 1st order piston theory (valid for M > 1.6)                      %
%   time-marching : newmark's constant average acceleration method                   %
%   boundaries    : simply supported at both ends[2]                                 %
%                 : or clamped flat plate[1,3]                                       %
%                                                                                    %
%   status        : free vibration problem: [ok]                                     %
%                 : with aerodynamics loads: unable to find the correct              %
%                  critical Mach number for case 2 (debug in progress ...)           %
%                                                                                    %
%   reference     : [1] Blom PhD Thesis, EPFL, 1998                                  %
%                 : [2] A-S. Monrouval, PhD INSA Rouen tel-00011632, 16 Feb. 2006    %
%                 : [3] S. Piperno, IJNMF, 25:1207-1226(1997)                        %
%                                                                                    %
%   todo-list     : compute the FFT at x=0.35                                        %
%                 : compute the eigenvalues and the mode shape of the free vibration %
%                 : compute the aeroelastic eigenvalues and mode shapes              %
%                                                                                    %
%   version       : 0.1a                                                             %
%   released      : 08 oct 2019 (angela.vincenti@upmc.fr)                            %
%                                                                                    %
% ================================================================================== %
   
    clear  all
    close  all
    
    format long
 
% ================================================================================== %
%                 Numerical parameters (to be DEFINED by the USER)
% ================================================================================== %

%  FEM and STRUCTURE
    Ne          = 10;           % number of elements 
    panel_model = 2;            % select the aeroelastic model 1 or 2

%  TIME INTEGRATION
    ntime       = 20000;        % number of time steps
    dt          = 4.E-6;        % size of the time step
    bi = 0.5;                   % parameter "beta"  for Newmark algorithm
    gi = 0.5;                   % parameter "gamma" for Newmark algorithm
    
%  -----------------------------
%   GRAPHIC OPTIONS
    NpptPlot    = 10;           % number of points to plot per element
    nshow       = 1.E3;         % 
 
%  --------------------------------------------------------------------------
%  AEROELASTIC MODEL #1: simply supported panel[2] (WARNING: not validated) !
%  --------------------------------------------------------------------------

    if panel_model == 1
%
%       fluid properties
        gam      = 1.4;
        rho_inf  = 0.4;         % freestream density           [kg/m3]
        prs_inf  = 13000;       % freestream pressure          [Pa]
        a_inf    = sqrt(gam*prs_inf/rho_inf);   %              [m/s]
        M_inf    = 2.1;                         %              [-]
        u_inf    = M_inf*a_inf;                 %              [m/s]

%       material properties of the plate
        rho_s    = 2710;        % density                      [kg/m3]
        nu       = 0.3;         % poisson coefficient          [-]
        length   = 0.5;         % length                       [m]
        thick    = 1.35*1.E-3;  % thickness                    [m]
        E        = 7.728*1.E10; % Young modulus                [N/m2]
        amp_init = 0.0001;      % Initial maximum amplitude    [m]
%  
        D        = E*thick^3/12/(1-nu^2);
        ms       = rho_s*thick;
        alpha    = rho_inf*u_inf^2/sqrt(M_inf^2-1.);
        beta     = rho_inf*u_inf*(M_inf^2-2.)/(M_inf^2-1.)^(3/2);
    end 
 
%  ------------------------------------------
%   AEROELASTIC MODEL #2: clamped panel[1,3]   
%                                              
%   Comment: the corresponding flutter         
%   conditions are obtained for M_inf ~ 2.25   
%  ------------------------------------------

    if panel_model == 2
%
%      fluid properties
        gam      = 1.4;
        rho_inf  = 0.4;         % freestream density          [kg/m3]
        prs_inf  = 25714;       % freestream pressure         [Pa]
        a_inf    = sqrt(gam*prs_inf/rho_inf);   %             [m/s]
        M_inf    = 1.4;                         %             [-]
        u_inf    = M_inf*a_inf;                 %             [m/s]
 
%       material properties of the plate
        rho_s    = 2710;        % density                     [kg/m3]
        nu       = 0.3;         % poisson coefficient         [-]
        length   = 0.5;         % length                      [m]
        thick    = 1.35*1.E-3;  % thickness                    [m]
        E        = 7.72810*1.E10; % Young modulus             [N/m2]
        d        = 1.;            % width of the plate          [m] 
        amp_init = 0.0001;        % Initial maximum amplitude [m]
%  
        D        = E*thick^3/12/(1-nu^2);
        ms       = rho_s*thick*d;
        alpha    = rho_inf*u_inf^2/sqrt(M_inf^2-1.);
        beta     = rho_inf*u_inf*(M_inf^2-2.)/(M_inf^2-1.)^(3/2); 
    end

%   debug only (no aerodynamics)
    alpha = 0;


% =========================================================================   
%                               Start the program
% =========================================================================   

%  ------------------
%  MESH : definition
%  ------------------

% Length of one element 
    Le = length / Ne;
    
% Total number of nodes
    Nn = Ne+1;
        
% Total number of d.d.l. (2-point Kirchhoff bending element : 2 ddl/node) 
    Nddl = 2*Nn;
   
% Table of nodal coordinates (Nn lines x 2 columns)
% number of the node / x of the node
    Nds = zeros(Nn,2); 
    for i=1:Nn
        Nds(i,1) = i;
        Nds(i,2) = Le * (i - 1);
    end 

% Table of connectivity (Ne lines x 3 columns)
% number of element / node 1 / node 2
    Connec = zeros(Ne,3);
    for e=1:Ne
        Connec(e,1) = e;
        Connec(e,2) = e;
        Connec(e,3) = e + 1;
    end

%  -----------------------------------------------
%  ELEMENTARY CONTRIBUTIONS TO MASS AND STIFFNESS 
%  -----------------------------------------------

% See theoretical developments : 2-point Kirchhoff circular bending element

%  MASS MATRIX (elementary: square, symmetric, 4x4) % see Mathematica calculation
    coeff_Me = ms * Le / 420;
    Me = coeff_Me * [[13/35     , 11*Le/210, 9/70      , -13*Le/420];
                     [11*Le/210 , Le^2/105 , 13*Le/420 , -Le^2/140 ];
                     [9/70      , 13*Le/420, 13/35     , -11*Le/210];
                     [-13*Le/420, -Le^2/140, -11*Le/210, Le^2/105  ]]*420; 

%  STIFFNESS MATRIX (elementary: square, symmetric, 4x4) % see Mathematica calculation
    coeff_Ke = D / (Le^3);
    Ke = coeff_Ke * [[ 12   ,  6*Le  , -12  ,  6*Le  ];
                     [  6*Le,  4*Le^2, -6*Le,  2*Le^2];
                     [-12   , -6*Le  , 12   , -6*Le  ];
                     [  6*Le,  2*Le^2, -6*Le,  4*Le^2]]; 

%  AEROELASTIC INFLUENCE MATRIX (elementary: square, anti-symmetric, 4x4)
    coeff_Afe = alpha;
    Afe = coeff_Afe * [[ -6/5 , -11*Le  /10,  6/5 ,   -Le  /10];
                       [-Le/10,  -2*Le^2/15, -6*Le,    Le^2/30];
                       [  6/5 ,     Le  /10, -6/5 , 11*Le  /10];
                       [-Le/10,     Le^2/30, Le/10, -2*Le^2/15]];
                 
% =========================================================================
%       ASSEMBLY : from the elementary contributions to the structure
% =========================================================================

% STIFFNESS MATRIX of the structure [K] (square, symmetric, Nddl x Nddl)
% MASS MATRIX of the structure [M] (square, symmetric, Nddl x Nddl)
% initialisation of "empty" matrices 
    K  = zeros(Nddl,Nddl);
    M  = zeros(Nddl,Nddl);
    Af = zeros(Nddl,Nddl);
 
% sum of elementary contributions
    for e=1:Ne
    
    % define the right values of d.o.f. indices in the global numbering 
    % corresponding to element "e" : edc is the vector of indices and indices
    % are functions of the element index e
        edc = [[2*e,2*e+1]; [2*e,2*e+1]];
    % sum up contribution of the e-th element in the global matrix         
        K(edc,edc)  = K(edc,edc)  + Ke; 
        M(edc,edc)  = M(edc,edc)  + Me; 
        Af(edc,edc) = Af(edc,edc) + Afe;
    end

%  -------------------------------------------------------------------------------------
%   BCs   - FREE VIBRATIONS     : [M]{Ü} + [K]{U} = {0} with fixed boundary conditions
%         - BOUNDARY CONDITIONS : (1) simply supported; (2) clamped-clamped
%  -------------------------------------------------------------------------------------

% Case "simply supported"
    if panel_model == 1
        %     vectors of d.o.f. indices in the global numbering corresponding to
        % fixed/imposed d.o.f (displacements/rotations)
        fixed_ddl = [1, Nddl-1]; 

        % free d.o.f (displacements/rotations)
        free_ddl = [2:Nddl-1];
        
        % number of free d.o.f.
        Nddl_red = Nddl - size(fixed_ddl,2);
    end
    
% -------------------------- %
        
% Case "clamped-clamped"
    if panel_model == 2
        %     vectors of d.o.f. indices in the global numbering corresponding to
        % fixed/imposed d.o.f (displacements/rotations)
        fixed_ddl = [1, 2, Nddl-1, Nddl]; 

        % free d.o.f (displacements/rotations)
        free_ddl = [3:Nddl-2];
        
        % number of free d.o.f.
        Nddl_red = Nddl - size(fixed_ddl,2);
    end

% REDUCED SYSTEM accounting for the given boundary conditions (this method
% is valid because imposed displacements are set to zero here!!!) :
% matrices [K_red] and [M_red]
    K_red  = K (free_ddl, free_ddl);
    M_red  = M (free_ddl, free_ddl);
    Af_red = Af(free_ddl, free_ddl);
  
%  -------------------------------------------------
%  MODAL ANALYSIS  : solving the eigenvalue problem 
%  -------------------------------------------------

%  Matrix of dynamic stiffness 
    DynStiff = K_red + Af_red;
    
%  Search for eigenvalues and eigenvectors 
%  lambda = matrix of eigenvalues (omega ^ 2)
%  X = matrix of eigenvectors or modal matrix
    [X_red, lambda] = eig(DynStiff);

%  Extracting vector Omega of free pulsations (non sorted vector)  
    MatrOmega = lambda ^ 0.5;
    Omega = zeros(Nddl_red,1);
    for i=1:Nddl_red
        Omega(i) = MatrOmega(i,i);
    end 
  
%  Vector of natural frequencies
    Freq = Omega / (2*pi);

% Mode 1 : frequency 'freq1', position 'pos1' and 
% associated eigenvector 'X1'
    [freq1,ind1] = min(Freq);
    X1_red = X_red(:,ind1);
    
% Mode 2 : frequency 'freq2' and position 'pos2'
    if ind1 == 1
        Freq_min1 = Freq(2:Nddl_red);
        X_red_min1 = X_red(:,2:Nddl_red);
    elseif ind1 == Nddl_red
        Freq_min1 = Freq(1:(ind1-1));
        X_red_min1 = X_red(:,1:(ind1-1));
    else
        act_ind = [1:(ind1-1) (ind1+1):Nddl_red];
        Freq_min1  = Freq(act_ind);
        X_red_min1 = X_red(:,act_ind);
    end     
    [freq2,ind2] = min(Freq_min1);
    X2_red = X_red_min1(:,ind2); 

% Reconstructing the whole vectors of modes X1 and X2  
    X1 = zeros(Nddl,1);
    X2 = zeros(Nddl,1);
    X1(free_ddl) = X1_red;
    X2(free_ddl) = X2_red;
  
%  ------------------------------------------------
%  POST-PROCESSING : plotting the modes n° 1 and 2
%  ------------------------------------------------

% Values of the base functions Ni(x) along the element

% Plot discretisation
    dx = Le / (NpptPlot - 1);
% Vector of x coordinates of points to plot along one single element
    NePlot = NpptPlot;
    xePlot = zeros(NpptPlot,1);
    for i=1:NpptPlot
        xePlot(i) = (i-1) * dx;
    end
    
% Vector of values of base functions Ni(x)
    N1 = zeros(NpptPlot,1);
    N2 = zeros(NpptPlot,1);
    N3 = zeros(NpptPlot,1);
    N4 = zeros(NpptPlot,1);
    for i=1:NpptPlot
        xsi = (xePlot(i) / Le);
        N1(i) = 1 - 3*xsi^2 + 2*xsi^3;
        N2(i) = Le*xsi*(1 - xsi)^2;
        N3(i) = 3*xsi^2 - 2*xsi^3;
        N4(i) = -1.*Le*(xsi^2 - xsi^3);
    end
    
% Plot of base functions Ni(x) over one single element
    figure;
    hold on;
    plot(xePlot, N1, '-ro');
    plot(xePlot, N3, '-bs');
    hold off;
  
    figure;
    hold on;
    plot(xePlot, N2, '-ro');
    plot(xePlot, N4, '-bs');
    hold off;

% Vector of x coordinates of points to plot along the structure
    NPlot = NpptPlot + (Ne - 1) * (NpptPlot - 1);
    xPlot = zeros(NPlot,1);
    for i=1:NPlot
        xPlot(i) = (i-1) * dx;
    end

% Vector of vertical displacements along the structure Vx
    Vx1 = zeros(NPlot,1);
    Vx2 = zeros(NPlot,1);

% Element 1
    for i=1:NePlot
        Vx1(i) = N1(i)*X1(1) + N2(i)*X1(2) + N3(i)*X1(3) + N4(i)*X1(4);
        Vx2(i) = N1(i)*X2(1) + N2(i)*X2(2) + N3(i)*X2(3) + N4(i)*X2(4);
    end
    
  figure;
  hold on;
  plot(xePlot, Vx1(1:NePlot), '-ro');
  hold off;
  
% Boucle sur les autres elements  
    for e=2:Ne
        i_start = NePlot + 1 + (e-2)*(NePlot - 1);
        i_fin = NePlot + (e-1)*(NePlot - 1);
        ind_1 = 2*e - 1;
        ind_2 = 2*e;
        ind_3 = 2*e + 1;
        ind_4 = 2*e + 2;
     
        for i=1:(NePlot-1)
            ii = i_start + i - 1; 
            Vx1(ii) = N1(i+1)*X1(ind_1) + N2(i+1)*X1(ind_2) + N3(i+1)*X1(ind_3) + N4(i+1)*X1(ind_4);
            Vx2(ii) = N1(i+1)*X2(ind_1) + N2(i+1)*X2(ind_2) + N3(i+1)*X2(ind_3) + N4(i+1)*X2(ind_4); 
        end
    end
    
    figure;
    hold on;
    title(['Mode n° 1, frequency ',num2str(freq1),' Hz']);
    plot(xPlot,Vx1,'-r');
    hold off;
 
    figure;
    hold on;
    title(['Mode n° 2, frequency ',num2str(freq2),' Hz']);
    plot(xPlot,Vx2,'-b');
    hold off;
  
% %  ----------------
% %  ANIMATION MODE 1
% %  ----------------
% 
%     figure;
%     set(gca,'nextplot','replacechildren');
%     for tt=1:100
%         coeff = (tt-1)/99;
%         plot(xPlot,coeff*Vx1,'-r')
%         axis equal
%         Film(tt) = getframe;
%     end
%     for tt=1:100
%         coeff = 1 - (tt-1)/99;
%         plot(xPlot,coeff*Vx1,'-r')
%         axis equal
%         Film(tt+100) = getframe;
%     end
%     for tt=1:100
%         coeff = -(tt-1)/99;
%         plot(xPlot,coeff*Vx1,'-r')
%         axis equal
%         Film(tt+200) = getframe;
%     end
%     for tt=1:100
%         coeff = -1 + (tt-1)/99;
%         plot(xPlot,coeff*Vx1,'-r')
%         axis equal
%         Film(tt+300) = getframe;
%     end

% %  ----------------
% %  ANIMATION MODE 2
% %  ----------------
% 
%     set(gca,'nextplot','replacechildren');
%     for tt=1:100
%         coeff = (tt-1)/99;
%         plot(xPlot,coeff*Vx2,'-b')
%         axis equal
%         Film(tt) = getframe;
%     end
%     for tt=1:100
%         coeff = 1 - (tt-1)/99;
%         plot(xPlot,coeff*Vx2,'-b')
%         axis equal
%         Film(tt) = getframe;
%     end
%     for tt=1:100
%         coeff = -(tt-1)/99;
%         plot(xPlot,coeff*Vx2,'-b')
%         axis equal
%         Film(tt) = getframe;
%     end
%     for tt=1:100
%         coeff = -1 + (tt-1)/99;
%         plot(xPlot,coeff*Vx2,'-b')
%         axis equal
%         Film(tt) = getframe;
%     end 

%%  ------------------------------------------------------------------------
%   Newmark algorithm for time integration of dynamic response: non damped system 
%   
%   Equations of reduced system (only free ddl): M_red Ü_red + K_red U_red = 0 
%   with initial displacements U0_red and initial speed V0_red

% Initial conditions (displacements and speeds)
U0_red = amp_init * X1_red;
V0_red = zeros(Nddl_red,1);

% Initialisation of accelerations
A0_red = inv(M_red) * K_red * U0_red;

% Initialisation of tables containing the history of mouvement: i-th column 
% is the result for time step (i-1)
% Reduced vectors: U_red = displacements; V_red = speeds; A_red = accelerations
U_red = zeros(Nddl_red, ntime + 1);
V_red = zeros(Nddl_red, ntime + 1);
A_red = zeros(Nddl_red, ntime + 1);

% Initialisation of column 1 = time step "0"
U_red(:,1) = U0_red;
V_red(:,1) = V0_red;
A_red(:,1) = A0_red;

% Initialisation of the auxiliary matrix S
S = M_red + bi * dt^2 * K_red;

% Vector of times
Tt = linspace(0, ntime*dt, ntime+1)';

% Iterations over time-steps (parameter "ntime"): span of time step "dt"
% At iteration "it", results are saved in the (i+1)th column of the history
% tables
for it=1:ntime
    % prediction of displacements at time-step "it"
    U_pred = zeros(Nddl_red,1);
    V_pred = zeros(Nddl_red,1);
    U_pred = U_red(:,it) + dt*V_red(:,it) + 0.5*dt^2*(1-2*bi)*A_red(:,it);
    V_pred = V_red(:,it) + dt*(1-gi)*A_red(:,it);
    
    % calculation of accelerations at time-step "it"
    A_red(:,it+1) = -1. * inv(S) * K_red * U_pred;
    
    % calculation of displacements and speeds at time-step "it"
    U_red(:,it+1) = U_pred + dt^2 * bi * A_red(:,it+1);
    V_red(:,it+1) = V_pred + dt * gi * A_red(:,it+1);
end

% Reconstruction of the history tables for the whole displacement, speed
% and acceleration vectors (all ddl, free and fixed)
U_tot = zeros(Nddl,ntime+1);
V_tot = zeros(Nddl,ntime+1);
A_tot = zeros(Nddl,ntime+1);

U_tot(free_ddl,:) = U_red;
V_tot(free_ddl,:) = V_red;
A_tot(free_ddl,:) = A_red;

% Extraction of displacement of the central point PC of the panel and plot of
% the mouvement U_PC(t)
% Number of node corresponding to PC, central point
n_PC = (Nn+1)/2;
% Index of ddl corresponding to its displacement
i_U_PC = 2 * n_PC - 1;
U_PC = U_red(i_U_PC,:)';

% Plot 
figure;
plot(Tt, U_PC)

figure;
plot(U_PC) 





















