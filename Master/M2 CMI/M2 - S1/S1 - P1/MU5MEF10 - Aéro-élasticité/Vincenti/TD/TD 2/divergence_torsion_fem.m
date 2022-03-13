%==================================================================================%
%                                                                                  %
%        Static analysis of the torsion divergence of a straight wing              %
%                             FINITE ELEMENT METHOD                                %
%        --------------------------------------------------------------            %
%                                                                                  %
%   structure    : Straight wing, torsion beam model, "2-point" torsion element    %
%   aerodynamics : 1st order linear theory                                         %
%   boundaries   : clamped-free                                                    %
%                                                                                  %
%   status       : static analysis under exterior moment due to the lift force     %
%                                                                                  %
%   reference    : Hodges, Introduction to aeroelasticity                          %
%                                                                                  %
%   todo-list    : convergence study                                               %
%                : parametric analysis of aeroelastic instability (divergence)     %
%                                                                                  %
%   version      : 0.1a                                                            %
%   released     : 29 oct 2021 (angela.vincenti@sorbonne-universite.fr)            %
%                                                                                  %
%==================================================================================%
   clear  all
   close  all
   
   format long
 
%  numerical parameters (to be DEFINED by the USER)
%  --------------------
%  FEM and STRUCTURE
   Ne          = 30 ;           % number of elements 
   
% =========================================================================
%  GRAPHIC OPTIONS
   NpptPlot    = 10 ;           % number of points to plot per element
   nshow       = 1.E3;          % 
   
%  --------------------------------------------------------------------------
%                                                                           !
%  AEROELASTIC MODEL : beam model of a wing, clamped-free
%                                                                           !
%  --------------------------------------------------------------------------
%
%     geometry
      l  = 0.7      ;        % span                        [m]
      c  = 0.0762   ;        % chord                       [m]
      ec = 0.25     ;        % eccentricity of aeroelastic center
      e  = ec * c   ;        % distance from aerodynamic to elastic center
      
      t  = 10^-3    ;        % thickness of the plate wing [m]
      Iz = c*t^3/12 ;        % quadratic moment of inertia along z
      Iy = c^3*t/12 ;        % quadratic moment of inertia along y
      
      J  = Iz + Iy  ;        % polar moment of inertia
      
      alpha_r_deg = 5 ;      % rigid angle of attack (degrees)
      alpha_r = alpha_r_deg*pi/180 ;  % rigid angle of attack (radians)

%     fluid properties
      rho_inf = 0.4 ;           % freestream density       [kg/m3]
      u_inf   = 150 ;           % air speed                [m/s]
      q = 0.5*rho_inf*u_inf^2 ; % dynamic pressure         [N/m2]
      a = 2 ;                   % angular lift coefficient [-] 
 
%     material properties of the plate
      nu      = 0.3;          % poisson coefficient        [-]
      E       = 7.728*1.E10;  % Young modulus              [N/m2]
      G = 0.5*E/(1 + nu) ;    % shear modulus              [N/m2]

% Loop over lineic forces
N = 10^5;
q_loop = zeros(1, N);

for i = 1:N
    q_loop(i) = (5*i / N) * 0.5*rho_inf*u_inf^2;
end

% FLUTTER
% Stock K matrix to then plot det(K) = f(q)
STOCK_det_K_mod = zeros(1, N);
STOCK_K_mod     = zeros(N, Ne+1 -1, Ne+1 -1); %-1 as we impose a BC

% ---------------------------------------------------
% loop over different loading to get divergence speed
% ---------------------------------------------------

% What is formula V_divergence and does it depend on loading ? else ?

for idx_loop = 1:N

      q = q_loop(idx_loop);
      
%     theoretical value of lambda (pulsation of torsional solution)
      lambda_theor = sqrt(e*q*c*a/G/J) ;

%     theoretical value of solution wave-length
      wl_theor = 2*pi/lambda_theor ;

%  Start the program
%  =================
% =========================================================================   
% MESH : definition
% =========================================================================   

% Length of one element 
   Le = l / Ne ;

% Total number of nodes
   Nn = Ne + 1;  

% Total number of d.d.l. (2-point Kirchhoff bending element : 2 ddl/node) 
   Nddl = Nn ;
   
% Table of nodal coordinates (Nn lines x 2 columns)
% number of the node / x of the node
   Nds = zeros(Nn, 2) ; 
   for i=1:Nn
      Nds(i,1) = i ;
      Nds(i,2) = Le * (i - 1) ;
   end    

% Table of connectivity (Ne lines x 3 columns)
% number of element / node 1 / node 2
   Connec = zeros(Ne,3) ;
   for e=1:Ne
      Connec(e,1) = e ;
      Connec(e,2) = e ;
      Connec(e,3) = e + 1 ;
   end

% =========================================================================   
% ELEMENTARY CONTRIBUTIONS
% =========================================================================   
% See theoretical developments : 2-point beam element, linear, torsion 
%  STIFFNESS MATRIX (elementary: square, symmetric, 2x2)
%  ===========
   coeff_Ke = G*J/Le ;
   Ke = coeff_Ke * [ 1 -1 ; 
                    -1  1 ] ; 
   
%  AEROELASTIC STIFFNESS MATRIX (elementary: square, anti-symmetric, 2x2)
%  ===========
   coeff_Ase = e*q*c*a*Le/6 ;
   Ase = coeff_Ase * [  2  1 ; 
                        1  2   ] ; 
                 
% NODAL VECTOR of MOMENTS from RIGID ANGLE CONTRIBUTION (elementary: 2x1)
   coeff_mom = e*q*c*a*alpha_r*Le/2 ;
   Fe = coeff_mom * [ 1 ; 1 ] ;
                    
% =========================================================================
% ASSEMBLY : from the elementary contributions to the structure
% =========================================================================   
% STIFFNESS MATRIX of the structure [K] (square, symmetric, Nddl x Nddl)
% AEROELASTIC STIFFNESS INFLUENCE MATRIX (square, symmetric, Nddl x Nddl)
% VECTOR of NODAL "FORCES" (moments) (size Nddlx1)
   K  = zeros(Nddl, Nddl) ;
   As = zeros(Nddl, Nddl) ;
   F  = zeros(Nddl, 1)    ;
   for e=1:Ne
      edc = [ e e+1 ] ;
      K(edc,edc)  = K(edc,edc)  + Ke ; 
      As(edc,edc) = As(edc,edc) + Ase ;
      F(edc)      = F(edc)      + Fe ;
   end

% =========================================================================
% =========================================================================
% STATIC ANALYSIS : solving the problem ([K]-[As]]{U}={F}
% =========================================================================
% BOUNDARY CONDITIONS : clamped-free (zero rotation angle at y=0)
% =========================================================================   
%     table of fixed displacements : number of node / value of fixed ddl
%     zero rotation angle at first node
      DepImp = [ 1  0 ] ;

%     indices of fixed ddl
      fixed_ddl = [1] ; 

%     indices of free ddl
      free_ddl = [2 : Nddl] ;
      Nddl_red = Nddl - 1 ;
  
% REDUCED SYSTEM accounting for the given boundary conditions (this method
% is valid because imposed displacements are set to zero here!!!) :
% matrices [K_red] and [M_red]
   K_red  = K(free_ddl,free_ddl) ;
   As_red = As(free_ddl,free_ddl) ;
   K_mod  = K_red - As_red ;
   F_red  = F(free_ddl) ;

   STOCK_det_K_mod(idx_loop) = det(K_mod);
   STOCK_K_mod(idx_loop,:,:) = K_mod;
   
end % loop over several lineic forces q

%% Post process

% FLUTTER
figure ; % fig det(K_mod) = f(q) to get flutter force for different modes
hold on;
plot(q_loop, STOCK_det_K_mod, '.');
hold off ;

norm_STOCK_det_K_mod = STOCK_det_K_mod/max(STOCK_det_K_mod);

idx_flutter = find(norm_STOCK_det_K_mod < 5/N & norm_STOCK_det_K_mod > -5/N);
q_flutter = q_loop(idx_flutter);

disp('Values of q where det(K)~=0 :')
disp(q_flutter)

% -------------------------------------------

% Structure (no aero influence)
disp('Values of det(K_struct) :')
disp(det(K_red))

%% Display - Displacement

% =========================================================================
% STATIC RESOLUTION  : solving the system [K_mod]{U}={F} (after reduction)
% =========================================================================   
   U_red = K_mod\F_red ;

% Reconstructing the whole vectors of modes X1 and X2  
   U = zeros(Nddl,1) ;
   U(free_ddl)  = U_red ;
   U(fixed_ddl) = DepImp(2); % BC zero angle at first node (lhs)
   
% =========================================================================
% =========================================================================
% POST-PROCESSING : plotting the modes n° 1 and 2
% =========================================================================   

% Values of the base functions Ni(x) along the element
% Plot discretisation
  dx = Le / (NpptPlot - 1) ;

% Vector of x coordinates of points to plot along one single element
  NePlot = NpptPlot ;
  xePlot = zeros(NpptPlot,1) ;
  for i=1:NpptPlot
      xePlot(i) = (i-1) * dx ;
  end
  
% Vector of values of base functions Ni(x)
  N1 = zeros(NpptPlot,1) ;
  N2 = zeros(NpptPlot,1) ;
  for i=1:NpptPlot
      csi = (xePlot(i) / Le) ;
      N1(i) = 1 - csi ;
      N2(i) = csi ;
  end
  
% Plot of base functions Ni(x) over one single element
  figure ;
  hold on;
  plot(xePlot, N1, '-ro');
  plot(xePlot, N2, '-bs');
  hold off ;
  
% Vector of x coordinates of points to plot along the structure
  NPlot = NpptPlot + (Ne - 1) * (NpptPlot - 1) ;
  xPlot = zeros(NPlot,1) ;
  for i=1:NPlot
      xPlot(i) = (i-1) * dx ;
  end

% Vector of torsional angles of rotation along the structure Vx
  theta_sol = zeros(NPlot, 1) ;

% Element 1
  for i=1:NePlot
      theta_sol(i) = N1(i)*U(1) + N2(i)*U(2) ;
  end 

%   figure ;
%   hold on;
%   plot(xePlot, Vx(1:NePlot), '-ro');
%   hold off ;
  
% Boucle sur les autres elements  
  for e=2:Ne
     i_start = NePlot + 1 + (e-2)*(NePlot - 1) ;
     i_fin = NePlot + (e-1)*(NePlot - 1) ;
     ind_1 = e ;
     ind_2 = e+1  ;
     
     for i=1:(NePlot-1)
         ii = i_start + i - 1 ; 
         theta_sol(ii) = N1(i+1)*U(ind_1) + N2(i+1)*U(ind_2) ;
     end
  end
  
  figure ;
  hold on;
  title(['Torsional angle along the beam (in radians)']) ;
  plot(xPlot, theta_sol, '-r');
  hold off ;

%% Modal displacement - FLUTTER

% =========================================================================
% STATIC RESOLUTION  : solving the system [K_mod]{U}={F} (after reduction)
% =========================================================================   

figure ;
hold on;
  
for idx_modes = idx_flutter
    A = reshape(STOCK_K_mod(idx_modes,:,:), [], size(STOCK_K_mod(idx_modes,:,:), 3), 1); %3D to 2D matrix due to bug
    U_red  = A \ F_red ;
    
% Reconstructing the whole vectors of modes X1 and X2  
   U = zeros(Nddl,1) ;
   U(free_ddl)  = U_red ;
   U(fixed_ddl) = DepImp(2); % BC zero angle at first node (lhs)

% =========================================================================
% =========================================================================
% POST-PROCESSING : plotting the modes n° 1 and 2
% =========================================================================   

% Values of the base functions Ni(x) along the element
% Plot discretisation
  dx = Le / (NpptPlot - 1) ;
  
% Vector of values of base functions Ni(x)
  N1 = zeros(NpptPlot,1) ;
  N2 = zeros(NpptPlot,1) ;
  for i=1:NpptPlot
      csi = (xePlot(i) / Le) ;
      N1(i) = 1 - csi ;
      N2(i) = csi ;
  end
  
% Vector of x coordinates of points to plot along the structure
  NPlot = NpptPlot + (Ne - 1) * (NpptPlot - 1) ;
  xPlot = zeros(NPlot,1) ;
  for i=1:NPlot
      xPlot(i) = (i-1) * dx ;
  end

% Vector of torsional angles of rotation along the structure Vx
  theta_sol = zeros(NPlot, 1) ;

% Element 1 %graph then reproduces a part of the displacement only
  for i=1:NePlot
      theta_sol(i) = N1(i)*U(1) + N2(i)*U(2);
  end 

% Boucle sur les autres elements  
  for e=2:Ne
     i_start = NePlot + 1 + (e-2)*(NePlot - 1) ;
     i_fin = NePlot + (e-1)*(NePlot - 1) ;
     ind_1 = e ;
     ind_2 = e+1  ;
     
     for i=1:(NePlot-1)
         ii = i_start + i - 1 ; 
         theta_sol(ii) = N1(i+1)*U(ind_1) + N2(i+1)*U(ind_2) ;
     end
  end
  
  title(['Torsional angle along the beam (in radians)']) ;
  plot(xPlot, theta_sol, '-');

end

legend('mode 1', 'mode 2a', 'mode 2b', 'Location', 'east');
hold off ;

%%
% ----------------
% convergence rate
% ----------------

% What to do ?