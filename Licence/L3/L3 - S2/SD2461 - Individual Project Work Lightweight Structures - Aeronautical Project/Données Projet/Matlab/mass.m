% mass estimation


L_tube=1.5;
S_frame=0.5;


% electronic
m_motor=0.095*2;    % 0.095*2
m_esc=0.1*2;      % 0.1*2
m_servo=0.0285*8;    % 0.0285*8
m_prop=0.1;     % 0.05*2
m_board=0.1;    % 0.05 + some cable
m_payload=0.5;   % payload including FPV system or other things
m_battery=0.17*2;
m_elec=m_motor+m_esc+m_servo+m_prop+m_board+m_payload+m_battery

% structure
dens_hull=1.4;  % kg/m^2
dens_tube=0.2;  % kg/m
dens_frame=1.5; % kg/m^2
m_hull=S_skin*dens_hull
m_frame=L_tube*dens_tube+S_frame*dens_frame
m_other=0.2     % other structures, i.e., connection parts or screws
m_struc=m_hull+m_frame+m_other

% total mass
m_tot=m_elec+m_struc




