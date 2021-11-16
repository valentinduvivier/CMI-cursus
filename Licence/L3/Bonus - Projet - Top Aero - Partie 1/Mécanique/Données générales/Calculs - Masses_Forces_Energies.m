%% Données générales - Dimension des différents composants

        g = 9.81;               % [m/s^2]
        Rho_air = 1.2;          % [kg*m^-3]

    % Electronique : moteurs, batteris, etc

        % Nb = Nombre
            Nb_moteurs = ??;            % []
            Nb_batterie = ??;           % []

            Nb_rasberry = ??            % []
            Nb_arduino = ??;            % []

            Nb_solaires = ??;           % []
            Nb_eolienne = ??;           % []

            Nb_capteurs = ??;           % []
            Nb_electro_aimant = ??;     % []
        
        % Masses
            Masse_moteur = ??;          % [kg]
            Masse_batterie = ??;        % [kg]

            Masse_raspberry = ??;       % [kg]
            Masse_arduino = ??;         % [kg]

            Masse_solaires = ??;        % [kg]
            Masse_eolienne = ??;        % [kg]

            Masse_fils = ??;            % [kg]
            Masse_capteurs = ??;        % [kg]
            Masse_electro_aimant = ??;  % [kg]  
        
        Masse_electronique = (Nb_moteurs*Masse_moteur + Nb_batterie*Masse_batterie + ...
            Nb_raspberry*Masse_raspberry + Nb_arduino*Masse_arduino + Nb_solaires*Masse_solaires + ...
            Nb_eolienne*Masse_eolienne + Masse_fils +  Nb_capteurs*Masse_capteurs + ...
            Nb_electro_aimant*Masse_electro_aimant);
    
        
    % Mecanique : ressorts, estimation poids coque, etc

        Nb_ressort_interne = ??;    % []
        Nb_ressort_externe = ??     % []
        
        Nb_rail = ??;               % []
        Nb_disque = ??;             % []
        
        
        Masse_rail = ??;            % [kg]
        Masse_disque = ??;          % [kg]

        Masse_detail = ??;          % [kg]  % Attaches, poids pour équilibre, etc
        Masse_coque_vide = ??;      % [kg]

        
        Masse_mecanique = Masse_coque_vide + Nb_ressort*Masse_ressort + Nb_rail*Masse_rail + ...
            Nb_disque*Masse_disque + Masse_detail;  % [kg] 
        
        
    % Physique : poids et dimensions parachutes, etc

        Rayon_parachute = ??;               % [m]
        Masse_parachute = ??;               % [kg]
        Volume_parachute = (2*pi*(Rayon_parachute^2))/3;   % [m^3]
        
        
        Masse_physique = Masse_parachute;   % [kg]
        
        
        
        Masse_totale = Masse_electronique + Masse_mecanique + Masse_physique;   % [kg]     

    
%% Phase 1 : Overture Parachute

    % Composantes mouvement
        Hauteur_ini = ??;         % [m]
        Vitesse_ini = 0;          % [m/s]

        Hauteur_avant_choc = ??;  % [m]
        Vitesse_avant_choc = ??;  % [m/s]
        Temps_choc = ??;          % [s]            

    % Forces
        Force_choc = ??;          % [N]  % Donné en fonction de la vitesse au moment de l'ouverture : Partie Physique

    % Energies
    
        % Energie du module avant le choc
            Ec1 = (Masse_totale *(Vitesse_avant_choc^2))/2;   % [J]
            Ep1 = (Masse_totale * g * Hauteur_avant_choc);    % [J]

        % Energie dissipée au sein du module (conservation de l'énergie
        % mécanique)
            Hauteur_apres_choc = ??;     % [m]
            Vitesse_apres_choc = ??;     % [m/s]

            Ec2 = (Masse_totale *(Vitesse_apres_choc^2))/2;      % [J]
            Ep2 = (Masse_totale * g * Hauteur_apres_choc);       % [J]
            
            E_dissipe = (Ec1 + Ep1) - (Ec2 + Ep2);      % [J]
    
%% Phase 2 Atterissage

    Vitesse_impacte = ??;     % [m/s]
    

    Force_impacte = ??;       % [N]
    
    
%% Centre de poussé et d'inertie


    % Centre d'inertie
        
        % Données de dimensionnement
        Hauteur_module = 0.2;           % [m]
        Rayon_parachute = 0.5;          % [m]
        Ecart_parachute_module = 0.2;   % [m]

        Volume_parachute = (2*pi*(Rayon_parachute^3))/3;     % [m^2]

        Masse_module = 1;                               % [kg]
        Masse_parachute = Rho_air*Volume_parachute;    % [kg]
        Masse_totale = Masse_module + Masse_parachute;  % [kg]

        % Origine = bas de la forme considérée
        OC_module = Hauteur_module/2;       % [m]
        OC_parachute = Rayon_parachute/2;   % [m]

        % Origine = bas du module
        Position_OC_total = (Masse_module*OC_module + Masse_parachute*(Hauteur_module + Ecart_parachute_module + OC_parachute))/Masse_totale;

    
    
    