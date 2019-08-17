function [SST_out,Budget] = BKT_MD_STP_2_MD_WOODEN_GRD_SIZ(true_SST,true_AT,e_air,...
    u_environment,s_environment,Cs,direct_ratio,deck_time,solar_shading,...
    zenith_angle,thickness)

    % Parameter of computation --------------------------------------------
    dt = 0.2;                  % unit: s

    % Some physical constant parameters -----------------------------------
    sigma         = 5.67e-8;   % S-B constant
    density       = 1023;      % density of water: kg/m^3
    cp_water      = 4200;      % specific heat of water: J/kg/K
    % s_environment = 7;       % ship speed: unit: m/s
    viscosity     = 1.5e-5;    % Terbulence viscosity: m^2/s
    thickness_water = 0.1;     % Thickness of water: unit: mm

    % Parameters of the bucket --------------------------------------------
    cp_bucket         = 1900;
    density_bucket    = 800;
    layers            = 5;
    % thickness         = 0.01;  % unit: m
    leakage_rate_haul = 0;  % unit: m/min
    leakage_rate_deck = 0;
    diameter          = 0.25;
    depth             = 0.2;
    albedo_bucket     = 0;
    cover_top         = 0;
    shading           = solar_shading;
    thc               = 0.3;   % Thermal Conductivity
    dx                = thickness / layers;
    alpha = thc/density_bucket/cp_bucket;

    % Parameter of measurement --------------------------------------------
    t_haul = 60;
    t_deck = deck_time;
    u_shield_haul = 0.6;
    u_shield_deck = 0.4;
    s_shield_haul = 1;
    s_shield_deck = 0.67;
    water_amount = 35;     % in the thermalmeter, unit: gram

    % Term Management -----------------------------------------------------
    Top_cooling_on   = 1;
    Top_sensible_on  = 1;
    Top_latent_on    = 1;
    Top_long_on      = 1;
    Top_solar_on     = 1;
    Wall_cooling_on  = 1;
    Wall_sensible_on = 1;
    Wall_latent_on   = 1;
    Wall_long_on     = 1;
    Wall_solar_on    = 1;
    Wall_heatflux_on = 1;
    Base_cooling_on  = 1;
    Base_sensible_on = 1;
    Base_latent_on   = 1;
    Base_long_on     = 1;
    Base_heatflux_on = 1;

    % reduce ambient wind speed -------------------------------------------
    u_reduced_haul = u_environment .* u_shield_haul;
    u_reduced_deck = u_environment .* u_shield_deck;

    % reduce ship speed ---------------------------------------------------
    s_reduced_haul = s_environment .* s_shield_haul;
    s_reduced_deck = s_environment .* s_shield_deck;

    % effective wind speed for different stage of measurement -------------
    u0_haul = sqrt(u_reduced_haul.^2 + s_reduced_haul.^2);
    u0_deck = sqrt(u_reduced_deck.^2 + s_reduced_deck.^2);

    % Determine the Rayleigh coefficient ----------------------------------
    Re_haul = u0_haul * diameter ./ viscosity;
    Re_deck = u0_deck * diameter ./ viscosity;

    % Bucket Size and Area coefficients -----------------------------------
    A_base = diameter.^2 ./ 4 .* pi;
    A_cycle = diameter .* pi .* depth;
    mass = A_base .* depth .* density;  % unit: kg

    % Sensiable heat flux coefficient -------------------------------------
    h_cycle_haul = 2.8 * (u0_haul./diameter).^0.5 .* (Re_haul < 1000) + ...
        4.3 * ((u0_haul).^(0.6))./((diameter).^(0.4)) .* (Re_haul >= 1000) ;
    h_cycle_deck = 2.8 * (u0_deck./diameter).^0.5 .* (Re_deck < 1000) + ...
        4.3 * ((u0_deck).^(0.6))./((diameter).^(0.4)) .* (Re_deck >= 1000);
    h_base_haul = 4.3 * (u0_haul./diameter).^0.5;
    h_base_deck = 4.3 * (u0_deck./diameter).^0.5;

    % ---------------------------------------------------------------------
    % Compute the temperature in the bucket -------------------------------
    time_step = (t_haul + t_deck) / dt;
    time_out  = (t_haul + t_deck) / 30 + 1;

    BT_wall     = nan([size(true_SST),layers]);
    BT_base     = nan([size(true_SST),layers]);

    % Set initial condition -----------------------------------------------
    SST = true_SST;
    for i = 1:layers
        BT_wall(:,:,:,:,i) = true_SST;
        BT_base(:,:,:,:,i) = true_SST;
    end
    mass_water = ones(size(true_SST)) * mass;
    SST_out    = nan([size(true_SST),time_out]);
    SST_out(:,:,:,:,1) = SST;
    ct_out = 1;

    mass_wall = A_cycle * density * thickness_water/1000 .* cp_water;
    mass_base = A_base  * density * thickness_water/1000 .* cp_water;

    % Compute the solar contaimination ------------------------------------
    clear('temp_cycle_direct','temp_cycle_base','temp_base')
    Cs_diffuse = Cs .* (1 - direct_ratio);
    Cs_direct  = Cs .* direct_ratio;
    temp_cycle_direct  = A_cycle ./pi .* (Cs_direct ./ cos(zenith_angle) .* sin(zenith_angle) .* (1-albedo_bucket));
    temp_cycle_diffuse = A_cycle .* (Cs_diffuse ./2 .* (1-albedo_bucket));
    temp_base  = A_base .* (Cs .* (1-cover_top));

    for t = 1:time_step

        e_bucket= 6.112 * exp(17.67 * (SST - 273.15)./(SST - 29.65)) .* 0.98;
        e_out_bucket_wall = 6.112 * exp(17.67 * (BT_wall(:,:,:,:,end) - 273.15)./(BT_wall(:,:,:,:,end) - 29.65)) .* 0.98;
        e_out_bucket_base = 6.112 * exp(17.67 * (BT_base(:,:,:,:,end) - 273.15)./(BT_base(:,:,:,:,end) - 29.65)) .* 0.98;

        % Compute each term of the rediative budget for the bucket top ----
        if((t * dt) <= t_haul) % During the process of hauling the bucket
            dQs_topdt = - (h_base_haul .* A_base) .* (SST - true_AT) .* Top_sensible_on;
            dQw_topdt = -1.7 * (h_base_haul .* A_base) .* (e_bucket - e_air) .* Top_latent_on;
        else            % During on deck stage
            dQs_topdt = - (h_base_deck .* A_base) .* (SST - true_AT) .* Top_sensible_on;
            dQw_topdt = -1.7 * (h_base_deck .* A_base) .* (e_bucket - e_air) .* Top_latent_on;
        end

        dQl_topdt = - (A_base) .* sigma .* (SST.^4 - true_AT.^4) .* Top_long_on;
        dQrdt_base  = (temp_base) * (1-shading) * Top_solar_on; % Solar Heating

        dQ_topdt = (dQs_topdt + dQw_topdt + dQl_topdt + dQrdt_base) .* Top_cooling_on;

        % Compute terms for bucket walls ----------------------------------
        if((t * dt) <= t_haul)
            dQs_walldt =  - (h_cycle_haul .* A_cycle) .* (BT_wall(:,:,:,:,end) - true_AT) .* Wall_sensible_on;
            dQw_walldt =  -1.7 * (h_cycle_haul .* A_cycle) .* (e_out_bucket_wall - e_air) .* Wall_latent_on;
        else
            dQs_walldt =  - (h_cycle_deck .* A_cycle) .* (BT_wall(:,:,:,:,end) - true_AT) .* Wall_sensible_on;
            dQw_walldt =  -1.7 * (h_cycle_deck .* A_cycle) .* (e_out_bucket_wall - e_air) .* Wall_latent_on;
        end
        dQl_walldt     = - (A_cycle) .* sigma .* (BT_wall(:,:,:,:,end).^4 - true_AT.^4) .* Wall_long_on;
        dQrdt_cycle    = (temp_cycle_direct + temp_cycle_diffuse) * (1-shading) .* Wall_solar_on;   % Solar Heating
        heat_flux_wall = A_cycle .* (BT_wall(:,:,:,:,end-1) - BT_wall(:,:,:,:,end)) ./ dx .* thc .* Wall_heatflux_on;

        dQ_walldt = dQs_walldt + dQw_walldt + dQl_walldt + dQrdt_cycle + heat_flux_wall;
        BT_wall(:,:,:,:,end) = BT_wall(:,:,:,:,end) + dQ_walldt ./ mass_wall .*dt;

        % Solve interitively for bucket walls -----------------------------
        for ly = layers-1:-1:2
            BT_wall(:,:,:,:,ly) = BT_wall(:,:,:,:,ly) + alpha ./ (dx.^2) .* (BT_wall(:,:,:,:,ly-1) - 2.*BT_wall(:,:,:,:,ly) + BT_wall(:,:,:,:,ly+1)) .*dt;
        end
        BT_wall(:,:,:,:,1) = BT_wall(:,:,:,:,1) + alpha ./ (dx.^2) .* (SST - 2.*BT_wall(:,:,:,:,1) + BT_wall(:,:,:,:,2)) .*dt;
        dQ_cycledt = - thc .* A_cycle .* (SST - BT_wall(:,:,:,:,1)) ./dx .* Wall_cooling_on;

        % Compute terms for bucket base -----------------------------------
        if((t * dt) <= t_haul)
            dQs_basedt =  - (h_base_haul .* A_base) .* (BT_base(:,:,:,:,end) - true_AT) .* Base_sensible_on;
            dQw_basedt =  -1.7 * (h_base_haul .* A_base) .* (e_out_bucket_base - e_air) .* Base_latent_on;
            dQl_basedt    = - (A_base) .* sigma .* (BT_base(:,:,:,:,end).^4 - true_AT.^4) .* Base_long_on;
        else
            dQs_basedt =  - (h_base_deck .* A_base) .* (BT_base(:,:,:,:,end) - true_AT) .* Base_sensible_on .* 0;
            dQw_basedt =  -1.7 * (h_base_deck .* A_base) .* (e_out_bucket_base - e_air) .* Base_latent_on .* 0;
            dQl_basedt    = - (A_base) .* sigma .* (BT_base(:,:,:,:,end).^4 - true_AT.^4) .* Base_long_on .* 0;
        end
        heat_flux_base = A_base .* (BT_base(:,:,:,:,end-1) - BT_base(:,:,:,:,end)) ./ dx .* thc .* Base_heatflux_on;

        dQ_base2dt = dQs_basedt + dQw_basedt + dQl_basedt + heat_flux_base;
        BT_base(:,:,:,:,end) = BT_base(:,:,:,:,end) + dQ_base2dt ./ mass_base .*dt;

        % Solve interitively for bucket base ------------------------------
        for ly = layers-1:-1:2
            BT_base(:,:,:,:,ly) = BT_base(:,:,:,:,ly) + alpha ./ (dx.^2) .* (BT_base(:,:,:,:,ly-1) - 2.*BT_base(:,:,:,:,ly) + BT_base(:,:,:,:,ly+1)) .*dt;
        end
        BT_base(:,:,:,:,1) = BT_base(:,:,:,:,1) + alpha ./ (dx.^2) .* (SST - 2.*BT_base(:,:,:,:,1) + BT_base(:,:,:,:,2)) .*dt;
        dQ_basedt = - thc .* A_base .* (SST - BT_base(:,:,:,:,1)) ./dx .* Base_cooling_on;

        % Solve the heat flux from Bucket wall and Base -------------------
        if ((t * dt) <= t_haul || (t * dt) >=t_haul + 30)
            dQdt = dQ_topdt + dQ_cycledt + dQ_basedt ;
        else
            dQdt = dQ_topdt + dQ_cycledt + dQ_basedt - 0.0012.* (SST - true_AT);
        end

        % Compute the leakage of mass -------------------------------------
        if (t * dt) <= t_haul,
            mass_water = mass_water - leakage_rate_haul./60 .* A_base .* density .*dt;
        else
            mass_water = mass_water - leakage_rate_deck./60 .* A_base .* density .*dt;
        end

        % Compute temperature tendency ------------------------------------
        if (t * dt) <= t_haul,
            dTdt = dQdt ./ (mass_water.*cp_water);
        else
            dTdt = dQdt ./ ((mass_water + water_amount/1000).*cp_water);
        end

        SST = SST + dTdt .*dt;

        if rem(t*dt,30) == 0,
            ct_out = ct_out + 1;
            SST_out(:,:,:,:,ct_out) = SST;
        end

        if rem(t*dt,10) == 0,
            disp([num2str(t*dt),' seconds finished'])
        end
    end
end
