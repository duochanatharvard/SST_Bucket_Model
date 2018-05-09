function [SST,Budget] = BKT_MD_STP_2_MD_WOODEN_SNG(true_SST,true_AT,e_air,u_environment,s_environment,Cs,direct_ratio,deck_time,solar_shading,zenith_angle)

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
    thickness         = 0.01;  % unit: m
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
    if Re_haul < 1000,
        h_cycle_haul = 2.8 * (u0_haul./diameter).^0.5;
    else
        h_cycle_haul = 4.3 * ((u0_haul).^(0.6))./((diameter).^(0.4));
    end
    if Re_deck < 1000,
        h_cycle_deck = 2.8 * (u0_deck./diameter).^0.5;
    else
        h_cycle_deck = 4.3 * ((u0_deck).^(0.6))./((diameter).^(0.4));
    end
    h_base_haul = 4.3 * (u0_haul./diameter).^0.5;
    h_base_deck = 4.3 * (u0_deck./diameter).^0.5;
    
    
    % ---------------------------------------------------------------------
    % Compute the temperature in the bucket -------------------------------
    time_step = (t_haul + t_deck) / dt;
    SST    = nan([size(true_SST),time_step]);
    BT_wall     = nan([size(true_SST),layers,time_step]);
    BT_base     = nan([size(true_SST),layers,time_step]);
    heat_flux_wall     = nan([size(true_SST),time_step]);
    heat_flux_base     = nan([size(true_SST),time_step]);
    dQs_topdt   = nan([size(true_SST),time_step]);
    dQs_walldt  = nan([size(true_SST),time_step]);
    dQs_basedt  = nan([size(true_SST),time_step]);
    dQw_topdt   = nan([size(true_SST),time_step]);
    dQw_walldt  = nan([size(true_SST),time_step]);
    dQw_basedt  = nan([size(true_SST),time_step]);
    dQl_topdt   = nan([size(true_SST),time_step]);
    dQl_walldt  = nan([size(true_SST),time_step]);
    dQl_basedt  = nan([size(true_SST),time_step]);
    dQ_topdt    = nan([size(true_SST),time_step]);
    dQ_cycledt  = nan([size(true_SST),time_step]);
    dQ_basedt   = nan([size(true_SST),time_step]);
    dQ_base2dt  = nan([size(true_SST),time_step]);
    dQ_walldt   = nan([size(true_SST),time_step]);
    dQdt   = nan([size(true_SST),time_step]);
    dTdt   = nan([size(true_SST),time_step]);
    dQrdt_cycle  = nan([size(true_SST),time_step]);
    dQrdt_base   = nan([size(true_SST),time_step]);
    mass_water = nan([size(true_SST),time_step]);
    e_bucket   = nan([size(true_SST),time_step]);
    e_out_bucket_wall   = nan([size(true_SST),time_step]);
    e_out_bucket_base   = nan([size(true_SST),time_step]);

    % Set initial condition -----------------------------------------------
    SST(:,:,1) = true_SST;
    for i = 1:layers
        BT_wall(:,:,i,1) = true_SST;
        BT_base(:,:,i,1) = true_SST;
    end
    mass_water(:,:,1) = mass;
    mass_wall = A_cycle * density * thickness_water/1000 .* cp_water;
    mass_base = A_base  * density * thickness_water/1000 .* cp_water;

    % Compute the solar contaimination ------------------------------------
    clear('temp_cycle_direct','temp_cycle_base','temp_base')
    Cs_diffuse = Cs .* (1 - direct_ratio);
    Cs_direct  = Cs .* direct_ratio;
    temp_cycle_direct  = A_cycle ./pi .* (Cs_direct ./ cos(zenith_angle) .* sin(zenith_angle) .* (1-albedo_bucket));
    temp_cycle_diffuse = A_cycle .* (Cs_diffuse ./2 .* (1-albedo_bucket));
    temp_base  = A_base .* (Cs .* (1-cover_top));

    for t = 1:time_step-1

        e_bucket(:,:,t) = 6.112 * exp(17.67 * (SST(:,:,t) - 273.15)./(SST(:,:,t) - 29.65)) .* 0.98;
        e_out_bucket_wall(:,:,t) = 6.112 * exp(17.67 * (BT_wall(:,:,end,t) - 273.15)./(BT_wall(:,:,end,t) - 29.65)) .* 0.98;
        e_out_bucket_base(:,:,t) = 6.112 * exp(17.67 * (BT_base(:,:,end,t) - 273.15)./(BT_base(:,:,end,t) - 29.65)) .* 0.98;

        % Compute each term of the rediative budget for the bucket top ----
        if((t * dt) <= t_haul) % During the process of hauling the bucket
            dQs_topdt(:,:,t) = - (h_base_haul .* A_base) .* (SST(:,:,t) - true_AT) .* Top_sensible_on;
            dQw_topdt(:,:,t) = -1.7 * (h_base_haul .* A_base) .* (e_bucket(:,:,t) - e_air) .* Top_latent_on;
        else            % During on deck stage
            dQs_topdt(:,:,t) = - (h_base_deck .* A_base) .* (SST(:,:,t) - true_AT) .* Top_sensible_on;
            dQw_topdt(:,:,t) = -1.7 * (h_base_deck .* A_base) .* (e_bucket(:,:,t) - e_air) .* Top_latent_on;
        end

        dQl_topdt(:,:,t) = - (A_base) .* sigma .* (SST(:,:,t).^4 - true_AT.^4) .* Top_long_on;
        dQrdt_base(:,:,t)  = (temp_base) * (1-shading) * Top_solar_on; % Solar Heating

        dQ_topdt(:,:,t) = (dQs_topdt(:,:,t) + dQw_topdt(:,:,t) + dQl_topdt(:,:,t) + dQrdt_base(:,:,t)) .* Top_cooling_on;

        % Compute terms for bucket walls ----------------------------------
        if((t * dt) <= t_haul)
            dQs_walldt(:,:,t) =  - (h_cycle_haul .* A_cycle) .* (BT_wall(:,:,end,t) - true_AT) .* Wall_sensible_on;
            dQw_walldt(:,:,t) =  -1.7 * (h_cycle_haul .* A_cycle) .* (e_out_bucket_wall(:,:,t) - e_air) .* Wall_latent_on;
        else
            dQs_walldt(:,:,t) =  - (h_cycle_deck .* A_cycle) .* (BT_wall(:,:,end,t) - true_AT) .* Wall_sensible_on;
            dQw_walldt(:,:,t) =  -1.7 * (h_cycle_deck .* A_cycle) .* (e_out_bucket_wall(:,:,t) - e_air) .* Wall_latent_on;
        end
        dQl_walldt(:,:,t)    = - (A_cycle) .* sigma .* (BT_wall(:,:,end,t).^4 - true_AT.^4) .* Wall_long_on;
        dQrdt_cycle(:,:,t)   = (temp_cycle_direct + temp_cycle_diffuse) * (1-shading) .* Wall_solar_on;   % Solar Heating
        heat_flux_wall(:,:,t) = A_cycle .* (BT_wall(:,:,end-1,t) - BT_wall(:,:,end,t)) ./ dx .* thc .* Wall_heatflux_on;

        dQ_walldt(:,:,t) = dQs_walldt(:,:,t) + dQw_walldt(:,:,t) + dQl_walldt(:,:,t) + dQrdt_cycle(:,:,t) + heat_flux_wall(:,:,t);
        BT_wall(:,:,end,t+1) = BT_wall(:,:,end,t) + dQ_walldt(:,:,t) ./ mass_wall .*dt;

        % Solve interitively for bucket walls -----------------------------
        for ly = layers-1:-1:2
            BT_wall(:,:,ly,t+1) = BT_wall(:,:,ly,t) + alpha ./ (dx.^2) .* (BT_wall(:,:,ly-1,t) - 2.*BT_wall(:,:,ly,t) + BT_wall(:,:,ly+1,t)) .*dt;
        end
        BT_wall(:,:,1,t+1) = BT_wall(:,:,1,t) + alpha ./ (dx.^2) .* (SST(:,:,t) - 2.*BT_wall(:,:,1,t) + BT_wall(:,:,2,t)) .*dt;
        dQ_cycledt(:,:,t) = - thc .* A_cycle .* (SST(:,:,t) - BT_wall(:,:,1,t)) ./dx .* Wall_cooling_on;

        % Compute terms for bucket base -----------------------------------
        if((t * dt) <= t_haul)
            dQs_basedt(:,:,t) =  - (h_base_haul .* A_base) .* (BT_base(:,:,end,t) - true_AT) .* Base_sensible_on;
            dQw_basedt(:,:,t) =  -1.7 * (h_base_haul .* A_base) .* (e_out_bucket_base(:,:,t) - e_air) .* Base_latent_on;
            dQl_basedt(:,:,t)    = - (A_base) .* sigma .* (BT_base(:,:,end,t).^4 - true_AT.^4) .* Base_long_on;
        else
            dQs_basedt(:,:,t) =  - (h_base_deck .* A_base) .* (BT_base(:,:,end,t) - true_AT) .* Base_sensible_on .* 0;
            dQw_basedt(:,:,t) =  -1.7 * (h_base_deck .* A_base) .* (e_out_bucket_base(:,:,t) - e_air) .* Base_latent_on .* 0;
            dQl_basedt(:,:,t)    = - (A_base) .* sigma .* (BT_base(:,:,end,t).^4 - true_AT.^4) .* Base_long_on .* 0;
        end
        heat_flux_base(:,:,t) = A_base .* (BT_base(:,:,end-1,t) - BT_base(:,:,end,t)) ./ dx .* thc .* Base_heatflux_on;

        dQ_base2dt(:,:,t) = dQs_basedt(:,:,t) + dQw_basedt(:,:,t) + dQl_basedt(:,:,t) + heat_flux_base(:,:,t);
        BT_base(:,:,end,t+1) = BT_base(:,:,end,t) + dQ_base2dt(:,:,t) ./ mass_base .*dt;

        % Solve interitively for bucket base ------------------------------
        for ly = layers-1:-1:2
            BT_base(:,:,ly,t+1) = BT_base(:,:,ly,t) + alpha ./ (dx.^2) .* (BT_base(:,:,ly-1,t) - 2.*BT_base(:,:,ly,t) + BT_base(:,:,ly+1,t)) .*dt;
        end
        BT_base(:,:,1,t+1) = BT_base(:,:,1,t) + alpha ./ (dx.^2) .* (SST(:,:,t) - 2.*BT_base(:,:,1,t) + BT_base(:,:,2,t)) .*dt;
        dQ_basedt(:,:,t) = - thc .* A_base .* (SST(:,:,t) - BT_base(:,:,1,t)) ./dx .* Base_cooling_on;

        % Solve the heat flux from Bucket wall and Base -------------------
        if ((t * dt) <= t_haul || (t * dt) >=t_haul + 30)
            dQdt(:,:,t) = dQ_topdt(:,:,t) + dQ_cycledt(:,:,t) + dQ_basedt(:,:,t) ;
        else
            dQdt(:,:,t) = dQ_topdt(:,:,t) + dQ_cycledt(:,:,t) + dQ_basedt(:,:,t)  - 0.0012.* (SST(:,:,t) - true_AT);
        end

        % Compute the leakage of mass -------------------------------------
        if (t * dt) <= t_haul,
            mass_water(:,:,t+1) = mass_water(:,:,t) - leakage_rate_haul./60 .* A_base .* density .*dt;
        else
            mass_water(:,:,t+1) = mass_water(:,:,t) - leakage_rate_deck./60 .* A_base .* density .*dt;
        end

        % Compute temperature tendency ------------------------------------
        if (t * dt) <= t_haul,
            dTdt(:,:,t) = dQdt(:,:,t) ./ (mass_water(:,:,t).*cp_water);
        else
            dTdt(:,:,t) = dQdt(:,:,t) ./ ((mass_water(:,:,t) + water_amount/1000).*cp_water);
        end

        SST(:,:,t+1) = SST(:,:,t) + dTdt(:,:,t) .*dt;
    end

    Budget = [];
end