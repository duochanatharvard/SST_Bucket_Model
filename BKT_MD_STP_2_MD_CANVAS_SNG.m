function [SST,Budget] = BKT_MD_STP_2_MD_CANVAS_SNG(true_SST,true_AT,e_air,u_environment,s_environment,Cs,direct_ratio,deck_time,solar_shading,zenith_angle)

    % Parameter of computation --------------------------------------------
    dt = 1;                    % unit: s

    % Some physical constant parameters -----------------------------------
    sigma         = 5.67e-8;   % S-B constant
    density       = 1023;      % density of water: kg/m^3
    cp_water      = 4200;      % specific heat of water: J/kg/K
    % s_environment = 7;       % unit: m/s
    viscosity     = 1.5e-5;    % Terbulence viscosity: m^2/s

    % Parameters of the bucket --------------------------------------------
    cp_bucket         = 3900;
    mass_bucket       = 0.2;
    leakage_rate_haul = 0.01;  % unit: m/min
    leakage_rate_deck = 0.005;
    diameter          = 0.163;
    depth             = 0.14;
    albedo_bucket     = 0.2;
    cover_top         = 0.5;
    shading           = solar_shading;

    % Parameter of measurement --------------------------------------------
    t_haul = 60;
    t_deck = deck_time;
    u_shield_haul = 0.6;
    u_shield_deck = 0.4;
    s_shield_haul = 1;
    s_shield_deck = 0.67;
    water_amount = 35;     % in the thermalmeter, unit: gram

    % Term Management -----------------------------------------------------
    sensible_on   = 1;
    latent_on     = 1;
    long_on       = 1;
    solar_top_on  = 1;
    solar_wall_on = 1;

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
    SST    = nan([size(true_SST),time_step ]);
    dQsdt  = nan([size(true_SST),time_step ]);
    dQwdt  = nan([size(true_SST),time_step ]);
    dQldt  = nan([size(true_SST),time_step ]);
    dQdt   = nan([size(true_SST),time_step ]);
    dTdt   = nan([size(true_SST),time_step ]);
    dQrdt_cycle  = nan([size(true_SST),time_step ]);
    dQrdt_base   = nan([size(true_SST),time_step ]);
    mass_water = nan([size(true_SST),time_step ]);
    e_bucket   = nan([size(true_SST),time_step ]);
    
    % Set initial condition -----------------------------------------------
    SST(:,:,1) = true_SST;
    mass_water(:,:,1) = mass;

    % Compute the solar contaimination ------------------------------------
    clear('temp_cycle_direct','temp_cycle_base','temp_base')
    Cs_diffuse = Cs .* (1 - direct_ratio);
    Cs_direct  = Cs .* direct_ratio;
    temp_cycle_direct  = A_cycle ./pi .* (Cs_direct ./ cos(zenith_angle) .* sin(zenith_angle) .* (1-albedo_bucket));
    temp_cycle_diffuse = A_cycle .* (Cs_diffuse ./2 .* (1-albedo_bucket));
    temp_base  = A_base .* (Cs .* (1-cover_top));

    for t = 1:time_step-1

        e_bucket(:,:,t) = 6.112 * exp(17.67 * (SST(:,:,t) - 273.15)./(SST(:,:,t) - 29.65)) .* 0.98;

        % Compute each term of the rediative budget -----------------------
        if((t * dt) <= t_haul) % During the process of hauling the bucket
            dQsdt(:,:,t) = - (h_cycle_haul .* A_cycle + h_base_haul .* A_base) .* (SST(:,:,t) - true_AT) .* sensible_on;
            dQwdt(:,:,t) = -1.7 * (h_cycle_haul .* A_cycle + h_base_haul .* A_base) .* (e_bucket(:,:,t) - e_air) .* latent_on;
        else            % During on deck stage
            dQsdt(:,:,t) = - (h_cycle_deck .* A_cycle + h_base_deck .* A_base) .* (SST(:,:,t) - true_AT) .* sensible_on;
            dQwdt(:,:,t) = -1.7 * (h_cycle_deck .* A_cycle + h_base_deck .* A_base) .* (e_bucket(:,:,t) - e_air) .* latent_on;
        end

        dQldt(:,:,t) = - (A_base + A_cycle) .* sigma .* (SST(:,:,t).^4 - true_AT.^4) .* long_on;
        dQrdt_cycle(:,:,t) = (temp_cycle_direct + temp_cycle_diffuse) .* (1 - shading) .* solar_wall_on;
        dQrdt_base(:,:,t)  = (temp_base) .* (1 - shading) .* solar_top_on;

        % Compute the total rediative budget ------------------------------
        if ((t * dt) <= t_haul || (t * dt) >=t_haul + 30)
            dQdt(:,:,t) = dQsdt(:,:,t) + dQwdt(:,:,t) + dQldt(:,:,t) + dQrdt_cycle(:,:,t) + dQrdt_base(:,:,t);
        else
            dQdt(:,:,t) = dQsdt(:,:,t) + dQwdt(:,:,t) + dQldt(:,:,t) + dQrdt_cycle(:,:,t) + dQrdt_base(:,:,t) - 0.0012.* (SST(:,:,t) - true_AT);
        end

        % Compute the leakage of mass -------------------------------------
        if (t * dt) <= t_haul,
            mass_water(:,:,t+1) = mass_water(:,:,t) - leakage_rate_haul./60.*A_base.*density .* dt;
        else
            mass_water(:,:,t+1) = mass_water(:,:,t) - leakage_rate_deck./60.*A_base.*density .* dt;
        end

        % Compute temperature tendency ------------------------------------
        if (t * dt) <= t_haul,
            dTdt(:,:,t) = dQdt(:,:,t) ./ (mass_bucket.*cp_bucket + mass_water(:,:,t).*cp_water);
        else
            dTdt(:,:,t) = dQdt(:,:,t) ./ (mass_bucket.*cp_bucket + (mass_water(:,:,t) + water_amount/1000).*cp_water);
        end

        SST(:,:,t+1) = SST(:,:,t) + dTdt(:,:,t) * dt;
    end
    
    Budget.dQsdt  = dQsdt; 
    Budget.dQwdt  = dQwdt;
    Budget.dQldt  = dQldt;
    Budget.dQrdt_cycle = dQrdt_cycle;
    Budget.dQrdt_base  = dQrdt_base;

end