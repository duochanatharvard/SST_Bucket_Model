% FP95 Canvas bucket model
% [SST_out] = BKT_MD_STP_2_MD_CANVAS_GRD_SIZ(true_SST,true_AT,e_air,...
%    u_environment,Cs,direct_ratio,zenith_angle,P)
% 
% P.deck_time,  P.solar_shading,   P.cover_top  P.s_environment
% P.diameter,   P.depth,           P.exp_id,    P.mass_bucket

function [SST_out] = BKT_MD_STP_2_MD_CANVAS_GRD_SIZ(true_SST,true_AT,e_air,...
    u_environment,Cs,direct_ratio,zenith_angle,P)
    % P.deck_time,    P.solar_shading,   P.cover_top
    % P.diameter,   P.depth,    P.exp_id,   P.mass_bucket
    % P.s_environment

    % Parameter of computation --------------------------------------------
    dt = 1;                    % unit: s

    % Some physical constant parameters -----------------------------------
    sigma         = 5.67e-8;   % S-B constant
    density       = 1023;      % density of water: kg/m^3
    cp_water      = 4200;      % specific heat of water: J/kg/K
    s_environment = P.s_environment;       % unit: m/s
    viscosity     = 1.5e-5;    % Terbulence viscosity: m^2/s

    % Parameters of the bucket --------------------------------------------
    cp_bucket         = 3900;
    mass_bucket       = P.mass_bucket;
    leakage_rate_haul = 0.01;  % unit: m/min
    leakage_rate_deck = 0.005;
    diameter          = P.diameter;
    depth             = P.depth;
    albedo_bucket     = 0.2;
    cover_top         = P.cover_top;
    shading           = P.solar_shading;

    % Parameter of measurement --------------------------------------------
    t_haul = 60;
    t_deck = P.deck_time;
    if P.exp_id == 1,
        disp('less exposure')
        u_shield_haul = 0.6;
        u_shield_deck = 0.25;
        s_shield_haul = 1;
        s_shield_deck = 0.42;
    else
        u_shield_haul = 0.6;
        u_shield_deck = 0.4;
        s_shield_haul = 1;
        s_shield_deck = 0.67;
    end
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
    
    % Set initial condition -----------------------------------------------
    SST = true_SST;
    mass_water = ones(size(true_SST)) * mass;
    SST_out    = nan([size(true_SST),time_out]);
    SST_out(:,:,:,:,1) = SST;
    ct_out = 1;

    % Compute the solar contaimination ------------------------------------
    clear('temp_cycle_direct','temp_cycle_base','temp_base')
    Cs_diffuse = Cs .* (1 - direct_ratio);
    Cs_direct  = Cs .* direct_ratio;
    temp_cycle_direct  = A_cycle ./pi .* (Cs_direct ./ cos(zenith_angle) .* sin(zenith_angle) .* (1-albedo_bucket));
    l          = abs(zenith_angle - pi/2) < 0.05   &  Cs_direct < 30;
    temp_cycle_direct(l) = 0;
    temp_cycle_diffuse = A_cycle .* (Cs_diffuse ./2 .* (1-albedo_bucket));
    temp_base  = A_base .* (Cs .* (1-cover_top));


    for t = 1:time_step

        e_bucket = 6.112 * exp(17.67 * (SST - 273.15)./(SST - 29.65)) .* 0.98;

        % Compute each term of the rediative budget -----------------------
        if((t * dt) <= t_haul) % During the process of hauling the bucket
            dQsdt = - (h_cycle_haul .* A_cycle + h_base_haul .* A_base) .* (SST - true_AT) .* sensible_on;
            dQwdt = -1.7 * (h_cycle_haul .* A_cycle + h_base_haul .* A_base) .* (e_bucket - e_air) .* latent_on;
        else            % During on deck stage
            dQsdt = - (h_cycle_deck .* A_cycle + h_base_deck .* A_base) .* (SST - true_AT) .* sensible_on;
            dQwdt = -1.7 * (h_cycle_deck .* A_cycle + h_base_deck .* A_base) .* (e_bucket - e_air) .* latent_on;
        end

        dQldt = - (A_base + A_cycle) .* sigma .* (SST.^4 - true_AT.^4) .* long_on;
        dQrdt_cycle = (temp_cycle_direct + temp_cycle_diffuse) .* (1 - shading) .* solar_wall_on;
        dQrdt_base  = (temp_base) .* (1 - shading) .* solar_top_on;

        % Compute the total rediative budget ------------------------------
        if ((t * dt) <= t_haul || (t * dt) >=t_haul + 30)
            dQdt = dQsdt + dQwdt + dQldt + dQrdt_cycle + dQrdt_base;
        else
            dQdt = dQsdt + dQwdt + dQldt + dQrdt_cycle + dQrdt_base - 0.0012.* (SST - true_AT);
        end

        % Compute the leakage of mass -------------------------------------
        if (t * dt) <= t_haul,
            mass_water = mass_water - leakage_rate_haul./60.*A_base.*density .* dt;
        else
            mass_water = mass_water - leakage_rate_deck./60.*A_base.*density .* dt;
        end

        % Compute temperature tendency ------------------------------------
        if (t * dt) <= t_haul,
            dTdt = dQdt ./ (mass_bucket.*cp_bucket + mass_water.*cp_water);
        else
            dTdt = dQdt ./ (mass_bucket.*cp_bucket + (mass_water + water_amount/1000).*cp_water);
        end

        SST = SST + dTdt * dt;
        
        if rem(t*dt,30) == 0,
            ct_out = ct_out + 1;
            SST_out(:,:,:,:,ct_out) = SST;
        end

        if rem(t*dt,30) == 0,
            disp([num2str(t*dt),' seconds finished'])
        end
    end
end