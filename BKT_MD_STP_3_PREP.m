% This function prepares for the enviromental variables to run the model 
% You can choose from two drivers:
% 1. ERA_interim: 1985-2014 renalysis
% 2. ICOADS3.0:   1950-1990 climatology on deck 
function [true_SST,true_AT,e_air,u_environment,Qs,direct_ratio,zenith_angle] = BKT_MD_STP_3_PREP(mode)

    dir_driver = BKT_OI('save_driver');

    % Prepare for the driver ----------------------------------------------
    if mode == 1,
        dew = load([dir_driver,'ERI-interim_5X5_d2m_1985_2014.mat']);
        air = load([dir_driver,'ERI-interim_5X5_t2m_1985_2014.mat']);
        wnd = load([dir_driver,'ERI-interim_5X5_Wnd_spd_1985_2014.mat']);
    else
        dew = load([dir_driver,'ICOADS_5X5_DPT_1950-1990.mat']);
        air = load([dir_driver,'ICOADS_5X5_AT_1950-1990.mat']);
        wnd = load([dir_driver,'ICOADS_5X5_WS_1950-1990.mat']);
        dew.clim_final = dew.clim_final + 273.41;
        air.clim_final = air.clim_final + 273.41;
    end
    ssrd = load([dir_driver,'ERI-interim_5X5_ssrd_1985_2014.mat']);
    ca   = load([dir_driver,'ICOADS_5X5_CA_1950-1990.mat']);
    sst  = load([dir_driver,'OI_SST_5X5_SST_1982-2014.mat']);

    dew = dew.clim_final;
    air = air.clim_final;
    wnd = wnd.clim_final;
    ssrd = ssrd.clim_final;
    ca = ca.clim_final;
    sst = sst.clim_final;

    mask_dew = mean(mean(dew,4),3);
    mask_air = mean(mean(air,4),3);
    mask_wnd = mean(mean(wnd,4),3);
    mask_ssrd = mean(mean(ssrd,4),3);
    mask_ca = mean(mean(ca,4),3);
    mask_sst = mean(mean(sst,4),3);

    mask = ~isnan(mask_dew) & ~isnan(mask_air) & ~isnan(mask_ssrd) &...
        ~isnan(mask_wnd) & ~isnan(mask_ca) & ~isnan(mask_sst);
    clear('mask_wnd','mask_dew','mask_air','mask_ssrd','mask_ca','mask_sst')

    mask = repmat(mask,1,1,24,12);
    dew(mask == 0) = nan;
    air(mask == 0) = nan;
    sst(mask == 0) = nan;
    ssrd(mask == 0) = nan;
    wnd(mask == 0) = nan;
    ca(mask == 0) = nan;
    clear('mask')

    % Prepare for the data fed into the model ---------------------------------
    true_SST = sst + 273.15;
    true_AT  = air;
    true_DT  = dew;
    u_environment = wnd;
    Qs       = ssrd;
    cld_tr   = ca + 1;
    clear('air','dew','sst','wnd','ca','ssrd')

    % Compute zenith angle ------------------------------------------------
    if 0,
        clear('zenith_angle')
        zenith_angle = nan(size(sst));
        for lon = 1
            for lat = 1:36
                for lcl = 1:24
                    for mon = 1:12
                        clear('location','time','sun')
                        location.longitude = lon*5 - 2.5;
                        location.latitude = lat*5 - 92.5;
                        location.altitude = 0;
                        time.year = 2010;
                        time.month = mon;
                        time.day = 15;
                        time.hour = lcl;
                        time.min = 0;
                        time.sec = 0;
                        time.UTC = location.longitude./15;
                        sun = sun_position(time, location);
                        zenith_angle(lon,lat,lcl,mon) = sun.zenith;
                    end
                end
            end
        end
        zenith_angle(2:72,:,:,:) = repmat(zenith_angle(1,:,:,:),71,1,1,1);
        save([dir_driver,'Zenith_angle.mat'],'zenith_angle','-v7.3');
    else
        load([dir_driver,'Zenith_angle.mat']);
    end
    logic_day = zenith_angle < 90;
    zenith_angle = zenith_angle /180 *pi;

    % Compute what ratio is diffused and what ratio is direct -------------
    clear_vertical_ratio = 0.95;
    Ld = 1 ./ ( -log(clear_vertical_ratio));
    L_true = cld_tr ./ cos(zenith_angle);
    L_true(logic_day == 0) = 999999999;
    direct_ratio = exp(- L_true./Ld);

    % Qs that are negative or during the night are 0 ----------------------
    Qs((logic_day == 0 | Qs < 0) & ~isnan(Qs)) = 0;

    % Compute the water vapor pressure in the air -------------------------
    e_air = 6.112 .* exp(17.67 .* (true_DT - 273.15)./(true_DT - 29.65));
end