% This function prepares for the enviromental variables to run the model
% You can choose from two drivers:
% 1. ERA_interim: 1985-2014 renalysis
% 2. ICOADS3.0:   1950-1990 climatology on deck (bucket records only)
% 3. ICOADS3.0:   1973-2002 climatology on deck (bucket records only)
% 4. NOCS:        1973-2002 climatology all at 10m
function [true_SST,true_AT,e_air,u_environment,Qs,Direct_ratio,zenith_angle] = BKT_MD_STP_3_PREP_2019(mode,env)

    if ~exist('env','var'), env = 1; end
    dir_driver = BKT_OI('save_driver',env);

    % Prepare for the driver ----------------------------------------------
    % diurnal cycle of SST estimated from buoy data
    if mode == 1,
        % diurnal cycle from ERA-interim
        dew = load([dir_driver,'ERI-interim_5X5_d2m_1985_2014.mat']);
        air = load([dir_driver,'ERI-interim_5X5_t2m_1985_2014.mat']);
        wnd = load([dir_driver,'ERI-interim_5X5_Wnd_spd_1985_2014.mat']);
        sst  = load([dir_driver,'OI_SST_5X5_SST_1982-2014.mat']);
    elseif mode == 2,
        % diurnal cycle of air variables estimated by ourselves from bucket records
        dew = load([dir_driver,'ICOADS_5X5_DPT_1950-1990.mat']);
        air = load([dir_driver,'ICOADS_5X5_AT_1950-1990.mat']);
        wnd = load([dir_driver,'ICOADS_5X5_WS_1950-1990.mat']);
        dew.clim_final = dew.clim_final + 273.15;
        air.clim_final = air.clim_final + 273.15;
        sst  = load([dir_driver,'OI_SST_5X5_SST_1982-2014.mat']);
    elseif mode == 3,
        % diurnal cycle of air variables estimated by ourselves from bucket records
        dew = load([dir_driver,'ICOADS_5X5_DPT_1973-2002.mat']);
        air = load([dir_driver,'ICOADS_5X5_AT_1973-2002.mat']);
        wnd = load([dir_driver,'ICOADS_5X5_WS_1973-2002.mat']);
        dew.clim_final = dew.clim_final + 273.15;
        air.clim_final = air.clim_final + 273.15;
        sst  = load([dir_driver,'OI_SST_5X5_SST_1982-2014.mat']);
    elseif mode == 4,
        % diurnal cycle of air variables estimated by ourselves from bucket records
        dew = load([dir_driver,'NOCS_5X5_DPT_1973-2002.mat']);
        air = load([dir_driver,'NOCS_5X5_AT_1973-2002.mat']);
        wnd = load([dir_driver,'NOCS_5X5_WS_1973-2002.mat']);
        dew.clim_final = dew.clim_final + 273.35;
        air.clim_final = air.clim_final + 273.35;
        sst  = load([dir_driver,'NOCS_5X5_SST_1973-2002.mat']);
    end
    ssrd = load([dir_driver,'ERI-interim_5X5_ssrd_1985_2014.mat']);


    dew = dew.clim_final;
    air = air.clim_final;
    wnd = wnd.clim_final;
    ssrd = ssrd.clim_final;
    sst = sst.clim_final;

    mask_dew = mean(mean(dew,4),3);
    mask_air = mean(mean(air,4),3);
    mask_wnd = mean(mean(wnd,4),3);
    mask_ssrd = mean(mean(ssrd,4),3);
    mask_sst = mean(mean(sst,4),3);

    mask = ~isnan(mask_dew) & ~isnan(mask_air) & ~isnan(mask_ssrd) &...
        ~isnan(mask_wnd) & ~isnan(mask_sst);
    clear('mask_wnd','mask_dew','mask_air','mask_ssrd','mask_ca','mask_sst')

    mask = repmat(mask,1,1,24,12);
    dew(mask == 0) = nan;
    air(mask == 0) = nan;
    sst(mask == 0) = nan;
    ssrd(mask == 0) = nan;
    wnd(mask == 0) = nan;
    clear('mask')

    % Prepare for the data fed into the model ---------------------------------
    if nanmean(sst(:)) < 200,
        true_SST = sst + 273.15;
    end
    true_AT  = air;
    true_DT  = dew;
    u_environment = wnd;
    Qs       = ssrd;
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

    % Compute ratio of ground insolation to toa insolation ----------------
    if 0,
        mon_list = [31 28 31 30 31 30 31 31 30 31 30 31];
        day_list = [0 cumsum(mon_list)];
        td       = (day_list(1:end-1) + day_list(2:end))/2;
        td       = repmat(reshape(td,1,1,1,12),72,36,24,1);
        SC       = 1367;
        
        lat = repmat(reshape(-87.5:5:87.5,1,36,1,1),72,1,24,12);
        th  = repmat(reshape(1:24,1,1,24,1),72,36,1,12);
        sin_delta = -sin(23.45/180*pi).*cos(2*pi*(td+10)/365);
        cos_delta = sqrt(1-sin_delta.^2);
        sin_beta = sin(lat/180*pi).*sin_delta + cos(lat/180*pi).*cos_delta.*cos(15*(th-12)/180*pi);
        
        S0       = SC .* (1 + 0.033 * cos(2 * pi * td / 365)) .* sin_beta; %.* cos(zenith_angle);
        S0(S0<0) = 0;
        Qs_vs_S0 = Qs./S0;
        direct_ratio = 2 * Qs_vs_S0 - 0.7;
        direct_ratio (direct_ratio > 0.95) = 0.95;
        direct_ratio (direct_ratio < 0) = 0;
        direct_ratio (Qs == 0) = 0;
        
        % Qs that are negative or during the night are 0 ------------------
        Qs((logic_day == 0 | Qs < 0) & ~isnan(Qs)) = 0;
        direct_ratio(Qs == 0) = 0;
        
        % Post-processing data --------------------------------------------
        Direct_ratio = direct_ratio;
        hr = 1:24;
        for id_lon = 1:72
            for id_lat = 1:36
                for id_mon = 1:12
                    a = squeeze(direct_ratio(id_lon,id_lat,:,id_mon))';
                    clear('l_ex')
                    l_ex = false(1,24);
                    for i = 14:1:22    if a(i+1) > a(i), l_ex(i+1) = 1;  end;  end
                    for i = 10:-1:2    if a(i-1) > a(i), l_ex(i-1) = 1;  end;  end
                    if ~all(isnan(a(~l_ex))),
                        b = interp1(hr(~l_ex),a(~l_ex),hr,'spline');
                        Direct_ratio(id_lon,id_lat,:,id_mon) = b;
                    end
                end
            end
        end
        Direct_ratio(Direct_ratio < 0) = 0;
        save([dir_driver,'Direct_ratio.mat'],'Direct_ratio','-v7.3');
    else
        load([dir_driver,'Direct_ratio.mat']);
    end
    

    % Compute the water vapor pressure in the air -------------------------
    e_air = 6.112 .* exp(17.67 .* (true_DT - 273.15)./(true_DT - 29.65));
end
