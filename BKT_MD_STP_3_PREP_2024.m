% This function prepares for the enviromental variables to run the model
% This version simply uses the most updated high-resolution driver (1x1) in 2024
function [true_SST,true_AT,e_air,u_environment,Qs,Direct_ratio,zenith_angle] = BKT_MD_STP_3_PREP_2024

    dir_driver = BKT_OI('save_driver');

    % Prepare for the driver ----------------------------------------------
    % diurnal cycle of SST estimated from buoy data
    % diurnal cycle of environmental variables estimated by ourselves
    % from ships taking bucket SSTs, 2024 updated version!
    disp('The most updated version of driver')
    dew  = load([dir_driver,'NOCS_1X1_DPT_1970_2014_new_version_2024.mat']);
    air  = load([dir_driver,'NOCS_1X1_AT_1970_2014_new_version_2024.mat']);
    wnd  = load([dir_driver,'NOCS_1X1_W_1970_2014_new_version_2024.mat']);
    dew.clim_final = dew.clim_final + 273.15;
    air.clim_final = air.clim_final + 273.15;
    sst  = load([dir_driver,'NOCS_1X1_SST_1970_2014_new_version_2024.mat']);
    ssrd = load([dir_driver,'ERI5_1X1_ssrd_1973_2006.mat']);

    dew  = double(dew.clim_final);
    air  = double(air.clim_final);
    wnd  = double(wnd.clim_final);
    ssrd = double(ssrd.clim_final);
    sst  = double(sst.clim_final);

    mask_dew  = mean(mean(dew,4),3);
    mask_air  = mean(mean(air,4),3);
    mask_wnd  = mean(mean(wnd,4),3);
    mask_ssrd = mean(mean(ssrd,4),3);
    mask_sst  = mean(mean(sst,4),3);

    mask      = ~isnan(mask_dew) & ~isnan(mask_air) & ~isnan(mask_ssrd) & ...
        ~isnan(mask_wnd) & ~isnan(mask_sst);
    clear('mask_wnd','mask_dew','mask_air','mask_ssrd','mask_ca','mask_sst')

    mask = repmat(mask,1,1,24,12);
    dew(mask == 0)  = nan;
    air(mask == 0)  = nan;
    sst(mask == 0)  = nan;
    ssrd(mask == 0) = nan;
    wnd(mask == 0)  = nan;
    clear('mask')

    % Prepare for the data fed into the model -----------------------------
    if nanmean(sst(:)) < 200
        true_SST = sst + 273.15;
    end
    true_AT       = air;
    true_DT       = dew;
    u_environment = wnd;
    Qs            = ssrd;
    clear('air','dew','sst','wnd','ca','ssrd')

    % Compute zenith angle ------------------------------------------------
    if ~isfile([dir_driver,'Zenith_angle_1X1.mat'])
        clear('zenith_angle')
        zenith_angle = nan(size(true_SST));
        for lon = 1
            for lat = 1:180
                for lcl = 1:24
                    for mon = 1:12
                        clear('location','time','sun')
                        location.longitude = lon*1 - 0.5;
                        location.latitude = lat*1 - 90.5;
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
        zenith_angle(2:360,:,:,:) = repmat(zenith_angle(1,:,:,:),359,1,1,1);
        zenith_angle = single(zenith_angle);
        save([dir_driver,'Zenith_angle_1X1.mat'],'zenith_angle','-v7.3');
    else
        load([dir_driver,'Zenith_angle_1X1.mat']);
        zenith_angle = double(zenith_angle);
    end
    logic_day = zenith_angle < 90;
    zenith_angle = zenith_angle /180 *pi;

    % Compute ratio of ground insolation to toa insolation ----------------
    if ~isfile([dir_driver,'Direct_ratio_1X1.mat'])
        reso     = 1;
        mon_list = [31 28 31 30 31 30 31 31 30 31 30 31];
        day_list = [0 cumsum(mon_list)];
        td       = (day_list(1:end-1) + day_list(2:end))/2;
        td       = repmat(reshape(td,1,1,1,12),360/reso,180/reso,24,1);
        SC       = 1367;
        
        lat = repmat(reshape(-89.5:1:89.5,1,180/reso,1,1),360/reso,1,24,12);
        th  = repmat(reshape(1:24,1,1,24,1),360/reso,180/reso,1,12);
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
        for id_lon = 1:360
            for id_lat = 1:180
                for id_mon = 1:12
                    a = squeeze(direct_ratio(id_lon,id_lat,:,id_mon))';
                    clear('l_ex')
                    l_ex = false(1,24);
                    for i = 14:1:22    if a(i+1) > a(i), l_ex(i+1) = 1;  end;  end
                    for i = 10:-1:2    if a(i-1) > a(i), l_ex(i-1) = 1;  end;  end
                    if ~all(isnan(a(~l_ex)))
                        b = interp1(hr(~l_ex),a(~l_ex),hr,'spline');
                        Direct_ratio(id_lon,id_lat,:,id_mon) = b;
                    end
                end
            end
        end
        Direct_ratio(Direct_ratio < 0) = 0;
        Direct_ratio = single(Direct_ratio);
        save([dir_driver,'Direct_ratio_1X1.mat'],'Direct_ratio','-v7.3');
    else
        load([dir_driver,'Direct_ratio_1X1.mat']);
        Direct_ratio = double(Direct_ratio);
    end
   
    % Compute the water vapor pressure in the air -------------------------
    e_air = 6.112 .* exp(17.67 .* (true_DT - 273.15)./(true_DT - 29.65));                      
end