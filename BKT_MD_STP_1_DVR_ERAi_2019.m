% This file prepares the hourly climatology of t2m, d2m, insolation, and
% 10m wind speed from 30yr ERA-interim 3-hourly data (1985-2014).
%
% The output of this file are not rearranged by hours, and is at the raw
% resolution (2x2) from the ERA-interim.
%
% The output units are:
% t2m(K), d2m(K), insolation(w/m2), wind speed(m/s)
%
% These outputs are further processed by "BKT_MD_STP_1_DVR_SUM.m"
% in order to get a 5x5 value with data rearranged by hours.


% #############################################################
%  Compute the climatology of t2m, d2m, and insolation       ##
% #############################################################
clear;
dir_load  = BKT_OI('load_ERA');
file_load = [dir_load,'bucket_model_environment_1985_2014.nc'];
var_list = {'ssrd','d2m','t2m'};

for var = 1:3

    % ************
    % read data **
    % ************
    varname = var_list{var};
    disp(varname);
    clear('temp','temp_clim','temp_clim2','temp_mon','temp_mon_clim','temp_hour_clim')
    for i = 1:88

        tic;
        disp(['i = ',num2str(i)]);
        if i == 88; len = 656; else len = 1000; end

        clear('temp_1')
        temp_1 = ncread(file_load,varname,[1 1 (i-1)*1000+1],[180 91 len]);
        temp(:,:,(i-1)*1000+1:min(i*1000,87656)) = temp_1;
        % temp(:,:,(i-1)*1000+1:min(i*1000,87656)) = temp_1(1:8:end,1:8:end,:);   % for debug
        toc
    end
    longitude = ncread(file_load,'longitude');
    latitude = ncread(file_load,'latitude');

    % ************************
    % remove the leap years **
    % ************************
    disp('remove leap years');
    logic = false(1,size(temp,3));
    for lpyr = 1:7
        logic((59 + (4*(lpyr) - 1)*365)*8 + [1:8]) = true;
    end
    temp(:,:,logic) = [];

    % ****************************
    % compute daily climatology **
    % ****************************
    disp('compute daily climatology');
    temp_clim = nan(size(temp,1),size(temp,2),365*8);
    for tz = 1:365*8
        temp_clim(:,:,tz) = nanmean(temp(:,:,tz:365*8:end),3);
    end

    % ************************************
    % for ssrd, convert to unit of w/m2 **
    % ************************************
    temp_clim_2 = nan(size(temp_clim));
    if(strcmp(varname,'ssrd'))
        temp_clim_2(:,:,1:4:end) = temp_clim(:,:,1:4:end);
        temp_clim_2(:,:,2:4:end) = temp_clim(:,:,2:4:end) - temp_clim(:,:,1:4:end);
        temp_clim_2(:,:,3:4:end) = temp_clim(:,:,3:4:end) - temp_clim(:,:,2:4:end);
        temp_clim_2(:,:,4:4:end) = temp_clim(:,:,4:4:end) - temp_clim(:,:,3:4:end);
        temp_clim_2 = temp_clim_2 / 3600/3;
    else
        temp_clim_2 = temp_clim;
    end
    clear('temp_clim','temp')

    % ***************************
    % compute the monthly mean **
    % ***************************
    disp('compute monthly climatology');
    temp_mon_clim = nan(size(temp_clim_2,1),size(temp_clim_2,2),96);
    mon_list = [31 28 31 30 31 30 31 31 30 31 30 31];
    mon_cum  = cumsum([0 mon_list])*8;
    for mon = 1:12
        clear('temp_mon')
        temp_mon = temp_clim_2(:,:,mon_cum(mon)+1:mon_cum(mon+1));
        for tz = 1:8
            temp_mon_clim(:,:,tz+(mon-1)*8) = nanmean(temp_mon(:,:,tz:8:end),3);
        end
    end

    % *******************************************
    % generate hourly data using interpolation **
    % *******************************************
    disp('compute hourly climatology');
    clear('temp_hour_clim')
    temp_hour_clim = nan(size(temp_mon_clim,1),size(temp_mon_clim,2),24,12);
    for lon = 1:size(temp_mon_clim,1)
        for lat = 1:size(temp_mon_clim,2)
            for mon = 1:12
                clear('tem')
                tem = squeeze(temp_mon_clim(lon,lat,[1:8]+8*(mon-1)));
                temp_hour_clim(lon,lat,[1:24],mon) = interp1([1.5:3:72],[tem;tem;tem],[1:24]+24,'spline');
                % temp_hour_clim(lon,lat,[1:24],mon) = interp1([2:3:72],[tem;tem;tem],[1:24]+24,'spline');
            end
        end
    end

    % ************
    % save data **
    % ************
    disp('saving data');
    dir_save = BKT_OI('save_driver',0);
    save([dir_save,'ERI-interim_2X2_',varname,'_1985_2014.mat'],'longitude','latitude','temp_hour_clim','-v7.3');
end

% #############################################################
%  Compute the climatology of 10 meter wind speed            ##
% #############################################################
clear;
dir_load  = BKT_OI('load_ERA');
file_load = [dir_load,'bucket_model_environment_1985_2014.nc'];
var_list = {'u10','v10'};

% ************
% read data **
% ************
disp('Wind speed');
for i = 1:88

    tic;
    disp(['i = ',num2str(i)]);
    if i == 88; len = 656; else len = 1000; end

    clear('temp_1','temp_2')
    var = 1;
    varname = var_list{var};
    temp_1 = ncread(file_load,varname,[1 1 (i-1)*1000+1],[180 91 len]);

    var = 2;
    varname = var_list{var};
    temp_2 = ncread(file_load,varname,[1 1 (i-1)*1000+1],[180 91 len]);

    temp(:,:,(i-1)*1000+1:min(i*1000,87656)) = sqrt(temp_1.^2 + temp_2.^2);
    clear('temp_1','temp_2')

    toc
end
longitude = ncread(file_load,'longitude');
latitude = ncread(file_load,'latitude');

% ************************
% remove the leap years **
% ************************
disp('remove leap years');
logic = false(1,size(temp,3));
for lpyr = 1:7
    logic((59 + (4*(lpyr) - 1)*365)*8 + [1:8]) = true;
end
temp(:,:,logic) = [];


% ****************************
% compute daily climatology **
% ****************************
disp('compute daily climatology');
temp_clim = nan(size(longitude,1),size(latitude,1),365*8);
for tz = 1:365*8
    temp_clim(:,:,tz) = nanmean(temp(:,:,tz:365*8:end),3);
end
temp_clim_2 = nan(size(temp_clim));
temp_clim_2 = temp_clim;
clear('temp_clim','temp')

% ***************************
% compute the monthly mean **
% ***************************
disp('compute monthly climatology');
temp_mon_clim = nan(size(longitude,1),size(latitude,1),96);
for mon = 1:12
    clear('temp_mon')
    temp_mon = temp_clim_2(:,:,240*(mon-1)+[1:240]);
    for tz = 1:8
        temp_mon_clim(:,:,tz+(mon-1)*8) = nanmean(temp_mon(:,:,tz:8:end),3);
    end
end

% *******************************************
% generate hourly data using interpolation **
% *******************************************
disp('compute hourly climatology');
clear('temp_hour_clim')
temp_hour_clim = nan(size(longitude,1),size(latitude,1),24,12);
for lon = 1:180
    for lat = 1:91
        for mon = 1:12
            clear('tem')
            tem = squeeze(temp_mon_clim(lon,lat,[1:8]+8*(mon-1)));
            temp_hour_clim(lon,lat,[1:24],mon) = interp1([2:3:72],[tem;tem;tem],[1:24]+24,'spline');
        end
    end
end

disp('saving data');
dir_save = BKT_OI('save_driver');
save([dir_save,'ERI-interim_2X2_Wnd_spd_1985_2014.mat'],'longitude','latitude','temp_hour_clim','-v7.3');
