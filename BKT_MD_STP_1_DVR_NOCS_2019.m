% This file prepares the hourly climatology of t10m, d10m, cloud cover, and 
% 10m wind speed from 30yr monthly NOCS dataset (1973-2002).
% 
% For AT, DPT and W, I add 45-yr (1970-2014) diurnal cycle estimates
% from ICOADS ships taking bucket SSTs to drive the bucket model.
% 
% Need to figure out what are 

% 
% The output units are:
% t10m(C), d10m(C), sst(C), and wind speed 10m (m/s)
%
% These outputs are further processed by "BKT_MD_STP_1_DVR_SUM.m"
% in order to get a 5x5 value with data rearranged by hours.


% ##########################################
%  Load monthly climatology from NOCSv2   ##
% ##########################################
clear;
dir_load  = BKT_OI('load_NOCS');

var_list = {'at','wspd','qair','sst'};

for var_id = 1:4
    
    varname = var_list{var_id};
    file_load = [dir_load,'NOCS_5X5_',varname,'_1973_2002.mat'];
    load(file_load);
    temp = clim_final;
    clear('clim_final')
    
    if var_id == 3,
        
        % for consistency, solve for the dew point temperature
        temp = temp / 1000;
        
        % Saturation vapor pressure 
        f_es_t = @(t) 6.112 .* exp(17.67 .* (t-273.15)./(t-29.65));
        
        % saturation specific humidity from temperature
        f_q_t = @(t,p) f_es_t(t) .* 0.622 ./ (p - 0.378 * f_es_t(t));
        
        % load the surface pressure
        file_load = [dir_load,'NOCS_5X5_slp_1973_2002.mat'];
        load(file_load);
        p = clim_final;
        clear('clim_final');
        
        % solve for dew-point temperature temp_x
        temp_x = zeros(size(temp)) + 273.15;

        err = 1;
        while err > 1e-6;
            q_x = f_q_t(temp_x,p);
            temp_x = temp_x - (q_x - temp) * 10;
            err =  max(abs(q_x(:) - temp(:)));
        end
        temp = temp_x - 273.15;
        
    end
    
    CLIM_DM(:,:,:,var_id) = temp;
end


% #######################################
% ICOADS3.0 5X5 diurnal amplitude (DA) ##
% #######################################
dir_ICOADS = BKT_OI('load_ICOADS');
var_list = {'AT','W','DPT','SST'};
CLIM_DA = nan(72,36,12,4);
SHP_DA  = nan(24,36,12,4);
for var_id = [1:4]
    
    clear('a')
    if var_id ~=4,
        file_load = [dir_ICOADS,'Diurnal_Shape_ship_bucket_',var_list{var_id},'.mat'];
        a = load(file_load);
        SHP_DA(:,:,:,var_id)  = a.Diurnal_Shape;

        file_load = [dir_ICOADS,'Diurnal_Amplitude_Ship_bucket_',var_list{var_id},'_1970_2014_climatology.mat'];
        a = load(file_load);
        CLIM_DA(:,:,:,var_id) = a.Diurnal_clim;
    else
        file_load = [dir_ICOADS,'Diurnal_Shape_buoy_SST.mat'];
        a = load(file_load);
        SHP_DA(:,:,:,var_id)  = a.Diurnal_Shape;
        
        file_load = [dir_ICOADS,'Diurnal_Amplitude_buoy_SST_1990_2014_climatology.mat'];
        a = load(file_load);
        CLIM_DA(:,:,:,var_id) = a.Diurnal_clim_buoy_1990_2014;
    end
end


% ############################################
% combine daily mean with the diurnal cycle ##
% ############################################
for var_id = [2:4]
    
    % **********************************
    % read data and the diurnal shape **
    % **********************************
    shape = SHP_DA(:,:,:,var_id);
    
    % ***************************************
    % Get the daily mean and diurnal cycle **
    % ***************************************
    clear('temp_dm','temp_da')
    temp_dm = CLIM_DM(:,:,:,var_id);
    temp_da = CLIM_DA(:,:,:,var_id);
    
    % **************************************
    % Combine the two components together **
    % **************************************
    clear('clim_final')
    clim_final = nan(72,36,24,12);
    for i = 1:72
        for j = 1:36
            for mon = 1:12
                clim_final(i,j,:,mon) = temp_dm(i,j,mon) + temp_da(i,j,mon) .* shape(:,j,mon);
            end
        end
    end
    
    % ************
    % save data **
    % ************
    file_save = [BKT_OI('save_driver',0),'NOCS_5X5_',var_list{var_id},'_1970_2014_new_version_2019.mat'];
    save(file_save,'clim_final','-v7.3');
end


% ################################################
% combine daily mean with the diurnal cycle -AT ##
% ################################################
% read data and the diurnal shape **
% **********************************
var_id = 1;
shape_at  = SHP_DA(:,:,:,1);
shape_sst = SHP_DA(:,:,:,4);

% ***************************************
% Get the daily mean and diurnal cycle **
% ***************************************
clear('temp_dm_at','temp_da_at')
temp_dm_at = CLIM_DM(:,:,:,1);
temp_da_at = CLIM_DA(:,:,:,1);
clear('temp_dm_sst','temp_da_sst')
temp_dm_sst = CLIM_DM(:,:,:,4);
temp_da_sst = CLIM_DA(:,:,:,4);

% **************************************
% Combine the two components together **
% **************************************
clear('clim_final')
clim_final = nan(72,36,24,12);
for i = 1:72
    for j = 1:36
        for mon = 1:12
            clear('ds_at','ds_sst','ds')
            ds_at  = temp_da_at(i,j,mon) .* shape_at(:,j,mon);
            ds_sst = temp_da_sst(i,j,mon) .* shape_sst(:,j,mon);
            ds = ds_at + nanmean(ds_sst([1:5 21:24]) - ds_at([1:5 21:24]));
            clim_final(i,j,:,mon) = temp_dm_at(i,j,mon) + ds;
        end
    end
end

% ************
% save data **
% ************
file_save = [BKT_OI('save_driver',0),'NOCS_5X5_',var_list{var_id},'_1970_2014_new_version_2019.mat'];
save(file_save,'clim_final','-v7.3');

