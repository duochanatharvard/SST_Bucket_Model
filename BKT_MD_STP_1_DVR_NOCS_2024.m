% This script is to calculate the climatology of at, cldc, qair, slp, sst,
% and wspd from the NOCSV2 dataset
%
% Step 1: Calculate Climatological Values from NOCS dataset at 1x1 reso
clear;
dir = '/Users/dc1e23/Downloads/nocs_v2/';
dir_save = ['/Users/dc1e23/Dropbox/',...
                '------------ Git_code ----------------/',...
                'ZZ_Completed/SST_Bucket_Model/NOCSV2_2024/'];

var_list = {'at','cldc','qair','slp','sst','wspd'};

for ct_var = 1:6

    clearvars -except dir dir_save var_list ct_var
    varname = var_list{ct_var};
    
    % Calculate the mean value directly -----------------------------------
    clim    = zeros(360,180,12);
    N       = 0;
    for yr  = 1973:2006
        disp(num2str(yr,'year: %6.0f'))
        fname = [dir,'nocs_v2_0_',varname,'_',num2str(yr),'.nc'];
        temp  = ncread(fname,varname);
        temp  = temp([181:360 1:180],:,:);
        clim  = (clim * N + temp) ./ (N+1);
        N     = N + 1;
    end
    
    % Do you need to do any smoothing? ------------------------------------
    for ct = 1:12
        clim_final(:,:,ct) = CDC_smooth2(clim(:,:,ct),3,0);
    end
    
    fsave = [dir_save,'NOCS_1X1_',varname,'_1973_2006.mat'];
    save(fsave,'clim','clim_final','-v7.3');
end

%% ######################################################################## 
% Step 2: Add diurnal cycles

clear;

% 2.1: Load in daily mean (DM) values
dir_load = ['/Users/dc1e23/Dropbox/',...
                '------------ Git_code ----------------/',...
                'ZZ_Completed/SST_Bucket_Model/NOCSV2_2024/'];

var_list = {'at','wspd','qair','sst'};

for var_id = 1:4

    varname = var_list{var_id};
    file_load = [dir_load,'NOCS_1x1_',varname,'_1973_2006.mat'];
    load(file_load);
    temp = clim_final;
    clear('clim_final')

    if var_id == 3

        % for consistency, solve for the dew point temperature
        temp = temp / 1000;

        % Saturation vapor pressure
        f_es_t = @(t) 6.112 .* exp(17.67 .* (t-273.15)./(t-29.65));

        % saturation specific humidity from temperature
        f_q_t = @(t,p) f_es_t(t) .* 0.622 ./ (p - 0.378 * f_es_t(t));

        % load the surface pressure
        file_load = [dir_load,'NOCS_1X1_slp_1973_2006.mat'];
        load(file_load);
        p = clim_final;
        clear('clim_final');

        % solve for dew-point temperature temp_x
        temp_x = zeros(size(temp)) + 273.15;

        err = 1;
        while err > 1e-6
            q_x = f_q_t(temp_x,p);
            temp_x = temp_x - (q_x - temp) * 10;
            err =  max(abs(q_x(:) - temp(:)));
        end
        temp = temp_x - 273.15;

    end

    CLIM_DM(:,:,:,var_id) = temp;  % DM stands for daily mean ...
end

% 2.2 Load in Diurnal Signals ---------------------------------------------
dir_ICOADS = BKT_OI('load_ICOADS');
var_list = {'AT','W','DPT','SST'};
CLIM_DA = nan(72,36,12,4);
SHP_DA  = nan(24,36,12,4);
for var_id = 1:4

    clear('a')
    if var_id ~=4
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

% Interpolate to 1 degree higher resolution
P.threshold = 1;
for ct_mon = 1:12
    for ct_var = 1:4
        temp                         = CDC_interp_high_reso(2.5:5:357.5,-87.5:5:87.5,CLIM_DA(:,:,ct_mon,ct_var),0.25:0.5:360,-89.75:0.5:90,'Ocean','linear');
        a                            = CDC_average_grid(0.25:0.5:360,-89.75:0.5:90,temp,0.5:1:359.5,-89.5:1:89.5,P);
        CLIM_DA_h(:,:,ct_mon,ct_var) = CDC_smooth2(a,3,1);
    end
end

SHP_DA_h = permute(interp1(-87.5:5:87.5,permute(SHP_DA,[2 1 3 4]),-89.5:1:89.5,'linear','extrap'),[2 1 3 4]);

% 2.3 combine daily mean with the diurnal cycle ---------------------------
for var_id = 2:4

    % **********************************
    % read data and the diurnal shape **
    % **********************************
    shape = SHP_DA_h(:,:,:,var_id);

    % ***************************************
    % Get the daily mean and diurnal cycle **
    % ***************************************
    clear('temp_dm','temp_da')
    temp_dm = CLIM_DM(:,:,:,var_id);
    temp_da = CLIM_DA_h(:,:,:,var_id);

    % **************************************
    % Combine the two components together **
    % **************************************
    clear('clim_final')
    clim_final = nan(360,180,24,12);
    for i = 1:360
        for j = 1:180
            for mon = 1:12
                clim_final(i,j,:,mon) = temp_dm(i,j,mon) + temp_da(i,j,mon) .* shape(:,j,mon);
            end
        end
    end

    % ************
    % save data **
    % ************
    clim_final = single(clim_final);
    file_save = [BKT_OI('save_driver'),'NOCS_1X1_',var_list{var_id},'_1970_2014_new_version_2024.mat'];
    save(file_save,'clim_final','-v7.3');
end

%2.4 combine daily mean with the diurnal cycle - AT needs some spetial treatment 
var_id = 1;
shape_at  = SHP_DA_h(:,:,:,1);
shape_sst = SHP_DA_h(:,:,:,4);

% ***************************************
% Get the daily mean and diurnal cycle **
% ***************************************
clear('temp_dm_at','temp_da_at')
temp_dm_at = CLIM_DM(:,:,:,1);
temp_da_at = CLIM_DA_h(:,:,:,1);
clear('temp_dm_sst','temp_da_sst')
temp_dm_sst = CLIM_DM(:,:,:,4);
temp_da_sst = CLIM_DA_h(:,:,:,4);

% **************************************
% Combine the two components together **
% **************************************
clear('clim_final')
clim_final = nan(360,180,24,12);
for i = 1:360
    for j = 1:180
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
clim_final = single(clim_final);
file_save = [BKT_OI('save_driver'),'NOCS_1X1_',var_list{var_id},'_1970_2014_new_version_2024.mat'];
save(file_save,'clim_final','-v7.3');
