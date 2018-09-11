% This file prepares the hourly climatology of t10m, d10m, cloud cover, and 
% 10m wind speed from 30yr monthly NOCS dataset (1973-2002).
% 
% We also add 30-yr diurnal cycle estimates
% from ICOADS (1973-2002) is added to the monthly value to drive the 
% bucket model. 
% 
% The output units are:
% t10m(C), d10m(C), sst(C), cld cover(octo), wind speed 10m (m/s)
%
% These outputs are further processed by "BKT_MD_STP_1_DVR_SUM.m"
% in order to get a 5x5 value with data rearranged by hours.


% ##########################################
%  Load monthly climatology from NOCSv2   ##
% ##########################################
clear;
dir_load  = BKT_OI('load_NOCS');

var_list = {'at','wspd','cldc','qair','slp','sst'};

for var_id = 1:6
    
    varname = var_list{var_id};
    file_load = [dir_load,'NOCS_5X5_',varname,'_1973_2002.mat'];
    load(file_load);
    temp = clim_final;
    clear('clim_final')
    
    if var_id == 3,
        temp = temp/12.5;
        
    elseif var_id == 4,
        % for consistency, solve for the dew point temperature
        temp = temp / 1000;
        
        % Saturation vapor pressure 
        f_es_t = @(t) 6.112 .* exp(17.67 .* (t-273.15)./(t-29.65));
        
        % saturation specific humidity from temperature
        f_q_t = @(t,p) f_es_t(t) .* 0.622 ./ (p - 0.378 * f_es_t(t));
        
        % load the surface pressure
        file_load = [dir_load,'NOCS_5X5_',var_list{5},'_1973_2002.mat'];
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
var_list = {'AT','WS','CA','DPT','SLP','SST'};
do_more_samples = 0;

clear('CLIM_DA')
CLIM_DA = nan(72,36,12,6);
for var_id = [1:4 6]
    
    % *************
    % read files **
    % *************
    clear('DA_WM','temp','temp_num')
    if var_id < 6
        file_load = [dir_ICOADS,'DA_',var_list{var_id},...
            '_Gridded_BCK_do_more_samples_',num2str(do_more_samples),'.mat'];
        load(file_load,'DA_WM','DA_NUM')
        temp_da = DA_WM(:,:,(1973-1850)*12+1:(2002-1849)*12);
        temp_num = DA_NUM(:,:,(1973-1850)*12+1:(2002-1849)*12);
    else
        file_load = [dir_ICOADS,'DA_',var_list{var_id},...
            '_Gridded_BUOY_do_more_samples_',num2str(do_more_samples),'.mat'];
        load(file_load,'DA_WM','DA_NUM')
        temp_da = DA_WM(:,:,(1990-1850)*12+1:(2014-1849)*12);
        temp_num = DA_NUM(:,:,(1990-1850)*12+1:(2014-1849)*12);
    end
    
    for mon = 1:12
        
        % *****************************
        % Compute the annual average **
        % *****************************
        clear('temp_clim','temp_mon_num','tem')
        tem = temp_da(:,:,mon:12:end);
        temp_num(abs(tem) > 5) = nan;
        tem(abs(tem) > 5) = nan;
        temp_clim = nansum(tem.*temp_num(:,:,mon:12:end),3)./nansum(temp_num(:,:,mon:12:end),3);
        temp_mon_num = nansum(temp_num(:,:,mon:12:end),3);
        temp_clim(temp_mon_num <= 5) = nan;
        
        % ****************************************
        % rule out outliers using zonal average **
        % ****************************************
        clear('tem')
        tem = nanmean(temp_clim,1);
        tem = smooth(tem);
        tem = repmat(tem',72,1) * 1.5;  % set a threshold, can be arbitrary
        tem(tem < 0) = 0;
        temp_clim(abs(temp_clim) > tem) = NaN;
        
        % *************************************
        % Interpolate to get global coverage **
        % *************************************
        [mask,topo,coast] = CDF_land_mask(5,1,5,10);
        MASK = ~mask';
        
        clear('temp_int')
        in_var = temp_clim;
        lon = 2.5: 5 :360;
        temp_int = nan(size(in_var));
        for i = 1:size(in_var,2)
            temp = in_var(:,i)';
            logic = ~isnan(temp);
            if nnz(logic),
                temp_int(:,i) = interp1([lon(logic)-360 lon(logic) lon(logic)+360],...
                    [temp(logic) temp(logic) temp(logic)],lon);
            end
        end
        temp_int(MASK == 0) = NaN;
        temp_int = CDC_smooth2(temp_int);
        
        % *******************************
        % Fill latitude with no values **
        % *******************************
        clear('tem','tem_2')
        tem = nanmean(temp_int,1);
        tem_lat = 1:36;
        tem_2 = interp1(tem_lat(~isnan(tem)),tem(~isnan(tem)),tem_lat);
        tem_2(1:find(~isnan(tem_2),1,'first')) = tem_2(find(~isnan(tem_2),1,'first'));
        tem_2(find(~isnan(tem_2),1,'last'):end) = tem_2(find(~isnan(tem_2),1,'last'));
        
        for i = 1:36
            if isnan(tem(i))
                temp_int(:,i) = tem_2(i);
            end
        end
        temp_int(MASK == 0) = nan;
        
        % ************
        % save data **
        % ************   
        CLIM_DA(:,:,mon,var_id) = temp_int;
    end
end
clear('MASK','mask','mon','tem','tem_2','tem_lat','temp_clim','temp_int','temp_mon_num','temp_num')

% ############################################
% combine daily mean with the diurnal cycle ##
% ############################################
for var_id = [1:4 6]
    
    % **********************************
    % read data and the diurnal shape **
    % **********************************
    file_load = [dir_ICOADS,'Diurnal_Shape_',var_list{var_id},'.mat'];
    clear('Diurnal_Shape','shape')
    load(file_load)
    
    if var_id < 6,
        shape = squeeze(Diurnal_Shape(:,:,1,:));
    else
        shape = squeeze(Diurnal_Shape(:,:,3,:));
    end
    
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
                clim_final(i,j,:,mon) = temp_dm(i,j,mon) + temp_da(i,j,mon) .* shape(j,:,mon);
            end
        end
    end
    
    % ************
    % save data **
    % ************
    file_save = [BKT_OI('save_driver'),'NOCS_5X5_',var_list{var_id},'_1973-2002.mat'];
    save(file_save,'clim_final','-v7.3');
end