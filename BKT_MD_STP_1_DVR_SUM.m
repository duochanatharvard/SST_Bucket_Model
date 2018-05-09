% Sum up all the input files such that the bucket model can take these 
% enviromental values and run them directly

clear;

% ###########################################
% ERA-interim 5X5 hourly rearranged fields ##
% ###########################################
dir_ERA = BKT_OI('save_driver');
var_list = {'t2m','d2m','Wnd_spd','ssrd'};

for varname = 1:4
        
    % *************
    % read files **
    % *************
    disp('read data...')
    file_load = [dir_ERA,'ERI-interim_2X2_',var_list{varname},'_1985_2014.mat'];
    clear('longitude','latitude','temp_hour_clim','clim_mediate')
    load(file_load);
    
    % **********************
    % rearrange the field **
    % **********************
    disp('sort by local time')
    clim_mediate = nan(size(temp_hour_clim));
    for i = 1:180
        for j = 1:91
            clear('temp','temp_2','lcl')
            temp = squeeze(temp_hour_clim(i,j,:,:));
            utc = 0.5:1:23.5;
            lcl = rem(utc + longitude(i)/15,24);
            
            % --------------------
            % fit the local time |
            % --------------------
            a = repmat(lcl,1,5);
            b = find(a(2:end)-a(1:end-1) < -20);
            for bb = 1:numel(b)
                a(b(bb)+1:end) = a(b(bb)+1:end)+24;
            end
            a = a - 48;
            
            % ------------------------------------
            % interpolation as a way of rearrange|
            % ------------------------------------
            for mon = 1:12
                clear('tem')
                tem = temp(:,mon);
                temp_2(:,mon) = interp1(a,repmat(tem(:)',1,5),[1:24],'spline');
            end
            
            clim_mediate(i,j,:,:)  = temp_2;
        end
    end
    
    % **********************
    % regrid to 1x1 grids **
    % **********************
    disp('regrid to 1x1 boxes')
    m_proj_nml(1,[1 1 1 1 1 1]);
    [altitude,lon,lat] = m_elev([1 359 -90 90]);
    clear('clim_mediate2')
    clim_mediate2 = nan(360,180,24,12);
    for i = 1:24
        for j = 1:12
            clear('temp')
            temp = interp2(repmat([longitude-360; longitude; longitude+360]',91,1),...
                repmat([latitude],1,540), [clim_mediate(:,:,i,j); ...
                clim_mediate(:,:,i,j); clim_mediate(:,:,i,j)]',lon,lat,'spline');
            
            temp(altitude > 0) = NaN;
            clim_mediate2(:,:,i,j) = temp';
        end
    end

    % **********************
    % regrid to 5x5 grids **
    % ********************** 
    disp('regrid to 5x5 boxes')
    clear('clim_final')
    clim_final = nan(72,36,24,12);
    lat_1 = [-89.5:1:89.5];
    for i = 1:72
        for j =1:36
            clear('temp')
            temp = clim_mediate2((i-1)*5+[1:5],(j-1)*5+[1:5],:,:);
            temp_lat = repmat(lat_1((j-1)*5+[1:5]),5,1);
            temp_weigh = cos(temp_lat.*pi./180);
            temp_weigh(isnan(temp(:,:,1))) = nan;
        
            for m = 1:12
                for hr = 1:24
                    tem = temp(:,:,hr,m);
                    clim_final(i,j,hr,m) = nansum(tem(:).*temp_weigh(:)) ./ nansum(temp_weigh(:));
                end
            end
        end
    end

    % ************
    % save data **
    % ************
    dir_save = BKT_OI('save_driver');
    file_save = [dir_save,'ERI-interim_5X5_',var_list{varname},'_1985_2014.mat'];
    save(file_save,'clim_final','-v7.3');
end

% #######################################
% ICOADS3.0 5X5 diurnal amplitude (DA) ##
% #######################################
dir_ICOADS = BKT_OI('load_ICOADS');
var_list = {'AT','WS','CA','DPT','SLP','SST'};
do_more_samples = 0;

clear('CLIM_DA')
CLIM_DA = nan(72,36,12,6);
for varname = 1:6
    
    % *************
    % read files **
    % *************
    clear('DA_WM','temp','temp_num')
    if varname < 6
        file_load = [dir_ICOADS,'DA_',var_list{varname},...
            '_Gridded_BCK_do_more_samples_',num2str(do_more_samples),'.mat'];
        load(file_load,'DA_WM','DA_NUM')
        temp_da = DA_WM(:,:,(1950-1850)*12+1:(1990-1849)*12);
        temp_num = DA_NUM(:,:,(1950-1850)*12+1:(1990-1849)*12);
    else
        file_load = [dir_ICOADS,'DA_',var_list{varname},...
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
        temp_int = smooth2CD(temp_int);
        
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
        CLIM_DA(:,:,mon,varname) = temp_int;
    end
end
clear('MASK','mask','mon','tem','tem_2','tem_lat','temp_clim','temp_int','temp_mon_num','temp_num')


% #######################################
% ICOADS3.0 5X5 daily mean climatology ##
% #######################################
clear('CLIM_DM')
CLIM_DM = nan(72,36,12,6);
for varname = 1:5
    
    % ************
    % read data **
    % ************
    file_load = [dir_ICOADS,'SHIP_GRID_C0_',var_list{varname},'.mat'];
    load(file_load,'DATA','DATA_NUM')
    
    % ************************
    % subset 1950-1990 data **
    % ************************
    clear('temp','temp_num')
    temp_dm = DATA(:,:,:,[1950:1990]-1799,:);
    temp_num = DATA_NUM(:,:,:,[1950:1990]-1799,:);

    % ****************
    % daily average **
    % **************** 
    temp_dm = squeeze(nansum(temp_dm.*temp_num,3) ./ nansum(temp_num,3));
    temp_num = squeeze(nansum(temp_num,3));

    % ********************
    % 1950-1990 average **
    % ********************    
    temp_dm = squeeze(nansum(temp_dm.*temp_num,3) ./ nansum(temp_num,3));
    temp_num = squeeze(nansum(temp_num,3));
    
    temp_dm(temp_num <= 5) = nan;
    
    for mon = 1:12
        
        % *************************************
        % Interpolate to get global coverage **
        % *************************************
        [mask,topo,coast] = CDF_land_mask(5,1,5,10);
        MASK = ~mask';
        
        clear('temp_int')
        in_var = temp_dm(:,:,mon);
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
        temp_int = smooth2CD(temp_int);
        
        CLIM_DM(:,:,mon,varname) = temp_int;
    end
end

% **************
% read OI-SST **
% **************
varname = 6;
load([BKT_OI('save_driver'),'OI_monthly_1982_2014.mat']);
CLIM_DM(:,:,:,varname) = OI_month;

% ############################################
% combine daily mean with the diurnal cycle ##
% ############################################
for varname = 1:6
    
    % **********************************
    % read data and the diurnal shape **
    % **********************************
    file_load = [dir_ICOADS,'Diurnal_Shape_',var_list{varname},'.mat'];
    clear('Diurnal_Shape','shape')
    load(file_load)
    
    if varname < 6,
        shape = squeeze(Diurnal_Shape(:,:,1,:));
    else
        shape = squeeze(Diurnal_Shape(:,:,3,:));
    end
    
    % ***************************************
    % Get the daily mean and diurnal cycle **
    % ***************************************
    clear('temp_dm','temp_da')
    temp_dm = CLIM_DM(:,:,:,varname);
    temp_da = CLIM_DA(:,:,:,varname);
    
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
    if varname < 6,
        file_save = [BKT_OI('save_driver'),'ICOADS_5X5_',var_list{varname},'_1950-1990.mat'];
    else
        file_save = [BKT_OI('save_driver'),'OI_SST_5X5_',var_list{varname},'_1982-2014.mat'];
    end
    save(file_save,'clim_final','-v7.3');
end