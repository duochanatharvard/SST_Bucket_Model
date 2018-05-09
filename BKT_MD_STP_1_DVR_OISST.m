% This file prepares the monthly 5X5 climatology from OI-SST (^oC)

% We take the input of 1982-2014 OI 11-day smoothed daily climatology 
% at 0.25X0.25 resolution, which is OI_clim_1982_2014.mat
%
% data: /Volumes/My Passport Pro/ICOADS_RE/Step_Miscellaneous_ICOADS_RE/
% code: generating this climatology is saved under the model folder
% but is not yet cleaned ...

clear;

% ************
% read data **
% ************
dir_load = BKT_OI('load_OI_SST');
load([dir_load,'OI_clim_1982_2014.mat']);

% ************************
% generate monthly mean **
% ************************
mon_list = [0 31 28 31 30 31 30 31 31 30 31 30 31];
mon_list = cumsum(mon_list);

clear('OI_temp','OI_month')
for mon = 1:12
    OI_temp(:,:,mon) = nanmean(OI_clim(:,:,mon_list(mon)+1:mon_list(mon+1)),3);
end
OI_temp(abs(OI_temp)>1000) = NaN;

% ****************************************
% regrid monthly mean to 5X5 resolution **
% ****************************************
lon_oi = 0.125:0.25:360;
lat_oi = -89.875:0.25:90;
clear('OI_month')
for i = 1:72
    for j = 1:36
        
        clear('temp','logic_lon','logic_lat','temp_lat','temp_weight')
        logic_lon = lon_oi >= (i-1)*5 & lon_oi < i*5;
        logic_lat = lat_oi >= (j-1)*5 - 90 & lat_oi < j*5 - 90;
        temp = OI_temp(logic_lon,logic_lat,:);
        temp_lat = repmat(lat_oi(logic_lat),20,1);
        temp_weigh = cos(temp_lat.*pi./180);
        temp_weigh(isnan(temp(:,:,1))) = nan;
        
        for m = 1:12
            tem = temp(:,:,m);
            OI_month(i,j,m) = nansum(tem(:).*temp_weigh(:)) ./ nansum(temp_weigh(:));
        end
    end
end

% ************
% save data **
% ************
dir_save = BKT_OI('save_driver');
save([dir_save,'OI_monthly_1982_2014.mat'],'OI_month','-v7.3')