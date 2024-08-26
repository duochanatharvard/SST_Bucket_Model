% This file prepares the hourly climatology of insolation (ssrd) 
% from 34yr ERA5 hourly data (1973--2006).

% These outputs are further processed by "BKT_MD_STP_1_DVR_SUM.m"
% in order to get a 5x5 value with data rearranged by hours.

% #############################################################
%  Compute the climatology of t2m, d2m, and insolation       ##
% #############################################################
clear;
dir_load  = '/n/holyscratch01/huybers_lab/dchan/ERA5/';
varname   = 'ssrd';

% read data ---------------------------------------------------------------
clim    = zeros(360,180,8760);
N       = 0;
for yr = 1973:2006
    disp(num2str(yr,'year: %6.0f'))
    fname = [dir_load,'ERA5_hourly_solar_',num2str(yr),'.nc'];
    temp  = ncread(fname,varname);
    temp  = temp(:,end:-1:1,:);
    temp  = (temp(:,1:end-1,:) + temp(:,2:end,:))/2;
    temp  = (temp + temp([2:end 1],:,:))/2;
    if size(temp,3) > 8760, temp = temp(:,:,[1:1416 1441:end]);   end
    clim  = (clim * N + temp) ./ (N+1);
    N     = N + 1;
end

clim      = clim / 3600;                      % convert to unit of w/m2
clim      = reshape(clim,360,180,24,365);
 
% compute the monthly mean ------------------------------------------------
disp('compute monthly climatology');
clim_mon = nan(360,180,24,12);
mon_list = [31 28 31 30 31 30 31 31 30 31 30 31];
mon_cum  = cumsum([0 mon_list]);
for mon = 1:12
    clim_mon(:,:,:,mon) = nanmean(clim(:,:,:,mon_cum(mon)+1:mon_cum(mon+1)),4);
end

% ERA-interim 5X5 hourly rearranged fields --------------------------------
disp('sort by local time')
clim_final = nan(size(clim_mon));
lon = 0.5:1:359.5;
for ct = 1:size(clim_mon,1)

    temp = permute(squeeze(clim_mon(ct,:,:,:)),[2 1 3]);
    lcl  = rem((-0.5:1:23) + lon(ct)/15,24);
    % Starting from -0.5 is important because the first time is 0, 
    % which should be the integral over 23-00h, 
    % so it actually centers on 23:30 rather than 0:00

    % fit the local time
    a = repmat(lcl,1,5);
    b = find(a(2:end)-a(1:end-1) < -20);
    for bb = 1:numel(b)
        a(b(bb)+1:end) = a(b(bb)+1:end)+24;
    end
    a = a - 48;

    temp2 = interp1(a,repmat(temp,5,1,1),(1:24),'spline');
    temp2 = permute(temp2,[2 1 3]);

    clim_final(ct,:,:,:) = temp2;
end

clim_final(clim_final<0) = 0;

% save data ---------------------------------------------------------------
clim_final = single(clim_final);
dir_save = '/n/holyscratch01/huybers_lab/dchan/ERA5/';
file_save = [dir_save,'ERI5_1X1_',varname,'_1973_2006.mat'];
save(file_save,'clim_final','-v7.3');
