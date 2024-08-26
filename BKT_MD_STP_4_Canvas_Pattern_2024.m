% A script for running the bucket model, more description below:
clear;
close all;

% #########################################################################
% Prepare for the enviromental driver 
% #########################################################################
[E.true_SST,E.true_AT,E.e_air,E.u_environment,E.Qs,E.direct_ratio,E.zenith_angle] = ...
                                              BKT_MD_STP_3_PREP_2024;

% ######################################
% Set parameters                      ##
% ######################################
P.deck_time     = 240;         % on deck insolation time
P.diameter      = 0.163;
P.depth         = 0.14;
P.mass_bucket   = 1.7;
P.cp_bucket     = 3900;
P.albedo_bucket = 0.2;
P.solar_shading = .5;          % shading to insolation
P.alpha         = 0;           % Residual air temperature in bucket water
P.s_environment = 4;           % Ship speed
P.exp_id        = 1;           % 2 - more exposure
P.leakage_rate_haul = 0.01;    % unit: m/min
P.leakage_rate_deck = 0.005;
% P.time_step     = 3.5;       % Plot results at which time_step
P.cover_top     = P.solar_shading;

% Initialize the model ----------------------------------------------------
E.init_SST      = E.true_AT * P.alpha + E.true_SST * (1-P.alpha);
% SST_raw         = BKT_canvas_model(E,P);

PP.sensible_on   = 1;
PP.long_on       = 1;
PP.latent_on     = 1;
PP.solar_top_on  = 1;
PP.solar_wall_on = 1;

SST_raw          = BKT_canvas_model(E,P,PP);
SST_bkt         = squeeze(nanmean(SST_raw,3));
SST_true        = squeeze(nanmean(E.true_SST,3));
Bias_bkt        = SST_bkt(:,:,:,end) - SST_true;

% Infill and infer regions without bias estimates ------------------------- 
clear('Bias_infill');
P.threshold = 1;
for ct_mon   = 1:12
    temp     = Bias_bkt(:,:,ct_mon);
    temp_ocn = CDC_interp_high_reso(0.5:1:360,-89.5:1:90,temp,0.25:0.5:360,-89.75:0.5:90,'Ocean','cubic');
    Bias_infill(:,:,ct_mon) = CDC_average_grid(0.25:0.5:360,-89.75:0.5:90,temp_ocn,0.5:1:359.5,-89.5:1:89.5,P);
end

% Find neighboring values -------------------------------------------------
Bucket_bias_pattern = Bias_infill;
for ct = 1:30
    Bucket_bias_pattern = get_value_from_neighbours(Bucket_bias_pattern);
end

% Smoothing ---------------------------------------------------------------
for ct = 1:10
    for ct_mon   = 1:12
        Bucket_bias_pattern(:,:,ct_mon) = CDC_smooth2(Bucket_bias_pattern(:,:,ct_mon),5,1);
    end
end

% Replace -----------------------------------------------------------------
a = Bias_infill;
a(isnan(a)) = Bucket_bias_pattern(isnan(a));
Bucket_bias_pattern = a;

% Mask out inner land -----------------------------------------------------
[mask,topo,coast] = CDF_land_mask(1,1,1,0);
coast = CDF_span_mask(coast');
for ct = 1:6
    coast = CDF_span_mask(coast);
end
l_rm = coast == 0 & isnan(Bias_infill(:,:,1)) & topo' > 0;
Bucket_bias_pattern(repmat(l_rm,1,1,12)) = nan;

% Show a plot -------------------------------------------------------------
ct_mon = 7;
lon = 0.5:1:360;
lat = -89.5:1:90;
subplot(2,1,1); cla;
CDF_pcolor(lon,lat,Bias_bkt(:,:,ct_mon)); caxis([-1 1]*.7); b2rCD(7,0);
CDF_boundaries;
subplot(2,1,2); cla;
CDF_pcolor(lon,lat,Bucket_bias_pattern(:,:,ct_mon)); caxis([-1 1]*.7); b2rCD(7,0);
CDF_boundaries;

% Data are stored according to month rather than season
% This change is implemented, i.e., the line commented, in Nov, 2022
% See "BKT_generate_patterns.m" for the first time of this changes
% That file is currently in "Git_Code/LME_intercompare/SAT_homogenization/backup_20231231/"
% Bucket_bias_pattern = [Bucket_bias_pattern(:,1:18,[7:12 1:6]) Bucket_bias_pattern(:,19:36,:)];

disp(['Get pattern of biases ...'])
dir_save = '/Users/dc1e23/Dropbox/------------ Git_code ----------------/ZZ_Completed/SST_Bucket_Model/';
file_pattern = [dir_save,'spatial_pattern_1x1_bucket_bias_solar_exposure_',num2str(1-P.solar_shading),'.mat'];
save(file_pattern,'Bucket_bias_pattern','-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_value_from_neighbours(in)

    clear('Neighbor')
    for ct1 = 1:size(in,1)
        for ct2 = 1:size(in,2)
            ii = [ct1-1:ct1+1];
            jj = [ct2-1:ct2+1];
            ii(ii<1) = ii(ii<1) + size(in,1);
            ii(ii>size(in,1)) = ii(ii>size(in,1)) - size(in,1);
            jj(jj<1) = [];
            jj(jj>size(in,2)) = [];
            temp = in(ii,jj,:,:);
            if ct2 == 1
                temp(2,1,:,:) = nan;              
            else
                temp(2,2,:,:) = nan;
            end
            temp = nanmean(reshape(temp,size(temp,1)*size(temp,2),size(temp,3),size(temp,4)),1);
            Neighbor(ct1,ct2,:,:) = temp;
        end
    end
    out = in;
    out(isnan(in)) = Neighbor(isnan(in));
end