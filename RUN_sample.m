% A sample script for running bucket models.
% It shows how to call: 
% >> "BKT_MD_STP_2_MD_CANVAS_GRD_SIZ.m" the canvas bucket model
% >> "BKT_MD_STP_2_MD_WOODEN_GRD_SIZ.m" the wooden bucket model
% 
% "GRD" means the model runs on individual grid boxes using vectorized
% computation in Matlab
% "SIZ" orginally meant we can tune bucket size (in an 2018 version), 
% but in the updated 2019 version, the script can take in a structure "P",
% which allows for playing with more model parameters.
%
% To check which parameters to play with, please type:
% >> help help BKT_MD_STP_2_MD_CANVAS_GRD_SIZ
% >> help BKT_MD_STP_2_MD_WOODEN_GRD_SIZ
% 
% Acknowlegement to Bowen Zhao for reminding me to update this sample
% script.  
% 
% [Note that] 
% "RUN_Japanese_bucket_cooling.m" and "RUN_national_relative_biases.m"
% are no longer compatible with the latest version of the bucket model.
% 
% References:
% [1] Folland, C. K., & Parker, D. E. (1995). Correction of instrumental 
% biases in historical sea surface temperature data. Quarterly Journal of 
% the Royal Meteorological Society, 121(522), 319-367.
% 
% [2] Chan, D., & Huybers, P. (2020). Systematic differences in bucket sea
% surface temperatures caused by misclassification of engine room intake
% measurements. Journal of Climate, 33(18), 7735-7753.
% 
% 
% Duo Chan
% Mar. 8, 2022
% WHOI


clear;

% #########################################################################
% Prepare for the enviromental driver 
% #########################################################################
driver_ERA = 5;                % The climatology from NOC [Recommended]
P.average_forcing  = 0;        % Run the model on 5x5 grids
                               % When equals to 1, the model will run using 
                               % zonal averaged environmental variables
[true_SST,true_AT,e_air,u_environment,Qs,direct_ratio,zenith_angle] = ...
                                    BKT_MD_STP_3_PREP_2019(driver_ERA,0,P);
init_SST = true_SST; % Initial SSTs are the same as unbiased actual SSTs


% #########################################################################
% Set the domain for the model [TODO]
% #########################################################################
do_regional = 0;                             % PLAY WITH THIS PARAMETER !!!
if do_regional == 1   % try a small patch
    cx = 68:69;
    cy = 23:24;
else                  % run globally at 5x5 resolution
    cx = 1:72;
    cy = 1:36;
end

% #########################################################################
% Call a canvas bucket model
% #########################################################################
clear('P')
P.deck_time     = 540;          % [s] - total intergration time 600s
P.s_environment = 7;            % [m/s]
P.solar_shading = .5;           % [fractional]
P.diameter      = 0.163;        % [m]
P.depth         = 0.14;         % [m]
P.mass_bucket   = 1.7;          % [kg]
P.exp_id        = 1;            % 1 - less exposure; 2 - more exposure
P.cover_top     = 1;            % whether to cover the top of the bucket

clear('SST_c')
SST_c  = BKT_MD_STP_2_MD_CANVAS_GRD_SIZ(...
    init_SST(cx,cy,:,:),true_AT(cx,cy,:,:),e_air(cx,cy,:,:),...
    u_environment(cx,cy,:,:),Qs(cx,cy,:,:),direct_ratio(cx,cy,:,:),...
    zenith_angle(cx,cy,:,:),P);


% #########################################################################
% Call a wooden bucket model
% #########################################################################
clear('P')
P.deck_time     = 540;                  % [s]
P.s_environment = 7;                    % [m/s]
P.solar_shading = 0.5;                  % [fractional] 
P.thickness     = 0.01;                 % [m]      

if do_regional == 1
    SST_w  = BKT_MD_STP_2_MD_WOODEN_GRD_SIZ(...
        init_SST(cx,cy,:,:),true_AT(cx,cy,:,:),e_air(cx,cy,:,:),...
        u_environment(cx,cy,:,:),Qs(cx,cy,:,:),direct_ratio(cx,cy,:,:),...
        zenith_angle(cx,cy,:,:),P);
end

% #########################################################################
% Display results
% #########################################################################
if do_regional == 1  
    
    SST_c = squeeze(nanmean(nanmean(SST_c,1),2));
    SST_w = squeeze(nanmean(nanmean(SST_w,1),2));
    
    figure(1); clf; 
    subplot(1,2,1); hold on;
    t = [1:size(SST_c,3)]*.5-.5;
    pic  = squeeze(nanmean(SST_c(:,[1],:),1));
    h(1) = plot(t,pic - pic(1),'b-');
    pic  = squeeze(nanmean(SST_c(:,[7],:),1));
    h(2) = plot(t,pic - pic(1),'r-');
    t = [1:size(SST_w,3)]*.5-.5;
    pic  = squeeze(nanmean(SST_w(:,[1],:),1));
    h(3) = plot(t,pic - pic(1),'b--');
    pic  = squeeze(nanmean(SST_w(:,[7],:),1));
    h(4) = plot(t,pic - pic(1),'r--');
    legend(h,{'canvas Jan.','canvas Jul.','wooden Jan.','wooden Jul.'},'location','southwest');
    xlabel('Time')
    ylabel('Bias [^oC]')
    
    subplot(1,2,2); hold on;
    Bias_c = SST_c - SST_c(:,:,1);
    contourf(Bias_c(:,:,end),'linest','none');
    colorbar;
    colormap(gca,b2r);
    caxis([-1 1]*.8);
    xlabel('Month')
    ylabel('Local Hour')
    title('Biases in a canvas bucket after 10 minutes','fontweight','normal')
 
else
    
    Bias_c = SST_c - SST_c(:,:,:,:,1);
    Bias_c = squeeze(nanmean(Bias_c(:,:,:,:,end),3));
    
    figure(2); clf;
    subplot(1,2,1); hold on;
    contourf(2.5:5:360,-87.5:5:90,Bias_c(:,:,1)','linest','none');
    colorbar;
    colormap(gca,b2r);
    caxis([-1 1]*1.2);
    title('Jan. canvas biases after 10 minutes','fontweight','normal')

    subplot(1,2,2); hold on;
    contourf(2.5:5:360,-87.5:5:90,Bias_c(:,:,7)','linest','none');
    colorbar;
    colormap(gca,b2r);
    caxis([-1 1]*1.2);
    title('Jul. canvas biases after 10 minutes','fontweight','normal')
    
end


function col = b2r
    col = [0.1774    0.0624    0.6376
    0.1456    0.0700    0.6945
    0.1073    0.0791    0.7499
    0.0902    0.1165    0.8032
    0.1038    0.1907    0.8541
    0.1368    0.2826    0.8856
    0.2022    0.3890    0.8846
    0.2664    0.4845    0.8849
    0.3287    0.5697    0.8871
    0.3888    0.6454    0.8914
    0.4466    0.7128    0.8982
    0.5019    0.7729    0.9073
    0.5553    0.8267    0.9184
    0.6070    0.8749    0.9312
    0.6578    0.9176    0.9448
    0.7085    0.9547    0.9586
    0.7600    0.9716    0.9582
    0.8131    0.9829    0.9588
    0.8688    0.9918    0.9646
    0.9274    0.9976    0.9765
    0.9845    0.9991    0.9259
    0.9809    0.9971    0.8634
    0.9862    0.9943    0.8018
    0.9908    0.9816    0.7408
    0.9867    0.9512    0.6804
    0.9821    0.9117    0.6206
    0.9769    0.8631    0.5613
    0.9710    0.8059    0.5027
    0.9643    0.7401    0.4450
    0.9564    0.6664    0.3883
    0.9473    0.5852    0.3330
    0.9365    0.4972    0.2792
    0.9239    0.4034    0.2274
    0.9092    0.3047    0.1777
    0.8921    0.2024    0.1302
    0.8409    0.1284    0.1170
    0.7674    0.1260    0.1665
    0.6994    0.1296    0.2105
    0.6377    0.1268    0.2397
    0.5811    0.1189    0.2576];
end
