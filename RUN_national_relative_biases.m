% A script for running the bucket model, more description below:
% This is to estimate the national relative offsets, and see if the
% difference can be explained by differences in canvas bucket parameters
% More specifically, here we consider two set of canceling parameters:
% 1. on deck time
% 2. bucket size

clear;
% ######################################
% Set parameters                      ##
% ######################################
driver_ERA = 1;
s_environment = 6;
deck_time  = 900;              % here it is directly set in the model input
solar_shading = 0.6;           % here it is directly set in the model input
cx = 1:72;                     % This is over the 60S-60N
cy = 7:30;

% ######################################
% Prepare for the enviromental driver ##
% ######################################
[true_SST,true_AT,e_air,u_environment,Qs,direct_ratio,zenith_angle] = BKT_MD_STP_3_PREP(driver_ERA);

% ##########################################################################
% Run the bucket models for 300s on deck using UK buckets as a reference  ##
% ##########################################################################
SST_ref  = BKT_MD_STP_2_MD_CANVAS_GRD_SIZ(true_SST(cx,cy,:,:),true_AT(cx,cy,:,:),...
    e_air(cx,cy,:,:),u_environment(cx,cy,:,:),s_environment,Qs(cx,cy,:,:),...
    direct_ratio(cx,cy,:,:),300,solar_shading,zenith_angle(cx,cy,:,:),0.163,0.14);
SST_ref = squeeze(nanmean(nanmean(SST_ref,1),2));
SST_ref = squeeze(nanmean(nanmean(SST_ref,1),2));

% ######################################################################
% Run the bucket models for 900s on deck using different bucket size  ##
% ######################################################################
ct = 0;
clear('SST_save')
for bucket_size = [0.5:0.1:2]

    ct = ct + 1;

    % ***************************
    % This is to run the model **
    % ***************************

    SST_exp  = BKT_MD_STP_2_MD_CANVAS_GRD_SIZ(true_SST(cx,cy,:,:),true_AT(cx,cy,:,:),...
        e_air(cx,cy,:,:),u_environment(cx,cy,:,:),s_environment,Qs(cx,cy,:,:),...
        direct_ratio(cx,cy,:,:),900,solar_shading,zenith_angle(cx,cy,:,:),0.163*bucket_size,0.14*bucket_size);
    SST_exp = squeeze(nanmean(nanmean(SST_exp,1),2));
    SST_exp = squeeze(nanmean(nanmean(SST_exp,1),2));

    SST_save(ct,:) = SST_exp;
end

dir_save = BKT_OI('save_output');
save([dir_save,'nation_relative_biases.mat'],'SST_save','SST_ref');
% ##########################################################
% Analysis model results and generate figures             ##
% ##########################################################
figure(1); clf; hold on;

pic = SST_save - SST_ref(end);

if 0,   % This is the colorful version
    contourf(0:1:15,0.5:0.1:2,pic,[-2:0.2:2],'w','linewi',2)
    caxis([-1 1]*2)
    color_shuang([.08 .03;.48 .6],[0.9 .2],[.7 .9],20,0);
    h = colorbar;
    ylabel(h,'Bias relative to 5 minute exposure of the reference bucket (^oC)')
else
    [C,h] = contour(0:1:15,0.5:0.1:2,pic,[-2:0.2:2],'k','linewi',1,'ShowText','on');
    clabel(C,h,'FontSize',15,'Color','k')
end

CDF_panel([0 15 0.5 2],'',{},'Exposure time (minutes)','Relative bucket size');
clear('h')
col = jet(7)*.92;
% nation_list = [0.0142 0.0165 0.0424 -0.0631 0.1054 -0.1051 -0.1348];
% nation_list = [0.01 0.02 0.03 -0.06 0.10 -0.10 -0.12];
nation_list = [0.0120   0.0234    0.0384   -0.0400    0.0860   -0.0954   -0.1224];
for i = 1:7
    [~,h(i)] = contour(0:1:15,0.5:0.1:2,pic,[1 1]*nation_list(i),'color',col(i,:),'linewi',3);
end

set(gca,'xtick',[0:2:14],'ytick',[0.5:0.25:2],'fontsize',12)
legend(h,{'DE','UK','US','JP','RU','NL','Deck 156'},'fontsize',12,...
    'fontweight','bold','location','northwest')
plot(5,1,'rp','markersize',15,'markerfacecolor','r')

CDF_save(1,'png',500,['/Users/zen/Dropbox/Research/SST_bucket_regional_bias/',...
    'SST_20180420_Method/Figure6/20180525_Bucket_model_case_2.png']);
