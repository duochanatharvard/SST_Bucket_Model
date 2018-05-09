% A script for running the bucket model, more description below:
% This is to estimate the Japanese bucket cooling, which is used in 
% Figure XX in Paper: 

clear;
% ######################################
% Set parameters                      ##
% ######################################
driver_ERA = 1;
s_environment = 6;
deck_time  = 900;              % here it is directly set in the model input
solar_shading = 1;             % here it is directly set in the model input
cx = 25:52;                    % This is over the North Pacific 
cy = 23:28;

% ######################################
% Prepare for the enviromental driver ##
% ######################################
[true_SST,true_AT,e_air,u_environment,Qs,direct_ratio,zenith_angle] = BKT_MD_STP_3_PREP(driver_ERA);

% ######################################
% Run the bucket models               ##
% ######################################
ct = 0;
for solar_shading = 0:0.1:1
    
    ct = ct + 1;

    % ***************************
    % This is to run the model **
    % ***************************
    
    SST_out_UK  = BKT_MD_STP_2_MD_CANVAS_GRD_SIZ(true_SST(cx,cy,:,:),true_AT(cx,cy,:,:),e_air(cx,cy,:,:),...
        u_environment(cx,cy,:,:),s_environment,Qs(cx,cy,:,:),direct_ratio(cx,cy,:,:),300,solar_shading,zenith_angle(cx,cy,:,:),0.163,0.14);

    SST_out_JP  = BKT_MD_STP_2_MD_CANVAS_GRD_SIZ(true_SST(cx,cy,:,:),true_AT(cx,cy,:,:),e_air(cx,cy,:,:),...
        u_environment(cx,cy,:,:),s_environment,Qs(cx,cy,:,:),direct_ratio(cx,cy,:,:),900,solar_shading,zenith_angle(cx,cy,:,:),0.2,0.4);

    % **************************************
    % This is to look at the mean cooling **
    % **************************************
    a_JP = SST_out_JP - repmat(SST_out_JP(:,:,:,:,1),1,1,1,1,size(SST_out_JP,5));
    a_JP = squeeze(nanmean(nanmean(a_JP,4),3));
    d_JP = squeeze(nanmean(a_JP,1));

    a_UK = SST_out_UK - repmat(SST_out_UK(:,:,:,:,1),1,1,1,1,size(SST_out_UK,5));
    a_UK = squeeze(nanmean(nanmean(a_UK,4),3));
    d_UK = squeeze(nanmean(a_UK,1));
    
    cooling(:,ct) = nanmean(d_JP,1) - nanmean(d_UK(:,6),1); 
    
    % ****************************************************
    % This is to look at the increase in seasonal cycle **
    % ****************************************************
    a_JP = SST_out_JP - repmat(SST_out_JP(:,:,:,:,1),1,1,1,1,size(SST_out_JP,5));
    b_JP = squeeze(nanmean(nanmean(a_JP(:,:,:,[6 7 8],:),4),3)) - ...
        squeeze(nanmean(nanmean(a_JP(:,:,:,[1 2 12],:),4),3));
    c_JP = squeeze(nanmean(b_JP,1));

    a_UK = SST_out_UK - repmat(SST_out_UK(:,:,:,:,1),1,1,1,1,size(SST_out_UK,5));
    b_UK = squeeze(nanmean(nanmean(a_UK(:,:,:,[6 7 8],:),4),3)) - ...
        squeeze(nanmean(nanmean(a_UK(:,:,:,[1 2 12],:),4),3));
    c_UK = squeeze(nanmean(b_UK,1));

    seasonal(:,ct) =nanmean(c_JP,1) - nanmean(c_UK(:,6),1);
    
    % ***************************************************
    % This is to look at the increase in diurnal cycle **
    % ***************************************************
    da = max(SST_out_JP,[],3) - min(SST_out_JP,[],3);
    a = da - repmat(da(:,:,:,:,1),1,1,1,1,size(da,5));
    a = squeeze(nanmean(a,4));
    a = squeeze(nanmean(a,1));
    a = squeeze(nanmean(a,1));
    diurnal(:,ct) = a;
end

% ##########################################################
% Analysis model results and generate figures             ##
% ##########################################################
driver_ERA = 1;
minutes = 15;

figure(2); 
clf; hold on;

col = color_shuang([.12 0;.6 .7],[.7 .3],[1 1],6,0);
plot([4 4]+1,[-1 1],'k-','linewi',2);
plot([0 minutes],[0 0],'k--','linewi',2);

for i = 1:11
    hold on;
    plot(0:minutes,cooling(:,i),'linewi',2,'color',col(12-i,:));
end

colormap(col(1:11,:));
h = colorbar;
caxis([-0.05 1.05])
ylabel(h,'Insolation ratio')
CDF_panel([0 minutes -0.5 0.5],'',{},'Exposure time (minutes)','Relative bias (^oC)');

% CDF_save(2,'png',500,'/Users/zen/Dropbox/Research/SST_Uniform_ETCW_Nature/Figures/20180124_FigS7.png');