% Linearized wooden bucket model

function [SST_bck_li,F_bck_li] = LR_Wooden_Model(P)

    dir = '/Volumes/Untitled/01_Research/03_DATA/Bucket_Model/';

    switch P.ct_reg,
        case 1,
            region_name = 'Trop';
        case 2,
            region_name = 'Trop_Atl';
        case 3,
            region_name = 'SubT';
        case 4,
            region_name = 'SubT_Atl';
    end

    file_load = [dir,'Linearized_wooden_bucket_model_',region_name,'.mat'];
    load(file_load,'LM')


    % -------------------------------------------------------------------------
    % Find Initial temperature and cooling rate
    % -------------------------------------------------------------------------
    clear('SST_init_li','SST_cool_li')
    F_init_li = (1 - P.r_at) .* LM.SST_true + P.r_at .* LM.AT_true;

    F_cool_li = P.shading .* LM.solar_dep + P.r_at .* LM.at_mix_dep + ...
                   P.bck_siz.^P.alpha .* LM.bck_siz_dep + LM.mean_cool;

    % -------------------------------------------------------------------------
    % Integrate the linearized model
    % -------------------------------------------------------------------------
    clear('SST_li')
    for t = 0:20
        F_li(t+1,:,:,:)  =  F_init_li + t * F_cool_li;
    end


    % -------------------------------------------------------------------------
    % Find and merge with ERI data
    % -------------------------------------------------------------------------
    LM.eri(1,:) = LM.SST_true(1,:) + P.bias_eri;
    F_eri       = repmat(reshape(LM.eri,[1,size(LM.eri)]),size(F_li,1),1,1);
    F_bck_li    = F_li .* (1 - P.r_eri) + F_eri .* P.r_eri;

              
    % -------------------------------------------------------------------------
    % Generate hourly outputs
    % -------------------------------------------------------------------------
    C0_LCL_hr = [1:24];
    omega = 2*pi/24;
    base_x2 = [ones(numel(C0_LCL_hr),1) sin(C0_LCL_hr'*omega) cos(C0_LCL_hr'*omega) ...
                                        sin(C0_LCL_hr'*omega*2) cos(C0_LCL_hr'*omega*2)];
    clear('SST_hr_li')
    for ct_t = 1:size(F_bck_li,1)
        for ct_sea = 1:size(F_bck_li,3)
            SST_bck_li(:,ct_sea,ct_t) = base_x2 * F_bck_li(ct_t,:,ct_sea)';
        end
    end
    
end