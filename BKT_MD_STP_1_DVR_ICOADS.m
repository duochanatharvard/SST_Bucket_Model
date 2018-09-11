% This file prepares the hourly climatology of air_temperature, dew point
% temperature, cloud cover, 10m wind speed, and sea level pressure from 
% ICOADS3.0 dataset, using only records whos SST measurement are valid and 
% are taken by buckets. 

% However, we did not compute these by ourselves here in this project, we
% use input from the preprocessing and diurnal project, that generate uses
% the following codes to generate the following results:

% preprocessing project
%       code: sst_function_step_0_pretreatment_step_5_MCLS_ICOADS_re.m
%       data: 5x5 monthly value of variables, rough version
%       source: SHIP_GRID_C0_XXXX.mat
%       memo: this function uses ICOADS QC flag <= 3, which constrains the 
%             data greatly to a designed climatology

% diurnal project 
%       code: sst_function_step_1_diurnal_step_6_sum_and_grid_ICOADS_re.m
%       data: gridded diurnal cycle of bucket and buoy
%       source: DA_XXXX_Gridded_BCK_do_more_samples_0.mat

% These files are saved in folder ICOADS

function BKT_MD_STP_1_DVR_ICOADS(num)

    dir_home = '/Volumes/My Passport Pro/ICOADS_RE/';

    % Generate 5 degree grids for each variable
    var_list = {'C0_SST','C0_OI_CLIM',''        ,'C0_AT'  ,'C1_AF' ,'C1_ANC',...
                'C0_WS' ,''          ,'C1_WNC'  ,'C0_DPT' ,'C1_RF' ,'C1_DNC',...
                'C0_CA' ,''          ,'C1_CNC'  ,'C0_SLP' ,'C1_PF' ,'C1_PNC'};

    varname = var_list(num*3-2:num*3);

    for i = 3:-1:1  if isempty(varname{i}), varname(i)=[];end; end

    varname
    
    sst_function_step_005_grid_ICOADS_re(dir_home,varname);
        
end

% -------------------------------------------------------------------------
function [DATA,DATA_ST,DATA_NUM] = sst_function_step_005_grid_ICOADS_re(dir_home,varname)

    %%
    reso_x = 5;
    reso_y = 5;

    clear('DATA','DATA_ST','DATA_NUM')
    DATA = nan(360/reso_x,180/reso_y,24,215,12);
    DATA_ST = nan(360/reso_x,180/reso_y,24,215,12);
    DATA_NUM = nan(360/reso_x,180/reso_y,24,215,12);
    
    %%
    for yr = 1950:2014
        if rem(yr,1) == 0, disp(['Starting Year: ',num2str(yr)]); end
        
        %%
        for mon = 1:12
            
            %
            clear('cmon');    cmon = '00';    cmon(end-size(num2str(mon),2)+1:end) = num2str(mon);
            file_load = [dir_home,'Step_0_pre_step_3_QCed_data_ICOADS_RE/IMMA1_R3.0.0_',num2str(yr),'-',cmon,'_QCed.mat'];
            
            fid = fopen(file_load);
            
            %
            if fid > 0,
                fclose(fid);
                
                %
                load(file_load,'C0_LON','C0_LAT','QC_FINAL','QC_NON_SST','C0_SI_4','C0_LCL');
                for i = 1:numel(varname)
                    load(file_load,varname{i});
                end
                
                % QC control for other meteorological data ----------------
                clear('logic')
                if strcmp(varname{1},'C0_SST')
                    logic = QC_FINAL;
                else
                    if size(varname,2) == 2,
                        eval(['logic = ',varname{2},' <= 3;']);
                    else
                        eval(['logic = ',varname{2},' <= 3 & ',varname{3},' <= 3;']);
                    end
                    logic = logic & QC_FINAL;
                end
                

                % Pickout ship measurement that have valid SST measurement 
                clear('logic_ship')
                logic_ship = logic & (C0_SI_4 >= 0 & C0_SI_4 <= 0.05);
                
                % Loop over local hours and grid the data -----------------
                clear('WM','ST','NUM')
                WM = nan(360/reso_x,180/reso_y,24);
                ST = nan(360/reso_x,180/reso_y,24);
                NUM = nan(360/reso_x,180/reso_y,24);
                for hr = 1:24
                    
                    % Pickout ship data to be gridded ---------------------
                    clear('logic','in_var','in_lon','in_lat','var_grd')
                    logic = logic_ship & abs(C0_LCL - hr) < 0.1;
                    if strcmp(varname{1},'C0_SST')
                        in_var = C0_SST(logic) - C0_OI_CLIM(logic);
                    else
                        eval(['in_var = ',varname{1},'(logic);']);
                    end
                    in_lon = C0_LON(logic);
                    in_lat = C0_LAT(logic);
                    % Gridding --------------------------------------------
                    [var_grd,~] = HM_function_pnt2grd_3d (in_lon,in_lat,[],in_var,[],reso_x,reso_y,[],2);
                    [WM(:,:,hr),ST(:,:,hr),NUM(:,:,hr)] = HM_function_grd_average(var_grd,[],[]);
                end
                DATA(:,:,:,yr-1949,mon) = WM;
                DATA_ST(:,:,:,yr-1949,mon) = ST;
                DATA_NUM(:,:,:,yr-1949,mon) = NUM;
                
            end
        end
    end
    dir_save = '/Users/zen/Research/Git_Code/SST_Bucket_Model/ICOADS/';
    save([dir_save,'BCK_GRID_',varname{1},'_1950_2014.mat'],'DATA','DATA_ST','DATA_NUM');
end