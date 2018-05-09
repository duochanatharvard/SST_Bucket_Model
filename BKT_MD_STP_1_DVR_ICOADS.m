% This file prepares the hourly climatology of air_temperature, dew point
% temperature, cloud cover, 10m wind speed, and sea level pressure from 
% ICOADS3.0 dataset.

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
