function dir = BKT_OI(input,env)

    if strcmp(input,'load_ERA'),
        dir = '/Volumes/My Passport Pro/ICOADS_3.0/';

    elseif strcmp(input,'save_driver')
        if env == 1,
            dir = '/n/home10/dchan/script/Peter/Hvd_SST/SST_Bucket_Model/Driver/';
        else
            dir = '/Users/zen/Research/Git_Code/SST_Bucket_Model/Driver/';
        end
    elseif strcmp(input,'load_OI_SST')
        dir = '/Volumes/My Passport Pro/ICOADS_RE/Step_Miscellaneous_ICOADS_RE/';

    elseif strcmp(input,'load_ICOADS')
        dir = '/Users/zen/Research/Git_Code/SST_Bucket_Model/ICOADS/';

    elseif strcmp(input,'save_output')
        if env == 1,
            dir = '/n/home10/dchan/script/Peter/Hvd_SST/SST_Bucket_Model/Outputs/';
        else
            dir = '/Users/zen/Research/Git_Code/SST_Bucket_Model/Outputs/';
        end

    elseif strcmp(input,'load_NOCS')
        dir = '/Users/zen/Research/Git_Code/SST_Bucket_Model/NOCSV2/';
    end
