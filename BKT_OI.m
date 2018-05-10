function dir = BKT_OI(input)

    if strcmp(input,'load_ERA'),
        dir = '/Volumes/My Passport Pro/ICOADS_3.0/';

    elseif strcmp(input,'save_driver')
        dir = '/Users/zen/Research/Git_Code/Bucket_Model/Driver/';

    elseif strcmp(input,'load_OI_SST')
        dir = '/Volumes/My Passport Pro/ICOADS_RE/Step_Miscellaneous_ICOADS_RE/';

    elseif strcmp(input,'load_ICOADS')
        dir = '/Users/zen/Research/Git_Code/Bucket_Model/ICOADS/';

    elseif strcmp(input,'save_output')
        dir = '/Users/zen/Research/Git_Code/Bucket_Model/';

    end
