function dir = BKT_OI(input)

    if strcmp(input,'home')
        dir = ['/Users/dc1e23/Dropbox/',...
                    '------------ Git_code ----------------/',...
                    'ZZ_Completed/SST_Bucket_Model/'];

    elseif strcmp(input,'load_ERA')
        dir = '/Volumes/My Passport Pro/ICOADS_3.0/';

    elseif strcmp(input,'save_driver')
        dir = [BKT_OI('home'),'Driver/'];

    elseif strcmp(input,'load_OI_SST')
        dir = '/Volumes/My Passport Pro/ICOADS_RE/Step_Miscellaneous_ICOADS_RE/';

    elseif strcmp(input,'load_ICOADS')
        dir = [BKT_OI('home'), 'ICOADS/'];

    elseif strcmp(input,'save_output')
        dir = [BKT_OI('home'),'Outputs/'];

    elseif strcmp(input,'load_NOCS')
        dir = [BKT_OI('home'),'NOCSV2/'];
    end
end