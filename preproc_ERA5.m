function [data_out] = preproc_ERA5(domain_def,RS_D)
    
    vars_in={       't2m', 'ssrd'};
    vars_conv={'- 273.15',     ''};
    % z: specifically referring to geopotential height at 500 hPa, original data is
    % geopotential height, need to /9.8 to convert it to the height
    %             1  2  3  4  5  6  7  8  9 10 11,12
    days_of_mon=[31,28,31,30,31,30,31,31,30,31,30,31];
    
    % Local climate data from ERA5
    global data_src;
    global chosen_prd;
    
    switch RS_D
        case 'RS_05'
            rs_dir='remap/0.5deg';
            rs_fstr='_0.5deg_remapycon';
        case 'RS_025'
            rs_dir='remap/0.25deg';
            rs_fstr='_0.25deg_remapycon';
        otherwise
            error('Invalid resolution!!');
    end
    
    [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, 'ERA5', RS_D);
    
    % period for complete annual data
    yr_s=chosen_prd(1);
    yr_e=chosen_prd(end);
        
    % check if the array dataset has created
    data_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, (yr_e-yr_s)+1, size(vars_in,2));

    for v=1:size(vars_in,2)
        for y=yr_s:yr_e
            % ---- Read data ---
            yi=y-yr_s+1;

            if strcmp(vars_in{v},'z')
                dir_in=sprintf('%s/ECMWF/ERA5/raw/atmosphere/mon/pressure_levels/%s/500/%s/',data_src,rs_dir,vars_in{v});
                fname=sprintf('%s_ECMWF-ERA5_rean_pressure_levels_500_mon_%d-%d%s.nc',vars_in{v},y,y,rs_fstr);
            else
                dir_in=sprintf('%s/ECMWF/ERA5/raw/surface/mon/single_level/%s/%s',data_src,rs_dir,vars_in{v});
                fname=sprintf('%s_ECMWF-ERA5_rean_single_level_mon_%d-%d%s.nc',vars_in{v},y,y,rs_fstr);
            end

            fprintf(sprintf('--> Processing file : %s\n',fname));

            fullname=sprintf('%s/%s',dir_in, fname);
            nc_var = ncgeodataset(fullname);

            % ---- Climate variables  ---
            % vars (month, lat, lon)
            dtmp=squeeze(double(nc_var.data(vars_in{v})));
            [s1 s2 s3]=size(dtmp);
            dtmp_uconv=nan(s2,s3);

            for m=1:12
                % Unit conversion
                dtmp_uconv=squeeze(eval(['dtmp(m,:,:)' vars_conv{v}]));
                
                data_out(:,:, m, yi, v)=squeeze(dtmp_uconv(y_s_e_g(1):y_s_e_g(2),x_s_e_g(1):x_s_e_g(2)))';
            end

            clearvars dtmp dir_in fname;
        end
    end
end