function [] = preproc_GLEAM(chosen_prd)
    % ---- Landcover ----
    global data_src;
    global domain_def;
    global x_s_e_g;
    global y_s_e_g;
    global dmn_lon_n_g;
    global dmn_lat_n_g;
    global DATA_Evap_05rs_out;
    global DATA_SMRoot_out;

    % Data source I: Europe domain (get: x_s_e, y_s_e, ind_dmn_mask)
    [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, 'GLEAM_E', 'RS_05');

    % Available data period
    data_sy=1980; % data start year (from Jan.)
    
    % period for complete annual data
    yr_s=chosen_prd(1);
    yr_e=chosen_prd(end);
    nyr=yr_e-yr_s+1;

    % check if the array dataset has created
    if size(DATA_Evap_05rs_out,2) <= 1
        % Create 
        %                        nlon         nlat          nmonth, nyear     # of Evap var.
        DATA_Evap_05rs_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, ((yr_e-yr_s)+1), 1);
    end
    if size(DATA_SMRoot_out,2) <= 1
        % Create 
        %                        nlon         nlat          nmonth, nyear     # of soil moisture var.
        DATA_SMRoot_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, ((yr_e-yr_s)+1)*12, 1);
    end
    
    % ET data
    % ---- Read data ---
    dir_in=sprintf('%s/GLEAM/3.5a/monthly/0.5deg',data_src);
    fname='E_1980-2020_GLEAM_v3.5a_MO_0.5deg_remapcon.nc';

    fprintf(sprintf('--> Processing file : %s\n',fname));

    fullname=sprintf('%s/%s',dir_in, fname);
    nc_var = ncgeodataset(fullname);

    % Var. dimension in netcdf file
    % (time (month x year), lon, lat)
    dtmp=double(squeeze(nc_var.data('E')));
    dtmp=permute(dtmp(((yr_s-data_sy)*12+1):((yr_e-data_sy+1)*12), x_s_e_g(1):x_s_e_g(2), y_s_e_g(1):y_s_e_g(2)),[2 3 1]);
    DATA_Evap_05rs_out=reshape(dtmp, [dmn_lon_n_g dmn_lat_n_g 12 ((yr_e-yr_s)+1)]);
    
    % SM data
    % ---- Read data ---
    dir_in=sprintf('%s/GLEAM/3.5a/monthly/0.5deg',data_src);
    fname='SMroot_1980-2020_GLEAM_v3.5a_MO_0.5deg_remapcon.nc';

    fprintf(sprintf('--> Processing file : %s\n',fname));

    fullname=sprintf('%s/%s',dir_in, fname);
    nc_var = ncgeodataset(fullname);

    % Var. dimension in netcdf file
    % (time (month x year), lon, lat)
    dtmp=double(squeeze(nc_var.data('SMroot')));
    dtmp=permute(dtmp(((yr_s-data_sy)*12+1):((yr_e-data_sy+1)*12), x_s_e_g(1):x_s_e_g(2), y_s_e_g(1):y_s_e_g(2)),[2 3 1]);
    DATA_SMRoot_out=reshape(dtmp, [dmn_lon_n_g dmn_lat_n_g 12 ((yr_e-yr_s)+1)]);
end
