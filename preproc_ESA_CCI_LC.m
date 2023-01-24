function [] = preproc_ESA_CCI_LC()
    % ---- Landcover ----
    global data_src;
    global domain_def;
    global x_s_e_g;
    global y_s_e_g;
    global dmn_lon_n_g;
    global dmn_lat_n_g;
    global DATA_ESA_CCI_LC_out;

    % Data source I: Europe domain (get: x_s_e, y_s_e, ind_dmn_mask)
    [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, 'ESA_CCI_LC', 'RS_05');

    % period for complete annual data
    yr_s=2015;
    yr_e=2015;

    % check if the array dataset has created
    if size(DATA_ESA_CCI_LC_out,2) <= 1
        % Create 
        %                   nlon         nlat          nyear     nigbp land class 
        DATA_ESA_CCI_LC_out = nan(dmn_lon_n_g, dmn_lat_n_g, (yr_e-yr_s)+1, 1);
    end
    
    var_name={
        'class_area_20'  % LCCS ID 20: cropland_irrigated, http://maps.elie.ucl.ac.be/CCI/viewer/download.php
    };

    for y=yr_s:yr_e

        % ---- Read data ---
        dir_in=sprintf('%s/ESA_CCI/LC/Aggregated/0.5deg',data_src);
        fname='ESACCI-LC-L4-LCCS-Map-300m-P1Y-aggregated-0.500000Deg-2015-v2.0.7.nc';
        
        fprintf(sprintf('--> Processing file : %s\n',fname));

        fullname=sprintf('%s/%s',dir_in, fname);
        nc_var = ncgeodataset(fullname);

        % ---- ESA_CCI land type  ---
        for v=1:size(var_name,1)
            dtmp=double(squeeze(nc_var.data(var_name{v}))');
            DATA_ESA_CCI_LC_out(:,:,(y-yr_s)+1,v)=fliplr(dtmp(x_s_e_g(1):x_s_e_g(2),y_s_e_g(1):y_s_e_g(2)));
        end
    end
end
