function [] = preproc_MODIS_NDVI_EVI_Global(pp_type)
    %% NDVI3g
    global data_src;
    global domain_def;
    global chosen_prd;
    global DATA_VI_05rs_out DATA_VI_05rs_sm_out;
    global nds nmd;
    
    [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, 'MODIS', pp_type);

    % period for complete annual data
    yr_s=chosen_prd(1);
    yr_e=chosen_prd(end);
    nyr=yr_e-yr_s+1;
    
    if strcmp(pp_type,'RS_05')
        
        dir_in=sprintf('%s/MODIS/NDVI-EVI/MOD13C2-v006/monthly/regrid/0.5deg/',data_src);
        
        % check if the array DATA_VI_05rs has created
        if size(DATA_VI_05rs_out,2) <= 1
            % Create 
            % Monthly data          nlon         nlat   nweek/nmonth  nyear   ndatasource
            DATA_VI_05rs_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, nyr,  nds);
            % Monthly data smoothed
            DATA_VI_05rs_sm_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, nyr,  nds, nmd);
        end
        
        for y=yr_s:yr_e
            for m=1:12
                % ---- Read data ---
                file_name=sprintf('MOD13C2_v006_NDVI-EVI_%d-%02d-01_0.5deg_remapbil.nc',y,m);
                nc_var = ncgeodataset(sprintf('%s/%s', dir_in, file_name));
                fprintf(sprintf('--> Processing file : %s\n',file_name));

                % ---- Vegetation indices  ---
                % NDVI(time, lat, lon)
                % Datasource ID
                ds_id=1;
                data_tmp=squeeze(double(nc_var.data('NDVI')));
                data_tmp(data_tmp<double(0.0))=nan;
                DATA_VI_05rs_out(:,:,m,(y-yr_s+1),ds_id)=squeeze(data_tmp(y_s_e_g(1):y_s_e_g(2),x_s_e_g(1):x_s_e_g(2)))';
                
                % EVI(time, lat, lon)
                % Datasource ID
                ds_id=2;
                data_tmp=squeeze(double(nc_var.data('EVI')));
                data_tmp(data_tmp<double(0.0))=nan;
                DATA_VI_05rs_out(:,:,m,(y-yr_s+1),ds_id)=squeeze(data_tmp(y_s_e_g(1):y_s_e_g(2),x_s_e_g(1):x_s_e_g(2)))';
                
                clearvars nc_var data_tmp;
            end
        end
        
        
        % Smoothing for a more stable estimate of growing season
        [s1 s2 s3 s4 s5]=size(DATA_VI_05rs_out);
        data_rs=reshape(DATA_VI_05rs_out, [s1 s2 s3*s4 s5]);

        % Smooth the data with different methods
        w=3; % windows size
        ts_sm1=smoothdata(data_rs,3,'sgolay',w, 'Degree',2);
        ts_sm2=smoothdata(data_rs,3,'movmean',w);
        ts_sm3=smoothdata(data_rs,3,'gaussian',w);

        DATA_VI_05rs_sm_out(:,:,:,:,:,1)=reshape(ts_sm1, [s1 s2 s3 s4 s5]);
        DATA_VI_05rs_sm_out(:,:,:,:,:,2)=reshape(ts_sm2, [s1 s2 s3 s4 s5]);
        DATA_VI_05rs_sm_out(:,:,:,:,:,3)=reshape(ts_sm3, [s1 s2 s3 s4 s5]);

end


