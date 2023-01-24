function [lons_out lats_out data_out] = preproc_SPEI_Global(pp_type, prd,speitype)
    % pp_type: RS_ORG, RS_05
    %% SPEI
    global data_src;
    global domain_def;
    
    [dmn_lon_n_g dmn_lat_n_g x_s_e_g y_s_e_g] = get_data_global_parm(domain_def, 'SPEI', 'RS_05');
    
    % period for complete annual data
    data_sy=1901; % data start year (from Jan.)
    data_ey=2020; % data start year (from Jan.)
    yr_s=prd(1);
    yr_e=prd(end);
    nyr=yr_e-yr_s+1;
    
    if strcmp(pp_type,'RS_05')
        
        % Create 
        % Monthly data          nlon         nlat   nweek/nmonth  nyear
        data_out = nan(dmn_lon_n_g, dmn_lat_n_g, 12, (yr_e-yr_s)+1);

        fprintf(('Arrays created! \n'));
        
        % ---- Read data ---
        dir_in=sprintf('%s/SPEI/',data_src);
        fname=dir(sprintf('%s/%s_1901-2020.nc',dir_in, speitype));

        if size(fname,1)==1

            fprintf(sprintf('--> Processing file : %s\n',fname.name));

            fullname=sprintf('%s/%s',dir_in, fname.name);
            nc_var = ncgeodataset(fullname);

            % ---- LON and LAT ---
            lons_tmp=nc_var.data('lon')';
            lats_tmp=nc_var.data('lat')';

            lon_n = length(lons_tmp);
            lat_n = length(lats_tmp);
            lons=repmat(lons_tmp', [1, lat_n]);
            lats=repmat(lats_tmp,  [lon_n, 1]);
            
            lons_out=lons(x_s_e_g(1):x_s_e_g(2), y_s_e_g(1):y_s_e_g(2));
            lats_out=lats(x_s_e_g(1):x_s_e_g(2), y_s_e_g(1):y_s_e_g(2));
            
            % spei(time, lat, lon)
            dtmp=squeeze(double(nc_var.data('spei')));            
            dtmp1=dtmp(((yr_s-data_sy)*12+1):(yr_e-data_sy+1)*12,y_s_e_g(1):y_s_e_g(2),x_s_e_g(1):x_s_e_g(2));
            clearvars dtmp;
            
            dtmp2=reshape(dtmp1,[12 nyr dmn_lat_n_g dmn_lon_n_g]);
            clearvars dtmp1;
            
            dtmp3=permute(dtmp2,[3 4 1 2]);
            clearvars dtmp2;
            
            % adjust to the right rotation
            for y=1:nyr
                for m=1:12
                    data_out(:,:,m,y) = rot90(fliplr((squeeze(dtmp3(:,:,m,y)))));
                end
            end
            clearvars dtmp3;
        end

        clearvars dtmp date_fmt dir_in fname;
    end
end