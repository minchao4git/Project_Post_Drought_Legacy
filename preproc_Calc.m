function [] = preproc_Calc(temp_int,pwi,use_gs, chosen_prd, chosen_prd_spei, ptype)

    global DATA_VI_05rs_out DATA_ERA5_05rs_out DATA_Evap_05rs_out DATA_SMRoot_out;
    global DATA_DGVM_05rs_out;
    global DATA_SPEI_out;
    global sos_map lgs_map;
    global mask_gs;
    
    % Output
    global DATA_05rs_out_gs DATA_05rs_dtds_out_gs;
    global DATA_05rs_dtdspw_out;
    global DATA_05rs_dtds_out_gs_inc;
    
    if strcmp(ptype,'SAT')
    
        VID=2; % use the EVI GS map

        % vegetation data and local climate data
        [s1 s2 nm ny s5]=size(DATA_VI_05rs_out);
        datain=nan(s1,s2,nm,ny,2); % dimension 5: 1:2, VI data, 3: SPEI data, 4: ET data

        datain(:,:,:,:,1:2)=squeeze(DATA_VI_05rs_out(:,:,:,:,1:2)); % 1=NDVI, 2=EVI
        datain(:,:,:,:,3)=DATA_SPEI_out(:,:,:,(chosen_prd(1)-chosen_prd_spei(1)+1):(chosen_prd(end)-chosen_prd_spei(1)+1));
        datain(:,:,:,:,4)=DATA_ERA5_05rs_out(:,:,:,:,1); % t2m
        datain(:,:,:,:,5)=DATA_ERA5_05rs_out(:,:,:,:,2); % ssrd
        datain(:,:,:,:,6)=DATA_SMRoot_out; % soil moisture
        datain(:,:,:,:,7)=DATA_Evap_05rs_out; % the first 4 dimension should be the same as the one for VI
        nv=size(datain,5);  % Number of variables to be processed
        datain=datain.*repmat(squeeze(mask_gs(:,:,1)),[1 1 nm ny nv]); % mask_gs is not used here.
        
    elseif strcmp(ptype,'DGVM')
        
        VID=1; % use the DGVMs GS map

        % vegetation data and local climate data
        [s1 s2 nm ny]=size(DATA_DGVM_05rs_out);
        datain=nan(s1,s2,nm,ny,2); % dimension 5: 1, DGVM, 2: SPEI data

        datain(:,:,:,:,1:(nv-1))=DATA_DGVM_05rs_out; %
        datain(:,:,:,:,nv)=DATA_SPEI_out(:,:,:,(chosen_prd(1)-chosen_prd_spei(1)+1):(chosen_prd(end)-chosen_prd_spei(1)+1));
        nv=size(datain,5);  % Number of variables to be processed
        datain=datain.*repmat(squeeze(mask_gs(:,:,1)),[1 1 nm ny nv]); % mask_gs is not used here.
        
    end
    
    a=reshape(datain,[s1 s2 nm*ny nv]);
    a=movmean(a,2,3); % Windows size: 2 months
    datain=reshape(a,[s1 s2 nm ny nv]);
    clearvars a;
    
    % ====     Get growing season absolute values         ====
    % ====  based on the gridwise-defined growing season  ====
    if use_gs
        [DATA_05rs_out_gs]=convert2gs(datain, sos_map(:,:,VID), lgs_map(:,:,VID));
    end

    % ====     Get growing season anomalies                  ====
    % ====  based on the gridwise-defined growing season     ====
    % ====  Detrend and deseasonalize growing season dataset ====
    dtmp_dtds=datain;
    if strcmp(temp_int,'monthly_dsclim')
        
        % Only process the selected variable
        ds_chosen=[1 2 4:7]; % not including SPEI
        for var = ds_chosen
            
            fprintf(sprintf('Processing var %d ...\n',var));
            
            % lon/lat lon/lat yr ds
            dtmp_movm=movmean(squeeze(nanmean(datain(:,:,:,:,var),3)),5, 3); % 5-year running mean of annual mean

            % Detrend the monthly data with the detrended annual mean
            atmp=repmat(dtmp_movm,[1 1 1 12]);
            dtmp_dtds(:,:,:,:,var)=dtmp_dtds(:,:,:,:,var)-permute(atmp,[1 2 4 3]);

            % manually deseason by substracting climatology seasonal cycle
            % lon/lat lon/lat m
            climssn=squeeze(nanmean(dtmp_dtds(:,:,:,:,var),4));
            dtmp_dtds(:,:,:,:,var)=dtmp_dtds(:,:,:,:,var)-repmat(climssn,[1 1 1 ny]);
            
            % Prewithening along the interannual direction
            nts=0;
            if pwi==1
                nhi=2;    % Highest AR order to consider
                k2=[1 2]; % set k2(1):1 to fit model from 1 to the defined nhi, 
                          % k2(2): 2 accept the lowest AIC

                for i=1:s1
                    for j=1:s2
                        for m=1:12
                            ts_tmp=squeeze(dtmp_dtds(i,j,m,:,var));
                            
                            if sum(isnan(ts_tmp))==0
                                nts=nts+1;
                                fprintf(sprintf('Prewhitening DS: %d, %d time series data ... \n',var, nts));
                                [DATA_05rs_dtdspw_out(i,j,m,:,var),dump1,dump2,dump3] = whit1(ts_tmp,nhi,k2);
                            end
                        end
                    end
                end
            end % Prewithening
        end
    end
    clearvars datain;
    
    % --- Calculate incremental anomalies ---
    if 1==2
    tmp_ts=reshape(dtmp_dtds, [s1 s2 nm*ny nv]);
    tmp_ts=tmp_ts(:,:,2:end,:)-tmp_ts(:,:,1:(end-1),:);
    
    dtmp_dtds_ts=nan(s1,s2,nm*ny,nv);
    dtmp_dtds_ts(:,:,2:end,:)=tmp_ts;
    clearvars tmp_ts;
    
    dtmp_dtds_ts_inc=reshape(dtmp_dtds_ts,[s1 s2 nm, ny nv]);
    clearvars dtmp_dtds_ts;
    
    [DATA_05rs_dtds_out_gs_inc]=convert2gs(dtmp_dtds_ts_inc, sos_map(:,:,VID), lgs_map(:,:,VID));
    clearvars dtmp_dtds_ts_inc;
    
    % Standardize VI anomalies with the annual mean
    atemp=repmat(squeeze(nanmax(nanmean(DATA_05rs_out_gs(:,:,:,:,1:(nv-1)),4),[],3)),[1 1 1 nm ny-1]);
    btemp=permute(atemp,[1 2 4 5 3]);
    DATA_05rs_dtds_out_gs_inc(:,:,:,:,1:(nv-1))=DATA_05rs_dtds_out_gs_inc(:,:,:,:,1:(nv-1))./btemp;
    clearvars atemp btemp;
    end
    
    % ---- Get growing season data based on the gridwise-defined growing season
    if use_gs
        fprintf('Standardizing data ...\n');
        
        [DATA_05rs_dtds_out_gs]=convert2gs(dtmp_dtds, sos_map(:,:,VID), lgs_map(:,:,VID));
        
        % VI: standardize anomalies by with the annual mean with the
        % maximum values of the climatology
        atemp=repmat(squeeze(nanmax(nanmean(DATA_05rs_out_gs,4),[],3)),[1 1 1 nm ny-1]);
        max_clm=permute(atemp,[1 2 4 5 3]);
        % convert to %
        ds_chosen_VI=1:2;
        DATA_05rs_dtds_out_gs(:,:,:,:,ds_chosen_VI)=DATA_05rs_dtds_out_gs(:,:,:,:,ds_chosen_VI)./max_clm(:,:,:,:,ds_chosen_VI)*100;
        clearvars atemp max_clm;
        
        % Others variables: standardized by standard deviation (e.g., z-score)
        atmp=reshape(DATA_05rs_dtds_out_gs,[s1 s2 nm*(ny-1) nv]);
        atmp_std=squeeze(std(atmp,0,3,'omitnan'));
        atmp_std=permute(repmat(atmp_std,[1 1 1 nm ny-1]),[1 2 4 5 3]);
        % convert to z-scores
        ds_chosen_clim=3:7;
        DATA_05rs_dtds_out_gs(:,:,:,:,ds_chosen_clim)=DATA_05rs_dtds_out_gs(:,:,:,:,ds_chosen_clim)./atmp_std(:,:,:,:,ds_chosen_clim);
        clearvars atemp atmp_std;
        
    end
    clearvars dtmp climssn dtmp_movm atmp;
    
end

